#!/usr/bin/env python3

##########################################################################
# parse_swarm_v2.py  –  Optimized SWARM parser                          #
# Anders Lanzen 2024 / optimized 2026                                    #
#                                                                        #
# Changes vs v1:                                                         #
#  1. INVERTED INDEX  isu → {sample: count}  (main speedup)             #
#     Eliminates O(N_samples) inner loop; now only touches samples       #
#     that actually contain each ISU.                                    #
#  2. Larger I/O buffers (64 MB) – helps on slow / NFS mounts.          #
#  3. Batch-flushed output – fewer write syscalls.                       #
#  4. Optional multiprocessing (N_WORKERS > 1, Linux/macOS only).       #
##########################################################################

import gc
import time
import multiprocessing as mp

# ── Configuration ─────────────────────────────────────────────────────────
MIN_AB      = 2
N_WORKERS   = 1          # set > 1 to enable multiprocessing (Linux/macOS)
CHUNK_LINES = 5_000      # swarm lines per worker chunk
IO_BUF      = 64 << 20  # 64 MB I/O buffer
FLUSH_EVERY = 20_000     # flush output every N rows
# ─────────────────────────────────────────────────────────────────────────

t0 = time.time()
def ts(): return f"[{time.time()-t0:6.0f}s]"


###########################################################################
# 1.  Load origins.tsv with INVERTED index: isu → {sample: count}
###########################################################################
# Original: samples[s][isu]  → must iterate ALL samples per ISU per swarm
#           O(N_samples × N_isus_per_swarm × N_swarms)
#
# New:      isu_samples[isu][s]  → only visits samples that have this ISU
#           O(avg_samples_per_ISU × N_isus_per_swarm × N_swarms)
#           typically 10–100× fewer operations for sparse eDNA data
###########################################################################

print(f"{ts()} Reading origins.tsv …")

isu_samples  = {}    # isu_name → {sample_name: abundance}
sample_order = []    # preserves column order for OTU table header
sample_seen  = set()

with open("origins.tsv", "r", buffering=IO_BUF) as fh:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        isu   = parts[0].split(";")[0]
        d     = {}
        for oc in parts[1:]:
            tok = oc.split(";")
            s   = tok[0]
            if s not in sample_seen:
                sample_seen.add(s)
                sample_order.append(s)
            d[s] = int(tok[1][5:])    # "size=N" → N
        isu_samples[isu] = d

all_samples = sample_order
n_samples   = len(all_samples)
sample_idx  = {s: i for i, s in enumerate(all_samples)}

print(f"{ts()} {len(isu_samples):,} ISUs loaded across {n_samples} samples.")
gc.collect()


###########################################################################
# 2.  Core per-swarm function (used by both single- and multi-threaded)
###########################################################################

def _process_swarm(i, line):
    """
    Process one line from SWARM.swarms.
    Returns (seed_isu, swarm_name, counts) or (seed_isu, None, None).
    Reads isu_samples / sample_idx / n_samples / MIN_AB as module globals
    — on Linux these are inherited by fork workers without pickling.
    """
    isus     = line.rstrip("\n").split(" ")
    main     = isus[0].split(";")
    seed_isu = main[0]
    main_ab  = int(main[1][5:])   # "size=N"

    # Same singleton / low-abundance filter as original
    if len(isus) < MIN_AB and main_ab < MIN_AB:
        return seed_isu, None, None

    swarm_name = f"SWARM_{i}"
    counts     = [0] * n_samples

    for tok in isus:
        isu_name = tok.split(";")[0]
        hit = isu_samples.get(isu_name)   # None if absent — no loop needed
        if hit:
            for s, cnt in hit.items():
                counts[sample_idx[s]] += cnt

    return seed_isu, swarm_name, counts


def _process_chunk(chunk):
    """Multiprocessing worker entry point. chunk = [(i, line), …]"""
    return [_process_swarm(i, line) for i, line in chunk]


###########################################################################
# 3.  Process SWARM.swarms → SWARM_table.tsv
###########################################################################
print(f"{ts()} Processing SWARM.swarms …")

otus      = {}   # seed_isu → swarm_name  (needed for FASTA step)
row_buf   = []
n_written = 0

with open("SWARM_table.tsv", "w", buffering=IO_BUF) as otu_out:
    otu_out.write("SWARM\t" + "\t".join(all_samples) + "\n")

    # ── single-threaded path ─────────────────────────────────────────────
    if N_WORKERS <= 1:
        with open("SWARM.swarms", "r", buffering=IO_BUF) as swarms:
            for i, line in enumerate(swarms, 1):
                seed_isu, swarm_name, counts = _process_swarm(i, line)

                if swarm_name is None:
                    continue

                otus[seed_isu] = swarm_name
                row_buf.append(
                    swarm_name + "\t" + "\t".join(map(str, counts)) + "\n")
                n_written += 1

                # Reclaim RAM as we go — same strategy as original
                for tok in line.rstrip("\n").split(" "):
                    isu_samples.pop(tok.split(";")[0], None)

                if len(row_buf) >= FLUSH_EVERY:
                    otu_out.writelines(row_buf)
                    row_buf.clear()

                if i % 500_000 == 0:
                    print(f"{ts()}   {i:,} swarms, {n_written:,} OTUs written")

    # ── multiprocessing path (Linux / macOS — fork semantics) ────────────
    # Workers inherit isu_samples via fork (copy-on-write).
    # We skip del/pop here to avoid COW pressure on every memory page;
    # this trades slightly higher peak RAM for parallel throughput.
    # NOTE: this will NOT work on Windows (spawn-based MP).
    else:
        print(f"{ts()}   {N_WORKERS} workers, chunk size {CHUNK_LINES:,} lines")

        def _chunks():
            with open("SWARM.swarms", "r", buffering=IO_BUF) as fh:
                buf = []
                for i, line in enumerate(fh, 1):
                    buf.append((i, line))
                    if len(buf) == CHUNK_LINES:
                        yield buf
                        buf = []
                if buf:
                    yield buf

        with mp.Pool(N_WORKERS) as pool:
            for chunk_n, results in enumerate(
                    pool.imap(_process_chunk, _chunks(), chunksize=1), 1):

                for seed_isu, swarm_name, counts in results:
                    if swarm_name is None:
                        continue
                    otus[seed_isu] = swarm_name
                    row_buf.append(
                        swarm_name + "\t" + "\t".join(map(str, counts)) + "\n")
                    n_written += 1

                if len(row_buf) >= FLUSH_EVERY:
                    otu_out.writelines(row_buf)
                    row_buf.clear()

                if chunk_n % 200 == 0:
                    print(f"{ts()}   chunk {chunk_n:,}, {n_written:,} OTUs")

    # Final flush
    if row_buf:
        otu_out.writelines(row_buf)
        row_buf.clear()

print(f"{ts()} OTU table complete — {n_written:,} OTUs written.")
del isu_samples
gc.collect()


###########################################################################
# 4.  Filter representative FASTA
###########################################################################
print(f"{ts()} Writing SWARM_OTUs_f.fasta …")

with open("SWARM_OTUs.fasta",   "r", buffering=IO_BUF) as fi, \
     open("SWARM_OTUs_f.fasta", "w", buffering=IO_BUF) as fo:
    included = False
    for line in fi:
        if line[0] == ">":
            isu = line[1:].split(";")[0]
            if isu in otus:
                included = True
                size = line.split(";")[1][5:].rstrip()   # strip trailing \n
                fo.write(f">{otus[isu]};size={size};\n")
            else:
                included = False
        elif included:
            fo.write(line)

print(f"{ts()} Done.")

