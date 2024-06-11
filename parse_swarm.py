#!/usr/bin/env python

#################### parse_swarm.py SWARM parser #########################
# Parses swarm output in its own format (SWARM.swarms) and using the     #
# FASTA output of SWARM (SWARM_OTUs.fasta) together with origin.tsv from #
# SLIM script. Replaces matrix_creation that did not work with SWARM v3  #
# Removes singletons by default but after SWARM parsing (not in origins) #
#                                                                        #
# Anders Lanzen 2024-06-04                                               #
##########################################################################

##
## Note: memory issues experienced with full MicroMon dataset. Origins being
## 2 Gb and SWARM.swarm 700M does not explain dict consumption (>200 Gb RAM)
## Inserted del dict[key] statements to save memory and some text output

import gc

minAb = 2

origins = open("origins.tsv","r")
swarms = open("SWARM.swarms","r")

### Read all data in origins.tsv to get sample distro of each unique sequence pre-SWARM

samples = {}
otus = {}

print("Reading origins.tsv data into memory to build OTU table..")

for line in origins:
    
    parts = line.rstrip().split("\t")
    isu = parts[0].split(";")[0]
    occurences = parts[1:]

    for oc in occurences:        
        #ISU_0;size=7118	S-Fer1;size=579	S-Fer2;size=636
        sample_reads=oc.split(";")
        s = sample_reads[0]
        if not s in samples:
            samples[s] = {}
        samples[s][isu]=int(sample_reads[1][5:]) #"size="


origins.close()

### Look up occurences for each ISU in swarm output and write to std out

# List of all samples
all_samples = list(samples.keys())


# All accumulated SWARM samples abundances, with first key being ISU name instead of sample as above
# (later allowing easier writing of OTU table)

print("Reading SWARM.swarms, building sample-wise abundances using data from origins.tsv and writing directly to OTU table...")

#swarm_isus = {}

### Write OTU table
otu_out=open("SWARM_table.tsv","w")
otu_out.write("SWARM\t")
otu_out.write("\t".join([s for s in all_samples]))
otu_out.write("\n")

print("%s samples:" % len(all_samples)) #DEBUG
print(", ".join([s for s in all_samples])) # DEBUG

i=0
for line in swarms:
    i=i+1
    isus = line.split(" ")

    main_reads = isus[0].split(";")
    swarm_isu=main_reads[0]
    main_ab=int(main_reads[1][5:]) #size=

    # Verify that this is not a singleton (or below min abundance)
    if (len(isus)>=minAb or main_ab>=minAb):

        #swarm_isus[swarm_isu]={s:0 for s in all_samples}
        swarm_isus = {s:0 for s in all_samples}

        # Make SWARM_n name and store
        swarm_name="SWARM_"+str(i)
        otus[swarm_isu]=swarm_name
        
        # iterate over all ISUs including first
        for isu in isus:
            isuName = isu.split(";")[0]
            # Iterate over all samples
            for s in all_samples:
                # Look up and set sReads as abundance in samples s
                if isuName in samples[s]:
                    sReads = samples[s][isuName]
                    # Depopulate samples accordingly dict to save up RAM
                    del samples[s][isuName]                                               
                    
                else:
                    sReads = 0
                
                # Add reads to SWARM for the present sample
                swarm_isus[s] = swarm_isus[s] + sReads
                
        # Write to OTU table
        otu_out.write(swarm_name+"\t")
        otu_out.write("\t".join([str(r) for r in swarm_isus.values()]))
        otu_out.write("\n")
        print("DEBUG: %s : %s, %s reads" % (swarm_isu, swarm_name, sum(swarm_isus.values())))
                
swarms.close()
otu_out.close()

print("Writing representative OTU sequences to SWARM_OTUs_f.fasta")

### Write OTU sequences in FASTA
fastaIn = open("SWARM_OTUs.fasta","r")
fastaOut = open("SWARM_OTUs_f.fasta","w")

for line in fastaIn:
    if line[0]==">":
        isu=line[1:].split(";")[0]
        if isu in otus:
            included=True
            # Singletons will not be included and are to be ignored
            #size=sum(swarm_isus[isu].values())
            size = line.split(";")[1][5:]
            name=otus[isu]
            fastaOut.write(">"+name+";size="+str(size)+";\n")
        else:
            included=False
    else:
        if included:
            fastaOut.write(line)

fastaIn.close()
fastaOut.close()
