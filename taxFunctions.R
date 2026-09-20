# Accepts a CREST4 assignments.txt file and makes a taxon_tableonomy data frame
## with rank specific info. Fixes some names with special symbols
## if domain_only==TRUE set all entries without domain to No hits
# v 0.5 updated with ARMS results (COI) to separate macrobenthos and zooplankton
# (has seen Metabridge, MicroMon and IndiRed 16S and 18S taxa, ARMS COI)

makeTaxonomy = function(crest4_assignment_file, domain_only=T) {
  taxon_table = read.table(crest4_assignment_file, sep="\t",row.names=1, quote = "\"")
  names(taxon_table) = "classification"
  
  taxon_table$classification = gsub("\\[","",taxon_table$classification)
  taxon_table$classification = gsub("\\]","",taxon_table$classification)
  
  ## Assumes OTU list in numerical order ;size=N tar present
  ## Remvove ";size=N" from row names if "size" in first row name
  if(grepl("size", row.names(taxon_table)[1])){
    row.names(taxon_table)[1:9] = substring(row.names(taxon_table)[1:9], 1,5)
    if(dim(taxon_table)[1]>9) {
      endRange = min(dim(taxon_table)[1], 99)
      row.names(taxon_table)[10:endRange] = substring(row.names(taxon_table)[10:endRange], 1,6)
    }
    if(dim(taxon_table)[1]>99) {
      endRange = min(dim(taxon_table)[1], 999)
      row.names(taxon_table)[100:endRange] = substring(row.names(taxon_table)[100:endRange], 1,7)
    }
    if(dim(taxon_table)[1]>999) {
      endRange = min(dim(taxon_table)[1], 9999)
      row.names(taxon_table)[1000:endRange] = substring(row.names(taxon_table)[1000:endRange], 1,8)
    }
    if(dim(taxon_table)[1]>9999) {
      endRange = min(dim(taxon_table)[1], 99999)
      row.names(taxon_table)[10000:endRange] = substring(row.names(taxon_table)[10000:endRange], 1,9)
    }
    if(dim(taxon_table)[1]>99999) {
      endRange = min(dim(taxon_table)[1], 999999)
      row.names(taxon_table)[100000:endRange] = substring(row.names(taxon_table)[100000:endRange], 1,10)
    }
  }
  
  ## Make rank specific or special taxonomy variables
  
  bestTx = array(dim=dim(taxon_table)[1])
  meta = array(dim=dim(taxon_table)[1])
  domain = array(dim=dim(taxon_table)[1])
  superkingdom = array(dim=dim(taxon_table)[1])
  kingdom = array(dim=dim(taxon_table)[1])
  phylum = array(dim=dim(taxon_table)[1])
  class = array(dim=dim(taxon_table)[1])
  orden = array(dim=dim(taxon_table)[1])
  fam = array(dim=dim(taxon_table)[1])
  genus = array(dim=dim(taxon_table)[1])
  sp = array(dim=dim(taxon_table)[1])
  for (i in 1:length(bestTx)) {
    
    # Make taxon_tablea path 
    taxon_tableapath = unlist(strsplit(as.character(taxon_table$classification[i]), split="; ", fixed=TRUE))
    
    bestTx[i] <- tail(taxon_tableapath,1)
    
    if (length(taxon_tableapath)>1) meta[i] = taxon_tableapath[2]
    else meta[i] = NA
    
    if (length(taxon_tableapath)>2) domain[i] = taxon_tableapath[3]
    else domain[i] = NA
    
    if (length(taxon_tableapath)>3) superkingdom[i] = taxon_tableapath[4]
    else superkingdom[i] = NA
    
    if (length(taxon_tableapath)>4) kingdom[i] = taxon_tableapath[5]
    else kingdom[i] = NA
    
    if (length(taxon_tableapath)>5) phylum[i] = taxon_tableapath[6]
    else phylum[i] = NA
    
    if (length(taxon_tableapath)>6) class[i] = taxon_tableapath[7]
    else class[i] = NA
    
    if (length(taxon_tableapath)>7) orden[i] = taxon_tableapath[8]
    else orden[i] = NA
    
    if (length(taxon_tableapath)>8) fam[i] = taxon_tableapath[9]
    else fam[i] = NA
    
    if (length(taxon_tableapath)>9) genus[i] = taxon_tableapath[10]
    else genus[i] = NA
    
    if (length(taxon_tableapath)>10) sp[i] = taxon_tableapath[11]
    else sp[i] = NA
    
  }
  
  if(domain_only) bestTx[is.na(domain)] = "No hits"
  
  ## Fix some entries - bestTx
  bestTx = gsub(" \\(superkingdom\\)","", bestTx)
  
  ## Fix some entries - domain
  domain = gsub(" \\(Chloroplast\\)","",domain)
  domain = gsub(" \\(Mitochondria\\)","",domain)
  
  ## Fix some entries - kingdom
  kingdom = gsub(" \\(kingdom\\)","",kingdom)
  kingdom = gsub(" kingdom incertae sedis","",kingdom)
  kingdom = gsub(" X$","",kingdom)
  
  ## Fix some entries - phylum
  phylum = gsub(" \\(Metazoa\\)","",phylum)
  phylum = gsub(" phylum incertae sedis","",phylum)
  phylum = gsub(" \\(phylum\\)","",phylum)
  phylum = gsub(" \\(Lophotrochozoa\\)","",phylum)
  phylum = gsub(" \\(Dikarya\\)","",phylum)
  phylum = gsub(" X$","",phylum)
  
  
  # sp. name before species
  phylum = gsub(" sp.","",phylum)
  class = gsub(" sp.","",class)
  orden = gsub(" sp.","",orden)
  fam = gsub(" sp.","",fam)
  genus = gsub(" sp.","",genus)
  
  # Factorize
  taxon_table$bestTx = as.factor(bestTx)
  taxon_table$meta = as.factor(meta)
  taxon_table$domain = as.factor(domain)
  taxon_table$superkingdom = as.factor(superkingdom)
  taxon_table$kingdom = as.factor(kingdom)
  taxon_table$phylum = as.factor(phylum)
  taxon_table$class = as.factor(class)
  taxon_table$order = as.factor(orden)
  taxon_table$fam = as.factor(fam)
  taxon_table$genus = as.factor(genus)
  taxon_table$species = as.factor(sp)
  
  ## Make category property and choose based on some heuristic criteria
  ## (experimental, and unable to distinguish e.g. phytoplankton from benthic autotrophs)
  fc = array(data=NA, dim=dim(taxon_table)[1])
  
  # ------- Fungi -----------
  fc[!is.na(kingdom) & kingdom=="Fungi"] = "Fungus-like"
  fc[!is.na(kingdom) & kingdom=="Aphelida"] = "Fungus-like"
  fc[!is.na(phylum) & phylum=="Hyphochytriomycota"] = "Fungus-like"
  fc[!is.na(phylum) & phylum=="Oomycota"] = "Fungus-like"
  fc[!is.na(phylum) & phylum=="Labyrinthulidia"] = "Fungus-like"
  
  
  # --- Meio-macrobenthos or zooplankton -------
  fc[!is.na(kingdom) & kingdom=="Metazoa"] = "Macrobenthos, zooplankton or biphasic"
  
  # Annelids are mainly benthic
  fc[!is.na(phylum) & phylum=="Annelida"] = "Macrobenthos"
  # ..with exceptions
  fc[!is.na(fam) & fam=="Tomopteridae"] = "Macrozooplankton" # Gossamer worms
  fc[!is.na(fam) & fam=="Alciopidae"] = "Macrozooplankton" # Camera eyed worms
  fc[!is.na(genus) & genus=="Swima"] = "Macrozooplankton" 
  fc[!is.na(genus) & genus=="Palola"] = "Macrozooplankton" 
  
  # Cnidaria are very complex, most benthic
  fc[!is.na(phylum) & phylum=="Cnidaria"] = "Macrobenthos"
  # ..but these weird anemones are not
  fc[!is.na(fam) & fam=="Boloceroididae"] = "Nektobenthos"
  
  # ..and jellyfish corals etc are biphasic
  fc[!is.na(class) & class=="Scyphozoa"] = "Biphasic"
  fc[!is.na(class) & class=="Cubozoa"] = "Biphasic"
  fc[!is.na(class) & class=="Hydrozoa"] = "Biphasic"
  # ..more weird exceptions
  fc[!is.na(genus) & genus=="Cassiopea"] = "Macrobenthos"  # Upside down jellies
  fc[!is.na(orden) & orden=="Siphonophorae"] = "Macrozooplankton" #Siphonopgores
  fc[!is.na(fam) & fam=="Hydridae"] = "Macrobenthos" # Hydra
  fc[!is.na(fam) & fam=="Tubulariidae"] = "Macrobenthos"
  fc[!is.na(genus) & genus=="Coryne"] = "Macrobenthos" 
  fc[!is.na(genus) & genus=="Nannocoryne"] = "Macrobenthos" 
  fc[!is.na(genus) & genus=="Bicorona"] = "Macrobenthos"
  fc[!is.na(genus) & genus=="Bimeria"] = "Macrobenthos"
  fc[!is.na(genus) & genus=="Pachycordyle"] = "Macrobenthos"
  fc[!is.na(genus) & genus=="Rhizorhagium"] = "Macrobenthos"
  fc[!is.na(genus) & genus=="Rhizogeton"] = "Macrobenthos"
  fc[!is.na(genus) & genus=="Merona"] = "Macrobenthos"
  fc[!is.na(genus) & genus=="Corydendrium"] = "Macrobenthos"
  fc[!is.na(fam) & fam=="Eudendriidae"] = "Macrobenthos" 
  fc[!is.na(fam) & fam=="Hydractiniidae"] = "Macrobenthos"
  fc[!is.na(fam) & fam=="Sertulariidae"] = "Macrobenthos" 
  fc[!is.na(fam) & fam=="Plumulariidae"] = "Macrobenthos" 
  fc[!is.na(fam) & fam=="Aglaopheniidae"] = "Macrobenthos" 
  fc[!is.na(fam) & fam=="Haleciidae"] = "Macrobenthos" 
  fc[!is.na(fam) & fam=="Halopterididae"] = "Macrobenthos" 
  fc[!is.na(fam) & fam=="Kirchenpaueriidae"] = "Macrobenthos"
  fc[!is.na(fam) & fam=="Plumulariidae"] = "Macrobenthos" 
  fc[!is.na(fam) & fam=="Sertulariidae"] = "Macrobenthos" 
  fc[!is.na(fam) & fam=="Kirchenpaueriidae"] = "Macrobenthos" 
  
  fc[!is.na(class) & class=="Staurozoa"] = "Macrobenthos"

  # Mollusks 
  fc[!is.na(phylum) & phylum=="Mollusca"] = "Macrozooplankton"
  
  fc[!is.na(phylum) & class=="Cephalopoda"] = "Nektobenthos" # most cephalopods are Nektobenthic
  fc[!is.na(genus) & genus=="Japetella"] = "Macrozooplankton" # Small planktonic octupii
  
  # Arhropoda are mainly mesozooplankton
  fc[!is.na(phylum) & phylum=="Arthropoda"] = "Mesozooplankton"
 
  fc[!is.na(orden) & orden=="Harpacticoida"] = "Meiobenthos"
  fc[!is.na(fam) & fam=="Miraciidae"] = "Hyperbenthos"
  fc[!is.na(fam) & fam=="Porcellidiidae"] = "Macrobenthos"
  fc[!is.na(fam) & fam=="Ectinosomatidae"] = "Hyperbenthos"
  fc[!is.na(fam) & fam=="Tisbidae"] = "Hyperbenthos"
  fc[!is.na(fam) & fam=="Tegastidae"] = "Hyperbenthos"
  
  fc[!is.na(fam) & fam=="Cyclopinidae"] = "Hyperbenthos"
  fc[!is.na(fam) & fam=="Schminkepinellidae"] = "Hyperbenthos"
  fc[!is.na(fam) & fam=="Cyclopinidae"] = "Hyperbenthos"
  
  fc[!is.na(genus) & genus=="Microsetella"] = "Mesozooplankton"
  
  # Malacostraca most are benthic and large Pycnogonida
  fc[!is.na(orden) & orden=="Euphausiacea"] = "Macrozooplankton" # krill
  fc[!is.na(fam) & fam=="Hyperiidae"] = "Macrozooplankton" # predatory or parasitic on jellies
  fc[!is.na(orden) & orden=="Mysida"] = "Hyperbenthos"
  fc[!is.na(orden) & orden=="Stomatopoda"] = "Biphasic" # Mantis shrimps
  
  fc[!is.na(orden) & orden=="Podocopida"] = "Mesozooplankton"
  fc[!is.na(orden) & orden=="Euphausiacea"] = "Macrozooplankton" # krill
  
  fc[!is.na(class) & class=="Pycnogonida"] = "Macrobenthos"
  fc[!is.na(class) & class=="Thecostraca"] = "Biphasic"
  
  
  fc[!is.na(phylum) & phylum=="Bryozoa"] = "Macrobenthos"
  
  # Tunicates
  fc[!is.na(class) & class=="Ascidiacea"] = "Macrobenthos"
  
  # Echinoderms
  fc[!is.na(phylum) & phylum=="Echinodermata"] = "Biphasic"
  fc[!is.na(genus) & genus=="Leptasterias"] = "Macrobenthos"
  
  fc[!is.na(phylum) & phylum=="Gastrotricha"] = "Meiobenthos"
  
  fc[!is.na(phylum) & phylum=="Nematoda"] = "Meiobenthos"
  fc[!is.na(fam) & fam=="Oncholaimidae"] = "Macrobenthos" 
  fc[!is.na(fam) & fam=="Enchelidiidae"] = "Macrobenthos" 
  
  fc[!is.na(phylum) & phylum=="Nemertea"] = "Meiobenthos"
  fc[!is.na(class) & class=="Pilidiophora"] = "Biphasic"
  fc[!is.na(genus) & genus=="Pelagonemertes"] = "Mesozooplankton"
  
  fc[!is.na(phylum) & phylum=="Phoronida"] = "Biphasic"
  
  fc[!is.na(phylum) & phylum=="Platyhelminthes"] = "Meiobenthos"
  
  fc[!is.na(phylum) & phylum=="Porifera"] = "Meiobenthos"
  
  fc[!is.na(phylum) & phylum=="Rotifera"] = "Meiobenthos"
  fc[!is.na(class) & class=="Monogononta"] = "Mesozooplankton"
  fc[!is.na(class) & class=="Eurotatoria"] = "Mesozooplankton"
  
  
  
  
  
  # ---- Nekton -----
  fc[!is.na(class) & class=="Actinopterygii"] = "Nekton"
  fc[!is.na(class) & class=="Actinopteri"] = "Nekton"
  fc[!is.na(class) & class=="Agnatha"] = "Nekton"
  fc[!is.na(class) & class=="Myxini"] = "Nekton"
  fc[!is.na(class) & class=="Chondrichthyes"] = "Nektobenthos" # Cjimera
  fc[!is.na(orden) & orden=="Pleuronectiformes"] = "Nektobenthos" # Flatfish
  
  # Some cephalopods are nekton
  fc[!is.na(genus) & genus=="Amphitretus"] = "Nekton" #  telescope octopus 
  fc[!is.na(genus) & genus=="Tremoctopus"] = "Nekton"
  fc[!is.na(genus) & genus=="Argonauta"] = "Nekton"
  fc[!is.na(genus) & genus=="Tremoctopus"] = "Nekton" # Blanket occi
  fc[!is.na(orden) & orden=="Teuthida"] = "Nekton" # Squid
  fc[!is.na(orden) & orden=="Oegopsida"] = "Nekton" # Squid (super-order)
  fc[!is.na(orden) & orden=="Spirulida"] = "Nekton" # Ram's Horn Squid
  fc[!is.na(fam) & fam=="Onychoteuthidae"] = "Nekton" # squid
  fc[!is.na(fam) & fam=="Nautilidae"] = "Nekton" # 
  fc[!is.na(fam) & fam=="Argonautidae"] = "Nekton" # 
  fc[!is.na(orden) & orden=="Teuthida"] = "Nekton" # True squid
  
   
  
  
  # Xeronekton that breath air and live partly on
  fc[!is.na(orden) & orden=="Testudines"] = "Xeronekton" # Turtles
  fc[!is.na(orden) & orden=="Odobenidae"] = "Xeronekton" # Walrus
  fc[!is.na(orden) & orden=="Otariidae"] = "Xeronekton" # Fur seals
  fc[!is.na(orden) & orden=="Phocidae"] = "Xeronekton" # True earless eals
  #fc[!is.na(orden) & orden=="Odobenidae"] = "Xeronekton" # Turtles
  
  
  # ------ Protists (plankton or microbenthos, not photosynthesizing) -------
  fc[!is.na(kingdom) & kingdom=="Amoebozoa"] = "Protists"
  fc[!is.na(kingdom) & kingdom=="Breviatae"] = "Protists" # Eukaryome
  fc[!is.na(kingdom) & kingdom=="Apusozoa"] = "Protists" # Surface associated
  fc[!is.na(kingdom) & kingdom=="Apusoozoa"] = "Protists" # Surface associated
  fc[!is.na(kingdom) & kingdom=="Centroheliozoa"] = "Protists" # Surface associated
  fc[!is.na(kingdom) & kingdom=="Choanoflagellida"] = "Protists"
  fc[!is.na(kingdom) & kingdom=="Choanoflagellata"] = "Protists"
  fc[!is.na(kingdom) & kingdom=="Choanoflagellozoa"] = "Protists" # Eukaryome
  fc[!is.na(kingdom) & kingdom=="Corallochytriozoa"] = "Protists" # Eukaryome
  fc[!is.na(orden) & orden=="Goniomonadales"] = "Protists"
  fc[!is.na(kingdom) & kingdom=="Discoba"] = "Protists"
  fc[!is.na(phylum) & phylum=="Katablepharidophyta"] = "Protists"
  fc[!is.na(phylum) & phylum=="Nucleariidea"] = "Protists" 
  fc[!is.na(phylum) & phylum=="Filasterea"] = "Protists"
  fc[!is.na(phylum) & phylum=="Filasteriae"] = "Protists" #Eukaryome
  fc[!is.na(phylum) & phylum=="Picozoa"] = "Protists" # Earlier thought to be phyto
  fc[!is.na(phylum) & phylum=="Protalveolata"] = "Protists (predatory)" # Colponema - Predators
  fc[!is.na(kingdom) & kingdom=="Telonemia"] = "Protists"
  fc[!is.na(kingdom) & kingdom=="Telonemae"] = "Protists"
  fc[!is.na(phylum) & phylum=="Cercozoa"] = "Protists"
  fc[!is.na(phylum) & phylum=="Pseudofungi"] = "Protists"
  fc[!is.na(phylum) & phylum=="Sagenista"] = "Protists"
  fc[!is.na(phylum) & phylum=="Opalozoa"] = "Protists"
  fc[!is.na(phylum) & phylum=="Ciliophora"] = "Protists"
  fc[!is.na(phylum) & phylum=="Radiolaria"] = "Protists"
  fc[!is.na(phylum) & phylum=="Bigyra"] = "Protists"
  fc[!is.na(phylum) & phylum=="Endomyxa"] = "Protists"
  fc[!is.na(phylum) & phylum=="Malawimonadidea"] = "Protists (anaerobic)" # bacteriovores
  fc[!is.na(phylum) & phylum=="Nucleariae"] = "Protists"
  fc[!is.na(kingdom) & kingdom=="Collodictyozoa"] = "Protists"
  #fc[!is.na(kingdom) & kingdom=="Mantamonada"] = "Protists" # some microbiome
  # fc[!is.na(kingdom) & kingdom=="Euglenozoa"] = "Protists" # some microbiome and parasites
  fc[!is.na(phylum) & phylum=="Petalomonadia"] = "Protists"
  fc[!is.na(kingdom) & kingdom=="Jakobae"] = "Protists"
  fc[!is.na(kingdom) & kingdom=="Hemimastigophora"] = "Protists"
  fc[!is.na(kingdom) & kingdom=="Meteorae"] = "Protists"
  fc[!is.na(class) & class=="Carpediemonadea"] = "Protists (anaerobic)"
  fc[!is.na(superkingdom) & superkingdom=="Provora"] = "Protists (predatory)"
  fc[!is.na(phylum) & phylum=="Petalomonadia"] = "Protists (phytobenthic)"
  fc[!is.na(phylum) & phylum=="Provora"] = "Protists (predatory)"
  fc[!is.na(phylum) & phylum=="Gromida"] = "Protists" # Very big saprotrophic amoeba (up to cm)
  fc[!is.na(phylum) & phylum=="Cantinia"] = "Protists (anaerobic)"
  fc[!is.na(genus) & genus=="Cyclidium"] = "Protists (anaerobic)"
  fc[!is.na(phylum) & phylum=="Developia"] = "Protists"
  fc[!is.na(phylum) & phylum=="Picophagophyta"] = "Protists"
  fc[!is.na(phylum) & phylum=="Pirsonida"] = "Protists"
  fc[!is.na(phylum) & phylum=="Placidida"] = "Protists"
  fc[!is.na(phylum) & phylum=="Solenicolida"] = "Protists"
  fc[!is.na(phylum) & phylum=="Placidida"] = "Protists"
  fc[!is.na(class) & class=="Aphagea"] = "Protists"
  fc[!is.na(genus) & genus=="Rapaza"] = "Protists (predatory)" #Kleptoparasitic phagotrophs that have their own plastids!
  
  
  
  
  #Additional from V4
  fc[!is.na(phylum) & phylum=="Ancyromonadida"] = "Protists"
  
  # ----- Mixoplankton (phyto) -------
  fc[!is.na(superkingdom) & superkingdom=="Archaeplastida"] = "Phytoplankton"
  fc[!is.na(superkingdom) & superkingdom=="CAM"] = "Phytoplankton"
  fc[!is.na(superkingdom) & kingdom=="Picozoa"] = "Protists"
  fc[!is.na(phylum) & phylum=="Dinoflagellata"] = "Phytoplankton"
  fc[!is.na(orden) & orden=="Cryptomonadales"] = "Phytoplankton"
  fc[!is.na(phylum) & phylum=="Cryptista ss"] = "Phytoplankton"
  fc[!is.na(orden) & orden=="Euglyphida"] = "Phytoplankton"
  fc[!is.na(phylum) & phylum=="Haptophyta"] = "Phytoplankton"
  fc[!is.na(class) & class=="Chlorarachniophyceae"] = "Phytoplankton"
  fc[!is.na(class) & class=="Eustigmatophyceae"] = "Phytoplankton"
  fc[!is.na(class) & class=="Fragilariophyceae"] = "Phytoplankton"
  fc[!is.na(phylum) & phylum=="Ochrophyta"] = "Phytoplankton"
  fc[!is.na(class) & class=="Raphidophyceae"] = "Phytoplankton"
  fc[!is.na(class) & class=="Mediophyceae"] = "Phytoplankton"
  fc[!is.na(class) & class=="Bacillariophyceae"] = "Phytoplankton"
  fc[!is.na(phylum) & phylum=="Bacillariophyta"] = "Phytoplankton"
  fc[!is.na(class) & class=="Bolidophyceae"] = "Phytoplankton"
  fc[!is.na(class) & class=="Coscinodiscophyceae"] = "Phytoplankton"
  fc[!is.na(class) & class=="Synurophyceae"] = "Phytoplankton"
  fc[!is.na(class) & class=="Xanthophyceae"] = "Phytoplankton"
  fc[!is.na(class) & class=="Chlorodendrophyceae"] = "Phytoplankton"
  fc[!is.na(orden) & orden=="Chlamydomonadales"] = "Phytoplankton"
  fc[!is.na(orden) & orden=="Sphaeropleales"] = "Phytoplankton"
  fc[!is.na(class) & class=="Chloropicophyceae"] = "Phytoplankton"
  fc[!is.na(class) & class=="Mamiellophyceae"] = "Phytoplankton"
  fc[!is.na(class) & class=="Nephroselmidophyceae"] = "Phytoplankton"
  fc[!is.na(class) & class=="Pyramimonadophyceae"] = "Phytoplankton"
  fc[!is.na(class) & class=="Coleochaetophyceae"] = "Phytoplankton (freshwater)" # Colony forming freshwater
  fc[!is.na(class) & class=="Prasinodermophyceae"] = "Phytoplankton (freshwater)" # Colony forming freshwater
  fc[!is.na(phylum) & phylum=="Zygnematophyta"] = "Phytoplankton"
  fc[!is.na(orden) & orden=="Pedospumella"] = "Phytoplankton"
  fc[!is.na(genus) & genus=="Spumella"] = "Phytoplankton"
  fc[!is.na(fam) & fam=="Dinobryaceae"] = "Phytoplankton"
  fc[!is.na(fam) & fam=="Paraphysomonadaceae"] = "Phytoplankton"
  fc[!is.na(fam) & fam=="Synuraceae"] = "Phytoplankton"
  fc[!is.na(phylum) & phylum=="Dictyochophyta"] = "Phytoplankton"
  fc[!is.na(phylum) & phylum=="Eustigmatophyta"] = "Phytoplankton"
  fc[!is.na(phylum) & phylum=="Pelagophyta"] = "Phytoplankton"
  fc[!is.na(phylum) & phylum=="Pinguiophyta"] = "Phytoplankton"
  fc[!is.na(phylum) & phylum=="Pelagophyta"] = "Phytoplankton"
  fc[!is.na(phylum) & phylum=="Synchromophyceae"] = "Phytoplankton"
  fc[!is.na(phylum) & phylum=="Synchromophyceae"] = "Phytoplankton"
  fc[!is.na(orden) & orden=="Chromulinales"] = "Phytoplankton"
  fc[!is.na(class) & class=="Chrysophyceae"] = "Phytoplankton"
  
  
  
  # -----  Macroalgae ------
  fc[!is.na(class) & class=="Florideophyceae"] = "Macroalgae"
  fc[!is.na(class) & class=="Bangiophyceae"] = "Macroalgae"
  fc[!is.na(class) & class=="Ulvophyceae"] = "Macroalgae"  # Seaweeds
  fc[!is.na(class) & class=="Phaeophyceae"] = "Macroalgae" # Brown algae
  fc[!is.na(fam) & fam=="Stylonematophyceae"] = "Macroalgae"
  fc[!is.na(fam) & fam=="Compsopogonophyceae"] = "Macroalgae"
  fc[!is.na(orden) & orden=="Chaetophorales"] = "Macroalgae (freshwater)" #Freshwater
  fc[!is.na(orden) & orden=="Oedogoniales"] = "Macroalgae (freshwater)" #Freshwater
  fc[!is.na(orden) & orden=="Alismatales"] = "Macroalgae" # Water plantains, including Zostera
  
  
  # ------- Parasites ----------
  
  fc[!is.na(phylum) & phylum=="Apicomplexa"] = "Parasites"
  fc[!is.na(phylum) & phylum=="Metamonada"] = "Parasites" # anaerobic
  fc[!is.na(phylum) & phylum=="Ichthyosporea"] = "Parasites"
  fc[!is.na(phylum) & phylum=="Ichthyosporidia"] = "Parasites" # EUKARYOME
  fc[!is.na(phylum) & phylum=="Mesomycetozoa"] = "Parasites"
  fc[!is.na(class) & class=="Endomyxa"] = "Parasites"
  fc[!is.na(class) & class=="Endomyxa-Phytomyxea"] = "Parasites"
  fc[!is.na(class) & class=="Haplosporida"] = "Parasites"
  fc[!is.na(class) & class=="Monogenea class incertae sedis"] = "Parasites"
  fc[!is.na(class) & class=="Monogenea"] = "Parasites"
  fc[!is.na(class) & class=="Trematoda class incertae sedis"] = "Parasites"
  fc[!is.na(class) & class=="Trematoda"] = "Parasites"
  fc[!is.na(class) & class=="Cestoda"] = "Parasites"
  fc[!is.na(phylum) & phylum=="Orthonectida"] = "Parasites"
  # fc[!is.na(class) & class=="Oomycota"] = "Parasites"
  # fc[!is.na(class) & class=="Oomycetes"] = "Parasites"
  fc[!is.na(class) & class=="Hyphochytriomyceta"] = "Parasites"
  # fc[!is.na(class) & class=="Labyrinthulomycetes"] = "Parasites"
  fc[!is.na(orden) & orden=="Thraustochytriales order incertae sedis"] = "Protists" # Exception
  fc[!is.na(phylum) & phylum=="Perkinsea"] = "Parasites"
  fc[!is.na(orden) & orden=="Syndiniales"] = "Parasites"
  fc[!is.na(fam) & fam=="Ichthyobodonidae"] = "Parasites"
  fc[!is.na(phylum) & phylum=="Perkinsia"] = "Parasites (of phytoplankton)"
  
  fc[!is.na(phylum) & phylum=="Planomonada"] = "Commensal"
  fc[!is.na(class) & class=="Oxymonadia"] = "Commensal"
  
    # V4
  fc[!is.na(genus) & genus=="Echinorhynchus"] = "Parasites" # Parasitic worms
  fc[!is.na(genus) & genus=="Nectonema"] = "Parasites" # Parasitic nematiodes
  fc[!is.na(phylum) & phylum=="Orthonectida"] = "Parasites" 
  
  # V5
  fc[!is.na(fam) & fam=="Pyramidellidae"] = "Parasites" # Parasitic micro-molluscs
  fc[!is.na(fam) & fam=="Ergasilidae"] = "Parasites" # Fish parasitic microcrustacea
  fc[!is.na(orden) & orden=="Poecilostomatoida"] = "Parasites" # Fish parasitic microcrustacea
  fc[!is.na(orden) & orden=="Siphonostomatoida"] = "Parasites" # Fish parasitic microcrustacea
  
  
  
  # ------- Terrestrial -------
  
  fc[!is.na(class) & class=="Embryophyta"] = "Terrestrial" # Plants, no longer a label
  fc[!is.na(class) & class=="Amphibia"] = "Terrestrial"
  fc[!is.na(class) & class=="Aves"] = "Terrestrial"
  fc[!is.na(class) & class=="Mammalia"] = "Terrestrial"
  fc[!is.na(class) & class=="Insecta"] = "Terrestrial"
  fc[!is.na(orden) & orden=="Araneae"] = "Terrestrial"
  fc[!is.na(class) & class=="Crocodylia class incertae sedis"] = "Terrestrial"
  fc[!is.na(class) & class=="Hirudinida"] = "Terrestrial" # Leeches
  fc[!is.na(orden) & orden=="Opiliones"] = "Terrestrial" # Daddy longlegs
  fc[!is.na(orden) & orden=="Ricinulei"] = "Terrestrial" # Spiders
  fc[!is.na(genus) & genus=="Penicillium"] = "Terrestrial" # Spiders
  fc[!is.na(orden) & orden=="Glomerellales"] = "Terrestrial"
  fc[!is.na(genus) & genus=="Catenulispora"] = "Terrestrial" # Bacteria from forest soil - forms mycelia
  fc[!is.na(phylum) & phylum=="Bryophyta"] = "Terrestrial" # Mosses (though some freshwater ones exist)
  fc[!is.na(class) & class=="Marchantiophyta"] = "Terrestrial" # Liverworts, non-vascular land plants
  fc[!is.na(orden) & orden=="Apiales"] = "Terrestrial"
  fc[!is.na(orden) & orden=="Asparagales"] = "Terrestrial"
  fc[!is.na(orden) & orden=="Cucurbitales"] = "Terrestrial"
  fc[!is.na(orden) & orden=="Ericales"] = "Terrestrial"
  fc[!is.na(orden) & orden=="Fabales"] = "Terrestrial"
  fc[!is.na(orden) & orden=="Fagales"] = "Terrestrial"
  fc[!is.na(orden) & orden=="Gentianales"] = "Terrestrial"
  fc[!is.na(orden) & orden=="Ericales"] = "Terrestrial"
  fc[!is.na(orden) & orden=="Lamiales"] = "Terrestrial"
  fc[!is.na(orden) & orden=="Laurales"] = "Terrestrial"
  fc[!is.na(orden) & orden=="Malpighiales"] = "Terrestrial"
  fc[!is.na(orden) & orden=="Malvales"] = "Terrestrial"
  fc[!is.na(orden) & orden=="Myrtales"] = "Terrestrial"
  fc[!is.na(orden) & orden=="Poales"] = "Terrestrial"
  fc[!is.na(orden) & orden=="Ranunculales"] = "Terrestrial"
  fc[!is.na(orden) & orden=="Sapindales"] = "Terrestrial"
  fc[!is.na(orden) & orden=="Solanales"] = "Terrestrial"
  fc[!is.na(class) & class=="Cupressopsida"] = "Terrestrial"
  fc[!is.na(class) & class=="Lycopodiopsida"] = "Terrestrial"
  fc[!is.na(class) & class=="Pinopsida"] = "Terrestrial"
  fc[!is.na(fam) & fam=="Dryopteridaceae"] = "Terrestrial"
  fc[!is.na(orden) & orden=="Poales"] = "Terrestrial"
  
  # v0.5 mollusks etc
  fc[!is.na(orden) & orden=="Stylommatophora"] = "Terrestrial"
  #fc[!is.na(orden) & orden=="Systellommatophora"] = "Terrestrial"
  #fc[!is.na(orden) & orden=="Stylommatophora"] = "Terrestrial"
  
  
  # Mammals come with exceptions (cetaceans etc)
  fc[!is.na(orden) & orden=="Sirenia"] = "Nekton" # Sea cows
  fc[!is.na(fam) & fam=="Balaenidae"] = "Nekton" # Whales
  fc[!is.na(fam) & fam=="Cetotheriidae"] = "Nekton" # Whales
  fc[!is.na(fam) & fam=="Balaenidae"] = "Nekton" # Whales
  fc[!is.na(fam) & fam=="Balaenopteridae"] = "Nekton" # Whales
  fc[!is.na(fam) & fam=="Eschrichtiidae"] = "Nekton" # Whales
  fc[!is.na(fam) & fam=="Platanistidae"] = "Nekton" # Whales
  fc[!is.na(fam) & fam=="Physeteridae"] = "Nekton" # Whales
  fc[!is.na(fam) & fam=="Kogiidae"] = "Nekton" # Whales
  fc[!is.na(fam) & fam=="Ziphiidae"] = "Nekton" # Whales
  fc[!is.na(fam) & fam=="Delphinida"] = "Nekton" # Dolphins
  
  
  # ------ PROKARYOTIC ---------
  
  # Organelles
  fc[!is.na(meta) & meta=="Chloroplast"] = "Chloroplasts" # Mostly phytoplankton
  fc[!is.na(meta) & meta=="Mitochondria"] = "Mitochondria" # Poor taxonomy
  
  # Set all bacteria as bacetrial heterotrophs by default
  fc[!is.na(domain) & domain=="Bacteria"] = "Heterotrophic aerobic prokaryotes" # Poor taxonomy
  # And some archaea
  fc[!is.na(phylum) & phylum=="Heimdallarchaeota"] = "Heterotrophic aerobic prokaryotes"
  #
  
  # # Aerobic bacterial heterotrophs - specialised
  # fc[!is.na(phylum) & phylum=="Planctomycetota"] = "Bact. heterotrophs spec."
  # Abditibacteriota, Acidobacteriota, Aminicenantia (though also anaerobic)
  # Actinobacteriota (incl Isobutyronitrile), Aerophobota, Aquificota, Armatimonadota, Marinifilaceae
  # Chitinophagaceae, Saprospiraceae, Amoebophilaceae, Cyclobacteriaceae, Flammeovirgaceae
  # Flavobacteriales, Sphingobacteriaceae, Ignavibacteria, Balneolales, Ktedonobacteria
  # Dadabacteria (though some are photohetero), Gemmatimonadota (some photohetero)
  # Myxococcota, 
  # Particle associated Planctomycetes: OM190, Phycisphaera, Gemmatales, Isosphaerales, Pirellulales, Planctomycetales
  # Proteobacteria: Acetobacteraceae, Caulobacteraceae, Hyphomonadaceae, Parvularculaceae, Emcibacteraceae
  # Kiloniellaceae (some photohetero), Kordiimonadales, Parvibaculales, Rhizobiales (some associated w diatoms)
  # Sphingomonadales, Thalassobaculales, Burkholderiales
  # Alcanivorax (basically oil degrading bacteria!)
  # Marinobacter (can degrade HCs)
  # Nitrincolaceae (probable phytoplankton assoicated)
  # Dasania (probable phytoplankton assoicated)
  # Pseudoxanthomonas taiwanensis (basically oil degrading bacteria!)
  # Verrucomicrobiota, Ferrovibrio, Dongiales, Acidibacter, Wenzhouxiangella
  
  # Cyanobacteria
  fc[!is.na(phylum) & phylum=="Cyanobacteria"] = "Photoautotrophic Cyanobacteria"
  
  # Green S
  fc[!is.na(phylum) & phylum=="Chlorobia"] = "Photoautotrophic green S bacteria"
  
  # Purple S
  fc[!is.na(fam) & fam=="Chromatiaceae"] = "Photoautotrophic Purple S" 
  fc[!is.na(fam) & fam=="Ectothiorhodospiraceae"] = "Photoautotrophic Purple S"
  
  # Prokaryotic symbionts or parasites with prokaryotic hosts
  fc[!is.na(kingdom) & kingdom=="DPANN"] = "Prok-prok symbionts or parasites" # obligate archaea symbionts w Chloroflexi hosts
  fc[!is.na(class) & class=="Bdellovibrionia"] = "Prok-prok symbionts or parasites"
  fc[!is.na(phylum) & phylum=="Gracilibacteria"] = "Prok-prok symbionts or parasites" 
  fc[!is.na(orden) & orden=="Micavibrionales"] = "Prok-prok symbionts or parasites" 
  fc[!is.na(orden) & orden=="Candidatus Endonucleariobacter order incertae sedis"] = "Prok-prok symbionts or parasites" 
  
  
  # Anaerobic heterotrophs
  fc[!is.na(phylum) & phylum=="Bathyarchaeota"] = "Heterotrophic anaerobic prok"
  fc[!is.na(class) & phylum=="Thermoplasmatota"] = "Heterotrophic anaerobic prok" # Thermofiles
  fc[!is.na(class) & class=="Thermococci"] = "Heterotrophic anaerobic prok" # Thermofiles
  fc[!is.na(fam) & fam=="Dysgonomonadaceae"] = "Heterotrophic anaerobic prok" 
  fc[!is.na(genus) & genus=="Paludibacter"] = "Heterotrophic anaerobic prok"
  fc[!is.na(fam) & fam=="Prolixibacteraceae"] = "Heterotrophic anaerobic prok"
  fc[!is.na(fam) & fam=="Tannerellaceae"] = "Heterotrophic anaerobic prok" 
  fc[!is.na(fam) & fam=="Lentimicrobiaceae"] = "Heterotrophic anaerobic prok" 
  fc[!is.na(phylum) & phylum=="Caldatribacteriota"] = "Heterotrophic anaerobic prok"
  fc[!is.na(phylum) & phylum=="Calditrichota"] = "Heterotrophic anaerobic prok"
  fc[!is.na(class) & class=="Anaerolineae"] = "Heterotrophic anaerobic prok"
  fc[!is.na(phylum) & phylum=="Deferrisoma"] = "Heterotrophic anaerobic prok"
  fc[!is.na(phylum) & phylum=="Deinococcota"] = "Heterotrophic anaerobic prok"
  fc[!is.na(phylum) & phylum=="Desulfobacterota"] = "Heterotrophic anaerobic prok"
  fc[!is.na(phylum) & phylum=="Dojkabacteria"] = "Heterotrophic anaerobic prok"
  fc[!is.na(phylum) & phylum=="Fermentibacterota"] = "Heterotrophic anaerobic prok"
  fc[!is.na(phylum) & phylum=="Fibrobacterota"] = "Heterotrophic anaerobic prok"
  fc[!is.na(genus) & genus=="Anaerobacillus"] = "Heterotrophic anaerobic prok"
  fc[!is.na(genus) & genus=="Anoxybacillus"] = "Heterotrophic anaerobic prok"
  fc[!is.na(orden) & orden=="Izemoplasmatales"] = "Heterotrophic anaerobic prok"
  fc[!is.na(orden) & orden=="Lactobacillales"] = "Heterotrophic anaerobic prok"
  fc[!is.na(class) & class=="Clostridia"] = "Heterotrophic anaerobic prok"
  fc[!is.na(genus) & genus=="Oceanirhabdus"] = "Heterotrophic anaerobic prok"
  fc[!is.na(genus) & genus=="Anaerofustis"] = "Heterotrophic anaerobic prok"
  fc[!is.na(orden) & orden=="Oscillospirales"] = "Heterotrophic anaerobic prok"
  fc[!is.na(orden) & orden=="Peptococcales"] = "Heterotrophic anaerobic prok"
  fc[!is.na(orden) & orden=="Peptostreptococcales-Tissierellales"] = "Heterotrophic anaerobic prok" # With exceptions below (host ass.)
  fc[!is.na(class) & class=="D8A-2"] = "Heterotrophic anaerobic prok" # syntrophic acetate-oxidizing
  fc[!is.na(class) & class=="Desulfotomaculia"] = "Heterotrophic anaerobic prok" # S reducers
  fc[!is.na(class) & class=="Dethiobacteria"] = "Heterotrophic anaerobic prok" # S reducers
  fc[!is.na(genus) & genus=="Hydrogenispora"] = "Heterotrophic anaerobic prok"
  fc[!is.na(class) & class=="Moorellia"] = "Heterotrophic anaerobic prok"
  fc[!is.na(class) & class=="Negativicutes"] = "Heterotrophic anaerobic prok"
  fc[!is.na(class) & class=="Thermoanaerobacteria"] = "Heterotrophic anaerobic prok"
  fc[!is.na(phylum) & phylum=="Fusobacteriota"] = "Heterotrophic anaerobic prok" # With exceptions below (host ass.)
  fc[!is.na(phylum) & phylum=="Hydrogenedentes"] = "Heterotrophic anaerobic prok" # With exceptions below (host ass.)
  fc[!is.na(phylum) & phylum=="Latescibacterota"] = "Heterotrophic anaerobic prok" # With exceptions below (host ass.)
  fc[!is.na(phylum) & phylum=="MBNT15"] = "Heterotrophic anaerobic prok" 
  fc[!is.na(phylum) & phylum=="Methylomirabilota"] = "Heterotrophic anaerobic prok"
  fc[!is.na(phylum) & phylum=="Modulibacteria"] = "Heterotrophic anaerobic prok" # Some are host associationed but there is no taxonomy
  fc[!is.na(fam) & fam=="Defluviicoccaceae"] = "Heterotrophic anaerobic prok"
  fc[!is.na(orden) & orden=="Micropepsales"] = "Heterotrophic anaerobic prok"
  fc[!is.na(phylum) & phylum=="Spirochaetota"] = "Heterotrophic anaerobic prok"
  fc[!is.na(phylum) & phylum=="Sumerlaeota"] = "Heterotrophic anaerobic prok"
  fc[!is.na(phylum) & phylum=="Synergistota"] = "Heterotrophic anaerobic prok"
  fc[!is.na(phylum) & phylum=="Thermotogota"] = "Heterotrophic anaerobic prok" # Hypertermophiles
  fc[!is.na(phylum) & phylum=="Caldisericota"] = "Heterotrophic anaerobic prok" # Hypertermophiles, S reducers
  fc[!is.na(class) & class=="Caldatribacteriia"] = "Heterotrophic anaerobic prok" # Hypertermophiles,
  fc[!is.na(class) & class=="Natranaerobiia"] = "Heterotrophic anaerobic prok"
  fc[!is.na(class) & class=="Thermacetogenia"] = "Heterotrophic anaerobic prok"  
  fc[!is.na(class) & class=="Holophagae"] = "Heterotrophic anaerobic prok"  # Marine Acidobacteria, anaerobic as opposed to class Acidoberiae
  fc[!is.na(phylum) & phylum=="Coprothermobacterota"] = "Heterotrophic anaerobic prok" # Hyperthermophilic w Archaea genetic elements
  fc[!is.na(phylum) & phylum=="Cloacimonadota"] = "Heterotrophic anaerobic prok"
  fc[!is.na(genus) & genus=="Actinomyces"] = "Heterotrophic anaerobic prok"
  fc[!is.na(class) & class=="Moorellia"] = "Heterotrophic anaerobic prok"
  
  
  # Host associated (eukaryotic host)
  fc[!is.na(orden) & orden=="Propionibacteriales"] = "Host associated prok"
  fc[!is.na(genus) & genus=="Prevotella"] = "Host associated prok"
  fc[!is.na(class) & class=="Vampirivibrionia"] = "Host associated prok"
  fc[!is.na(phylum) & phylum=="Dependentiae"] = "Host associated prok"
  fc[!is.na(phylum) & phylum=="Elusimicrobiota"] = "Host associated prok"
  fc[!is.na(phylum) & phylum=="Entotheonellaeota"] = "Host associated prok"
  fc[!is.na(orden) & orden=="Entomoplasmatales"] = "Host associated prok"
  fc[!is.na(orden) & orden=="Erysipelotrichales"] = "Host (mammal) associated prok"
  fc[!is.na(genus) & genus=="Enterococcus"] = "Host associated prok" # Mammal / human
  fc[!is.na(genus) & genus=="Ligilactobacillus"] = "Host associated prok"
  fc[!is.na(genus) & genus=="Streptococcus"] = "Host associated prok"
  fc[!is.na(fam) & fam=="Mycoplasmataceae"] = "Host associated prok"
  fc[!is.na(genus) & genus=="Staphylococcus"] = "Host associated prok"
  fc[!is.na(genus) & genus=="Christensenella"] = "Host associated prok" # Mammal / human
  fc[!is.na(orden) & orden=="Lachnospirales"] = "Host associated prok"
  fc[!is.na(genus) & genus=="Ruminococcus"] = "Host associated prok"
  fc[!is.na(genus) & genus=="Anaerovoracaceae"] = "Host associated prok" # Peptostreptococcales
  fc[!is.na(genus) & genus=="Anaerococcus"] = "Host associated prok" # Peptostreptococcales
  fc[!is.na(genus) & genus=="Peptoniphilus"] = "Host associated prok" # Peptostreptococcales
  fc[!is.na(genus) & genus=="Romboutsia"] = "Host associated prok" # Peptostreptococcales
  fc[!is.na(genus) & genus=="Psychrilyobacter"] = "Host associated prok" # Fusobacteriota
  fc[!is.na(fam) & fam=="Leptotrichiaceae"] = "Host associated prok" # Fusobacteriota
  fc[!is.na(phylum) & phylum=="Margulisbacteria"] = "Host associated prok" 
  fc[!is.na(phylum) & phylum=="NB1-j"] = "Host associated prok" # Parasites or predators of algae
  fc[!is.na(phylum) & phylum=="Parcubacteria"] = "Host associated prok" # Parasites or predators of algae
  fc[!is.na(phylum) & phylum=="Poribacteria"] = "Host associated prok" # Sponge symbionts
  fc[!is.na(orden) & orden=="Holosporaceae"] = "Host associated prok" # Ciliate host
  fc[!is.na(orden) & orden=="Paracaedibacterales"] = "Host associated prok" # Amoeba host
  fc[!is.na(genus) & genus=="Halocynthiibacter"] = "Host associated prok" # Isolated from sea squirt / pineapple
  fc[!is.na(genus) & genus=="Sedimentitalea"] = "Host associated prok"
  fc[!is.na(genus) & genus=="Candidatus Riegeria"] = "Host associated prok"
  fc[!is.na(orden) & orden=="Rickettsiales"] = "Host associated prok"
  fc[!is.na(orden) & orden=="Candidatus Berkiella order incertae sedis"] = "Host associated prok"
  fc[!is.na(orden) & orden=="Candidatus Ovatusbacter order incertae sedis"] = "Host associated prok"
  fc[!is.na(orden) & orden=="Coxiellales"] = "Host associated prok"
  fc[!is.na(genus) & genus=="Escherichia-Shigella"] = "Host associated prok"
  fc[!is.na(genus) & genus=="Moritella"] = "Host associated prok"
  fc[!is.na(genus) & genus=="Aliivibrio"] = "Host associated prok" # bioluminescent squid symbionts
  fc[!is.na(genus) & genus=="Photobacterium"] = "Host associated prok" # bioluminescent symbionts (fish etc)
  fc[!is.na(orden) & orden=="Francisellales"] = "Host associated prok" # weird and scary pathogens
  fc[!is.na(genus) & genus=="Legionella"] = "Host associated prok"
  fc[!is.na(fam) & fam=="Piscirickettsiaceae"] = "Host associated prok"
  fc[!is.na(fam) & fam=="Thiomicrospiraceae"] = "Host associated prok" # Ciliate host
  fc[!is.na(phylum) & phylum=="Saccharimonadia"] = "Host associated prok" 
  fc[!is.na(orden) & orden=="Bifidobacteriales"] = "Host associated prok" 
  fc[!is.na(orden) & orden=="Brevinematales"] = "Host associated prok" 
  fc[!is.na(genus) & genus=="Arcanobacterium"] = "Host associated prok" 
  fc[!is.na(genus) & genus=="Trueperella"] = "Host associated prok"
  fc[!is.na(genus) & genus=="Varibaculum"] = "Host associated prok" 
  fc[!is.na(orden) & orden=="Euzebyales"] = "Host associated prok" 
  #fc[!is.na(genus) & genus==""] = "Host associated prok" 
  
  # Methanogens (anaerobic hetero also but easy to identify)
  fc[!is.na(class) & class=="Methanobacteria"] = "Methanogens"
  fc[!is.na(class) & class=="Methanocellia"] = "Methanogens"
  fc[!is.na(class) & class=="Methanomicrobia"] = "Methanogens"
  fc[!is.na(class) & class=="Methanonatronarchaeia"] = "Methanogens"
  fc[!is.na(class) & class=="Methanosarcinia"] = "Methanogens"
  fc[!is.na(class) & class=="Methanopyri"] = "Methanogens"
  fc[!is.na(orden) & orden=="Methanomassiliicoccales"] = "Methanogens"
  fc[!is.na(orden) & orden=="Methanofastidiosales"] = "Methanogens"
  fc[!is.na(class) & class=="Syntrophomonadia"] = "Methanogens" # In reality obligately synthropic w methanogens
  fc[!is.na(genus) & genus=="Gelria"] = "Methanogens" # In reality obligately synthropic w methanogens
  fc[!is.na(class) & class=="Syntrophorhabdia"] = "Methanogens" # In reality obligately synthropic w methanogens. Can degrade phenol
  
  # Methanotrophs
  fc[!is.na(fam) & fam=="Beijerinckiaceae"] = "Methanotrophs"
  fc[!is.na(fam) & fam=="Methylocystaceae"] = "Methanotrophs"
  fc[!is.na(fam) & fam=="Methylophilaceae"] = "Methanotrophs"
  fc[!is.na(genus) & genus=="Methyloversatilis"] = "Methanotrophs"
  fc[!is.na(orden) & orden=="Methylococcales"] = "Methanotrophs"
  fc[!is.na(fam) & fam=="Methylophagaceae"] = "Methanotrophs"
  fc[!is.na(fam) & startsWith(fam,"ANME")] = "Methanotrophs" # Anaerobic
  fc[!is.na(fam) & startsWith(genus,"ANME")] = "Methanotrophs" # Anaerobic
  
  
  # Anaerobic chemoautotrophs
  fc[!is.na(phylum) & phylum=="Hadarchaeota"] = "Chemoautotrophic anaerobic prok"
  fc[!is.na(phylum) & phylum=="Hydrothermarchaeota"] = "Chemoautotrophic anaerobic prok"
  fc[!is.na(phylum) & phylum=="Deferribacterota"] = "Chemoautotrophic anaerobic prok"
  fc[!is.na(phylum) & phylum=="Brocadiae"] = "Chemoautotrophic anaerobic prok" # Annamox
  fc[!is.na(class) & class=="Symbiobacteriia"] = "Chemoautotrophic anaerobic prok" # Anaerobic propane oxidation coupled to denotrofocation
  fc[!is.na(class) & class=="Thermaerobacteria"] = "Chemoautotrophic anaerobic prok" # Anaerobic propane oxidation coupled to denotrofocation
  fc[!is.na(class) & class=="Thermodesulfobacteria"] = "Chemoautotrophic anaerobic prok" # Fe reducers
  fc[!is.na(class) & class=="Thermolithobacterales class incertae sedis"] = "Chemoautotrophic anaerobic prok" # Fe reducers
  fc[!is.na(orden) & orden=="Thermincolales"] = "Chemoautotrophic anaerobic prok" # Fe reducers
  fc[!is.na(fam) & fam=="Desulforudaceae"] = "Chemoautotrophic anaerobic prok" # Survives at 3 km sub-sediment on radioactive decay
  
  # Aerobic chemoautotrophs
  fc[!is.na(phylum) & phylum=="Thaumarchaeota"] = "Chemoautotrophic aerobic prok" 
  fc[!is.na(phylum) & phylum=="Campylobacterota"] = "Chemoautotrophic aerobic prok"
  fc[!is.na(phylum) & phylum=="Nitrospinota"] = "Chemoautotrophic aerobic prok" # Nitrifying bact.
  fc[!is.na(fam) & fam=="Nitrosomonadaceae"] = "Chemoautotrophic aerobic prok" # AOB
  fc[!is.na(genus) & genus=="Nitrosococcus"] = "Chemoautotrophic aerobic prok" # AOB
  fc[!is.na(genus) & genus=="Nitrobacter"] = "Chemoautotrophic aerobic prok" # Nitrifying bact.
  fc[!is.na(phylum) & phylum=="Nitrospirota"] = "Chemoautotrophic aerobic prok" # Nitrifying bact. w genera like Leptospirillum that is Fe oxidising
  fc[!is.na(orden) & orden=="Thermodesulfovibrionales"] = "Heterotrophic aerobic prok" # Exception
  fc[!is.na(genus) & genus=="Magnetospira"] = "Chemoautotrophic aerobic prok" 
  fc[!is.na(orden) & orden=="Acidiferrobacterales"] = "Chemoautotrophic aerobic prok" 
  fc[!is.na(orden) & orden=="Arenicellales"] = "Chemoautotrophic aerobic prok"
  fc[!is.na(orden) & orden=="Beggiatoales"] = "Chemoautotrophic aerobic prok"
  fc[!is.na(orden) & orden=="Gallionellaceae"] = "Chemoautotrophic aerobic prok" # Fe Ox
  fc[!is.na(genus) & genus=="Sulfuritalea"] = "Chemoautotrophic aerobic prok" # Fe Ox
  fc[!is.na(genus) & genus=="Thiohalobacter"] = "Chemoautotrophic aerobic prok" # S
  fc[!is.na(fam) & fam=="Sedimenticolaceae"] = "Chemoautotrophic aerobic prok" # S
  fc[!is.na(fam) & fam=="Nitrosococcaceae"] = "Chemoautotrophic aerobic prok"
  fc[!is.na(fam) & fam=="Sulfuriflexus"] = "Chemoautotrophic aerobic prok"
  fc[!is.na(fam) & fam=="Tenderiaceae"] = "Chemoautotrophic aerobic prok"
  fc[!is.na(fam) & fam=="Thioprofundaceae"] = "Chemoautotrophic aerobic prok"
  fc[!is.na(genus) & genus=="Leucothrix"] = "Chemoautotrophic aerobic prok"
  fc[!is.na(genus) & genus=="Thiothrix"] = "Chemoautotrophic aerobic prok"
  fc[!is.na(class) & class=="Magnetococcia"] = "Chemoautotrophic aerobic prok"
  fc[!is.na(class) & class=="Zetaproteobacteria"] = "Chemoautotrophic aerobic prok"
  fc[!is.na(phylum) & phylum=="Sva0485"] = "Chemoautotrophic aerobic prok"
  fc[!is.na(fam) & fam=="Halothiobacillaceae"] = "Chemoautotrophic aerobic prok"  # S Ox
  fc[!is.na(class) & class=="Acidimicrobiia"] = "Chemoautotrophic aerobic prok" # Fe ox
  
  
  # Photoheterotrophs
  fc[!is.na(class) & class=="Chloroflexia"] = "Photoheterophic prok"
  fc[!is.na(genus) & genus=="Pelagibius"] = "Photoheterophic prok"
  fc[!is.na(orden) & orden=="Pelagibacterales"] = "Photoheterophic prok"
  fc[!is.na(orden) & orden=="Puniceispirillales"] = "Photoheterophic prok"
  fc[!is.na(genus) & genus=="Amylibacter"] = "Photoheterophic prok"
  fc[!is.na(genus) & genus=="Dinoroseobacter"] = "Photoheterophic prok"
  fc[!is.na(genus) & genus=="Donghicola"] = "Photoheterophic prok"
  fc[!is.na(genus) & genus=="Litoreibacter"] = "Photoheterophic prok"
  fc[!is.na(genus) & genus=="Pseudophaeobacter"] = "Photoheterophic prok"
  fc[!is.na(genus) & genus=="Dokdonia"] = "Photoheterophic prok"
  fc[!is.na(genus) & genus=="Pseudophaeobacter"] = "Photoheterophic prok"
  fc[!is.na(genus) & genus=="Erythrobacter"] = "Photoheterophic prok"
  fc[!is.na(genus) & genus=="Acidiphilium"] = "Photoheterophic prok"
  fc[!is.na(genus) & genus=="Rhodobacter"] = "Photoheterophic prok"
  fc[!is.na(fam) & fam=="Nitrobacteraceae"] = "Photoheterophic prok"
  fc[!is.na(genus) & genus=="Rhodomicrobium"] = "Photoheterophic prok"
  fc[!is.na(fam) & fam=="Heliobacteriaceae"] = "Photoheterophic prok"
  fc[!is.na(genus) & genus=="Roseovarius"] = "Photoheterophic prok"
  fc[!is.na(genus) & genus=="Rubellimicrobium"] = "Photoheterophic prok"
  fc[!is.na(genus) & genus=="Shimia"] = "Photoheterophic prok"
  fc[!is.na(genus) & genus=="Luminiphilus"] = "Photoheterophic prok"
  fc[!is.na(class) & class=="Halobacteria"] = "Photoheterophic prok"
  
  taxon_table$category = as.factor(fc)
  
  # Return
  return(taxon_table)
}

## Return a list of OTU names that are potential mocks
identifyMocks = function(taxon_table, presentIn, minAbundance=0) {
  
  present = presentIn>minAbundance
  potMocks = vector(length=dim(taxon_table)[1])
  
  # 18S - the six chordata are unresolved
  #potMocks[!is.na(taxon_table$phylum) & taxon_table$phylum=="Chordata"] <- TRUE
  potMocks[grep("Chordata", taxon_table$classification)] <- TRUE
  
  
  # Earthworm
  #potMocks[!is.na(taxon_table$fam) & taxon_table$fam=="Lumbricidae"] <- TRUE
  potMocks[grep("Lumbricidae", taxon_table$classification)] <- TRUE
  
  # Leech
  # potMocks[!is.na(taxon_table$class) & taxon_table$class=="Hirudinida"] <- TRUE
  # potMocks[!is.na(taxon_table$order) & taxon_table$order=="Hirudinida"] <- TRUE
  potMocks[grep("Hirudinida", taxon_table$classification)] <- TRUE
  
  # Slug or snail
  # potMocks[!is.na(taxon_table$genus) & taxon_table$genus=="Euthyneura"] <- TRUE
  # potMocks[!is.na(taxon_table$order) & taxon_table$order=="Stylommatophora"] <- TRUE
  potMocks[grep("Euthyneura", taxon_table$classification)] <- TRUE
  potMocks[grep("Stylommatophora", taxon_table$classification)] <- TRUE
  
  # Insects (honey bee, red ant, houseflyy, leaf beetle, orangetip, dragonfly, libellulla)
  potMocks[!is.na(taxon_table$class) & taxon_table$class=="Insecta"] <- TRUE
  
  # Spiders
  potMocks[!is.na(taxon_table$class) & taxon_table$class=="Arachnida"] <- TRUE
  
  # Woodlouse - none found in V1V2
  potMocks[!is.na(taxon_table$genus) & taxon_table$genus=="Porcellio"] <- TRUE # COI - change to isopoda?
  potMocks[!is.na(taxon_table$fam) & taxon_table$fam=="Isopoda sp. BOLD"] <- TRUE # COI - change to isopoda?
  
  ## 16 S has:
  # Pseudomonas aeruginosa
  potMocks[!is.na(taxon_table$genus) & taxon_table$genus=="Pseudomonas"] <- TRUE
  
  # Escherichia coli
  potMocks[!is.na(taxon_table$genus) & taxon_table$genus=="Escherichia-Shigella"] <- TRUE
  
  # Salmonella enterica
  potMocks[!is.na(taxon_table$order) & taxon_table$order=="Enterobacterales"] <- TRUE 
  
  # Lactobacillus fermentum
  potMocks[!is.na(taxon_table$fam) & taxon_table$fam=="Lactobacillaceae"] <- TRUE 
  
  # Enterococus faecalis
  potMocks[!is.na(taxon_table$genus) & taxon_table$genus=="Enterococcus"] <- TRUE
  
  # Staphylococcus aureus
  potMocks[!is.na(taxon_table$genus) & taxon_table$genus=="Staphylococcus"] <- TRUE
  
  # Listeria monocytogenes
  potMocks[!is.na(taxon_table$genus) & taxon_table$genus=="Listeria"] <- TRUE 
  
  # Bacillus subtilis
  potMocks[!is.na(taxon_table$genus) & taxon_table$genus=="Bacillus"] <- TRUE 
  
  # Saccharomyces cerevisiae, Cryptococcus neoformans (fungi)
  
  return(row.names(taxon_table)[potMocks & present])
}


