library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(tidyr)
library(BSgenome.Celegans.UCSC.ce11)
library(ggpubr)
library(tidyr)

ce11seqinfo<-seqinfo(Celegans)

#####################-
## Chromatin states Evans et al. (2016) (Ahinger lab)------
#####################-
# https://www.pnas.org/content/pnas/113/45/E7020.full.pdf

# Dataset S1. Coordinates of EE and L3 chromatin states
# File of chromosome, start position, end position, and state number. Coordinates are in WS220 and follow BED conventions (start positions are in zero-based coordinates and end positions in one-based coordinates).
outPath="."
dir.create("./publicData")
remakeFiles=T
if(remakeFiles | !file.exists(paste0(outPath,"/publicData/chromDomains_L3_Evans2016_ce11.bed"))){
  chromStatesURL<-"https://www.pnas.org/highwire/filestream/623778/field_highwire_adjunct_files/0/pnas.1608162113.sd01.xlsx"
  download.file(chromStatesURL,paste0(outPath,"/publicData/",basename(chromStatesURL)))
  ce10toCe11url<-"http://hgdownload.soe.ucsc.edu/goldenPath/ce10/liftOver/ce10ToCe11.over.chain.gz"
  ce10toCe11<-"ce10Toce11.over.chain"
  download.file(ce10toCe11url,paste0(outPath,"/publicData/",ce10toCe11,".gz"))
  system(paste0("gunzip ",outPath,"/publicData/",ce10toCe11,".gz"))
  file.remove(paste0(outPath,"/publicData/",ce10toCe11,".gz"))

  chrAstates<-readxl::read_excel(paste0(outPath,"/publicData/pnas.1608162113.sd01.xlsx"),sheet="L3 autosome states",col_names=c("chr","start","end","state"))
  chrXstates<-readxl::read_excel(paste0(outPath,"/publicData/pnas.1608162113.sd01.xlsx"),sheet="L3 chr X states",col_names=c("chr","start","end","state"))


  chrStates<-rbind(chrAstates,chrXstates)
  chrStatesGR<-GRanges(paste0("chr",chrStates$chr,":",
                              (chrStates$start+1),"-",
                              chrStates$end))
  chrStatesGR$score<-c(chrAstates$state,chrXstates$state)
  chainCe10toCe11<-import.chain(paste0(outPath,"/publicData/",ce10toCe11))
  chrStatesGR_ce11<-unlist(liftOver(chrStatesGR,chain=chainCe10toCe11))

  seqinfo(chrStatesGR_ce11)<-ce11seqinfo
  export(chrStatesGR_ce11,
         con=paste0(outPath,"/publicData/chromStates_L3_Evans2016_ce11.bed"),
         format="bed")
  file.remove(paste0(outPath,"/publicData/",basename(chromStatesURL)))

  # Dataset S2. Coordinates of EE and L3 domains
  # Excel file of chromosome, start position, end position of EE and L3 active domains, border regions, and regulated domains (each in separate tab). Additionally, border regions have strand information in column six to indicate if active domain is on the left (−) or on the right (+). Coordinates are in WS220 and follow BED conventions.
  chromDomainsURL<-"https://www.pnas.org/highwire/filestream/623778/field_highwire_adjunct_files/1/pnas.1608162113.sd02.xlsx"
  download.file(chromDomainsURL,paste0(outPath,"/publicData/",basename(chromDomainsURL)))

  l3active<-readxl::read_excel(paste0(outPath,"/publicData/",
                                      basename(chromDomainsURL)),
                               sheet="L3 active domains",
                               col_names=c("chr","start","end"))
  l3active$name<-"active"
  l3active$score<-"."
  l3active$strand<-"*"
  l3regulated<-readxl::read_excel(paste0(outPath,"/publicData/",
                                         basename(chromDomainsURL)),
                                  sheet="L3 regulated domains",
                                  col_names=c("chr","start","end"))
  l3regulated$name<-"regulated"
  l3regulated$score<-"."
  l3regulated$strand<-"*"
  l3borders<-readxl::read_excel(paste0(outPath,"/publicData/",basename(chromDomainsURL)),sheet="L3 borders",col_names=c("chr","start","end","name","score","strand"))
  l3borders$name<-"border"


  chrDomains<-rbind(l3active,l3regulated,l3borders)
  chrDomainsGR<-GRanges(paste0("chr",chrDomains$chr,":",(chrDomains$start+1),"-",
                               chrDomains$end,":",chrDomains$strand))
  mcols(chrDomainsGR)<-chrDomains[,c("name","score")]
  chainCe10toCe11<-import.chain(paste0(outPath,"/publicData/",ce10toCe11))
  chrDomainsGR_ce11<-unlist(liftOver(chrDomainsGR,chain=chainCe10toCe11))

  seqinfo(chrDomainsGR_ce11)<-ce11seqinfo
  chrDomainsGR_ce11$score<-NULL
  export(chrDomainsGR_ce11,
         con=paste0(outPath,"/publicData/chromDomains_L3_Evans2016_ce11.bed"),
         format="bed")
  file.remove(paste0(outPath,"/publicData/",basename(chromDomainsURL)))
}



states<-import.bed("./publicData/chromStates_L3_Evans2016_ce11.bed")
domains<-import.bed("./publicData/chromDomains_L3_Evans2016_ce11.bed")
seqlevels(domains)<-seqlevels(states)

stateClrs<-c("#fe0003","#f59745","#008100","#74943a",
             "#c4d69c","#05ff7f","#ceff65","#fd0082",
             "#ff70d0","#ffb6c4","#7f00ff","#1900fe",
             "#528dd4","#814006","#b8cde4","#808080",
             "#dcdac3","#c4bc97","#938b54","#141414")

domainClrs<-c("#cc0000","#5b5b5b","#bcbcbc")

if(!file.exists("./publicData/chromStates_L3_Evans2016_ce11_rgb.bed")){
  stateRGB<-apply(col2rgb(stateClrs),2,paste,collapse=",")
  states2bed<-states
  md<-data.frame(name=states2bed$score, score=states2bed$score,
                 strand=".",
                 thickStart=start(states), thickEnd=end(states),
                 itemRGB=stateRGB[states$score],blockCount="1",
                 blockSizes=width(states),blockStarts=start(states))
  mcols(states2bed)<-md
  trackLine<-'track name="chromStates" description="Evans et al. 2016" visibility=1 itemRgb="On"\n'
  export.bed(states2bed,"./publicData/chromStates_L3_Evans2016_ce11_rgb.bed")
  states2bed<-read.delim("./publicData/chromStates_L3_Evans2016_ce11_rgb.bed",header=F)
  states2bed<-cbind(states2bed,md[,c("itemRGB")])
  write.table(states2bed,
              file="./publicData/chromStates_L3_Evans2016_ce11_rgb.bed",
              sep="\t",row.names=F,col.names=F,quote=F)
  # add trackline (but not necessary)
  fConn <- file("./publicData/chromStates_L3_Evans2016_ce11_rgb.bed", 'r+')
  Lines <- readLines(fConn)
  writeLines(c(trackLine, Lines), con = fConn)
  close(fConn)
}

if(!file.exists("./publicData/chromDomains_L3_Evans2016_ce11_rgb.bed")){
  domainRGB<-apply(col2rgb(domainClrs),2,paste,collapse=",")
  domains2bed<-domains
  domains2bed$score<-as.numeric(factor(domains$name,levels=c("active","regulated","border")))
  md<-data.frame(name=domains2bed$name, score=domains2bed$score,
                 strand=".",
                 thickStart=start(domains), thickEnd=end(domains),
                 itemRGB=domainRGB[domains2bed$score],blockCount="1",
                 blockSizes=width(domains),blockStarts=start(domains))
  mcols(domains2bed)<-md
  trackLine<-'track name="chromDomains" description="Evans et al. 2016" visibility=1 itemRgb="On"\n'
  export.bed(domains2bed,"./publicData/chromDomains_L3_Evans2016_ce11_rgb.bed")
  domains2bed<-read.delim("./publicData/chromDomains_L3_Evans2016_ce11_rgb.bed",header=F)
  domains2bed<-cbind(domains2bed,md[,c("itemRGB")])
  write.table(domains2bed,
              file="./publicData/chromDomains_L3_Evans2016_ce11_rgb.bed",
              sep="\t",row.names=F,col.names=F,quote=F)
  # add trackline (but not necessary)
  fConn <- file("./publicData/chromDomains_L3_Evans2016_ce11_rgb.bed", 'r+')
  Lines <- readLines(fConn)
  writeLines(c(trackLine, Lines), con = fConn)
  close(fConn)
}

#split domains
if(!file.exists("./publicData/chromDomains_active_L3_Evans2016_ce11.bed")){
  for (domainType in unique(domains$name)){
    forBed<-domains[domains$name==domainType]
    seqinfo(forBed)<-seqinfo(Celegans)
    export.bed(forBed,con = paste0("./publicData/chromDomains_",
                                   domainType,"_L3_Evans2016_ce11.bed"))
  }
}



#split states
if(!file.exists("./publicData/chromStates_state1_L3_Evans2016_ce11.bed")){
  for (stateType in unique(states$score)){
    forBed<-states[states$score==stateType]
    seqinfo(forBed)<-seqinfo(Celegans)
    export.bed(forBed,con = paste0("./publicData/chromStates_state",
                                   stateType,"_L3_Evans2016_ce11.bed"))
  }
}

