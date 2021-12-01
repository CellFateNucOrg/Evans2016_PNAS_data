library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)

urls<-read.csv("EvansL3urls.csv")

options(timeout = max(600, getOption("timeout")))
workDir="."
if(!dir.exists(paste0(workDir,""))){
  dir.create(paste0(workDir,""))
}

WBtoUCSCseqlevels<-function(gr){
  if(any(grepl("CHROMOSOME_",seqlevels(gr)))){ # wrombase CHROMOSOME_ style
    seqlevels(gr)<-gsub("CHROMOSOME_","chr",seqlevels(gr))
    seqlevels(gr)<-gsub("chrMtDNA","chrM",seqlevels(gr))
  } else if(!any(grepl("chr",seqlevels(gr)))){ # wormbase I II..MtDNA style
    seqlevels(gr)<-paste0("chr",seqlevels(gr))
    seqlevels(gr)<-gsub("chrMtDNA","chrM",seqlevels(gr))
  }
  return(gr)
}


# # get ce10
# ce10seqinfo<-seqinfo(BSgenome.Celegans.UCSC.ce10::Celegans)
# ws220seqinfo<-ce10seqinfo
# seqnames(ws220seqinfo)<-gsub("chr","CHROMOSOME_",seqnames(ws220seqinfo))
# seqnames(ws220seqinfo)<-gsub("_M","_MtDNA",seqnames(ws220seqinfo))
# genome(ws220seqinfo)<-"ws220"

# get ce6
ce6seqinfo<-seqinfo(BSgenome.Celegans.UCSC.ce6::Celegans)

# get ce6 chain file
ce6toCe11url<-"https://hgdownload.soe.ucsc.edu/goldenPath/ce6/liftOver/ce6ToCe11.over.chain.gz"

ce6toCe11<-"ce6ToCe11.over.chain"
download.file(ce6toCe11url,paste0(workDir,"/",ce6toCe11,".gz"))
system(paste0("gunzip ",workDir,"/",ce6toCe11,".gz"))
file.remove(paste0(workDir,"/",ce6toCe11,".gz"))
chainCe6toCe11<-import.chain(paste0(workDir,"/",ce6toCe11))

# get WS180
ws180chromsizes<-read.delim("c_elegans.WS180.genomic.chrom.sizes",
                            header=F)
ws180seqinfo<-ce6seqinfo
seqnames(ws180seqinfo)<-ws180chromsizes$V1
seqlengths(ws180seqinfo)<-ws180chromsizes$V2
genome(ws180seqinfo)<-"WS180"
ws180seqinfo<-WBtoUCSCseqlevels(ws180seqinfo)
chainWS180toWS235<-import.chain("WS180toWS235_ucscFormat.chain")

chainWS180toWS235


# get WS190
ws190chromsizes<-read.delim("c_elegans.WS190.genomic.chrom.sizes",
                            header=F)
ws190seqinfo<-ce6seqinfo
seqnames(ws190seqinfo)<-ws190chromsizes$V1
seqlengths(ws190seqinfo)<-ws190chromsizes$V2
genome(ws190seqinfo)<-"WS190"
ws190seqinfo<-WBtoUCSCseqlevels(ws190seqinfo)
chainWS190toWS235<-import.chain("WS190toWS235_ucscFormat.chain")


# get WS235
ws235chromsizes<-read.delim("c_elegans.WS235.genomic.chrom.sizes",
                            header=F)
ws235seqinfo<-ce6seqinfo
seqnames(ws235seqinfo)<-ws235chromsizes$V1
seqlengths(ws235seqinfo)<-ws235chromsizes$V2
genome(ws235seqinfo)<-"WS235"
ws235seqinfo<-WBtoUCSCseqlevels(ws235seqinfo)






evansDownloadAndLiftover<-function(url,genomeVer,workDir=".",chainFile,
                                   oldseqinfo, newseqinfo){
  if(!file.exists(paste0(workDir,"/",basename(url))) &
     !file.exists(gsub("\\.gz","",paste0(workDir,"/",basename(url))))){
    download.file(url,paste0(workDir,"/",basename(url)))
  }
  if(!file.exists(gsub("\\.gz","",paste0(workDir,"/",basename(url))))){
    system(paste0("gunzip ",paste0(workDir,"/",basename(url))))
  }
  wg<-rtracklayer::import.wig(gsub("\\.gz","",paste0(workDir,"/",basename(url))))
  # if it is a GRangeslist
  if(length(wg)<10){
    print("converting grl to gr")
    names(wg)
    names(wg)<-gsub("^.*_chr","chr",names(wg))
    wg<-unlist(wg)
    rownames(wg)<-NULL
    wg<-GenomicRanges::resize(wg,width=1,fix="start")
    # ol<-findOverlaps(wg)
    # ol<-ol[queryHits(ol)!=subjectHits(ol),]
    # if(length(ol)!=0){
    #   "remove overlaps from grl"
    #   wg<-wg[-queryHits(ol)]
    # }
  }
  #table(width(disjoin(wg)))
  wg<-GenomicRanges::resize(wg,width=1,fix="start")
  ol<-GenomicRanges::findOverlaps(wg,drop.self=T,drop.redundant=T)
  length(ol)
  wg<-WBtoUCSCseqlevels(wg)
  #wg<-GenomicRanges::trim(wg)
  if(genomeVer=="WS180"){
    GenomeInfoDb::seqlevels(wg)<-GenomeInfoDb::seqlevels(ws180seqinfo)
    GenomeInfoDb::seqinfo(wg)<-ws180seqinfo
    wgce11<-unlist(rtracklayer::liftOver(wg,chain=chainFile))
    seqinfo(wgce11)<-newseqinfo
  } else if(genomeVer=="WS190"){
    GenomeInfoDb::seqlevels(wg)<-GenomeInfoDb::seqlevels(ws190seqinfo)
    GenomeInfoDb::seqinfo(wg)<-ws190seqinfo
    wgce11<-unlist(rtracklayer::liftOver(wg,chain=chainFile))
    seqinfo(wgce11)<-newseqinfo
  } else if(genomeVer=="ce6"){
    GenomeInfoDb::seqlevels(wg)<-GenomeInfoDb::seqlevels(ce6seqinfo)
    GenomeInfoDb::seqinfo(wg)<-ce6seqinfo
    wgce11<-unlist(rtracklayer::liftOver(wg,chain=chainFile))
    seqinfo(wgce11)<-newseqinfo
  }
  # extend ranges to cover gaps
  prcd<-precede(wgce11,select="all")
  dd<-distance(wgce11[queryHits(prcd)],wgce11[subjectHits(prcd)])
  dd1<-sapply(dd,min,50)
  wgce11<-resize(wgce11[queryHits(prcd)],width=width(wgce11[queryHits(prcd)])+dd1,fix="start")
  rtracklayer::export.bw(wgce11,paste0(workDir,"/",gsub("\\.wig.gz","_ce11.bw",basename(url))))
  file.remove(gsub("\\.gz","",paste0(workDir,"/",basename(url))))
}

i=1
url=urls$url[i]
genomeVer=urls$genomeVer[i]
if(genomeVer=="ce6"){
  chainFile=chainCe6toCe11
  oldseqinfo=ce6seqinfo
} else if(genomeVer=="WS180"){
  chainFile=chainWS180toWS235
  oldseqinfo=ws180seqinfo
} else if(genomeVer=="WS190"){
  chainFile=chainWS190toWS235
  oldseqifo=ws190seqinfo
}
newseqinfo=seqinfo(BSgenome.Celegans.UCSC.ce11::Celegans)
#newseqinfo=ws235seqinfo
for(i in 1:nrow(urls)){
  print(i)
  print(urls$url[i])
  print(urls$genomeVer[i])
  evansDownloadAndLiftover(urls$url[i],urls$genomeVer[i],
                           chainFile=chainFile,
                           oldseqinfo=oldseqinfo,
                           newseqinfo=newseqinfo)
}
url=urls$url[i]
genomeVer=urls$genomeVer[i]
