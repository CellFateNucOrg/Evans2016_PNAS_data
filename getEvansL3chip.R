library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)

urls<-read.csv("EvansL3urls.csv")

options(timeout = max(600, getOption("timeout")))
workDir="."
if(!dir.exists(paste0(workDir,""))){
  dir.create(paste0(workDir,""))
}

# # get ce10
# ce10seqinfo<-seqinfo(BSgenome.Celegans.UCSC.ce10::Celegans)
# ws220seqinfo<-ce10seqinfo
# seqnames(ws220seqinfo)<-gsub("chr","CHROMOSOME_",seqnames(ws220seqinfo))
# seqnames(ws220seqinfo)<-gsub("_M","_MtDNA",seqnames(ws220seqinfo))
# genome(ws220seqinfo)<-"ws220"

# get ce6
ce6seqinfo<-seqinfo(BSgenome.Celegans.UCSC.ce6::Celegans)
ws190seqinfo<-ce6seqinfo
seqnames(ws190seqinfo)<-gsub("chr","",seqnames(ce6seqinfo))
seqnames(ws190seqinfo)<-gsub("M","MtDNA",seqnames(ws190seqinfo))
genome(ws190seqinfo)<-"ws190"

# get ce4



# # get ce10 chain file
# ce10toCe11url<-"http://hgdownload.soe.ucsc.edu/goldenPath/ce10/liftOver/ce10ToCe11.over.chain.gz"
# ce10toCe11<-"ce10Toce11.over.chain"
# download.file(ce10toCe11url,paste0(workDir,"/",ce10toCe11,".gz"))
# system(paste0("gunzip ",workDir,"/",ce10toCe11,".gz"))
# file.remove(paste0(workDir,"/",ce10toCe11,".gz"))
# chainCe10toCe11<-import.chain(paste0(workDir,"/",ce10toCe11))


# get ce6 chain file
ce6toCe11url<-"https://hgdownload.soe.ucsc.edu/goldenPath/ce6/liftOver/ce6ToCe11.over.chain.gz"
ce6toCe11<-"ce6ToCe11.over.chain"
download.file(ce6toCe11url,paste0(workDir,"/",ce6toCe11,".gz"))
system(paste0("gunzip ",workDir,"/",ce6toCe11,".gz"))
file.remove(paste0(workDir,"/",ce6toCe11,".gz"))
chainCe6toCe11<-import.chain(paste0(workDir,"/",ce6toCe11))

# new2old<-c(CHROMOSOME_I="chrI",
#            CHROMOSOME_II="chrII",
#            CHROMOSOME_III="chrIII",
#            CHROMOSOME_IV="chrIV",
#            CHROMOSOME_V="chrV",
#            CHROMOSOME_X="chrX",
#            CHROMOSOME_MtDNA="chrM")

ce6DownloadAndLiftover<-function(url,genomeVer,workDir=".",chainFile, ws190seqinfo){
  if(!file.exists(paste0(workDir,"/",basename(url))) & !file.exists(gsub("\\.gz","",paste0(workDir,"/",basename(url))))){
  download.file(url,paste0(workDir,"/",basename(url)))
  }
  if(!file.exists(gsub("\\.gz","",paste0(workDir,"/",basename(url))))){
    system(paste0("gunzip ",paste0(workDir,"/",basename(url))))
  }
  wg<-rtracklayer::import.wig(gsub("\\.gz","",paste0(workDir,"/",basename(url))))
  if(grepl("WS",genomeVer) & sum(grepl("chrIII",seqnames(wg)))<1){
    GenomeInfoDb::seqlevels(wg)<-GenomeInfoDb::seqlevels(ws190seqinfo)
    GenomeInfoDb::seqlevelsStyle(wg)<-"ucsc"
  }
  GenomeInfoDb::seqinfo(wg)<-GenomeInfoDb::seqinfo(BSgenome.Celegans.UCSC.ce6::Celegans)
  wg<-GenomicRanges::trim(wg)
  wg<-GenomicRanges::resize(wg,width=1,fix="center")
  wgce11<-unlist(rtracklayer::liftOver(wg,chain=chainFile))
  GenomeInfoDb::seqlevels(wgce11)<-GenomeInfoDb::seqlevels(BSgenome.Celegans.UCSC.ce11::Celegans)
  GenomeInfoDb::seqinfo(wgce11)<-GenomeInfoDb::seqinfo(BSgenome.Celegans.UCSC.ce11::Celegans)
  rtracklayer::export.bw(wgce11,paste0(workDir,"/",gsub("\\.wig.gz","_ce11.bw",basename(url))))
  file.remove(gsub("\\.gz","",paste0(workDir,"/",basename(url))))
}

i=31
url=urls$url[i]
genomeVer=urls$genomeVer[i]
chainFile=chainCe6toCe11
for(i in 31:nrow(urls)){
  ce6DownloadAndLiftover(urls$url[i],urls$genomeVer[i],chainFile=chainCe6toCe11,
                         ws190seqinfo=ws190seqinfo)
}
