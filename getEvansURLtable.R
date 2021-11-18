library(tidyr)
library(stringr)

df<-read.csv("TableS3_Evans2016.csv")
names(df)<-c("target","stage","GEO","used")

l3<-df[df$stage=="L3" & df$used=="x",]
l3$rep1<-NA
l3$rep2<-NA
l3$genomeVer<-NA
i=6
for(i in 17:length(l3$GEO)){
  matrixFile<-paste0(l3$GEO[i],"_series_matrix.txt")
  if(!file.exists(paste0(matrixFile,".gz")) & !file.exists(matrixFile)){
    download.file(paste0("ftp://ftp.ncbi.nlm.nih.gov/geo/series/",substr(l3$GEO[i],1,5),"nnn/",l3$GEO[i],"/matrix/",matrixFile,".gz"),destfile=paste0(matrixFile,".gz"))
  }
  if(!file.exists(matrixFile)){
    system(paste0("gunzip ",matrixFile,".gz"))
  }
  genomeVer<-system(paste0("grep 'WS' ",matrixFile),intern=T)
  l3$genomeVer[i]<-str_extract(genomeVer, "WS[[:digit:]]{3}")[1]
  if(is.na(l3$genomeVer[i])){
    gv<-system(paste0("grep 'ce[1-9]' ",matrixFile),intern=T)
    l3$genomeVer[i]<-str_extract(gv,"ce[:digit:]{1,2}")[1]
  }
  line<-system(paste0("grep ftp.*\\.wig\\.gz ", matrixFile),intern=T)
  urls<-unlist(strsplit(line,'\t\"'))
  urls<-urls[!grepl("Input",urls)]
  idxs<-grep("\\.wig\\.gz",urls)
  l3$rep1[i]<-gsub('\"','',urls[idxs[1]])
  l3$rep2[i]<-gsub('\"','',urls[idxs[2]])
  file.remove(matrixFile)
}

l3wide<-l3 %>% pivot_longer(cols=c("rep1","rep2"),names_to="replicate",values_to="url")

l3wide$used<-NULL

write.table(l3wide,file="EvansL3urls.csv",sep=",",row.names=F)
