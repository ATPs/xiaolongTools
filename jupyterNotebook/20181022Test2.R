library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)

folder_files <- '/home/xcao/w/20180905Junonia_coenia/20190131New/structure/results109/'
result_files <- list.files(folder_files,full.names = FALSE)
basetypes <- c("zChrCDS","CDS","zChr")

f.getbasetype <- function(x){
  for (keyword in basetypes){
    if (grepl(keyword,x)) return(keyword)
  }
}

f.getSplitgroup <- function(x){
  if (grepl('.split',x)) {
    return(as.integer(gsub("(.*split\\d*\\.)|(\\.STRUCTUREresult.*)","",x)))
  } else return(0)
}

f.getRandint <- function(x){
  #x="Calephelis266.zChr.split5.4.STRUCTUREresult.s1000.pop8_f" return 1000
  return(as.integer(gsub('(.*STRUCTUREresult.s)|(.pop.*_f)','',x)))
}

f.getPopulation <- function(x) {
  return(as.integer(gsub('(.*STRUCTUREresult.s\\d*.pop)|(_f)','',x)))
}

df_files <- data.frame(baseType=sapply(result_files, f.getbasetype),
                       splitID=sapply(result_files, f.getSplitgroup),
                       seed=sapply(result_files, f.getRandint),
                       population=sapply(result_files,f.getPopulation),
                       filename=file.path(folder_files,result_files), stringsAsFactors = FALSE)

df_files <- df_files[order(df_files$baseType, df_files$population,df_files$seed,df_files$splitID),]

F.renameHead <- function(x){
  x[1] = 'sample'
  x[2:length(x)] = paste0('pop',1:(length(x)-1))
  return(x)
}

filename2 <- '/home/xcao/w/20180905Junonia_coenia/20180919Info/20190130Junonia-all-sample-summary.xlsx'
xdf_info <- read_excel(filename2,sheet = 'Sheet1')
# xdf_info <- xdf_info[!is.na(xdf_info$order_266mitoTree),]
# v_id2Treename <- xdf_info$Treename
# names(v_id2Treename) <- sapply(xdf_info$ID, function(x)substring(x,2))
v_id2species <- xdf_info$`Taxon name`
names(v_id2species) <- sapply(xdf_info$ID, function(x)substring(x,2))
species_in_order <- xdf_info[order(xdf_info$species),]$ID
species_in_order <- substr(species_in_order,2,length(species_in_order))
species_in_order <- species_in_order[species_in_order %in% xdf_s$sample]
species_in_order <- c(species_in_order,'18098H12F', '18126E01F', '18126E02F', '18126E03F', '18126E04F', '18126E05F', '18126E06F', '18126E07F', '18126E08F', '18126E09F')

oneline = df_files[1,]
xdf_s <- read.csv(oneline$filename,header = FALSE,sep = '\t', stringsAsFactors = FALSE)
colnames(xdf_s) <- F.renameHead(colnames(xdf_s))
xdf_s$species <- v_id2species[xdf_s$sample]
#xdf_s$treename <- v_id2Treename[xdf_s$sample]
xdf_s.clear <- gather(xdf_s,key = 'pop',value = 'fraction',-sample,-species)
xdf_s.clear <- xdf_s.clear[order(xdf_s.clear$species,xdf_s.clear$pop,-xdf_s.clear$fraction),]

p <- ggplot(xdf_s.clear,aes(x=sample,y=fraction,fill=pop)) + geom_bar(stat = 'identity') + scale_x_discrete(limits=species_in_order)



F.plotOneline <- function(oneline,xtext=FALSE){
  xdf_s <- read.csv(oneline$filename,header = FALSE,sep = '\t', stringsAsFactors = FALSE)
  colnames(xdf_s) <- F.renameHead(colnames(xdf_s))
  xdf_s$species <- v_id2species[xdf_s$sample]
  #xdf_s$treename <- v_id2Treename[xdf_s$sample]
  xdf_s.clear <- gather(xdf_s,key = 'pop',value = 'fraction',-sample,-species)
  xdf_s.clear <- xdf_s.clear[order(xdf_s.clear$species,xdf_s.clear$pop,-xdf_s.clear$fraction),]
  
  p <- ggplot(xdf_s.clear,aes(x=sample,y=fraction,fill=pop)) + geom_bar(stat = 'identity') + scale_x_discrete(limits=species_in_order[species_in_order %in% xdf_s$sample]) + ylab(paste0('seed',oneline$seed,'_',oneline$baseType,'_K=',oneline$population," ",oneline$splitID))
  if (xtext) {
    p <- p + geom_text(aes(x=sample,y=0,label=paste(sample,species)),angle=90,hjust=0)+ theme(legend.position="none", axis.text.x= element_blank(), axis.title.x = element_blank())
    #p <- p + theme(legend.position="none", axis.text.x= element_text(angle = 90, hjust = 1), axis.title.x = element_blank())
  }
  else {
    p <-  p + theme(legend.position="none", axis.text.x= element_blank(), axis.title.x = element_blank())
  }
  return(p)
}

F.plotMultipleLines <- function(filelines,xtext=FALSE){
  p <- list()
  for (i in 1:nrow(filelines)){
    p[[i]] <- F.plotOneline(filelines[i,],xtext = xtext)
  }
  #p[[i+1]] <- F.plotOneline(filelines[i+1,], xtext = TRUE)
  return(grid.arrange(grobs=p,ncol=1))
}

pdf('/home/xcao/w/20180905Junonia_coenia/20190131New/structure/Junonia109structure2.pdf',width = 24, height = 48)
filelines = df_files[df_files$baseType == 'CDS',]
F.plotMultipleLines(filelines = filelines,xtext = TRUE)
filelines = df_files[df_files$baseType == 'zChrCDS',]
F.plotMultipleLines(filelines = filelines,xtext = TRUE)
filelines = df_files[df_files$baseType == 'zChr',]
F.plotMultipleLines(filelines = filelines,xtext = TRUE)
dev.off()