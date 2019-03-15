library(strataG)


Sys.setenv(http_proxy='http://proxy.swmed.edu:3128')
Sys.setenv(https_proxy='https://proxy.swmed.edu:3128')
Sys.setenv(PATH="/usr/local/sbin:/usr/local/bin:/usr/bin:/usr/sbin:/sbin:/bin:/home/xcao/p/bin")
Sys.setenv(PATH="/usr/local/sbin:/usr/local/bin:/usr/bin:/usr/sbin:/sbin:/bin:/home/xcao/p/anaconda3_5.2.0/bin")

source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("treeio")
biocLite("ggtree")


library(ggplot2)

filename = '/home/xcao/w/20180905Junonia_coenia/20181010Trees/20181026MergeSample/fraglen10000_tree.sum.csv'
xdf = read.csv(filename,sep = '\t',stringsAsFactors = FALSE)
xdf$outgroupTogether <- xdf$outgroupTogether =='True'
xdf$Boots75 <-  xdf$minBoots >= 75

xdf$minBootsRange <- cut(xdf$minBoots,seq(0,100,20))
#ggplot(xdf,aes(x=minBootsRange)) + geom_boxplot(aes(y=outgroup_relative_distance))
ggplot(xdf) + geom_histogram(aes(x=outgroup_relative_distance))
ggplot(xdf) + geom_histogram(aes(x=outgroup_relative_distance)) + facet_grid(outgroupTogether~.)
ggplot(xdf) + geom_histogram(aes(x=outgroup_relative_distance)) + facet_grid(Boots75~.)
ggplot(xdf) + geom_histogram(aes(x=outgroup_relative_distance)) + facet_grid(outgroupTogether~minBootsRange)
