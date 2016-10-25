library(spls)
library(doParallel)
library(foreach)
library(randomForest)
source("layer2script.R")
########
args = commandArgs(trailingOnly=TRUE)
args[1] = "../data/DREAM5/Network3/input data/net3_expression_data.tsv"
args[2] = "../data/DREAM5/Network3/input data/net3_transcription_factors.tsv"
args[3] = "./DREAM_results/Network3/net3_spls_top20.txt"
args[4] = 20

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
#######


all.expre <- read.table(args[1],header = T,sep = "\t")#All expression data
#all.expre <- all.expre[1:67,]
tf.list <- read.table(args[2],header=T) # TF list
tf.exp <<- subset(all.expre, select = tf.list$name)
myvars <- names(all.expre) %in% c("unknow","Experiment","Perturbations","PerturbationLevels","Treatment","DeletedGenes","OverexpressedGenes","Time","Repeat","Name","Experiment_Description")
gene.exp <- all.expre[!myvars]
myvars <- names(gene.exp) %in% tf.list$name
gene.exp <- all.expre[!myvars]

gene.ids <- data.frame(gene = 1 : dim(gene.exp)[2])

#ls <- func(1,20)# test is it works

cl<-makeCluster(8)
registerDoParallel(cl)
ls <-NULL
#dim(gene.ids)[1]
start <- Sys.time()
ls<-foreach(i=1:dim(gene.ids)[1],.combine =rbind,.packages=c("spls","randomForest")) %dopar% {
 func(i,args[4])
}
write.table(ls,args[3],row.names = FALSE,sep = "\t")
end <- Sys.time()
print(end - start)
stopCluster(cl)

