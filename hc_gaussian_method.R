library(bnlearn)
library(spls)
tf.input <- read.csv("../data/arabidopsis_128MA/At_Stem_TF1622TF_uniq_exp_HW.txt",sep = "\t")
anno <-read.csv("../data/arabidopsis_128MA/annotations.csv",sep = ",")
colnames(anno)[1]<-"Gene"
tf.input<-as.data.frame(merge(tf.input,anno,by="Gene"))
tf.input <- unique(tf.input)
tf.input <- within(tf.input, rm(Gene))
genes <- tf.input$sym
colnames(tf.input)
tf.input <- as.data.frame(t(tf.input[,-129]))
colnames(tf.input) <- genes


pw.input <- read.csv("../data/arabidopsis_128MA/lignin_7_pwg.csv",sep = ",")
pw.input<-as.data.frame(merge(pw.input,anno,by="Gene"))
pw.input <- unique(pw.input)
pw.input <- within(pw.input, rm(Gene))
genes<-pw.input$sym
pw.input <- as.data.frame(t(pw.input[,-129]))
colnames(pw.input) <- genes

topTFs = 20
combined_df <- data.frame()
wl <-data.frame()
for (i in 1:dim(pw.input)[2]){
  yVal <- pw.input[,i]
  yName <- names(pw.input)[i]
  temp <- NULL
  cv <-cv.spls(tf.input,pw.input[,i],eta = seq(0.1,0.9,0.1),K = c(5:10), kappa=0.5, select="pls2", fit="simpls",scale.x=TRUE, scale.y=TRUE, plot.it=F)
  f <- spls(tf.input,pw.input[,i],eta = cv$eta.opt,K = cv$K.opt)
  coef.f <- coef(f)
  temp <- as.data.frame(coef.f)
  temp$xName <- rownames(temp)
  temp <- temp[temp['V1'] > 0,]
  temp <- temp[ order(-temp[,1]),]
  #topTFs<-round((max(max(temp$V1)) - ((max(temp$V1) - min(temp$V1)) * 0.2))*100)
  temp['TO'] <- rep(yName, times = dim(temp)[1])
  temp <- temp[1:topTFs,]
  combined_df<-rbind(combined_df,temp)
  wl<-rbind(wl,data.frame(V1=combined_df$xName,V2=combined_df$TO))
  wl<-rbind(wl,as.data.frame(t(combn(temp$xName,2))))
  
  
  
}
#create node dataset
top20x7<-as.data.frame(combined_df$xName)
selected<-unique(top20x7)
colnames(selected)[1]<-"TF"
tf.input<-as.data.frame(t(tf.input))
TF<-rownames(tf.input)
tf.input$TF<-TF
selected<-merge(selected,tf.input,by = "TF")
selected.tf<-selected$TF
selected<-as.data.frame(t(selected))

#create all edge list
selected<-selected[-1,]
colnames(selected)<-selected.tf
network.nodes<-cbind(selected,pw.input)
al<-data.frame(t(combn(colnames(network.nodes),2)))
al$id <-paste(al$X1,al$X2,sep="_")
al$id2<-paste(al$X2,al$X1,sep="_")


#prepare whitelist
wl$id <-paste(wl$V1,wl$V2,sep="_")
wl$id2<-paste(wl$V2,wl$V1,sep="_")


filter1<-al[ !(al$id %in% wl$id), ]
bl<-filter1[ !(filter1$id2 %in% wl$id2), ]
bl <- within(bl, rm(id,id2))

df<-data.frame()
for (i in 1:length(colnames(pw.input))){
  from <-rep(colnames(pw.input)[i],length(colnames(network.nodes)))
  to <-colnames(network.nodes)
  df<-rbind(df,data.frame(from,to))
}

networkNodes<-read.csv("./test.csv")

colnames(bl)<-c("from","to")
bl<-rbind(bl,df)
# hc come up with a network based on the values.
bn5<-mmhc(networkNodes,whitelist = NULL, blacklist = bl, score = "loglik-g",debug = FALSE, restart = 0, perturb = 20, max.iter = Inf, optimized = TRUE)
graphviz.plot(bn5,layout = "dot",shape = "circle", main = NULL, sub = NULL,highlight = list(nodes=c("SND1","SND2","SND3","MYB103","MYB43","MYB46","MYB63","MYB85","NST1","NST2","NAC042"),fill="red") )

