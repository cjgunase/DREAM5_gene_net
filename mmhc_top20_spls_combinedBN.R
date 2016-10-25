library(bnlearn)
library(spls)
tf.input <- read.csv("../data/arabidopsis_128MA/At_Stem_TF1622TF_uniq_exp_HW.txt",sep = "\t")
anno <-read.csv("../data/arabidopsis_128MA/annotations.csv",sep = ",")
colnames(anno)[1]<-"Gene"
tf.input<-as.data.frame(merge(tf.input,anno,by="Gene"))
tf.input <- unique(tf.input)
tf.input <- within(tf.input, rm(Gene))
genes <- tf.input$sym
tf.input <- as.data.frame(t(tf.input[,-129]))
colnames(tf.input) <- genes

write.csv(tf.input,"tf_input.csv")
pw.input <- read.csv("../data/arabidopsis_128MA/lignin_7_pwg.csv",sep = ",")
pw.input<-as.data.frame(merge(pw.input,anno,by="Gene"))
pw.input <- unique(pw.input)
pw.input <- within(pw.input, rm(Gene))
genes<-pw.input$sym
pw.input <- as.data.frame(t(pw.input[,-129]))
colnames(pw.input) <- genes
# 
# write.csv(pw.input,"pw_input.csv")

#tf.input<-read.csv("../data/synthetic_wingping/yeast1_tf_data_200.csv")
#pw.input<-read.csv("../data/synthetic_wingping/yeast1_pathway_data_200.csv")

########################################################################
#load TF and pathway from lists
tf.list <- read.csv("../data/human_data/selected/TF_list_2179.txt")
pw.list1<- read.csv("../data/human_data/selected/PW_gene_list_1.txt")

setdiffVar <-setdiff(tf.list$TF,pw.list1$PW)
tf.list <- data.frame(TF = as.vector(setdiffVar))


human.exp <- read.csv("../data/human_data/selected/human_stem_cell_expression_116_sample.txt",sep="\t")

tf.exp<-as.data.frame(merge(tf.list,human.exp,by.x="TF", by.y="X.Gene", all.x=TRUE))
tf.list<-as.vector(tf.exp$TF)
temp<-as.data.frame(t(tf.exp[,-1]))
colnames(temp) <-tf.list
tf.input <- temp

pw.exp <-as.data.frame(merge(pw.list1,human.exp,by.x="PW", by.y="X.Gene", all.x=TRUE))
pw.list<-as.vector(pw.exp$PW)
temp<-as.data.frame(t(pw.exp[,-1]))
colnames(temp) <-pw.list
pw.input <- temp
#########################################################################

topTFs = 20 #not to do cut off
combined_df <- data.frame()
wl <-data.frame()
topscore <-data.frame()
for (i in 1:dim(pw.input)[2]){
  yVal <- pw.input[,i]
  yName <- names(pw.input)[i]
  temp <- NULL
  cv <-cv.spls(tf.input,yVal,eta = seq(0.1,0.9,0.1),K = c(5:10), kappa=0.5, select="pls2", fit="simpls",scale.x=TRUE, scale.y=TRUE, plot.it=F)
  f <- spls(tf.input,yVal,eta = cv$eta.opt,K = cv$K.opt)
  coef.f <- coef(f)
  temp <- as.data.frame(coef.f)
  temp$xName <- rownames(temp)
  temp <- temp[temp['V1'] > 0,]
  temp <- temp[ order(-temp[,1]),]
  #topTFs<-round((max(max(temp$V1)) - ((max(temp$V1) - min(temp$V1)) * 0.2))*100)
  temp['TO'] <- rep(yName, times = dim(temp)[1])
  temp <- temp[1:topTFs,]
  temp<-na.omit(temp)
  topscore<-rbind(topscore,temp)
  #write.csv(temp,"CYP98A3_varImp.csv")
  network.nodes<-cbind(tf.input[temp$xName],pw.input[yName])
  
  df<-data.frame()
  from <-rep(colnames(pw.input)[i],length(colnames(network.nodes)))
  to <-colnames(network.nodes)
  df<-rbind(df,data.frame(from,to))
  
  bn5<-mmhc(network.nodes,whitelist = NULL, blacklist = df, score = "loglik-g",debug = FALSE, restart = 0, perturb = 20, max.iter = Inf, optimized = TRUE)
  strength.1<-arc.strength(bn5, network.nodes,  criterion = "loglik-g", debug = FALSE)
  combined_df<-rbind(combined_df,strength.1)

  
}

write.csv(topscore,"./top20TFscore_each_pathway_file2.csv")
write.csv(combined_df,"./edge_file_top20TF_pathway_file2.csv")
#strength.plot(bn5,strength.1)


# TP <-read.table("../data/synthetic_wingping/yeast1_positive_tf.txt")
# foundedPTF <- intersect(unique(inferredtop100$from,inferredtop100$to),TP$V1)
# 
# dt1<-synthtic_gold_100[is.element(synthtic_gold_100$from,foundedPTF ),]
# 
# dt2<-inferredtop100[is.element(inferredtop100$from,foundedPTF ),]
# 
# sum(dt2$from %in% dt1$from & dt2$to %in% dt1$to)


