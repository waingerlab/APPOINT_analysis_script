##### APPOINT analysis script ####
#### Analysis of calcium imaging data collected using ImageXpress Micro Confocal high content imager with automated fluidics and analyzed using MetaXpress software custom analysis pipeline

### General set up ####
# Necessary libraries
#install.packages(c("tidyverse","rstudioapi","cowplot","matrixStats","data.table","factoextra","caret","AppliedPredictiveModeling","randomForest"))
library(tidyverse) #v1.3.0, includes ggplot2_3.3.2, tibble_3.0.1, tidyr_1.0.2, readr_1.3.1, purrr_0.3.4, dplyr_0.8.5, stringr_1.4.0, forcats_0.5.0
library(rstudioapi) #v0.11
library(cowplot) #v1.0.0
library(matrixStats) #v0.56.0 
library(data.table) #v1.12.8
library(factoextra) #v1.0.7, includes lattice_0.20-35
library(caret) #v6.0-86
library(AppliedPredictiveModeling) #v1.1-7
library(randomForest) #v4.6-14
library(DescTools) #v0.99.25
theme_set(theme_cowplot())

### Custom functions ####
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), ncol = cols, nrow = ceiling(numPlots/cols)) }
  if (numPlots==1) { print(plots[[1]]) } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,layout.pos.col = matchidx$col)) } } }
normalize<-function(input,base_begin,base_end) {
  (input-mean(input[base_begin:base_end]))/(mean(input[base_begin:base_end])) }
Peak_median3<-function(input,time_of_peak) {
  t2<-time_of_peak
  t1<-t2-1
  t3<-t2+1
  median(c(input[t1],input[t2],input[t3])) }
trigger_count<-function(input,time_start,time_end,thresh) {
  if(max(input[time_start:time_end])>thresh) {
    return("Yes") } else {
      return("No") } }
which.median<-function(df) {
  which(df==median(df))}

### Initial setup of data and script options ####
# Set options, load data
pdf.options(useDingbats=F, width=10.5, height=8)
setwd("") #insert path to folder containing Data and Metadata folders to be analyzed 
files<-list.files(path="Data/",pattern=".csv",full.names=T)
data<-lapply(files,fread)
# Load csv of metadata information
Meta<-list.files("Metadata/",pattern=".csv",full.names=T)
Meta<-fread(Meta)
# Create directory for saving analysis results and go there
if(!(dir.exists("Analysis"))) {dir.create("Analysis")}
setwd('Analysis')


### General data processing ####
# Data pre-processing and isolation of intensity values
names<-str_split(files,pattern="_",simplify=T) %>% .[,c(-1,-4)] %>%
  str_c(.[,2],sep="_") %>% .[c(1:length(files))]
names(data)<-names
for( i in seq_along(data)){
  data[[i]]$key <- rep(names[i],nrow(data[[i]]))}
Data<-as_tibble(do.call(rbind,data))
Data<-Data[,c(1,4,5,6,7)]
colnames(Data)<-c("Plane","Object","Area","Intensity","key")
Intensity<-subset(Data,Object>0,select=c(Plane,Object,Intensity,key))
# Transform Intensity from long to wide format and normalize values to initial baseline intensity
Intensity<-unite(Intensity,col="name", key, Object,sep="_")
Int<-spread(Intensity,key=name, value=Intensity)
IntNorm<-cbind(Int[,1],map_df(Int[,-1], ~normalize(.x,base_begin=1,base_end=10)))
##plot traces of IntNorm data directly into pdf; each page as well
INP<-IntNorm %>% gather(key="name",value="Int",-Plane)
INP$key<-str_c(INP$name %>% str_split(pattern="_",simplify=T) %>% .[,1],
               INP$name %>% str_split(pattern="_",simplify=T) %>% .[,2],
               sep="_")
wells<-INP$key %>% unique()
pdf("Example traces- all cells, paged by well.pdf")
for(i in seq_along(wells)) { #i<-1
  pINP<-ggplot(Meta[1,]) +
    labs(tag=wells[i],
         caption=(Meta[str_detect(Meta$key,wells[i]),"Condition"])) +
    ylim(-0.5,5) + xlim(-1,NA) + ylab("") + xlab("") +
    geom_segment(aes(x=-0.5,xend=-0.5,y=-0.3,yend=0.2),size=1.2) +
    geom_text(aes(x=-0.8,y=-0.3,label="0.5 dF/F"),angle=90,vjust=0,hjust=0) +
    geom_segment(aes(x=-0.5,xend=4.5,y=-0.3,yend=-0.3),size=1.2) +
    geom_text(aes(x=-0.5,y=-0.32,label="5 s"),vjust=1,hjust=0) +
    geom_line(data=(INP %>% filter(key==wells[i])),
                    aes(x=Plane,y=Int,group=name),size=0.9,alpha=0.3) +
    theme(axis.line=element_blank(),axis.ticks=element_blank(),
          axis.text=element_blank()) + 
    geom_vline(xintercept=c(5,25,45),linetype=2,alpha=0.4,color="blue") 
  print(pINP)
  rm(pINP) }
dev.off()
# Calculate first derivative of smoothed intensity profiles
IntNormDiff<-map_df(IntNorm[,-1], ~diff(.x))
IntNormDiff<-cbind(Plane=IntNorm[-61,1],IntNormDiff)


### Quantify response properties for each stimulus response ####
Response<-str_split(colnames(IntNorm)[-1],pattern="_",n=3,simplify=T) %>%
  as_tibble()
colnames(Response)<-c("Plate","Well","Object")
Response$key<-str_c(Response$Plate,Response$Well,sep="_")
Meta<-Meta %>% filter(key %in% (Response$key %>% unique()))
Response$name<-colnames(IntNorm)[-1]
Response<-left_join(Response,Meta,by=c("key"))
Response$Base1<-Int[1:5,-1] %>% map_dbl(mean)
Response$Base2<-Int[21:25,-1] %>% map_dbl(mean)
Response$Base3<-Int[41:45,-1] %>% map_dbl(mean)
Response$ToP1 <- IntNorm[5:20,-1] %>% map_dbl(~which(.x==max(.x))[1]+4)
Response$ToP2 <- IntNorm[25:40,-1] %>% map_dbl(~which(.x==max(.x))[1]+24)
Response$ToP3 <- IntNorm[45:60,-1] %>% map_dbl(~which(.x==max(.x))[1]+44)
Response$Peak1<-map2_dbl(Int[,-1],Response$ToP1,~Peak_median3(.x,.y))
Response$Peak2<-map2_dbl(Int[,-1],Response$ToP2,~Peak_median3(.x,.y))
Response$Peak3<-map2_dbl(Int[,-1],Response$ToP3,~Peak_median3(.x,.y))
Response$Amp1<-(Response$Peak1-Response$Base1)/Response$Base1
Response$Amp2<-(Response$Peak2-Response$Base2)/Response$Base2
Response$Amp3<-(Response$Peak3-Response$Base3)/Response$Base3
Response$Rise1<-IntNormDiff[5:20,-1] %>% map_dbl(max)
Response$Rise2<-IntNormDiff[25:40,-1] %>% map_dbl(max)
Response$Rise3<-IntNormDiff[45:60,-1] %>% map_dbl(max)

## Split responses to each stimulus within cell; allows for individual classification ###
base<-Response %>% select(name,Base1,Base2,Base3) %>%
  gather(key="stim",value='base',-name)
base$stim<-base$stim %>% str_replace("Base","stim")
top<-Response %>% select(name,ToP1,ToP2,ToP3) %>%
  gather(key="stim",value='top',-name)
top$stim<-top$stim %>% str_replace("ToP","stim")
peak<-Response %>% select(name,Peak1,Peak2,Peak3) %>%
  gather(key="stim",value='peak',-name)
peak$stim<-peak$stim %>% str_replace("Peak","stim")
amp<-Response %>% select(name,Amp1,Amp2,Amp3) %>%
  gather(key="stim",value='amp',-name)
amp$stim<-amp$stim %>% str_replace("Amp","stim")
rise<-Response %>% select(name,Rise1,Rise2,Rise3) %>%
  gather(key="stim",value='rise',-name)
rise$stim<-rise$stim %>% str_replace("Rise","stim")
response2<-inner_join(base,top,by=c("name","stim")) %>%
  inner_join(peak,by=c("name","stim")) %>% 
  inner_join(amp,by=c("name","stim")) %>%
  inner_join(rise,by=c("name","stim"))
rm(base,peak,amp,rise,top)

### Make random forest classifier to identify responses vs non-responses ####
training<-rbind(response2 %>% filter(stim=='stim1',amp<0.3)%>%
                  sample_n(100,replace=T),
                response2 %>% filter(stim=='stim3',amp>0.1)%>%
                  sample_n(100,replace=T))
training$resp<-c(rep("No",100),rep("Yes",100))
fit<-train(resp~., data=training %>% select(resp,base,peak,amp,rise),
           method='rf',seed=5067,selectionFunction="oneSE")
### Validate rf model using 100 new subsets of responses
model_testing<-data.frame(it=c(1:100),false_neg=0,false_pos=0,
                          true_neg=0,true_pos=0)
for(i in c(1:100)) { #i<-1
  test_neg<-response2 %>% filter(stim=="stim1") %>%
    sample_n(100,replace=T)
  test_neg$model<-predict(fit,test_neg)
  model_testing$false_pos[i]<-test_neg %>% filter(model=="Yes") %>%
    nrow()
  model_testing$true_neg[i]<-test_neg %>% filter(model=="No") %>%
    nrow()
  test_pos<-response2 %>% filter(stim=="stim3",amp>0.5) %>%
    sample_n(100,replace=T)
  test_pos$model<-predict(fit,test_pos)
  test_pos$prob<-predict(fit,test_pos,type='prob')$Yes
  model_testing$true_pos[i]<-test_pos %>% filter(model=="Yes") %>%
    nrow()
  model_testing$false_neg[i]<-test_pos %>% filter(model=="No") %>%
    nrow()
  rm(test_neg,test_pos) }
mod_results<-model_testing %>% summarize(
  fp_mean=mean(false_pos) %>% round(1),fp_sd=sd(false_pos) %>% round(1),
  tp_mean=mean(true_pos) %>% round(1),tp_sd=sd(true_pos) %>% round(1),
  tn_mean=mean(true_neg) %>% round(1),tn_sd=sd(true_neg) %>% round(1),
  fn_mean=mean(false_neg) %>% round(1),fn_sd=sd(false_neg) %>% round(1))
rm(training)
saveRDS(fit,"RF model- no,stim1,under0.3- yes,stim3,over0.1.rds")
## Predict response v non-response for each stimulus
response2$model<-predict(fit,response2) %>% as.character()
response2$prob<-predict(fit,response2,type="prob")$Yes

### Use HClust to identify positive and negative responses ####
rlow<-0.5;rhigh<-1.5;alow<-0.075;ahigh<-0.225
res2<-response2 %>% filter((rise>=rlow &rise<=rhigh &amp>=alow &amp<=ahigh))
hc_res2<-res2 %>% select(amp,rise) %>% as.data.frame() %>% 
  scale(center=T,scale=T)%>%dist(method="euclidean")%>%
  hclust(method="ward.D2")
res2$cluster<-cutree(hc_res2,3)
clusters<-res2 %>% group_by(cluster) %>%
  summarize_at(vars("amp","rise"),mean)
clus_no<-c(clusters$cluster[clusters$amp %>% which.min()],
           clusters$cluster[clusters$rise %>% which.min()]) %>% unique()
res2$res<-"Yes"
res2$res[which(res2$cluster%in%clus_no)]<-"No"
response2<-left_join(response2,res2 %>% select(name,stim,res),
                     by=c("name","stim"))
response2$res[is.na(response2$res)]<-"depends"
thr_r<-res2 %>% filter(res=='No') %>% top_n(50,rise) %>%
  pull(rise) %>% median()
thr_a<-res2 %>% filter(res=='No') %>% top_n(50,amp) %>% 
  pull(amp) %>% median()
response2$res[which(response2$res=="depends"&
                      response2$rise<thr_r & response2$amp<thr_a)]<-"No"
response2$res[which(response2$res=="depends")]<-"Yes"
response2$res[is.na(response2$res)]<-"depends"
rm(thr_r,thr_a,clusters,hc_res2,res2)


### Compile results of random forest and hc response analyses within individual cells ####
r2p<-response2 %>% select(name,stim,prob) %>%
  spread(key="stim",value="prob")
colnames(r2p)<-c("name","s1_rf_prob","s2_rf_prob","s3_rf_prob")
r2f<-left_join(response2 %>% select(name,stim,model) %>%
                 spread(key="stim",value="model"),
               response2 %>% select(name,stim,res) %>% 
                 spread(key="stim",value="res"),
               by="name")
colnames(r2f)<-c("name","s1_rf","s2_rf","s3_rf","s1_hc","s2_hc","s3_hc")
r2f$rf_patt<-str_c(r2f$s1_rf,r2f$s2_rf,r2f$s3_rf,sep="_")
r2f$hc_patt<-str_c(r2f$s1_hc,r2f$s2_hc,r2f$s3_hc,sep="_")
Response<-left_join(Response,r2p,by="name")
Response<-left_join(Response,r2f,by="name")
##Make consensus call (if either says yes, call yes)
Response$stim1<-"No";Response$stim2<-"No";Response$stim3<-"No"
Response$stim1[which(Response$s1_rf=="Yes" | Response$s1_hc=="Yes")]<-"Yes"
Response$stim2[which(Response$s2_rf=="Yes" | Response$s2_hc=="Yes")]<-"Yes"
Response$stim3[which(Response$s3_rf=="Yes" | Response$s3_hc=="Yes")]<-"Yes"
Response$pattern<-str_c(Response$stim1,Response$stim2,Response$stim3,
                        sep="_")
cells_resp_any<-Response %>% filter(stim1=="Yes"|stim2=="Yes"|
                                      stim3=="Yes") %>% pull(name)
Resp_filt<-Response

### Count cells based on Yes/No response calls ####
counts<-left_join(Resp_filt %>% group_by(Plate,Well,key,Condition) %>%
                    summarize(total=length(name)),
                  Resp_filt %>% filter(stim1=="Yes") %>% group_by(key) %>%
                    summarize(noise=length(name)),
                  by="key") %>% 
  left_join(Resp_filt %>%filter(stim1=="No"&(stim2=="Yes"|stim3=="Yes"))%>% 
              group_by(key) %>% summarize(live=length(name)),
            by="key") %>% 
  left_join(Resp_filt %>%filter(stim1=="No"&(stim3=="Yes"))%>% 
              group_by(key) %>% summarize(pc=length(name)),
            by="key") %>% 
  left_join(Resp_filt %>% filter(stim1=="No"&stim2=="Yes") %>% 
              group_by(key) %>% summarize(exp=length(name)),
            by="key") %>%
  left_join(Resp_filt %>%filter(stim1=="No"&(stim2=="Yes"&stim3=="Yes"))%>%
              group_by(key) %>% summarize(dual=length(name)),
            by="key") %>% ungroup()
counts[is.na(counts)]<-0
counts$pct_nt<-(counts$noise / counts$total *100) %>% round(2)
counts$pct_el<-(counts$exp / counts$live *100) %>% round(2)
counts$pct_de<-counts$dual / counts$exp *100 %>% round(2)
counts$pct_dp<-counts$dual / counts$pc *100 %>% round(2)
counts$pct_dl<-counts$dual / counts$live *100 %>% round(2)
counts[is.na(counts)]<-0

### Output analysis results to csv files ####
write_csv(counts,str_c("APPOINT analysis results- plate #",
                       str_split(wells[1],pattern="_",simplify=T)[1],
                       "- final cell counts by well.csv"))
write_csv(Response,str_c("APPOINT analysis results- plate #",
                        str_split(wells[1],pattern="_",simplify=T)[1],
                        "- response metrics by cell.csv"))
write_csv(response2,str_c("APPOINT analysis results- plate #",
                         str_split(wells[1],pattern="_",simplify=T)[1],
                         "- response metrics by stimulus.csv"))
write_csv(model_testing,str_c("APPOINT analysis results- plate #",
                             str_split(wells[1],pattern="_",simplify=T)[1],
                             "- results of RF model testing.csv"))


### Organize final results for plotting- specific to each experiment ####
Response$stimulus<-Response$Condition %>% factor(levels=c("Vehicle","L/P","Ki","L","Ki/L"),
                                                 labels=c("saline","LPA/PGE2","LPA/PGE2","LPA/PGE2","LPA/PGE2"))
Response$treatment<-Response$Condition %>% factor(levels=c("Vehicle","L/P","Ki","L","Ki/L"),
                                                  labels=c("vehicle","vehicle","Ki-16425","L-798,106","Ki+L"))
counts$stimulus<-counts$Condition %>% factor(levels=c("Vehicle","L/P","Ki","L","Ki/L"),
                                            labels=c("saline","LPA/PGE2","LPA/PGE2","LPA/PGE2","LPA/PGE2"))
counts$treatment<-counts$Condition %>% factor(levels=c("Vehicle","L/P","Ki","L","Ki/L"),
                                             labels=c("vehicle","vehicle","Ki-16425","L-798,106","Ki+L"))
counts<-counts %>% filter(treatment!="Ki+L")
counts$baseline<-counts %>% filter(stimulus=="LPA/PGE2",treatment=="vehicle") %>% pull(pct_el) %>% mean()
counts$pct_change<-(counts$pct_el - counts$baseline) / counts$baseline *100


### Plot final results into pdf file- specific to each experiment ####
pdf(str_c("APPOINT analysis results- plate #",
          str_split(wells[1],pattern="_",simplify=T)[1],
          "- final plot results.pdf"))
ggplot(counts %>% filter(treatment=="vehicle"),
       aes(x=stimulus,y=pct_el)) +
  scale_y_continuous(limits=c(0,100),expand=c(0,0),
                     name="Activated cells (% live)") +
  xlab("Stimulus") +
  stat_summary(fun.data=mean_se,geom="crossbar",fill='grey50',alpha=0.2) +
  geom_jitter(height=0,width=0.1,size=4,alpha=0.7) +
  theme(axis.line.x=element_line(lineend='round',size=1.2),
        axis.line.y=element_line(lineend='round',size=1.2),
        axis.ticks.y=element_line(lineend='round',size=1.2),
        axis.text.x=element_text(size=15,angle=35,hjust=1),
        axis.text.y=element_text(size=15,hjust=1),
        axis.title=element_text(size=15),
        axis.ticks.x=element_blank())
ggplot(counts %>% filter(treatment%in%c("Ki-16425","L-798,106")),
       aes(x=treatment,y=pct_change)) +
  scale_y_continuous(limits=c(-100,100),expand=c(0,0),
                     name="Drug effect (% change)") +
  scale_x_discrete(name="Treatment") +
  geom_hline(yintercept=0,size=1.2) +
  stat_summary(fun.data=mean_se,geom="crossbar",fill='grey50',alpha=0.2) +
  geom_jitter(height=0,width=0.1,size=4,alpha=0.7) +
  theme(axis.line.x=element_blank(),
        axis.line.y=element_line(lineend='round',size=1.2),
        axis.ticks.y=element_line(lineend='round',size=1.2),
        axis.text.x=element_text(size=15,angle=35,hjust=1),
        axis.text.y=element_text(size=15,hjust=1),
        axis.title=element_text(size=15),
        axis.ticks.x=element_blank())
ggplot(Response %>% filter(treatment!="Ki+L"),
       aes(x=Amp2,color=interaction(stimulus,treatment,sep="_"))) +
  stat_ecdf(pad=F,size=1.2,alpha=0.5) +
  scale_color_manual(values=c("grey50","green3","darkblue","cyan3"),
                     name="Stimulus_Treatment") +
  scale_x_continuous(limits=c(NA,10),expand=c(0,0),
                     name="Peak response amplitude (dF/F)") +
  scale_y_continuous(limits=c(0,1),expand=c(0,0),breaks=c(0,0.25,0.5,0.75,1),
                    labels=c(0,25,50,75,100),name="Responses (%)") +
  theme(axis.line.x=element_line(lineend='round',size=1.2),
        axis.line.y=element_line(lineend='round',size=1.2),
        axis.ticks.y=element_line(lineend='round',size=1.2),
        axis.text.x=element_text(size=15,angle=35,hjust=1),
        axis.text.y=element_text(size=15,hjust=1),
        axis.title=element_text(size=15),
        axis.ticks.x=element_line(lineend='round',size=1.2))
ggplot(Response %>% filter(treatment!="Ki+L"),
       aes(x=Amp3,color=treatment)) +
  stat_ecdf(pad=F,size=1.2,alpha=0.5) +
  scale_color_manual(values=c("grey50","darkblue","cyan3"),
                     name="Treatment") +
  scale_x_continuous(limits=c(NA,7),expand=c(0,0),breaks=seq(from=0,to=7,by=1),
                     name="Peak KCl response amplitude (dF/F)") +
  scale_y_continuous(limits=c(0,1),expand=c(0,0),breaks=c(0,0.25,0.5,0.75,1),
                     labels=c(0,25,50,75,100),name="Responses (%)") +
  theme(axis.line.x=element_line(lineend='round',size=1.2),
        axis.line.y=element_line(lineend='round',size=1.2),
        axis.ticks.y=element_line(lineend='round',size=1.2),
        axis.text.x=element_text(size=15,angle=35,hjust=1),
        axis.text.y=element_text(size=15,hjust=1),
        axis.title=element_text(size=15),
        axis.ticks.x=element_line(lineend='round',size=1.2))
ggplot(counts %>% filter(treatment!="Ki+L")) +
  scale_x_continuous(expand=c(0,0),limits=c(0,600),
                     name="Cells per well") +
  scale_y_continuous(expand=c(0,0),name="") +
  stat_density(aes(x=total,y=..scaled..),fill='grey50',adjust=0.8,alpha=0.3) +
  stat_density(aes(x=live,y=..scaled..),fill='blue3',adjust=0.8,alpha=0.3) +
  stat_density(aes(x=noise,y=..scaled..),fill='red3',adjust=0.8,alpha=0.3) +
  geom_text(data=counts[1,],aes(x=600,y=1,label="Total"),color='grey50',alpha=0.8,size=6,
            hjust=1,vjust=1) +
  geom_text(data=counts[1,],aes(x=600,y=0.95,label="Live"),color='blue3',alpha=0.8,size=6,
            hjust=1,vjust=1) +
  geom_text(data=counts[1,],aes(x=600,y=0.9,label="Saline-responsive"),color='red3',alpha=0.8,size=6,
            hjust=1,vjust=1) +
  theme(axis.line.x=element_line(lineend='round',size=1.2),
        axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=15,angle=35,hjust=1),
        axis.text.y=element_blank(),
        axis.title=element_text(size=15),
        axis.ticks.x=element_line(lineend='round',size=1.2))
dev.off()