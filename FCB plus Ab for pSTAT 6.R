setwd("C:/Users/artat/Documents/Flow data")
rm(list=ls())
library(flowCore)
library(flowClust)
library(flowViz)
library(flowWorkspace)
library(ggcyto)
library(flowType)
library(ggplot2)

filenames <- c("10-02-18_pSTAT3 nonstimulated_025")

###read the data
BCD_stain <- read.FCS(paste(filenames[1], ".fcs", sep=""), transformation=FALSE, alter.names = T)
ggcyto(BCD_stain, aes(x="FSC.A", y="SSC.A")) + geom_hex(bins=128) + xlim(c(0,2.5e5)) + ylim(c(0,2e4))
####Let us remove debris
nodebris <- rectangleGate(filterId="no Debris", "FSC.A"=c(20000,1.2e5), "SSC.A"=c(0,Inf))
BCD_stain_nodebris <- Subset(BCD_stain, filter(BCD_stain,nodebris))
ggcyto(BCD_stain_nodebris, aes(x="FSC.A", y="SSC.A")) + geom_hex(bins=128) + xlim(c(2e4,1.2e5)) + ylim(c(0,1e4))
Lympho_mono <- flowClust(BCD_stain_nodebris, varNames = c("FSC.A", "SSC.A"), K=2, B=100)
plot(Lympho_mono, data=BCD_stain_nodebris, level=0.8, z.cutoff=0, xlim=c(2e4,1.2e5), ylim=c(0,1e4), 
     pch=c(16,16), col = c("red", "blue"), 
     cex=1, cex.axis=2, xlab="FSC-A", ylab="SSC-A", cex.lab=2, las=1)
legend("topright", pch = c(16,16), 
       col = c("red", "blue"), legend = Lympho_mono@prior$order)
###Extract clusters as lympho and mono objects
FCB_pop1 <- split(BCD_stain_nodebris, Lympho_mono, population=list(p1=1,p2=2))
G1 <- FCB_pop1$p1
G2 <- FCB_pop1$p2
g1 <- summary(G1)
g2 <- summary(G2)
g1a <- g1[6,4]
g2a <- g2[6,4]
Max_SSC <- c(g1a, g2a)
get_order1 <- Max_SSC[order(Max_SSC, decreasing = TRUE)]
Maximum_SSC <- head(get_order1)[1]
Minimum_SSC <- head(get_order1)[2]
if (g1a == Maximum_SSC) {
  (Mono <- G1)
} else {
  (Mono <- G2)
}
if (g1a == Minimum_SSC) {
  (Lympho <- G1)
} else {
  (Lympho <- G2)
}
###Let us transform the data
biexp  <- biexponentialTransform("myTransform",w=0)
after.2 <- transform(Lympho, transformList(c("FSC.A", "FSC.H", "FSC.W", "SSC.A", "SSC.H", "SSC.W", 
                                                 "PerCP.Cy5.5.A","PE.A", "PE.Cy7.A","FITC.A","APC.H7.A", "BUV395.A",  
                                                 "Pacific.Orange.A", "Time"), biexp))
tf <- transformList(from=colnames(Lympho), tfun=biexp)
BCD_stain_trans<-tf %on% Lympho
after.3 <- transform(Mono, transformList(c("FSC.A", "FSC.H", "FSC.W", "SSC.A", "SSC.H", "SSC.W",
                                             "PerCP.Cy5.5.A","PE.A", "PE.Cy7.A","FITC.A","APC.H7.A", "BUV395.A",  
                                             "Pacific.Orange.A", "Time"), biexp))
tf2 <- transformList(from=colnames(Mono), tfun=biexp)
BCD_mono_trans <-tf2 %on% Mono

######Let us cluster for FCB dyes
BCD_gate <- rectangleGate(filterId = "FluoRegion1","Pacific.Orange.A"=c(5,12), "BUV395.A"=c(0,12))
BCD_filter = Subset(BCD_stain_trans, filter(BCD_stain_trans,BCD_gate))
ggcyto(BCD_filter, aes(x="Pacific.Orange.A", y="BUV395.A")) + geom_hex(bins=128)
FCB_clust <- flowClust(BCD_filter, varNames = c("Pacific.Orange.A","BUV395.A"), K=9, B=100)

png(width=2000, height=2000, res=200, file=paste(filenames[1], "_clust.png", sep=""))
plot(FCB_clust, data=BCD_filter, level=0.8, z.cutoff=0, xlim=c(5,12), ylim=c(0,10), 
     pch=c(16,16,16,16,16,16,16,16,16), col = c("red", "green", "blue", "cyan", "darkorange1", "gold", "hotpink", "deeppink2","deepskyblue2"), 
     cex=1, cex.axis=2, xlab="Pacific Orange", ylab="DyLight 350", cex.lab=2, las=1)
legend("topright", pch = c(16,16,16,16,16,16,16,16,16), 
       col = c("red", "green", "blue", "cyan", "darkorange1", "gold", "hotpink", "deeppink2","deepskyblue2"), legend = FCB_clust@prior$order)
dev.off()

###Extract clusters as separate objects
FCB_pop <- split(BCD_filter, FCB_clust, population=list(p1=1,p2=2,p3=3,p4=4,p5=5,p6=6,p7=7,p8=8,p9=9))
pop1 <- FCB_pop$p1
pop2 <- FCB_pop$p2
pop3 <- FCB_pop$p3
pop4 <- FCB_pop$p4
pop5 <- FCB_pop$p5
pop6 <- FCB_pop$p6
pop7 <- FCB_pop$p7
pop8 <- FCB_pop$p8
pop9 <- FCB_pop$p9
FCB_fs <- c(pop1, pop2, pop3, pop4, pop5, pop6, pop7, pop8, pop9) ###Combine object

results_CD3CD20 = list()
results_CD4CD8 = list()
plotsCD3CD20 = list()
plotsCD4CD8 = list()
lst=list()
lst2=list()

for (ii in 1:9) {
  FCB_clust2 <- flowClust(FCB_fs[[ii]], varNames = c("APC.H7.A","PerCP.Cy5.5.A"), K=3, B=200)
  plot1 = plot(FCB_clust2, data=FCB_fs[[ii]], level=0.8, z.cutoff=0, xlim=c(4,8), ylim=c(2,9), 
       pch=c(16,16,16), col = c("red", "green", "blue"), 
       cex=1, cex.axis=2, ylab="CD3", xlab="CD20", cex.lab=2, las=1)
  legend("bottomleft", pch = c(16,16,16), title="Cluster    Percent    ",
         col = c("red", "green", "blue"), legend = c(FCB_clust2@prior$order, (FCB_clust2@w*100)), ncol = 2)
  plot1
  plotsCD3CD20[[ii]] = plot1
  dev.copy(jpeg,filename=paste("CD3CD20_plot_", ii, ".jpeg", sep=""))
  dev.off ()
  FCB_pop2 <- split(FCB_fs[[ii]], FCB_clust2, population=list(d1=1,d2=2,d3=3))
  gate1 <- FCB_pop2$d1
  gate2 <- FCB_pop2$d2
  gate3 <- FCB_pop2$d3
  a <- summary(gate1)
  b <- summary(gate2)
  c <- summary(gate3)
  a1 <- a[6,8]
  b1 <- b[6,8]
  c1 <- c[6,8]
  Max_values <- c(a1, b1, c1)
  get_order <- Max_values[order(Max_values, decreasing = TRUE)]
  Maximum <- head(get_order)[1]
  Maximum
  if (a1==Maximum) {
    (CD3_cells <- gate1)
  } else if (b1 == Maximum) {
    (CD3_cells <- gate2)
  } else {
    (CD3_cells <- gate3)
  }
  FCB_clust3 <- flowClust(CD3_cells, varNames = c("FITC.A","PE.Cy7.A"), K=3, B=200)
  d1 = plot(FCB_clust3, data=CD3_cells, level=0.8, z.cutoff=0, xlim=c(4,12), ylim=c(3,9), 
       pch=c(16,16,16), col = c("red", "green", "blue"), 
       cex=1, cex.axis=2, xlab="CD8", ylab="CD4", cex.lab=2, las=1)
  legend("bottomleft", pch = c(16,16,16), 
         col = c("red", "green", "blue"), legend = c(FCB_clust3@prior$order, (FCB_clust3@w*100)), ncol=2)
  d1
  plotsCD4CD8[[ii]] = d1
  dev.copy(jpeg,filename=paste("CD4CD8_plot_", ii, ".jpeg", sep=""))
  dev.off ()
  
  percent_CD3CD20 <- (FCB_clust2@w *100)
  percent_CD4CD8 <- (FCB_clust3@w *100)
  results_CD3CD20[[ii]] = percent_CD3CD20
  results_CD4CD8[[ii]] = percent_CD4CD8
 
  results_CD3CD20[[ii]] <- data.frame("Population"= ii, "Proportions"= c(FCB_clust2@w *100), check.rows = T, check.names = T)
  results_CD4CD8[[ii]] <- data.frame("Population"= ii, "Proportions"= c(FCB_clust3@w *100), check.rows = T, check.names = T)
  
}
lst <- do.call(cbind, results_CD3CD20) ###combine data in one table

lst2 <- do.call(cbind, results_CD4CD8) ###combine data in one table


######Let us cluster for FCB dyes on Monocytes
BCD_gate <- rectangleGate(filterId = "FluoRegion1","Pacific.Orange.A"=c(5,14), "BUV395.A"=c(0,12))
BCD2_filter = Subset(BCD_mono_trans, filter(BCD_mono_trans,BCD_gate))
ggcyto(BCD2_filter, aes(x="Pacific.Orange.A", y="BUV395.A")) + geom_hex(bins=128)
FCB_clust4 <- flowClust(BCD2_filter, varNames = c("Pacific.Orange.A","BUV395.A"), K=9, B=100)

png(width=2000, height=2000, res=200, file=paste(filenames[1], "_clust_mono.png", sep=""))
plot(FCB_clust4, data=BCD2_filter, level=0.8, z.cutoff=0, xlim=c(5,14), ylim=c(0,12), 
     pch=c(16,16,16,16,16,16,16,16,16), col = c("red", "green", "blue", "cyan", "darkorange1", "gold", "hotpink", "deeppink2","deepskyblue2"), 
     cex=1, cex.axis=2, xlab="Pacific Orange", ylab="DyLight 350", cex.lab=2, las=1)
legend("topright", pch = c(16,16,16,16,16,16,16,16,16), 
       col = c("red", "green", "blue", "cyan", "darkorange1", "gold", "hotpink", "deeppink2","deepskyblue2"), legend = FCB_clust4@prior$order)
dev.off()
FCB_pop4 <- split(BCD2_filter, FCB_clust4, population=list(m1=1,m2=2,m3=3,m4=4,m5=5,m6=6,m7=7,m8=8,m9=9))
M1 <- FCB_pop4$m1
M2 <- FCB_pop4$m2
M3 <- FCB_pop4$m3
M4 <- FCB_pop4$m4
M5 <- FCB_pop4$m5
M6 <- FCB_pop4$m6
M7 <- FCB_pop4$m7
M8 <- FCB_pop4$m8
M9 <- FCB_pop4$m9
FCB_mono_fs <- c(M1, M2, M3, M4, M5, M6, M7, M8, M9) ###Combine object
results_CD14= list()
plotsCD14 = list()
lst3=list()
for (kk in 1:9) {
  data_M <- as.matrix(summary(FCB_mono_fs[[kk]]))
  data_M
  CD14_gate <- rectangleGate("SSC.A"=c(data_M["1st Qu.", "SSC.A"], data_M["Max.", "SSC.A"]), 
                           "PE.A"=c(data_M["Min.", "PE.A"], data_M["Max.", "PE.A"]))
  BCD.CD14.filter <- filter(FCB_mono_fs[[kk]], CD14_gate)
  percent_CD14 <- summary(BCD.CD14.filter)
  percent_CD14_b <- (percent_CD14@p *100)
  results_CD14[[kk]] = percent_CD14_b
}

lst3 <- do.call(cbind, results_CD14) ###combine data in one table

library(openxlsx)
write.xlsx(lst,"C:/Users/artat/Documents/Flow data/Results_percentages_CD3CD20.xlsx", row.names = TRUE, col.names = TRUE)
write.xlsx(lst2,"C:/Users/artat/Documents/Flow data/Results_percentages_CD4CD8.xlsx", row.names = TRUE, col.names = TRUE)
write.xlsx(lst3,"C:/Users/artat/Documents/Flow data/Results_percentages_CD14.xlsx", row.names = TRUE, col.names = TRUE)


