# clear workspace
rm(list = ls())

# using library
library(coda);
library(plotrix)

# Setting Path
rfDistPath = '/Users/esaymuh/Desktop/Academia/DomainDLRS_data/Running/Results/Simulations2/';

combinePlot = T;
if(combinePlot) {

# Reading domainDLRS results:
dt <-read.table( paste(rfDistPath,"ZNF91/domainDLRS_MCMC_ZNF.map.rfd"  , "", sep=""),header=F,sep="\t");
data91 <- data.matrix(dt);

dt <-read.table( paste(rfDistPath,"ZNF468/domainDLRS_MCMC_ZNF.map.rfd"  , "", sep=""),header=F,sep="\t");
data468 <- data.matrix(dt);

dt <-read.table( paste(rfDistPath,"ZNF679/domainDLRS_MCMC_ZNF.map.rfd"  , "", sep=""),header=F,sep="\t");
data679 <- data.matrix(dt);


# Reading MrBayes Results:
dt <-read.table( paste(rfDistPath,"ZNF91/simZNF91_MrBayes_ZNF.map.rfd" , "", sep=""),header=F,sep="\t");
mbData91 <- data.matrix(dt);

dt <-read.table( paste(rfDistPath,"ZNF468/simZNF468_MrBayes_ZNF.map.rfd" , "", sep=""),header=F,sep="\t");
mbData468 <- data.matrix(dt);

dt <-read.table( paste(rfDistPath,"ZNF679/simZNF679_MrBayes_ZNF.map.rfd" , "", sep=""),header=F,sep="\t");
mbData679 <- data.matrix(dt);


# Prepare the data frame
method= c( rep("DomainDLRS", length(data91)),  rep("MrBayes",  length(mbData91)), 
           rep("DomainDLRS", length(data468)), rep("MrBayes",  length(mbData468)),
           rep("DomainDLRS", length(data679)), rep("MrBayes",  length(mbData679)) ) ;
         
family= c( rep("SZNF91", 2*length(data91)),rep("SZNF468", 2*length(data91)),rep("SZNF679", 2*length(data91)));           

rfdistances = c(data91,mbData91,data468,mbData468,data679,mbData679);

plotData <- data.frame( Method=method,DomainFamily=family,RF_Distance=rfdistances)


#position="identity"
#position="dodge"
# + theme(plot.title = element_text(lineheight=.8, face="bold"))

require(ggplot2)

pl = ggplot(plotData) + theme_bw() +
     geom_histogram(aes(x=RF_Distance, fill=Method), colour="grey50", alpha=0.6, position="identity") +
     theme(text = element_text(size=16), axis.text.x = element_text(angle=0, vjust=1)) +
     facet_grid(DomainFamily ~ .) + xlab("RF Distance") + ylab("Frequency")

#+ ggtitle("Robinson Floud distance of ZNF domain trees inferred\n by domainDLRS and MrBayes from ground-truth") 

pl

#Saving bar chart
outPath="/Users/esaymuh/Desktop/Academia/DomainDLRS_data/Running/Results/Simulations2/rfDistances/figures/";
pdf( paste0(outPath,"rfDistancePlot_SZNF_DomainDLRS_MrBayes.pdf", sep=""));
plot(pl)
dev.off()

}

######################################################################################################
# Please note this will produce individual plots for Zinc-Finger domain tree by domainDLRS and MrBayes 
# for ZNF91, ZNF468, ZNF679 ...
######################################################################################################
individualPlot = F;
if(individualPlot) {

# Setting Path
rfDistPath = "/Users/auwnm/Desktop/Running/Results/Simulations2/ZNF679/";

# Data Reading
dt1 <-read.table( paste(rfDistPath,"domainDLRS_MCMC_ZNF.map.rfd"  , "", sep=""),header=F,sep="\t");
data1 <- data.matrix(dt1);

dt2 <-read.table( paste(rfDistPath,"simZNF679_MrBayes_ZNF.map.rfd" , "", sep=""),header=F,sep="\t");
data2 <- data.matrix(dt2);

plotData <- data.frame( Method =c(rep("domainDLRS On Zinc-Finger domain", length(data1)), rep("MrBayes On Zinc-Finger domain", length(data2)) ), RF_Distance=c(data1,data2) )


pdf( paste (path,"RFDistances_Plots/individualPlot_ZNF679.pdf", sep=""));

require(ggplot2)
pl = ggplot(plotData,aes(RF_Distance,fill=Method),binwidth=1) + geom_bar() + facet_grid(Method ~ .) +  xlab("RF Distances") 


plot(pl)
dev.off();

}


