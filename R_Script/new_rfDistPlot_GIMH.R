# clear workspace
rm(list = ls())

# using library
library(coda);
library(plotrix)

# Setting Path
rfDistPath = '/Users/esaymuh/Desktop/Academia/DomainDLRS_data/Running/Results/Simulations2/';


# Reading domainDLRS results:
dt <-read.table( paste(rfDistPath,"ZNF91/domainDLRS_HC_ZNF.map.rfd" , "", sep=""),header=F,sep="\t");
data_HC <- data.matrix(dt);

dt <-read.table( paste(rfDistPath,"ZNF91/domainDLRS_MCMC_Sampling_ZNF.map.rfd" , "", sep=""),header=F,sep="\t");
data_GIMH <- data.matrix(dt);

dt <-read.table( paste(rfDistPath,"ZNF91/domainDLRS_MCMC_ZNF.map.rfd" , "", sep=""),header=F,sep="\t");
data_hGIMH <- data.matrix(dt);


# Prepare the data frame
technique= c(   rep("Hill Climbing", length(data_HC)),  
                rep("GIMH"         , length(data_GIMH)), 
                rep("hGIMH"        , length(data_hGIMH)) ) ;
       
rfdistances = c(data_HC, data_GIMH, data_hGIMH);

plotData <- data.frame(Technique=technique,RF_Distance=rfdistances)


require(ggplot2)

pl = ggplot(plotData) + theme_bw() +
     geom_histogram(aes(x=RF_Distance), fill="#FF9999", colour = "black", position="identity", bins = 30) + ylim(0,15) + 
     theme(text = element_text(size=16), axis.text.x = element_text(angle=0, vjust=1) ,strip.text.y = element_text(size = 18, angle = 270) ) +
     facet_grid(Technique ~ .) + xlab("RF Distance") + ylab("Frequency") + theme(legend.position='none') 
     

#Saving bar chart
outPath="/Users/esaymuh/Desktop/Academia/DomainDLRS_data/Running/Results/Simulations2/rfDistances/figures/";
pdf( paste0(outPath,"rfDistancePlot_SZNF91_GIMH.pdf", sep=""));
plot(pl)
dev.off()
