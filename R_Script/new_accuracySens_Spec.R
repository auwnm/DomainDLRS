# clear workspace
rm(list = ls())

library(plotrix)
library(ggplot2)
library(dplyr)
library(reshape2)
library(plyr)

# Setting Path
path = '/Users/esaymuh/Desktop/Academia/DomainDLRS_data/Running/Results/Simulations2/';
setwd(path);

fam <- c("ZNF91", "ZNF468","ZNF679")

data <- list()
N <- 3

for (i in 1:N) {
 
   print(fam[i])

  dd_file <- paste0(path,fam[i],"/sensitivity_specificity_sim",fam[i],"_domainDLRS_ZNF.csv");
  mb_file <- paste0(path,fam[i],"/sensitivity_specificity_sim",fam[i],"_MrBayes_ZNF.csv");
  
  fam_name <- paste0('S',fam[i]);
  df_dd <- cbind( read.csv(file= dd_file, sep= '\t', header = FALSE), rep("DomainDLRS",100), rep(fam_name,100));
  df_mb <- cbind( read.csv(file= mb_file, sep= '\t', header = FALSE), rep("MrBayes-MPR",100), rep(fam_name,100));

  colnames(df_dd) <- c("Sensitivity", "Specificity", "method","family");
  colnames(df_mb) <- c("Sensitivity", "Specificity", "method", "family");

  data[[i]]<- rbind(df_dd,df_mb);

}

sim_znf<-data.frame(do.call(rbind, data))
melt_df<-melt(data=sim_znf, id.vars = c("method","family"), measure.vars = c("Sensitivity", "Specificity"))

colnames(melt_df) <- c("Method", "Family", "Measure","Percentage");
plot.df <- ddply(melt_df, c("Method", "Family", "Measure"), summarize, Mean = mean(Percentage), SD = sd(Percentage))

pl = ggplot(plot.df, aes(Family, Mean, fill= Method)) + 
  geom_bar(stat="identity",  position=position_dodge(width=0.9)) + facet_grid(. ~ Measure ) +
  theme(text = element_text(size=14), axis.text.x = element_text(angle=90, vjust=1)) +
  ylim(0.0, 1.0) + geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, inherit.aes = TRUE,  position=position_dodge(width=0.9))

pl = ggplot(melt_df, aes(x=Family, y=Percentage, color=Method, fill=Method)) + theme_bw() +
  theme(text = element_text(size=16), axis.text.x = element_text(angle=90, vjust=1)) + 
  xlab('Family') + ylab('Percentage') +
  geom_boxplot(fill="white",outlier.colour = NA, position = position_dodge(width=0.9))+
  geom_point(position=position_jitterdodge(dodge.width=0.9),size = 0.7) +
  facet_grid(. ~ Measure) 
pl


#Saving bar chart
outPath="/Users/esaymuh/Desktop/Academia/DomainDLRS_data/Running/Results/Simulations2/rfDistances/figures/";
pdf( paste0(outPath,"accuracyPlot_SZNF_Domain_sensitivity_specificity.pdf", sep=""));
plot(pl)
dev.off()


