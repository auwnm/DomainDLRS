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
 
  fam_name <- paste0('S',fam[i]);
  
  dd_sens_file <- paste0(path,fam[i],"/accuracy_sensitivity_sim",fam[i],"_domainDLRS.csv");
  dd_spec_file <- paste0(path,fam[i],"/accuracy_specificity_sim",fam[i],"_domainDLRS.csv");
  
  dd_sens_df <- cbind(read.csv(file= dd_sens_file, sep= '\t', header = TRUE), Method=rep("DomainDLRS",100), Family=rep(fam_name,100), Measure=rep("Sensitivity",100))
  dd_sens_df <- melt(data=dd_sens_df, id.vars = c("Method","Family", "Measure"), measure.vars = c("Ancient", "Recent", "All" )) ;
  
  dd_spec_df <- cbind(read.csv(file= dd_spec_file, sep= '\t', header = TRUE), Method=rep("DomainDLRS",100), Family=rep(fam_name,100), Measure=rep("Specificity",100))
  dd_spec_df <- melt(data=dd_spec_df, id.vars = c("Method","Family", "Measure"), measure.vars = c("Ancient", "Recent", "All" )) ;
  
  df_dd <- rbind(dd_sens_df,dd_spec_df);
  
  mb_sens_file <- paste0(path,fam[i],"/accuracy_sensitivity_sim",fam[i],"_MrBayes.csv");
  mb_spec_file <- paste0(path,fam[i],"/accuracy_specificity_sim",fam[i],"_MrBayes.csv");
  
  mb_sens_df <- cbind(read.csv(file= mb_sens_file, sep= '\t', header = TRUE), Method=rep("MrBayes-MPR",100), Family=rep(fam_name,100),Measure=rep("Sensitivity",100))
  mb_sens_df <- melt(data=mb_sens_df, id.vars = c("Method","Family", "Measure"), measure.vars = c("Ancient", "Recent", "All" )) ;
  
  mb_spec_df <- cbind(read.csv(file= mb_spec_file, sep= '\t', header = TRUE), Method=rep("MrBayes-MPR",100), Family=rep(fam_name,100),Measure=rep("Specificity",100))
  mb_spec_df <- melt(data=mb_spec_df, id.vars = c("Method","Family", "Measure"), measure.vars = c("Ancient", "Recent", "All" )) ;
  
  df_mb <- rbind(mb_sens_df,mb_spec_df);
  
  data[[i]]<- rbind(df_dd,df_mb);

}

sim_znf<-data.frame(do.call(rbind, data))
colnames(sim_znf) <- c("Method", "Family", "Measure","Variable","Percentage");


# Plotting individual box plots for each event
events = c("Ancient", "Recent", "All" );
for (i in 1:length(events)) {
  
  
  plot.df<-subset(sim_znf, Variable==events[i])
  pl_event = ggplot(plot.df, aes(x=Family, y=Percentage, color=Method, fill=Method)) + theme_bw() +
    theme(text = element_text(size=16), axis.text.x = element_text(angle=90, vjust=1)) + 
    xlab('Family') + ylab('Percentage') +
    geom_boxplot(fill="white",outlier.colour = NA, position = position_dodge(width=0.9))+
    geom_point(position=position_jitterdodge(dodge.width=0.9),size = 0.7) +
    facet_grid(. ~ Measure) 
  
  
  #Saving bar chart
  outPath="/Users/esaymuh/Desktop/Academia/DomainDLRS_data/Running/Results/Simulations2/rfDistances/figures/";
  pdf( paste0(outPath,"accuracy_boxPlot_SZNF_sensitivity_specificity_event=",events[i],".pdf", sep=""));
  plot(pl_event)
  dev.off()
}


plot.df <- ddply(sim_znf, c("Method", "Family", "Measure","Variable"), summarize, Mean = mean(Percentage), SD = sd(Percentage))
pl_bar = ggplot(plot.df, aes(Family, Mean, fill= Method)) + 
  geom_bar(stat="identity",  position=position_dodge(width=0.9)) + facet_grid(Variable ~ Measure ) +
  theme(text = element_text(size=14), axis.text.x = element_text(angle=90, vjust=1)) +
  ylim(0.0, 1.0) + geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, inherit.aes = TRUE,  position=position_dodge(width=0.9))

pl_bar

#Saving bar chart
outPath="/Users/esaymuh/Desktop/Academia/DomainDLRS_data/Running/Results/Simulations2/rfDistances/figures/";
pdf( paste0(outPath,"accuracy_barPlot_SZNF_sensitivity_specificity_grid.pdf", sep=""));
plot(pl_bar)
dev.off()

pl_box = ggplot(sim_znf, aes(x=Family, y=Percentage, color=Method, fill=Method)) + theme_bw() +
  theme(text = element_text(size=16), axis.text.x = element_text(angle=90, vjust=1)) + 
  xlab('Family') + ylab('Percentage') +
  geom_boxplot(fill="white",outlier.colour = NA, position = position_dodge(width=0.9))+
  geom_point(position=position_jitterdodge(dodge.width=0.9),size = 0.7) +
  facet_grid(Variable ~ Measure) 
pl_box


#Saving bar chart
outPath="/Users/esaymuh/Desktop/Academia/DomainDLRS_data/Running/Results/Simulations2/rfDistances/figures/";
pdf( paste0(outPath,"accuracy_boxPlot_SZNF_sensitivity_specificity_grid.pdf", sep=""));
plot(pl_box)
dev.off()


