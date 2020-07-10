# clear workspace
rm(list = ls())

library(plotrix)
library(ggplot2)
library(dplyr)
library(reshape2)
library(plyr)

# Setting Path
path = '/Users/esaymuh/Desktop/Academia/DomainDLRS_data/Running/Results/';
setwd(path);

domainDLRS_path = "/Users/esaymuh/Desktop/Academia/DomainDLRS_data/Running/Results/DomainDLRS/new_results_2020/"
MrBayes_path = "/Users/esaymuh/Desktop/Academia/DomainDLRS_data/Running/Results/MrBayes/new_results_2020/"

fam <- c("ZNF91", "ZNF468","ZNF679","ZNF558","ZNF611","ZNF764")

for (i in 1:length(fam)) {
  
  i <- 1;
  
  dd_post_file <- paste0(domainDLRS_path,"post_domainDLRS_",fam[i]);
  mb_post_file <- paste0(MrBayes_path,"post_MrBayes_",fam[i]);
  
  dd_post_df<- read.csv(file= dd_post_file, sep= '\t', header = TRUE)[,1]
  mb_post_df<- read.csv(file= mb_post_file, sep= '\t', header = TRUE)[,1]
  
  dd_post_df <- cbind( dd_post_df, Method=rep("DomainDLRS",length(dd_post_df)))
  mb_post_df <- cbind( mb_post_df, Method=rep("MrBayes",length(mb_post_df)))
  
  df_post <- rbind(dd_post_df,mb_post_df);
  colnames(df_post) <- c("Posterior", "Method");
  
  #ggplot(df_post, aes(x=Posterior, fill=Method, color=Method)) +
  #  geom_histogram(position="identity")
  
  
  ggplot(data=df_post, aes(log10(df_post ))) + geom_histogram() +
    xlab('Product of clade posterior probability for each tree (log)') + ylab('Frequency') 
  
}
  