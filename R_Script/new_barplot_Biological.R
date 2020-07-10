# clear workspace
rm(list = ls())


library(plotrix)

# Setting Path
#path = '/Users/esaymuh/Desktop/Academia/DomainDLRS_data/Running/Results/BiologicalResults/';
#setwd(path);


families = c("ZNF91","ZNF558","ZNF468","ZNF611","ZNF679","ZNF764");
N = c(122,66,86,138,62,44);
families_sizes = 2 * N -1;

#Domain Duplications
domainDup_DomainDLRS = c(38,9,17,23,10,6);
domainDup_MrBayes    = c(70,29,46,68,31,18);

#Domain Bifurcation
domainBifur_DomainDLRS = c(83,56,68,114,51,37);
domainBifur_MrBayes    = c(51,36,39,69,30,25);

#Domain Losses
domainLosses_DomainDLRS = c(2,0,3,19,1,2);
domainLosses_MrBayes    = c(109,68,113,198,78,44);

#Pre-root Duplications
#domainPreDup_DomainDLRS = c(22,8,11,17,7,5,4);
#domainPreDup_MrBayes    = c(58,19,36,35,20,7,13); 
#events =   c(rep("Domain Duplications",14),rep("Domain Bifurcations",14), rep("Domain Losses",14));
#domainPreDup_DomainDLRS / families_sizes, domainPreDup_MrBayes / families_sizes 


events =   c(rep("Domain Duplications",12),rep("Domain Bifurcations",12), rep("Domain Losses",12));


method =   c(rep("DomainDLRS",6),rep("MrBayes-MPR",6) );


percentage = c( domainDup_DomainDLRS    / families_sizes, domainDup_MrBayes   / families_sizes,
                domainBifur_DomainDLRS  / families_sizes, domainBifur_MrBayes / families_sizes,
                domainLosses_DomainDLRS / families_sizes, domainLosses_MrBayes / families_sizes ) ;

             
                
df <- data.frame( Events = events, Families = c(rep(families,6)), Method=c(rep(method,3)), Percentage = percentage );


require(ggplot2)
#pl = qplot(Events,Precision,data=data1, geom="bar", fill=Events, ylim=c(0, 1),stat="identity" ) 

#pl = qplot( Events, data = data1, geom = "bar", fill = Events, weight=intensity , position = "dodge",ylim=c(0, 100),facets = . ~ Factor) + theme(text = element_text(size=12),
#        axis.text.x = element_text(angle=90, vjust=1))


pl = ggplot(df, aes(Families, Percentage, fill= Method)) + theme_bw() + 
    geom_bar(stat="identity", position="dodge") + ylim(0,1.0) + facet_grid(. ~ Events ) +
    theme(strip.text.x = element_text(size = 10),text = element_text(size=14), axis.text.x = element_text(angle=90, vjust=1)) +
    geom_text(aes(label = sprintf("%.02f",Percentage)), position = position_dodge(width = .8), vjust = -0.5, size=2) +
    scale_fill_manual(values=c('#F8766D','#00BFC4'))

pl

#coord_flip()
#facet_grid(Events ~ .) 
#facet_wrap(~ Events)

#Saving bar chart
outPath="/Users/esaymuh/Desktop/Academia/DomainDLRS_data/Running/Results/Simulations2/rfDistances/figures/";
pdf( paste0(outPath,"realData_DomainEventPlot.pdf", sep=""));
plot(pl)
dev.off()







DLPlot = F;
if(DLPlot){
s = c(0.49,0.48,0.35,0.63,0.61,0.54, 0.12,0.23,0.06,0.73,0.87,0.78);
familyVect = c("ZNF91","ZNF468","ZNF679");
data1 <- data.frame( Family=c(rep(familyVect,2)), Method=c(rep("domainDLRS",3),rep("MrBayes MPR",3)),Measure=c(rep("duplication Score",6),rep("Loss Score",6)), Score=s)


require(ggplot2)
pl = ggplot(data1, aes(Family, Score, fill= Method)) + geom_bar(stat="identity", position="dodge") + facet_grid(Measure ~ .) + theme(text = element_text(size=12), axis.text.x = element_text(angle=90, vjust=1)) + ylim(0,1) + scale_fill_manual(values=c("#CC6666", "#9999CC", "#66CC99"))

#Saving bar chart
pdf( paste (path,"dupLossPlot1.pdf", sep=""));
plot(pl)
dev.off()

}
