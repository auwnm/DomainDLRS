# clear workspace
rm(list = ls())

# using library
library(coda);
library(plotrix)


# Setting Path
rfDistPath = '/Users/esaymuh/Desktop/Academia/DomainDLRS_data/Running/Results/Simulations2/';

#MrBayes data Reading
dt <-read.table( paste(rfDistPath,"ZNF91/simZNF91_MrBayes_Tandem.map.rfd" , "", sep=""),header=F,sep="\t");
znf_91_tandem <- data.matrix(dt);

dt <-read.table( paste(rfDistPath,"ZNF91/simZNF91_MrBayes_Random.map.rfd" , "", sep=""),header=F,sep="\t");
znf_91_random <- data.matrix(dt);


dt <-read.table( paste(rfDistPath,"ZNF679/simZNF679_MrBayes_Tandem.map.rfd" , "", sep=""),header=F,sep="\t");
znf_679_tandem <- data.matrix(dt);

dt <-read.table( paste(rfDistPath,"ZNF679/simZNF679_MrBayes_Random.map.rfd" , "", sep=""),header=F,sep="\t");
znf_679_random <- data.matrix(dt);


dt <-read.table( paste(rfDistPath,"ZNF468/simZNF468_MrBayes_Tandem.map.rfd" , "", sep=""),header=F,sep="\t");
znf_468_tandem <- data.matrix(dt);

dt <-read.table( paste(rfDistPath,"ZNF468/simZNF468_MrBayes_Random.map.rfd" , "", sep=""),header=F,sep="\t");
znf_468_random <- data.matrix(dt);


#DomainDLRS data reading
dt <-read.table( paste(rfDistPath,"ZNF91/domainDLRS_MCMC_Gene.map.rfd" , "", sep=""),header=F,sep="\t");
znf_91_DomainDLRS <- data.matrix(dt);

dt <-read.table( paste(rfDistPath,"ZNF679/domainDLRS_MCMC_Gene.map.rfd" , "", sep=""),header=F,sep="\t");
znf_679_DomainDLRS <- data.matrix(dt);

dt <-read.table( paste(rfDistPath,"ZNF468/domainDLRS_MCMC_Gene.map.rfd" , "", sep=""),header=F,sep="\t");
znf_468_DomainDLRS <- data.matrix(dt);


# Preparing data frame
tools <- c( rep("DomainDLRS",100), rep("MrBayes (Tandem Model)", 100), rep("MrBayes (Random Model)", 100) );
families <- c(rep("SZNF91",300),rep("SZNF679",300),rep("SZNF468",300))

znf91_data  <- c(znf_91_DomainDLRS,  znf_91_tandem, znf_91_random);
znf679_data <- c(znf_679_DomainDLRS, znf_679_tandem, znf_679_random);
znf468_data <- c(znf_468_DomainDLRS, znf_468_tandem, znf_468_random);

plotData <- data.frame( Tools=c(rep(tools,3) ), Families = families, RF_Distance=c(znf91_data, znf679_data, znf468_data) );

pl = ggplot(plotData, aes(factor(RF_Distance),fill=factor(Tools, levels=c("DomainDLRS","MrBayes (Tandem Model)","MrBayes (Random Model)"),ordered=TRUE) ),stat="identity") + scale_x_discrete(drop = FALSE) +
       geom_bar(position = position_dodge(width = .8), width = 0.7)  +
       facet_grid(Families ~ .)  + scale_fill_discrete(name="Method") + xlab("RF Distance") + ylab("Frequency") 


# New plotting starts from here
plotData$RF_Distance = as.factor(plotData$RF_Distance)
plotData$Families = as.factor(plotData$Families)
plotData$Tools = factor(plotData$Tools, levels=c("DomainDLRS","MrBayes (Tandem Model)","MrBayes (Random Model)"),ordered=TRUE) 

df1 <- ddply(plotData, .(Families,Tools), summarise, 
            types = as.factor(names(table(RF_Distance))), 
            counts = as.numeric(table(RF_Distance)))

pl = ggplot(df1, aes(x = types, y = counts, fill = Tools)) + theme_bw() +
  geom_bar(stat='identity',colour="black", position = position_dodge(width = .8), width = 0.7) +  
  scale_fill_manual(values=c('#F8766D','#00BFC4','#619CFF')) + 
  theme(text = element_text(size=16), axis.text.x = element_text(angle=0, vjust=1)) +
  geom_text(aes(label = counts), position = position_dodge(width = .8), vjust = -0.5, size=2) +
  facet_grid(Families ~ .) + xlab("RF Distance") + ylab("Percentage")
  
# scale_fill_brewer(palette = "Set2")
# facet_wrap(~ Families) 

#Saving bar chart
outPath="/Users/esaymuh/Desktop/Academia/DomainDLRS_data/Running/Results/Simulations2/rfDistances/figures/";
pdf( paste0(outPath,"rfDistancePlot_SZNF_Random_Tandem_Model.pdf", sep=""));
plot(pl)
dev.off()




#dt1 <-read.table( paste("/Users/auwnm/Documents/Jworkspace/JPrIME_runs/GenPhyloData/Simulator/sim3/" , "log.sum.state.out", sep=""),header=F,sep="\t");
#data1 <- data.matrix(dt1);

#dt3 <-read.table( paste(rfDistPath,"simZNF679_MrBayes_Tandem.map.rfd" , "", sep=""),header=F,sep="\t");
#data3 <- data.matrix(dt3);

#dt3 <-read.table( paste(rfDistPath,"simZNF679_MrBayes_Tandem.map.rfd" , "", sep=""),header=F,sep="\t");
#data3 <- data.matrix(dt3);

#dt4 <-read.table( paste("simTest2_PhyML.rfd" , "", sep=""),header=F,sep="\t");
#data4 <- data.matrix(dt4);

# Data Selection
#dist <- data2[,1];






#plotData <- data.frame( Method =c(rep("Hill Climbing on state", length(data1)), rep("Sampling on state", length(data2)),rep("Hill Climbing", length(data3)) ), RF_Distance=c(data1,data2,data3) )
#plotData <- data.frame( Method =c(rep("domainDLRS", length(data1)), rep("MrBayes (Random Model)", length(data2)),rep("MrBayes (Tandem Model)", length(data3)) ), RF_Distance=c(data1,data2,data3) )

#pl = qplot(RF_Distance,data=plotData, geom = "histogram" ,binwidth = 1,fill = software );
#pl = qplot(RF_Distance, data = plotData, geom = "bar",binwidth = 2,stat="bin", xlim = c(0, 10)) + facet_grid(tool ~ .)
#pl = ggplot(plotData, aes(RF_Distance, fill=tool),stat="identity",geom = "bar") + geom_bar(position="dodge")
#pl = ggplot(plotData, aes(RF_Distance, fill=tool,width=5),stat="identity",geom="bar",binwidth = 2) + geom_bar(position="dodge")
#pl = qplot(factor(RF_Distance), data= plotData, geom="bar",stat="bin", fill=Tool) + facet_grid(Tool ~ .) + ylab("Frequency") +   xlab("RF Distances") 
#pl = qplot(factor(RF_Distance), data= plotData, geom="bar", binwidth = .01,stat="bin",fill=Method) + facet_grid(Method ~ .) + ylab("Frequency") +   xlab("RF Distances") 

#pl = ggplot(plotData, aes(RF_Distance, fill= Method)) + geom_bar(position="dodge")
#pl = ggplot(plotData,aes(RF_Distance),binwidth=1) + geom_bar() + facet_wrap(~ Method) 
#pl = qplot(factor(RF_Distance), data= plotData, geom="bar",stat="bin", binwidth=0.2,fill=Method) + facet_grid(Method ~ .) + ylab("Frequency") +   xlab("RF Distances") 
#pl = qplot(x=tool, y=RF_Distance, fill= RF_Distance,binwidth = 1, data=plotData, geom="bar", stat="bin", position="dodge") 




# Histogram analysis
if(FALSE) {
h1 = hist(dist,breaks=seq(min(dist), max(dist),by= 1));
pdf("rfdistGene.pdf");
plot(h1,col='gray'); #freq = F,xlim = c(0,30)
dev.off();
}


# Pi chart analysis
X = dist;
pieChart = FALSE;
if(pieChart) {
	grp1=X[(X==0)];
	grp2=X[(X==2)];
	grp3=X[(X==4)];
	grp4=X[(X==6)];
	grp5=X[(X > 6)];
	
	#Step 3: calculating %ages for slice
	slice<-c(length(grp1)/length(X),length(grp2)/length(X),length(grp3)/length(X),length(grp4)/length(X),length(grp5)/length(X))*100;
	
	#Step 4: Making Labels
	leg <- c("d=0","d=2","d=4","d=6","d>=8");
	colors <- c("blue","lightblue","violet","Yellow","Red");
	pct <- round(slice);
	lbls <- paste(pct,"%",sep=""); 
	
	#Step 5: Preparing the output File
	pdf("rfdist.pdf");
	
	#Step 5: Plotting
	#pie(slice,labels= lbls,col=rainbow(4),main="Pie chart of Graph Degree");
	pie3D(slice,labels= lbls,col=colors,explode=0.1,shade=0.8,main="RF distances between original tree topology and MCMC chain\n for gene tree");
	legend("top",leg,fill=colors,cex=0.87,bty = "n",horiz = TRUE);
	
	#Step 6: Closing
	dev.off();
}


#hist(dist,range=0)
# which(birthRate == min(birthRate))


# Simple Pie Chart
#slices <- h1
#pie(dist);

#product <- log(-1.0*data2[,6]);
#product <- data2[,5];




#sum <- log2(-1*data1[,1]);
#productMin <- min(product);
#productMax <- max(product);



#diff = (productMax - productMin) / 50.0;
#h = hist( product, breaks=seq(productMin, productMax,by= diff) );
#plot(h);

#max(product);
#plot(sum,ylim = c(sumMin, sumMax),type = "l");



# Plotting
#chain1 = mcmc(sum);
#minX = min(sum);
#maxX = max(sum);

#pdf(paste("/Users/auwnm/Documents/Jworkspace/JPrIME_runs/GenPhyloData/Simulator/sim3/sumConverge",".pdf",sep=""));

#summary(chain1)
#plot(chain1)
#plot(sum, type = "l");

#dev.off();

