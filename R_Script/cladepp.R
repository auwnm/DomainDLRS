library("plotrix")

#counts
transpose<-function(x, prog){
  row.names(x)=x$DomainFamily
  x$DomainFamily = NULL
#  x=as.data.frame(apply(x,1,function(y) y=c(y[1],y[-1]-y[-10])))
  x=as.data.frame(apply(x,1,function(y) y=y)) 
  names(x)=paste(prog, names(x), sep=".")
  return(x)
}


#setwd("/Users/bengts/SciLifeLab/auwn/articles/DomainDLRS/R-files")


## pre.ddlrs=read.table("cladepp_pre_ddlrs.txt",header=T)
## pre.mb=read.table("cladepp_pre_mb.txt",header=T)
## post.ddlrs=read.table("cladepp_post_ddlrs.txt",header=T)
## post.mb=read.table("cladepp_post_mb.txt",header=T)
all = readLines("cladepp_all.txt")
pre.ddlrs=read.table(textConnection(all[2:8]), header=T)
print(pre.ddlrs)
post.ddlrs=read.table(textConnection(all[11:17]), header=T)
print(post.ddlrs)
pre.mb=read.table(textConnection(all[20:26]), header=T)
print(pre.mb)
post.mb=read.table(textConnection(all[29:35]), header=T)
print(pre.ddlrs)
all.ddlrs=data.frame(DomainFamily = pre.ddlrs$DomainFamily)
all.mb=list(DomainFamily = pre.mb$DomainFamily)
for(i in paste0("lt_",c('0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'))){
  all.ddlrs[[i]] = pre.ddlrs[[i]]+post.ddlrs[[i]]
  all.mb[[i]] = pre.mb[[i]]+post.mb[[i]]
}
all.ddlrs=data.frame(all.ddlrs)
all.mb=data.frame(all.mb)

res=cbind(transpose(pre.ddlrs, "pre.ddlrs"),transpose(pre.mb,"pre.mb"),transpose(post.ddlrs, "post.ddlrs"),transpose(post.mb,"post.mb"))
res$lims=as.numeric(sub("lt_","",row.names(res)))
print(res)

vars=c("ZNF91","ZNF558","ZNF468","ZNF611","ZNF679","ZNF764")
library(plotrix)
#pdf(file="cladepp_counts.pdf")
par(mfrow=c(2,3))
for(v in vars){
#  vals=as.vector(outer(paste(c("pre","post"),c("ddlrs","mb"), sep="."),v, FUN=paste,sep="."))
  vals=paste(as.vector(t(outer(c("pre","post"), c("ddlrs","mb"),paste,sep="."))),v, sep=".")
  print(vals)
  yran=range(res[,names(res) %in% vals])
  #yran=range(c(res[,paste("ddlrs",v)],res[,paste("mb",v)]))
  print(yran)
  yby=ceiling(yran[2]/100)*10
  print(yby)
  ymin=floor(min(res[10,names(res) %in% vals[3:4]])/10)*10
  print(ymin)
  #  ymin=floor(min(res[10,paste("ddlrs",v)],res[10,paste("mb",v)])/10)*10
  nymax=ceiling(max(res[-10,names(res) %in% vals],res[10,names(res) %in% vals[1:2]])/10)*10
  print(nymax)
  lytic=seq(0,nymax,by=yby)
  print(lytic)
  uytic=seq(ymin,ceiling(yran[2]/10)*10, by=yby)
  print(uytic)
  ytic=c(lytic,uytic)
  print(ytic)
  gap=c(nymax*1.1, ymin*0.9)
  print(gap)
  col = list("red","blue","orange","darkgreen")
  names(col)=vals
  gap.plot(res$lims,res[[vals[1]]], gap=gap, gap.axis="y",ylab="raw counts",xlab="post clade prob",main=v,type='p',pch=19,lwd=1,col=col[[vals[1]]], ytics=ytic, yticlab=format(ytic,digit=1),ylim=yran)
  str(res)
  if(v==vars[1]){
    legend("topleft", legend=c("ddlrs pp ancient","mb pp ancient", "ddlrs pp recent","mb pp recent"), col=unlist(col),pch=19, cex=0.7)
  }
  for(i in vals[2:4]){
    print(i)
    print(res[,i])
    print(res[[i]])
    gap.plot(res$lims,res[[i]], gap=gap, gap.axis="y",ylab=v,type='p',pch=19,lwd=1,col=col[[i]], ylim=yran,add=TRUE)
  }
#  gap.plot(res$lims,res[[paste("post ddlrs",v)]], gap=gap, gap.axis="y",ylab=v,type='o',lwd=1,col="green", ylim=yran,add=TRUE)
#  gap.plot(res$lims,res[[paste("pos mb",v)]], gap=gap, gap.axis="y",ylab=v,type='o',lwd=1,col="orange", ylim=yran,add=TRUE)
}
#dev.off()  


#freq


transpose<-function(x, prog){
  row.names(x)=x$DomainFamily
  x$DomainFamily = NULL
  x=as.data.frame(apply(x,1,function(y) y=y/y[length(y)])) #c(y[1],y[-1]-y[-10])))
  names(x)=paste(prog, names(x), sep=".")
  return(x)
}


#setwd("/Users/bengts/SciLifeLab/auwn/articles/DomainDLRS/R-files")



#pre.ddlrs=read.table("cladepp_pre_ddlrs.txt",header=T)
#pre.mb=read.table("cladepp_pre_mb.txt",header=T)
#post.ddlrs=read.table("cladepp_post_ddlrs.txt",header=T)
#post.mb=read.table("cladepp_post_mb.txt",header=T)
res=cbind(transpose(pre.ddlrs, "pre.ddlrs"),transpose(pre.mb,"pre.mb"),
          transpose(post.ddlrs, "post.ddlrs"),transpose(post.mb,"post.mb"),
          transpose(all.ddlrs, "all.ddlrs"),transpose(all.mb,"all.mb")
          )
res$lims=as.numeric(sub("lt_","",row.names(res)))
print(res)

vars=c("ZNF91","ZNF558","ZNF468","ZNF611","ZNF679","ZNF764")


library(plotrix)
#pdf(file="cladepp_freq.pdf",height=3.5 )
par(mfrow=c(1,2))
for(v in c("ddlrs","mb")){
  vals=as.vector(outer(paste(c("pre","post", "all"),v, sep="."),vars, FUN=paste,sep="."))
  print(vals)
  yran=range(res[,names(res) %in% vals])
  ytic=seq(0,1,by=0.1)#c(lytic,uytic)
  gap=c(nymax*1.1, ymin*0.9)
  col = as.list(rep(c("red","green", "black"),6))
  names(col)=vals
  plot(res$lims,res[[vals[1]]],ylim=yran,ylab="freq",xlab="post clade prob",main=v,type='l',pch=19,lwd=2,col=col[[vals[1]]]) #,ylim=yran)
  if(v=="ddlrs"){ #vars[1]){
    legend("topleft", legend=c("pp ancient","pp recent","pp_all"), col=unlist(col),pch=19, cex=0.7)
  }
  for(i in vals[2:18]){
    lines(res$lims,res[[i]],type='l',pch=19,lwd=2,col=col[[i]])
  }
}
#dev.off()  

le0_8 = list()
for(v in as.vector(outer(c("pre","post"),c("ddlrs","mb"), paste, sep="."))){ 
  vals=as.vector(outer(v,vars, FUN=paste,sep="."))
  le0_8[[v]] = mean(unlist(res[8,vals]))
  print(paste("le0_8[[",v,"]] = ",le0_8[[v]]))
}

res2=res
res2["lt_0.0", ] = rep(0,ncol(res2))
res2$name = row.names(res2)
row.names(res2) = as.character(c(1,2,3,4,5,6,7,8,9,10,0))
for(v in names(res2)[!names(res2)%in%c('name', 'lims')]){
  for(j in 10:1){
    i = as.character(j)
    k = as.character(j-1)
    print(paste(i,v))
    res2[i, v] = res2[i,v]-res2[k,v]
    print(res2[i,v])
  }
}
res2 = res2[row.names(res2) != '0',]
row.names(res2) =res2$name
res2$name=NULL



library(plotrix)
#pdf(file="cladepp_freq.pdf",height=3.5 )
par(mfrow=c(1,2))
for(v in c("ddlrs","mb")){
  vals=as.vector(outer(paste(c("pre","post", "all"),v, sep="."),vars, FUN=paste,sep="."))
  print(vals)
  yran=range(res2[,names(res2) %in% vals])
  ytic=seq(0,1,by=0.1)#c(lytic,uytic)
  gap=c(nymax*1.1, ymin*0.9)
  col = as.list(rep(c("red","blue", "green"),6))
  names(col)=vals
  plot(res2$lims,res2[[vals[1]]], 
       ylim=yran, ylab="freq", xlab="post clade prob", main=v, 
       type='l', pch=19, lwd=2, col=col[[vals[1]]]) #,ylim=yran)
  if(v=="ddlrs"){ #vars[1]){
    legend("topleft", legend=c("pp ancient","pp recent","pp_all"), col=unlist(col),pch=19, cex=0.7)
  }
  for(i in vals[2:18]){
    lines(res2$lims,res2[[i]],type='l',pch=19,lwd=2,col=col[[i]])
  }
}
#dev.off()  