library("pegas")
library("ape")

cluster_file<-"Community/Haplotype_cluster/Af_Am_level2_level3_cluster_5cM_nonredundant_v2_ethnicity.txt"
cluster<-read.table(cluster_file,sep="\t",head=F,colClasses=c("character","character","character",rep("numeric",16),"character"),skip=1)
old_chrom=""
Result=list()
index=0
#for(i in 1:length(cluster[,1])){
for(i in 15:16){
  chrom=cluster[i,1]
  print(c(chrom,i,"start"))
  chr_id=as.integer(substr(chrom,4,nchar(chrom)))
  cluster_name=cluster[i,3]
  if(chrom!=old_chrom){
    if(old_chrom!=""){
      run_haplonet(haplotype_cluster,marker_set,old_chrom,start,end)
    }
    if(chr_id<10){
      marker_set <- read.table(paste("Community/Haplotype_cluster/chr0",chr_id,".bim",sep=""), as.is=T, head=F, sep=" ")
    } else{
      marker_set <- read.table(paste("Community/Haplotype_cluster/",chrom,".bim",sep=""), as.is=T, head=F, sep=" ")
    }
    start=min(cluster[which(cluster$V1==chrom),11])+1
    end=max(cluster[which(cluster$V1==chrom),12])+1
    hap=strsplit(cluster[i,20],"")
    hap=as.integer(hap[[1]])
    hap=t(as.matrix(hap))
    haplotype_cluster<-rep("n",end-start+1)
    haplotype_cluster[(cluster[i,11]-start+2):(cluster[i,12]-start+2)]=hap
    Af_Am=cluster[i,18]
    non_Af_Am=cluster[i,17]-Af_Am
    haplotype_cluster=t(replicate(cluster[i,17],haplotype_cluster))
    rownames(haplotype_cluster)=c(rep("Af_Am",times=Af_Am),rep("non_Af_Am",times=non_Af_Am))
    old_chrom=chrom
  }else{
    hap=strsplit(cluster[i,20],"")
    hap=as.integer(hap[[1]])
    hap=t(as.matrix(hap))
    hap_cluster<-rep("n",end-start+1)
    hap_cluster[(cluster[i,11]-start+2):(cluster[i,12]-start+2)]=hap
    Af_Am=cluster[i,18]
    non_Af_Am=cluster[i,17]-Af_Am
    hap_cluster=t(replicate(cluster[i,17],hap_cluster))
    rownames(hap_cluster)=c(rep("Af_Am",times=Af_Am),rep("non_Af_Am",times=non_Af_Am))
    haplotype_cluster=rbind(haplotype_cluster,hap_cluster)
  }
}
net=run_haplonet(haplotype_cluster,marker_set,old_chrom,start,end)
index=index+1
Result[[index]]=net

run_haplonet<-function(haplotyep_cluster,marker_set,chrom,start,end){
  haplotype_cluster[haplotype_cluster==0]="a"
  for(pos in start:end){
    for(ind in 1:dim(haplotype_cluster)[1]){
      if(haplotype_cluster[ind,pos-start+1]==1){
        haplotype_cluster[ind,pos-start+1]=tolower(marker_set[pos,6])
      } 
    }
  }
  print("prepare input for haplonet")
  DNAbin_cluster<-as.DNAbin(haplotype_cluster)
  h <- haplotype(DNAbin_cluster,d=NULL)
  (net <- haploNet(h))
  ind.hap<-with(
    stack(setNames(attr(h, "index"), rownames(h))), 
    table(hap=ind, pop=rownames(haplotype_cluster)[values])
  )
  print("start drawing")
  pdf(paste("Community/Haplotype_cluster/Af_Am_",chrom,"_haplonet_after_count.pdf",sep=""),10,7)  
  plot(net,size = attr(net,"freq"), col = "black", bg = "white",col.link = "black", lwd = 1, lty = 1, pie = ind.hap,labels = TRUE, font = 1, cex = 1, scale.ratio = 6,asp = 1, legend = FALSE, fast = FALSE, show.mutation = 2,main=paste("Haplotype Network enriched in Af Am ",chrom,sep=""))
  legend("topleft", colnames(ind.hap), col=rainbow(ncol(ind.hap)), pch=20)
  dev.off()
  return(net)
}

plot(net,size = attr(net,"freq"), col = "black", bg = "white",col.link = "black", lwd = 1, lty = 1, pie = ind.hap,labels = TRUE, font = 1, cex = 1, scale.ratio = 6,asp = 1, legend = FALSE, fast = FALSE, show.mutation = 0,main=paste("Haplotype Network enriched in Af Am ",chrom,sep=""))
legend("topleft", colnames(ind.hap), col=rainbow(ncol(ind.hap)), pch=20)
#text(100,1000,"13",col="blue",cex=0.8)
#text(-400,100,"19",col="blue",cex=0.8)
#text(1800,2500,"28",col="blue",cex=0.8)
#text(800,1400,"30",col="blue",cex=0.8)

slices <- c(688, 150)
lbls <- c("Af_Am", "non_Af_Am")
pie(slices, labels = lbls,col=rainbow(2))
slices <- c(396, 54)
lbls <- c("Af_Am", "non_Af_Am")
pie(slices, labels = lbls,col=rainbow(2))
