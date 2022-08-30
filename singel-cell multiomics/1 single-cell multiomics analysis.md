#Singel-cell multiomics analysis

##1 Tert gene expression level

###1.1 Data download
```Shell
wget -c https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE136714&format=file&file=GSE136714%5Fraw%2Ecounts%2Efor%2Egeo%2Exlsx
```

###1.2 Analysis process
```R
library(data.table)
library(DESeq2)

library(readxl)
dat<-read_xlsx("GSE136714_raw.counts.for.geo.xlsx")
dat1<-data.frame(t(dat))
colnames(dat1)<-dat1[1,]
dat2<-dat1[-1,]
dat2$cell<-gsub("split1","split2",rownames(dat2))
dat2$cell_raw<-rownames(dat2)
dat3<-data.frame(t(dat2))[1:37405,]
for (i in 1:401) {
  dat3[,i]<-as.numeric(dat3[,i])
}

condition<-factor(dat2$cell)
colData<-data.frame(row.names=colnames(dat3),condition)
dds<-DESeqDataSetFromMatrix(countData = dat3,colData =colData,design = ~condition)

dds$sample<-gsub("split1","split2",colnames(dds))

ddsColl<-collapseReplicates(dds,dds$sample)

final<-data.frame(ddsColl@assays@data@listData[["counts"]])
final<-final["ENSMUSG00000021611",]
final1<-data.frame(t(final))
final1$Tert<-log(1+final1$ENSMUSG00000021611,base=2)
final1$cell<-sapply(strsplit(rownames(final1),"_"),"[",1)
final2<-subset(final1,cell=="X8cell"| cell=="X32cell")
table(final2$cell)
final2$cell_raw<-rownames(final2)
final2$cell_raw<-gsub("_split1","",final2$cell_raw)
final2$cell_raw<-gsub("_split2","",final2$cell_raw)

cellid<-read.csv("cellID.csv")
cellid$cell_raw<-paste0("X",cellid$cell_raw)
cellid$cell_raw<-gsub(" ","",cellid$cell_raw)
final3<-merge(final2,cellid,by="cell_raw")
write.csv(final3,"8cell_32cell_Tert_counts.csv",quote=F,row.names=F)
```


##2 Tert gene methylation level

###2.1 Data download and transformation
```Shell
wget -c https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE136715&format=file
tar xvf GSE136715_RAW.tar
bigWigToBedGraph ${sample}.WCG.bw ${sample}.WCG.bedGraph
```

###2.2 Analysis process
```R
library(data.table)
a<-list.files(pattern="*.WCG.bedGraph",full.names=T)
for(i in 1:length(a)){
	dat = read.table(a[i], header = F)
	ID<-gsub("./","",a[i])
	ID1<-paste0("/public/user/zj2020/meth_single/meth_res3/",ID)
	dat1<-subset(dat,V1=="chr13")
	dat2<-subset(dat1,V2>=73625001 & V3<=73649041) #73627101
	write.table(dat2,ID1,row.names=F,quote=F,col.names=F)
}


library(data.table)
a<-list.files(pattern="*.WCG.bedGraph",full.names=T)
dat = read.table(a[1], header = F)
library(data.table)
a<-list.files(pattern="*.WCG.bedGraph",full.names=T)
dat = read.table(a[1], header = F)
ID1<-sapply( strsplit(a[1] , "_") , "[" , 1 )
ID2<-gsub("./","",ID1)

final<-NULL
final$IID<-ID2
final$meth_sum_c<-sum(dat$V4)
final$pos_c<-nrow(dat)
final<-data.frame(final)


for (i in 2:length(a)){
	dat = fread(a[i])
	ID1<-sapply( strsplit(a[i] , "_") , "[" , 1 )
	ID2<-gsub("./","",ID1)
	tmp<-NULL
	tmp$IID<-ID2
	tmp$meth_sum_c<-sum(dat$V4)
	tmp$pos_c<-nrow(dat)
	tmp<-data.frame(tmp)
	final<-rbind(final,tmp)
}

write.csv(final,"GEO_meth_level.csv",row.names=F,quote=F)
```

##3 Multiomics analysis
```R
meth=read.csv("GEO_meth_level.csv")
rna=read.csv("8cell_32cell_Tert_counts.csv")
tert<-merge(rna,meth,by.x="ID",by.y="IID")
tert1<-subset(tert,pos_c!=0)
wilcox.test(tert1$Tert~tert1$cell)

#whole gene methylaion level
meth2<-read.csv("meth_level_res.csv")
tert1<-merge(tert1,meth2,by.x="ID",by.y="IID")
wilcox.test(tert1$level~tert1$cell)


tert1$meth_level_group<-ifelse(tert1$level<=0.5,0,1)
wilcox.test(tert1$Tert~tert1$meth_level_group)
```