library(ggplot2)
library(reshape2)
library(gridExtra)
library(cowplot)
library(gplots)
library(scales)
library(annotables)


grch <- as.matrix(grch37)
### load marker genes
markers <- read.csv("~/review_paper/pH_homeostasis_responders.csv")
markers <- as.matrix(markers)
mm.names <- names(table(markers[,1]))


### Urea staff
urea <- c("CPS1", "OTC", "ASS1", "ASL", "ARG1")
polyA <- c("ODC1", "AGMAT", "ARG1", "ARG2")
urea.id <- getids(urea)
ployA.id <- getids(polyA)
urea.fig <- Figure.draw(urea.id)
ployA.fig <- Figure.draw(ployA.id)
setwd("~/Downloads/projects/2018Spring/Review_Article/data/files/Heatmaps")
pdf("./polyA_enzyme.pdf", height=7, width=10, onefile=T)
multi.fig(ployA.fig)
graphics.off()
### correlation
cancers <- names(stages.data)
cancers <- names(All.DE)
cancersN <- paste("~/Downloads/projects/GDC/R_data_space/", cancers, "_FPKM_N.RData", sep="")
cancersT <- paste("~/Downloads/projects/GDC/R_data_space/", cancers, "_FPKM_T.RData", sep="")


setwd("~/Downloads/projects/2018Spring/Review_Article/data/files/Heatmaps")
bcaa <- c("SLC16A1", "SLC16A2", "SLC16A4", "SLC16A5", "SLC16A6", "SLC16A7",
"BCAT1", "SLC7A5", "SLC7A6", "SLC7A7", "SLC7A8", "SLC43A1", "SLC43A2",
"SLC43A3", "BCAT2", "BCKDHA", "BCKDHB", "DBT")
bcaa.id <- getids(bcaa)
glutamine <- c("SLC1A5", "SLC6A14", "SLC6A19", "SLC7A5", "SLC7A6", "SLC7A7", "SLC7A8", "SLC7A9",
"SLC38A1", "SLC38A2", "SLC38A3", "SLC38A5", "SLC38A7", "SLC38A8", "GLS", "GLS2","GLUD1",
"GLUD2", "GPT", "GOT1", "GOT2", "OGDH", "DLD", "DLST", "SUCLG1", "SUCLG2", "SUCLA2",
"SDHA", "SDHB", "SDHC", "SDHD", "FH", "MDH1", "MDH2", "CS", "ME1", "ME2", "ME3", "ACLY",
"FASN")
glutamine.id <- getids(glutamine)
serine <- c("PHGDH", "PSAT1", "PSPH")
serine.id <- getids(serine)

pdf("./arginine_correlation_heatmap.pdf", height=7, width=10, onefile=T)
par(cex.main=0.9)
for(i in 1:length(cancersT))
{
	print(i)
	load(cancersT[i])
	data.m <- cor_analysis(all.id, all.id, data_t)
	heatmap.2(data.m, dendrogram="none", col=bluered(128), Rowv=F, Colv=F, 
		breaks=seq(-0.8, 0.8, length.out=129), 
		scale="none",density.info="none", trace="none", cexCol=0.9, cexRow=0.9, 
		key=T, keysize=1, main=cancers[i])
}
dev.off()


###
multi.fig <- function(p)
{
	for(i in 1:length(p))
	{
		pp <- p[[i]]
		print(i)
		print(plot_grid(pp[[1]], pp[[2]], pp[[3]], pp[[4]], ncol=2))
	}
}

### Function Correlated SQ
cor_analysis <- function(x, y, data)
{
	cor.m <- matrix(0, ncol=nrow(x), nrow=nrow(y))
	cor.t <- matrix(0, ncol=nrow(x), nrow=nrow(y))
	x.d <- data[x[,1], ]
	y.d <- data[y[,1], ]
	for(i in 1:nrow(x.d))
	{
		for(j in 1:nrow(y.d))
		{
			cor.m[j,i] <- cor(x.d[i, ], y.d[j, ])
			cor.t[j,i] <- cor.test(x.d[i, ], y.d[j, ])$p.value
		}
	}
	cor.m[cor.t > 0.05] <- 0
	colnames(cor.m) <- rownames(x)
	rownames(cor.m) <- rownames(y)
	return(cor.m)
}


setwd("~/GDC/GDC/R_data_space/")
load("~/GDC/GDC/TCGA_TN_DE_Result.RData")
cancers <- names(All.DE)
cancersN <- paste("~/GDC/GDC/R_data_space/", cancers, "_FPKM_N.RData", sep="")
cancersT <- paste("~/GDC/GDC/R_data_space/", cancers, "_FPKM_T.RData", sep="")


load(cancersT[1])
rown <- unlist(lapply(rownames(data_t), function(x) strsplit(x, split=".", fixed=T)[[1]][1]))

### All genes
genes <- read.table("~/SQ/genes.txt", sep="\t", header=T)
genes <- as.matrix(genes)
genes.id <- getids(genes[,2])

for(i in 1:length(stages.data))
{
	print(i)
	data <- matrix(0, ncol=10, nrow=nrow(genes.id))
	dd <- stages.data[[i]]
	dd1 <- as.matrix(dd[[1]])
	dd2 <- as.matrix(dd[[2]])
	dd3 <- as.matrix(dd[[3]])
	dd4 <- as.matrix(dd[[4]])
	dd5 <- as.matrix(dd[[5]])
	if(i != 12)
	{
		data[,1] <- apply(dd1[genes.id[,1], ], 1, mean)
		data[,2] <- apply(dd2[genes.id[,1], ], 1, mean)
		data[,3] <- apply(dd3[genes.id[,1], ], 1, mean)
		data[,4] <- apply(dd4[genes.id[,1], ], 1, mean)
		data[,5] <- apply(dd5[genes.id[,1], ], 1, mean)
		data[,6] <- apply(dd1[genes.id[,1], ], 1, median)
		data[,7] <- apply(dd2[genes.id[,1], ], 1, median)
		data[,8] <- apply(dd3[genes.id[,1], ], 1, median)
		data[,9] <- apply(dd4[genes.id[,1], ], 1, median)
		data[,10] <- apply(dd5[genes.id[,1], ], 1, median)
	}
	else
	{
		data[,1] <- dd1[genes.id[,1], 1]
		data[,2] <- apply(dd2[genes.id[,1], ], 1, mean)
		data[,3] <- apply(dd3[genes.id[,1], ], 1, mean)
		data[,4] <- apply(dd4[genes.id[,1], ], 1, mean)
		data[,5] <- apply(dd5[genes.id[,1], ], 1, mean)
		data[,6] <- dd1[genes.id[,1], 1]
		data[,7] <- apply(dd2[genes.id[,1], ], 1, median)
		data[,8] <- apply(dd3[genes.id[,1], ], 1, median)
		data[,9] <- apply(dd4[genes.id[,1], ], 1, median)
		data[,10] <- apply(dd5[genes.id[,1], ], 1, median)
	}
	colnames(data) <- c("Norm_mean", "St1_mean", "St2_mean", "St3_mean", "St4_mean", "Norm_median", 
		"St1_median", "St2_median", "St3_median", "St4_median")
	rownames(data) <- rownames(genes.id)
	write.csv2(data, file=paste(names(stages.data)[i], "_expression.csv", sep=""), quote=F)
}

fatty.absorb <- markers[markers[,1]==mm.names[5], 2]
fatty.absorb <- getids(fatty.absorb)

### define function to get gene ids
getids <- function(genes)
{
	index <- match(genes, grch[,3])
	gg.en <- grch[index[!is.na(index)], 1]
	ggs <- genes[!is.na(index)]
	ggs <- ggs[!is.na(match(gg.en, rown))]
	gg.id <- match(gg.en, rown)[!is.na(match(gg.en, rown))]
	data <- matrix(0, ncol=1, nrow=length(gg.id))
	data[,1] <- gg.id
	rownames(data) <- ggs
	return(data)
}


### Clinical information with stages
cancers <- tolower(unlist(lapply(cancers, function(x) strsplit(x, split="-", fixed=T)[[1]][2])))
clinics <- paste("nationwidechildrens.org_clinical_patient_", cancers, ".txt.csv", sep="")
##Only BLCA & UCEC has grade information
if("neoplasm_histologic_grade" %in% colnames(clinical)) 
{print(clinics[i])
	print(table(clinical[ ,"neoplasm_histologic_grade"]))}

### Clinical data
setwd("~/GDC/GDC/OldVersion/Clinical/")
clinics <- list.files()
clinics <- clinics[grep(".txt.csv", clinics)]
cancer.clinic <- unlist(lapply(clinics, function(x) strsplit(x, split="_", fixed=T)[[1]][4]))
cancer.clinic <- unlist(lapply(cancer.clinic, function(x) strsplit(x, split=".txt")[[1]][1]))
cancer.clinic <- toupper(cancer.clinic)
cancersN <- paste("~/GDC/GDC/R_data_space/TCGA-", cancer.clinic, "_FPKM_N.RData", sep="")
cancersT <- paste("~/GDC/GDC/R_data_space/TCGA-", cancer.clinic, "_FPKM_T.RData", sep="")


nums <- c()
stages.data <- list()
j <- 1
for(i in 1:length(clinics))
{
	clinical <- read.csv(clinics[i])
	clinical <- as.matrix(clinical)
#	siaGALNAC.fig[[i]] <- list()
#	siaGAL.fig[[i]] <- list()
#	siaSIA.fig[[i]] <- list()
	if("ajcc_pathologic_tumor_stage" %in% colnames(clinical))
	{
	print(clinics[i])
	print(clinics[i])
	nums <- c(nums, i)
	st1 <- clinical[clinical[,"ajcc_pathologic_tumor_stage"] %in% c("Stage Tis", "Stage X", "Stage I", "Stage IA", "Stage IB", "Stage IC"), 2]
	st2 <- clinical[clinical[,"ajcc_pathologic_tumor_stage"] %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC"), 2]
	st3 <- clinical[clinical[,"ajcc_pathologic_tumor_stage"] %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"), 2]
	st4 <- clinical[clinical[,"ajcc_pathologic_tumor_stage"] %in% c("Stage IV", "Stage IVA", "Stage IVB", "Stage IVC"), 2]
	load(cancersT[i])
	load(cancersN[i])
	coln <- colnames(data_t)
	coln <- unlist(lapply(coln, function(x) paste(strsplit(x, split="-")[[1]][1:3], collapse="-")))
#	tn.coln <- match(tn.samples, coln)
	st1.index <- match(st1, coln)
	st2.index <- match(st2, coln)
	st3.index <- match(st3, coln)
	st4.index <- match(st4, coln)
	st1.data <- data_t[, st1.index[!is.na(st1.index)]]
	st2.data <- data_t[, st2.index[!is.na(st2.index)]]
	st3.data <- data_t[, st3.index[!is.na(st3.index)]]
	st4.data <- data_t[, st4.index[!is.na(st4.index)]]
	stages.data[[j]] <- list()
	stages.data[[j]][[1]] <- data_n
	stages.data[[j]][[2]] <- st1.data
	stages.data[[j]][[3]] <- st2.data
	stages.data[[j]][[4]] <- st3.data
	stages.data[[j]][[5]] <- st4.data
	names(stages.data[[j]]) <- c("Normal", "Stage I", "Stage II", "Stage III", "Stage IV")
	j <- j + 1
	}
}
names(stages.data) <- cancer.clinic[nums]


### define function to draw figures
Figure.draw <- function(genes)
{
	fig <- list()
	for(i in 1:length(stages.data))
	{
		fig[[i]] <- list()
		print(i)
		n <- stages.data[[i]][[1]][genes[,1], ]
		n <- as.matrix(n)
		t1 <- stages.data[[i]][[2]][genes[,1], ]
		t2 <- stages.data[[i]][[3]][genes[,1], ]
		t3 <- stages.data[[i]][[4]][genes[,1], ]
		t4 <- stages.data[[i]][[5]][genes[,1], ]
		for(j in 1:nrow(genes))
		{
			data <- c(n[j,], t1[j,], t2[j,], t3[j,], t4[j,])
			dg <- log(data+1)
			dg <- data.frame(y=dg, x=c(rep("Normal", ncol(n)), rep("Stage I", ncol(t1)), 
				rep("Stage II", ncol(t2)), rep("Stage III", ncol(t3)), rep("Stage IV", ncol(t4))))
			fig[[i]][[j]] <- ggplot(dg, aes(x=x, y=y, fill=factor(x))) + geom_boxplot() + labs(y="Log FPKM Expression", x="") 
			fig[[i]][[j]] <- fig[[i]][[j]] + theme(legend.text=element_text(size=5), 
				axis.text=element_blank(), axis.title=element_text(size=8), legend.position="bottom")
			fig[[i]][[j]] <- fig[[i]][[j]] + scale_fill_discrete("") 
			fig[[i]][[j]] <- fig[[i]][[j]] + ggtitle(paste(names(stages.data)[i], rownames(genes)[j], sep="."))
		}
		names(fig[[i]]) <- rownames(genes)
	}
	names(fig) <- names(stages.data)
	return(fig)
}

### run figures
fatty.absorb.fig <- Figure.draw(fatty.absorb)
setwd("~/review_paper/")
pdf("./Fatty_acids_aborption_gene_changes.pdf", height=8, width=11, onefile=T)
multi.fig(fatty.absorb.fig)
graphics.off()


###TN information
	tn.st1 <- intersect(tn.coln, st1.index)
	tn.st2 <- intersect(tn.coln, st2.index)
	tn.st3 <- intersect(tn.coln, st3.index)
	tn.st4 <- intersect(tn.coln, st4.index)
	st1.data <- data_t[, tn.st1[!is.na(tn.st1)]]
	st2.data <- data_t[, tn.st2[!is.na(tn.st2)]]
	st3.data <- data_t[, tn.st3[!is.na(tn.st3)]]
	st4.data <- data_t[, tn.st4[!is.na(tn.st4)]]
	brca.tn <- list()
	brca.tn[[1]] <- data_n
	brca.tn[[2]] <- st1.data
	brca.tn[[3]] <- st2.data
	brca.tn[[4]] <- st3.data
	brca.tn[[5]] <- st4.data
	names(brca.tn) <- c("Normal", "Stage I", "Stage II", "Stage III", "Stage IV")
####


### Cytoscape
setwd("~/Downloads/projects/2018Spring/Review_Article/data/files/fig1/")
fig <- read.delim("~/Downloads/projects/2018Spring/Review_Article/data/Fig1.txt", header=F)
fig <- as.matrix(fig)
fig <- fig[,c(1:3)]
for(i in 1:length(All.DE))
{
	write.table(fig, file=paste(names(All.DE)[i], "_fig1.txt", sep=""), sep="\t", 
		row.names=F, col.names=F, quote=F)
}

library("RCy3")
load()
PRMT <- c("SLC22A2", "PRMT1", "PRMT2", "PRMT3", "PRMT5", "PRMT6", "PRMT7", "PRMT8", "PRMT9")
DDAH <- c("DDAH1", "DDAH2", "PYCR1")
ASS <- c("ASS1", "ASL", "AZIN2", "SRM", "SMS", "RAI1", "SRMS")
chain <- c("ALDH18A1", "OAT", "ODC1", "SLC7A1", "SAT1", "SLC26A1", "SMO", "SMOX", 
	"SLC3A2")

all <- c(PRMT, DDAH, ASS, chain)
all.id <- getids(all)

setwd("~/Downloads/projects/2018Spring/Review_Article/data/data_cancers/")
#cancers <- paste("TCGA-", names(stages.data), sep="")
for(i in 1:length(All.DE))
{
	data <- All.DE[[i]]
	data <- data[all.id[,1], ]
	colors <- rep("black", nrow(data))
	colors[data[,1]<=0.05 & data[,2]>1] <- "red"
	colors[data[,1]<=0.05 & data[,2]<1] <- "blue"
	data <- cbind(rownames(all.id), data, colors)
	colnames(data)[1] <- "genes" 
	write.table(data, file=paste("All_", names(All.DE)[i], ".txt", sep=""), quote=F,
		row.names=F, col.names=T, sep="\t")
}

Cytoscape
1. import network fig1
2. import file edge
3. stlye -- sample1
4. line type dash
5. Label colors
6. Stroke Color ()


### Linux change
cat All_TCGA-BLCA.txt|cut -f4|sed 's/0/black/g; s/-1/blue/g; s/1/red/g; s/s/s_n/g' |
paste -d "\t" All_TCGA-BLCA.txt - > ./All_TCGA-BLCA.txt

for f in *.txt; do
	cat $f |cut -f4|sed 's/0/black/g; s/-1/blue/g; s/1/red/g; s/s/s_n/g' | 
	paste -d "\t" $f - > $(echo new_$f)
done

for f in *.txt; do
	echo new_$f
done

for i in *.doc ; do mv "$i" $(echo $i | sed s/doc/txt/) ; done






