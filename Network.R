###Cell.Cycle
##load library
library(tidyverse)
library(stringr)
library(scales)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(cowplot)
library(gplots)
library(scales)
library(annotables)
library(ggrepel)
library(ggraph)
library(widyr)
suppressPackageStartupMessages(library(igraph))

##load data
load("~/Cancers_Express.RData")
load("~/SQ/grch37.RData")
setwd("~/GDC/GDC/R_data_space/")
load("~/GDC/GDC/TCGA_TN_DE_Result.RData")
cancers <- names(All.DE)
cancersN <- paste("~/GDC/GDC/R_data_space/", cancers, "_FPKM_N.RData", sep="")
cancersT <- paste("~/GDC/GDC/R_data_space/", cancers, "_FPKM_T.RData", sep="")
load(cancersT[1])
rown <- unlist(lapply(rownames(data_t), function(x) strsplit(x, split=".", fixed=T)[[1]][1]))


load(cancersT[1])
rown <- unlist(lapply(rownames(data_t), function(x) strsplit(x, split=".", fixed=T)[[1]][1]))

##functions
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


### CellCycles genes
cAMP <- c("CREB1", "ATF1", "CREBBP", "CREBZF", "ATF4", "DDIT3", "FOXO1", "ATF2", "CREM")
cSIG <- c("AKT1", "AKT2", "AKT3", "PIK3CA", "HRAS", "KRAS", "NRAS", "MAPK3", "PRKCA", "MTOR")
cReg <- c("CDC25A", "CDC25C", "CDKN1B", "CTTN", "ARPC1A", "ARPC1B", "ARPC2", "ARPC3", "ARPC4", 
	"ARPC5")
ccg <- c(cAMP, cSIG, cReg)
### gene ids
cAMP.id <- getids(cAMP)
cSIG.id <- getids(cSIG)
cREG.id <- getids(cReg)
ccg.id <- c(cAMP.id, cSIG.id, cREG.id)
ccg.id <- as.matrix(ccg.id)
rownames(ccg.id) <- ccg

### Calculate correlation
options(repr.plot.width = 13, repr.plot.height = 8)
for(i in 1:length(cancersT))
{
	load(cancersT[i])
	data.m <- cor_analysis(ccg.id, ccg.id, data_t)
	data.m[lower.tri(data.m, diag=T)] <- NA
	cors_data <- melt(data.m, na.rm=T)
	cors_data <- as.matrix(cors_data)
	cors_data <- tbl_df(cors_data)
	cors_data <- cors_data %>% mutate(value=as.numeric(value))
	cors_data <- cors_data %>% filter(abs(value) > .05)
	cors_data <- cors_data %>% mutate(col=ifelse(value>=0, 1, -1))
	vertices <- All.DE[[i]][ccg.id[,1], ]
	vertices <- as.data.frame(vertices)
	vertices$nodes <- ccg
	vertices$class <- c(rep("cAMP", 9), rep("Signaling", 10), rep("Actin", 10))
	vertices <- vertices[,c(3,1,2,4)]
	vertices <- tbl_df(vertices)
	set.seed(2017)
	cors_data %>% graph_from_data_frame(vertices=vertices) %>%
    ggraph(layout = "fr") +
    geom_edge_fan(aes(edge_alpha=value, alpha=0.25, colour=col), show.legend=F, check_overlap=T) +
    geom_edge_density() +
    geom_node_point(aes(size=2, shape=class, color=pmin(2, fold.change))) +
    scale_size_continuous(range = c(.01, 5), breaks = seq(.05, .2, .05)) +
    geom_node_text(aes(label = name), repel = TRUE, size = 4)  +                
	theme_void(base_size = 16) + 
    scale_color_gradient2(low="#512E5F", high="#78281F", trans="log", breaks=2^seq(-2, 2, 0.5), labels=rep("", length(seq(-2, 2, 0.5)))) +
    labs(title=names(All.DE)[i],
        subtitle="Expression Network of cAMP, signaling and Actin", 
        size="", color="Fold-Change") +
   expand_limits(x = -1)
}

# + facet_nodes(~class)
# + facet_edges(~col)
# + facet_graph(col ~ class, margin=T)
