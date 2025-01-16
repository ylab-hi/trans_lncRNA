library("edgeR")
library("limma")


# args <- commandArgs(trailingOnly = TRUE);
# input.matrix = args[1];

src.data <- read.delim("expression.txt", row.names="GeneID", check.names=FALSE);
data <- src.data[,seq(2, dim(src.data)[2])];

contract <- read.delim("./expression.contract", check.names=FALSE, sep="\t", header=T);

# normal: 52; primary: 500; su2c:101
#batch <- c(rep(1, 552), rep(2, 101))

group <- as.factor(contract$group)


y <- DGEList(counts=data, group=group);
y$genes <- data.frame(Length=src.data$Length);
rownames(y$genes) <- rownames(y$counts);

#y$genes$Symbol <- data.frame(Length=src.data$Symbol);
#rownames(y$symbol) <- rownames(y$counts);

# filtering out low expressed genes
#  the minimum number of samples in each group is three, over here.
keep <- rowSums(cpm(y)>1) >= 52
table(keep);
# FALSE  TRUE
# 7223 10431
# 4527 16416
y <- y[keep, keep.lib.sizes=FALSE];


design <- model.matrix(~0+group);
rownames(design) <- colnames(y)
design
# the filtering does not bias the subsequent differential expression analysis.
# The TMM normalization is applied to account for the compositional biases:


# normalization by the library sizes
y <- calcNormFactors(y);
y$samples;

data.rpkm <- rpkm(y);
write.table(data.frame("GeneID"=rownames(data.rpkm), data.rpkm, check.names=FALSE), file=paste0('rpkm','.count.txt'), sep='\t', row.names=FALSE,quote = FALSE);
## write.table(rpkm(y), file=paste0('rpkm','.count.txt'), sep='\t', row.names = TRUE, quote = FALSE);

barplot(y$samples$lib.size,names=colnames(y),las=2)
# Add a title to the plot
title("Barplot of library sizes")

# plotMDS(y);

# normalise the read counts with 'voom' function
v <- voom(y,design,plot = TRUE)
# extract the normalised read counts
counts.voom <- v$E

boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(v$E),col="blue")

# save normalised expression data into output dir
#write.table(counts.voom,file="counts.voom.txt",row.names=T,quote=F,sep="\t");
write.table(data.frame("GeneID"=rownames(counts.voom), counts.voom, check.names=FALSE),file="counts.voom.txt",row.names=F,quote=F,sep="\t");

######################################################################################
######################################################################################
# fit linear model for each gene given a series of libraries
fit <- lmFit(v, design)


############## Mets VS Normal #################
# construct the contrast matrix corresponding to specified contrasts of a set of parameters
matrix.3vs1 <- makeContrasts(group3-group1,levels=design)
# compute estimated coefficients and standard errors for a given set of contrasts
fit.3vs1 <- contrasts.fit(fit, matrix.3vs1)

# compute moderated t-statistics of differential expression by empirical Bayes moderation of the standard errors
fit.3vs1 <- eBayes(fit.3vs1)
options(digits=3)

#colnames(fit$coefficients)

summary(decideTests(fit.3vs1, p.value=0.05,lfc=1))

num = length(fit.3vs1$genes$Length)
degs.3vs1 <- topTable(fit.3vs1, coef="group3 - group1", confint=TRUE, number = num)
#write.table(degs.3vs1, file=paste0('Mets.vs.Normal','.degs.txt'), sep='\t',row.names = TRUE, quote = FALSE);
write.table(data.frame("GeneID"=rownames(degs.3vs1), degs.3vs1, check.names=FALSE), file=paste0('Mets.vs.Normal','.degs.txt'), sep='\t',row.names = FALSE, quote = FALSE);
