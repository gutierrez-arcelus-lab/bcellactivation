library(DESeq2)
library(vsn)
library(BiocParallel)

count.table <- read.delim("./results/bwa/merged_library/macs2/narrow_peak/consensus/consensus_peaks.mLb.clN.featureCounts.txt", header = TRUE, skip = 1)
colnames(count.table) <- gsub("\\.mLb\\.\\clN\\.sorted\\.bam", "", colnames(count.table))
outprefix <- "consensus_peaks.mLb.clN"
opt <- data.frame(outsuffix = '', outprefix = outprefix, 
		  outdir = "./results_deseq2", cores = 4, 
		  vst = FALSE)

rownames(count.table) <- count.table$Geneid
interval.table <- count.table[,1:6]
count.table <- count.table[,7:ncol(count.table),drop=FALSE]

## RUN DESEQ2

if (file.exists(opt$outdir) == FALSE) {
    dir.create(opt$outdir,recursive=TRUE)
}

setwd(opt$outdir)

samples.vec <- sort(colnames(count.table))
groups <- sub("_[^_]+$", "", samples.vec)
print(unique(groups))

if (length(unique(groups)) == 1) {
    quit(save = "no", status = 0, runLast = FALSE)
}

DDSFile <- paste(opt$outprefix,".dds.rld.RData",sep="")

if ( !file.exists(DDSFile) ) {
    counts <- count.table[,samples.vec,drop=FALSE]
    coldata <- data.frame(row.names=colnames(counts),condition=groups)
    dds <- DESeqDataSetFromMatrix(countData = round(counts), colData = coldata, design = ~ condition)
    dds <- DESeq(dds, parallel=TRUE, BPPARAM=MulticoreParam(opt$cores))
    if (!opt$vst) {
        rld <- rlog(dds)
    } else {
        rld <- vst(dds)
    }
    save(dds,rld,file=DDSFile)
}

## PCA
pca.data <- DESeq2::plotPCA(rld, intgroup = "condition", ntop = 5000, returnData=TRUE)
saveRDS(pca.data, paste0(opt$outprefix, "_pcadata_5000peaks.rds"))

SizeFactorsDir <- "sizeFactors"
if ( !file.exists(SizeFactorsDir) ) {
    dir.create(SizeFactorsDir, recursive=TRUE)
}

NormFactorsFile <- paste(SizeFactorsDir, opt$outprefix, ".sizeFactors.RData", sep = "")
if ( !file.exists(NormFactorsFile) ) {
    normFactors <- sizeFactors(dds)
    save(normFactors, file = NormFactorsFile)

    for (name in names(sizeFactors(dds))) {
        sizeFactorFile <- paste(SizeFactorsDir, name, opt$outsuffix, ".sizeFactor.txt", sep = "")
        if (file.exists(sizeFactorFile) == FALSE) {
            write(as.numeric(sizeFactors(dds)[name]), file = sizeFactorFile)
        }
    }
}

## WRITE LOG FILE

LogFile <- paste(opt$outprefix, ".log", sep = "")
if (file.exists(LogFile) == FALSE) {
    cat("\nSamples =", samples.vec, "\n\n", file = LogFile, append = TRUE, sep = ', ')
    cat("Groups =", groups, "\n\n", file = LogFile, append = TRUE, sep = ', ')
    cat("Dimensions of count matrix =", dim(counts), "\n\n", file = LogFile, append = FALSE, sep = ' ')
    cat("\n", file = LogFile, append = TRUE, sep = '')
}

## LOOP THROUGH COMPARISONS                   

ResultsFile <- paste(opt$outprefix, ".results.txt",sep="")
if (file.exists(ResultsFile) == FALSE) {

    raw.counts <- counts(dds, normalized = FALSE)
    colnames(raw.counts) <- paste(colnames(raw.counts), 'raw', sep = '.')
    pseudo.counts <- counts(dds, normalized = TRUE)
    colnames(pseudo.counts) <- paste(colnames(pseudo.counts), 'pseudo', sep = '.')

    deseq2_results_list <- list()
    #comparisons <- combn(unique(groups),2)
    comparisons <- combn(c("unst_0", "unst_24", unique(groups)[1:8]), 2)
    for (idx in 1:ncol(comparisons)) {

        control.group <- comparisons[1,idx]
        treat.group <- comparisons[2,idx]
        CompPrefix <- file.path("./results_deseq2", paste(control.group, treat.group, sep = "vs"))
        #CompPrefix <- paste(control.group,treat.group,sep="vs")
        cat("Saving results for ", CompPrefix, " ...\n", sep = "")

        CompOutDir <- paste(CompPrefix, '/', sep = "")
        if (file.exists(CompOutDir) == FALSE) {
            dir.create(CompOutDir, recursive = TRUE)
        }

        control.samples <- samples.vec[which(groups == control.group)]
        treat.samples <- samples.vec[which(groups == treat.group)]
        comp.samples <- c(control.samples, treat.samples)

        comp.results <- results(dds, contrast = c("condition", c(control.group, treat.group)))
        comp.df <- as.data.frame(comp.results)
        comp.table <- cbind(interval.table, as.data.frame(comp.df), raw.counts[ ,paste(comp.samples, 'raw', sep = '.')], pseudo.counts[ ,paste(comp.samples, 'pseudo', sep = '.')])

        ## WRITE RESULTS FILE
        CompResultsFile <- paste(CompOutDir, CompPrefix, opt$outsuffix, ".deseq2.results.txt", sep = "")
        write.table(comp.table, file = CompResultsFile, col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)

        ## FILTER RESULTS BY FDR & LOGFC AND WRITE RESULTS FILE
        #pdf(file=paste(CompOutDir,CompPrefix,opt$outsuffix,".deseq2.plots.pdf",sep=""),width=10,height=8)
        if (length(comp.samples) > 2) {
            for (MIN_FDR in c(0.01,0.05)) {

                ## SUBSET RESULTS BY FDR
                pass.fdr.table <- subset(comp.table, padj < MIN_FDR)
                pass.fdr.up.table <- subset(comp.table, padj < MIN_FDR & log2FoldChange > 0)
                pass.fdr.down.table <- subset(comp.table, padj < MIN_FDR & log2FoldChange < 0)

                ## SUBSET RESULTS BY FDR AND LOGFC
                pass.fdr.logFC.table <- subset(comp.table, padj < MIN_FDR & abs(log2FoldChange) >= 1)
                pass.fdr.logFC.up.table <- subset(comp.table, padj < MIN_FDR & abs(log2FoldChange) >= 1 & log2FoldChange > 0)
                pass.fdr.logFC.down.table <- subset(comp.table, padj < MIN_FDR & abs(log2FoldChange) >= 1 & log2FoldChange < 0)

                ## WRITE RESULTS FILE
                CompResultsFile <- paste(CompOutDir, CompPrefix, opt$outsuffix, ".deseq2.FDR",MIN_FDR, ".results.txt", sep = "")
                CompBEDFile <- paste(CompOutDir, CompPrefix, opt$outsuffix, ".deseq2.FDR", MIN_FDR, ".results.bed", sep = "")
                write.table(pass.fdr.table, file = CompResultsFile, col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)
                write.table(pass.fdr.table[, c("Chr", "Start", "End", "Geneid", "log2FoldChange", "Strand")], file = CompBEDFile, col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)

                ## ADD COUNTS TO LOGFILE
                cat(CompPrefix, " genes with FDR <= ", MIN_FDR, ": ", nrow(pass.fdr.table), " (up=",nrow(pass.fdr.up.table),", down=", nrow(pass.fdr.down.table),")","\n", file = LogFile, append = TRUE, sep = "")
                cat(CompPrefix, " genes with FDR <= ", MIN_FDR, " & FC > 2: ", nrow(pass.fdr.logFC.table), " (up=",nrow(pass.fdr.logFC.up.table),", down=",nrow(pass.fdr.logFC.down.table),")","\n", file = LogFile,append = TRUE, sep = "")

            }
            cat("\n", file = LogFile, append = TRUE, sep = "")
        }

        colnames(comp.df) <- paste(CompPrefix, ".", colnames(comp.df), sep = "")
        deseq2_results_list[[idx]] <- comp.df

    }

    ## WRITE RESULTS FROM ALL COMPARISONS TO FILE
    deseq2_results_table <- cbind(interval.table, do.call(cbind, deseq2_results_list), raw.counts, pseudo.counts)
    write.table(deseq2_results_table, file = ResultsFile, col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)

}

## R SESSION INFO

RLogFile <- "R_sessionInfo.log"
if (file.exists(RLogFile) == FALSE) {
    sink(RLogFile)
    a <- sessionInfo()
    print(a)
    sink()
}
