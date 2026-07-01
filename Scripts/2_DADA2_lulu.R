# =============================================================================
# 2_DADA2_lulu.R
# DADA2 amplicon sequence variant (ASV) inference, LULU curation,
# negative-control decontamination, and phyloseq object construction
#
# Input:  primer-trimmed .fq.gz files per plate from 1_demultiplexing.R
# Output: phyloseq .rds objects (rg2 PCR-replicate filter; nopool/pool/pspool)
#
# Run inside R (same conda environment used for demultiplexing)
#   conda activate my_r_env
#   R
# =============================================================================

# ---- LIBRARIES --------------------------------------------------------------
library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

library(stringr)
library(abind)
library(tidyverse)
library(dplyr)
library(data.table)
library(magrittr)

library(devtools)
library(lulu)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
library(phyloseq)
library(seqinr)
library(ggplot2)

# =============================================================================
# SECTION 1 — DADA2 READ PROCESSING (P1 + P2 processed together)
# =============================================================================

## Making filepath based on where the trimmed files are.
## Follow the DADA2 tutorial while trying these steps:
## https://benjjneb.github.io/dada2/tutorial.html

setwd("/home/daniel/Ptarmigan/trimmed/alltrim/")
path <- "/home/daniel/Ptarmigan/trimmed/alltrim/"
list.files(path)

# Both plates (P1 + P2) from the same sequencing run — load all files together
fnFs <- sort(list.files(path, pattern="\\.trim1\\.fq\\.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="\\.trim2\\.fq\\.gz", full.names = TRUE))

# ---- Tidy up sample names and replicate labels ------------------------------
# Extract sample name (everything before "_concat__") and assign rep labels.
# Files sort as rep1/rep2 alternating per sample, so rep() cycles correctly.
sample.names <- sub("_concat__.*", "", basename(fnFs))
reps <- rep(c("r1", "r2"), times = length(fnFs) / 2)

sample.names <- paste0(sample.names, reps)

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

# ---- Quality inspection -----------------------------------------------------
plotQualityProfile(fnFs[1:20])  
plotQualityProfile(fnRs[1:20]) 
# Quality looks fine - no truncation needed

# ---- Quality filtering ------------------------------------------------------
# Did not truncate reads (already short, no tail quality problems)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(2,2), minLen=50, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows: multithread=FALSE
head(out, n=20)

out_df <- as.data.frame(out)
out_df$pct_kept <- round(out_df$reads.out / out_df$reads.in * 100, 1)
out_df[grepl("^(S_|W_)", rownames(out_df)), ]
head(out_df)
# ---- Identify samples that failed filtering ---------------------------------
df.fe <- data.frame(theref  = file.exists(filtFs),
                    therer  = file.exists(filtRs),
                    filef   = filtRs,
                    filer   = filtRs)
subset(df.fe, theref == "FALSE") 
subset(df.fe, therer == "FALSE")

# If any samples failed filtering (check subset() output above), remove by:
# failed <- !file.exists(filtFs)
#filtFs <- filtFs[!failed]
#filtRs <- filtRs[!failed]
#sample.names <- sample.names[!failed]


# ---- Error modelling (standard approach) ------------------------------------
errF <- learnErrors(filtFs, multithread=TRUE, MAX_CONSIST=20)
errR <- learnErrors(filtRs, multithread=TRUE, MAX_CONSIST=20)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# ---- NovaSeq quality binning fix for error models ---------------------------
# NovaSeq bins quality scores, which violates DADA2 error model assumptions.
# Implementing the solution by JacobRPrice:
# https://github.com/benjjneb/dada2/issues/1307
# (alter loess arguments + enforce monotonicity)

loessErrfun_mod <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # Guillem Salazar's solution
        # https://github.com/benjjneb/dada2/issues/938
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot), span = 2)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      }
    }
  }
  
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # Enforce monotonicity  
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  return(err)
}

errF <- learnErrors(filtFs, multithread=TRUE, errorEstimationFunction=loessErrfun_mod, verbose=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE, errorEstimationFunction=loessErrfun_mod, verbose=TRUE)

# Check new error plots - should show monotonic decrease and better fit to black points
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# ---- Dereplicate ------------------------------------------------------------
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

# ---- DADA2 core denoising: three pooling strategies -------------------------
# Checkpointing: each dada() result is saved immediately after completion.
# If R crashes, re-running this block reloads completed results and skips them.
ckpt <- "/home/daniel/Ptarmigan/trimmed/dada2_checkpoints"
dir.create(ckpt, showWarnings=FALSE)

dada_or_load <- function(path, expr) {
  if (file.exists(path)) { message("Loading checkpoint: ", path); return(readRDS(path)) }
  res <- expr; saveRDS(res, path); res
}

# multithread=TRUE may only use 1 thread if RcppParallel defaults are capped by the
# environment (OMP_NUM_THREADS etc). Set explicitly to use 50 of the available CPUs.
RcppParallel::setThreadOptions(numThreads=50)

# No pooling (default)
dadaFs   <- dada_or_load(file.path(ckpt,"dadaFs.rds"),
              dada(derepFs, err=errF, multithread=50))
dadaRs   <- dada_or_load(file.path(ckpt,"dadaRs.rds"),
              dada(derepRs, err=errR, multithread=50))

# True pooling
dadaFPPs <- dada_or_load(file.path(ckpt,"dadaFPPs.rds"),
              dada(derepFs, err=errF, multithread=50, pool=TRUE))
dadaRPPs <- dada_or_load(file.path(ckpt,"dadaRPPs.rds"),
              dada(derepRs, err=errR, multithread=50, pool=TRUE))

# Pseudo-pooling — disable SSE/vectorised alignment to avoid segfault seen at
# selfConsist step 21 on samples with >50k unique sequences. Results are identical.
setDadaOpt(SSE=0, VECTORIZED_ALIGNMENT=FALSE)
dadaFpsPPs <- dada_or_load(file.path(ckpt,"dadaFpsPPs.rds"),
               dada(derepFs, err=errF, multithread=50, pool="pseudo"))
dadaRpsPPs <- dada_or_load(file.path(ckpt,"dadaRpsPPs.rds"),
               dada(derepRs, err=errR, multithread=50, pool="pseudo"))
setDadaOpt(SSE=2, VECTORIZED_ALIGNMENT=TRUE)  # restore defaults

# ---- Merge paired ends ------------------------------------------------------
mergers      <- mergePairs(dadaFs,    filtFs, dadaRs,    filtRs, verbose=TRUE)
mergersPP    <- mergePairs(dadaFPPs,  filtFs, dadaRPPs,  filtRs, verbose=TRUE)
mergers_psPP <- mergePairs(dadaFpsPPs,filtFs, dadaRpsPPs,filtRs, verbose=TRUE)

# ---- Sequence tables --------------------------------------------------------
seqtab    <- makeSequenceTable(mergers)
seqtabPP  <- makeSequenceTable(mergersPP)
seqtabpsPP <- makeSequenceTable(mergers_psPP)

dim(seqtab); dim(seqtabPP); dim(seqtabpsPP)
table(nchar(getSequences(seqtab)))
table(nchar(getSequences(seqtabPP)))
table(nchar(getSequences(seqtabpsPP)))

# Save combined sequence table (both plates already included)
dir.create("/home/daniel/Ptarmigan/trimmed/mergedPlates", showWarnings = FALSE)
saveRDS(seqtab,    "/home/daniel/Ptarmigan/trimmed/mergedPlates/seqtab.rds")
saveRDS(seqtabPP,  "/home/daniel/Ptarmigan/trimmed/mergedPlates/seqtabPP.rds")
saveRDS(seqtabpsPP,"/home/daniel/Ptarmigan/trimmed/mergedPlates/seqtabpsPP.rds")

# =============================================================================
# SECTION 2 — CHIMERA REMOVAL
# =============================================================================

setwd("/home/daniel/Ptarmigan/trimmed/mergedPlates/")

# Both plates were processed together in Section 1 — no merging needed
seqtab    <- readRDS("seqtab.rds")
seqtabPP  <- readRDS("seqtabPP.rds")
seqtabpsPP <- readRDS("seqtabpsPP.rds")

merged_seqtab    <- seqtab
merged_seqtabPP  <- seqtabPP
merged_seqtabpsPP <- seqtabpsPP

# Remove chimeras
merged_seqtab.nochim    <- removeBimeraDenovo(merged_seqtab,    method="consensus", multithread=TRUE, verbose=TRUE)
merged_seqtabPP.nochim  <- removeBimeraDenovo(merged_seqtabPP,  method="consensus", multithread=TRUE, verbose=TRUE)
merged_seqtabpsPP.nochim <- removeBimeraDenovo(merged_seqtabpsPP,method="consensus", multithread=TRUE, verbose=TRUE)

dim(merged_seqtab.nochim)
dim(merged_seqtabPP.nochim)
dim(merged_seqtabpsPP.nochim)

# Sequence length summaries
summary(nchar(getSequences(merged_seqtab.nochim)))
summary(nchar(getSequences(merged_seqtabPP.nochim)))

# =============================================================================
# SECTION 3 — LULU CURATION
# Requires BLAST installed in PATH (or conda environment):
#   conda install -c bioconda blast
# =============================================================================

# Export FASTA files for BLAST
uniquesToFasta(merged_seqtab.nochim,    fout="merged_seqtab.nochim.fasta",    ids=paste0("OTU", seq(length(getSequences(merged_seqtab.nochim)))))
uniquesToFasta(merged_seqtabPP.nochim,  fout="merged_seqtabPP.nochim.fasta",  ids=paste0("OTU", seq(length(getSequences(merged_seqtabPP.nochim)))))
uniquesToFasta(merged_seqtabpsPP.nochim,fout="merged_seqtabpsPP.nochim.fasta",ids=paste0("OTU", seq(length(getSequences(merged_seqtabpsPP.nochim)))))

# Make LULU OTU tables (OTUs: rows, samples: columns)
npool.lulu <- merged_seqtab.nochim
colnames(npool.lulu) <- paste0("OTU", seq(length(getSequences(merged_seqtab.nochim))))
npool.lulu <- t(npool.lulu)

pool.lulu <- merged_seqtabPP.nochim
colnames(pool.lulu) <- paste0("OTU", seq(length(getSequences(merged_seqtabPP.nochim))))
pool.lulu <- t(pool.lulu)

pspool.lulu <- merged_seqtabpsPP.nochim
colnames(pspool.lulu) <- paste0("OTU", seq(length(getSequences(merged_seqtabpsPP.nochim))))
pspool.lulu <- t(pspool.lulu)

save.image(file = "ITS.RData")

# ---- BLAST (run these in the bash terminal, not in R) -----------------------
# makeblastdb -in merged_seqtab.nochim.fasta    -parse_seqids -dbtype nucl
# makeblastdb -in merged_seqtabPP.nochim.fasta  -parse_seqids -dbtype nucl
# makeblastdb -in merged_seqtabpsPP.nochim.fasta -parse_seqids -dbtype nucl
#
# blastn -db merged_seqtab.nochim.fasta    -outfmt '6 qseqid sseqid pident' -out NoPool_match_list.txt  -qcov_hsp_perc 80 -perc_identity 84 -query merged_seqtab.nochim.fasta
# blastn -db merged_seqtabPP.nochim.fasta  -outfmt '6 qseqid sseqid pident' -out Pool_match_list.txt    -qcov_hsp_perc 80 -perc_identity 84 -query merged_seqtabPP.nochim.fasta
# blastn -db merged_seqtabpsPP.nochim.fasta -outfmt '6 qseqid sseqid pident' -out psPool_match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query merged_seqtabpsPP.nochim.fasta

# ---- Run LULU algorithm in R ------------------------------------------------
setwd("/home/daniel/Ptarmigan/trimmed/mergedPlates/")

NoPool_match_list.txt  <- read.table("NoPool_match_list.txt")
Pool_match_list.txt    <- read.table("Pool_match_list.txt")
psPool_match_list.txt  <- read.table("psPool_match_list.txt")

nopool.nochim.curated_result <- lulu(as.data.frame(npool.lulu),  NoPool_match_list.txt)
pool.nochim.curated_result   <- lulu(as.data.frame(pool.lulu),   Pool_match_list.txt)
pspool.nochim.curated_result <- lulu(as.data.frame(pspool.lulu), psPool_match_list.txt)

# Summary of OTUs before/after LULU
print(paste0("Not Pooled: OTUs after Lulu: ",  nopool.nochim.curated_result$curated_count,
             " --- OTUs before Lulu: ",         nrow(nopool.nochim.curated_result$original_table)))
print(paste0("Pooled: OTUs after Lulu: ",       pool.nochim.curated_result$curated_count,
             " --- OTUs before Lulu: ",         nrow(pool.nochim.curated_result$original_table)))
print(paste0("PseudoPooled: OTUs after Lulu: ", pspool.nochim.curated_result$curated_count,
             " --- OTUs before Lulu: ",         nrow(pspool.nochim.curated_result$original_table)))

# Recover kept OTU indices (column numbers in original merged table)
nopool.kept.otus  <- as.numeric(gsub("OTU","", rownames(nopool.nochim.curated_result$curated_table),  perl=TRUE))
pool.kept.otus    <- as.numeric(gsub("OTU","", rownames(pool.nochim.curated_result$curated_table),    perl=TRUE))
pspool.kept.otus  <- as.numeric(gsub("OTU","", rownames(pspool.nochim.curated_result$curated_table),  perl=TRUE))

# Restore sequence names as column names
nopool.lulu  <- t(nopool.nochim.curated_result$curated_table)
colnames(nopool.lulu)  <- colnames(merged_seqtab.nochim[,  nopool.kept.otus])

pool.lulu    <- t(pool.nochim.curated_result$curated_table)
colnames(pool.lulu)    <- colnames(merged_seqtabPP.nochim[, pool.kept.otus])

pspool.lulu  <- t(pspool.nochim.curated_result$curated_table)
colnames(pspool.lulu)  <- colnames(merged_seqtabpsPP.nochim[,pspool.kept.otus])

# ---- Tidy up sample name formatting -----------------------------------------
name.change <- function(x) {
  rownames(x) <- gsub("_P4_r", "_P4r", rownames(x), perl=TRUE)
  rownames(x) <- gsub("_P3_r", "_P3r", rownames(x), perl=TRUE)
  rownames(x) <- gsub("_P2_r", "_P2r", rownames(x), perl=TRUE)
  rownames(x) <- gsub("_P1_r", "_P1r", rownames(x), perl=TRUE)
  return(x)
}

nopool.lulu  <- name.change(nopool.lulu)
pool.lulu    <- name.change(pool.lulu)
pspool.lulu  <- name.change(pspool.lulu)

print(rownames(nopool.lulu))
print(rownames(pool.lulu))
print(rownames(pspool.lulu))

# ---- Per-sample index summary -----------------------------------------------
index.info <- function(x){
  y <- data.frame(matrix(NA, nrow=nrow(x), ncol=0))
  y$rep     <- str_sub(rownames(x), start=-1)
  y$sample  <- str_sub(rownames(x), 1, nchar(rownames(x))-2)
  y$full    <- rownames(x)
  y$totseq  <- rowSums(x)
  y$OTUs    <- rowSums(x > 0)
  return(y)
}

nopool.lulu.index  <- index.info(nopool.lulu)
pool.lulu.index    <- index.info(pool.lulu)
pspool.lulu.index  <- index.info(pspool.lulu)

save.image(file="merged_ITS.RData")

# =============================================================================
# SECTION 4 — NEGATIVE CONTROL DECONTAMINATION
# Subtract maximum NTC / extraction blank reads per OTU
# =============================================================================

# IMPORTANT: prepare this metadata file before running Section 4.
# Required columns: sample, extract_blank, Fieldcontrol, sampletype.
# Source: /home/daniel/Ptarmigan/metadata/ITS_metadata_ptarmigan_clean.xlsx
# Verify column names match the functions below before proceeding.
exdata <- readxl::read_excel("/home/daniel/Ptarmigan/metadata/ITS_metadata_ptarmigan_clean.xlsx")
exdata <- as.data.frame(exdata)

# Align column names to what Section 4 functions expect
exdata <- dplyr::rename(exdata,
  sample        = SampleID,
  extract_blank = extraction_blank,
  sampletype    = SampleType
)
row.names(exdata) <- exdata$sample

# Extraction blanks (sampletype E) have NA in extract_blank in the metadata;
# assign them their own batch label so blank.change() can find them inside
# each split group (e.g. B13_P1_2A → extract_blank = "B13").
exdata$extract_blank[exdata$sampletype == "E"] <-
  sub("_.*", "", exdata$sample[exdata$sampletype == "E"])

# Data entry fix: S_1_8_P1_8F has extract_batch=13 (=B14) but extract_blank was omitted.
exdata$extract_blank[exdata$sample == "S_1_8_P1_8F"] <- "B14"

# ---- Subtract max NTC reads per OTU from all samples ----------------------
# NTC_ = no-template controls (primers + water); sets the per-OTU noise floor.
ntc.change <- function(x) {
  mind <- apply(x[grep("^NTC", rownames(x)), , drop = FALSE], 2,
                function(y) max(y, na.rm = TRUE))
  pmax(sweep(x, 2, mind), 0)
}

nopool.lulu.ntc <- ntc.change(nopool.lulu)
pool.lulu.ntc   <- ntc.change(pool.lulu)
pspool.lulu.ntc <- ntc.change(pspool.lulu)

# ---- 4 Remove PCR blanks (NTCs) ------------------------------------------
build_sample_info <- function(index_df) {
  sample.info <- as.data.frame(index_df)
  row.names(sample.info) <- sample.info$full
  sample.info
}

ntc.to.blankcontrol <-  function(ntc_mat, sample_index, exdata_df){
  z <- ntc_mat
  sample.info <- build_sample_info(sample_index)
  z = merge(z, sample.info, by =  'row.names', all.x=TRUE)
  z = merge(z, exdata_df, by.x="sample", by.y="sample", all.x=TRUE)
  rownames(z) <- z$Row.names
  z = split(z, z$extract_blank)
  dropnames <- colnames(z[[2]][, c(which(nchar(colnames(z[[2]]))< 30))])
  z <- lapply(z, function(x) x[!(names(x) %in% dropnames)])
  z = lapply(z, as.matrix)
  return(z)
}

nopool.lulu.ntc.blankC <- ntc.to.blankcontrol(nopool.lulu.ntc, nopool.lulu.index, exdata)
pool.lulu.ntc.blankC <- ntc.to.blankcontrol(pool.lulu.ntc, pool.lulu.index, exdata)
pspool.lulu.ntc.blankC <- ntc.to.blankcontrol(pspool.lulu.ntc, pspool.lulu.index, exdata)

## getting rid of experimental controls with no blank samples (NA extract_blank)
nopool.lulu.ntc.blankC1 <- nopool.lulu.ntc.blankC[!is.na(names(nopool.lulu.ntc.blankC))]
pool.lulu.ntc.blankC1   <- pool.lulu.ntc.blankC[!is.na(names(pool.lulu.ntc.blankC))]
pspool.lulu.ntc.blankC1 <- pspool.lulu.ntc.blankC[!is.na(names(pspool.lulu.ntc.blankC))]

unique(exdata$extract_blank) # copy in blank extract names
blank.change <- function(x){
  mind <- apply(x[grep("^B", rownames(x)), , drop = FALSE], 2, function(y) max(y, na.rm = TRUE))
  x1 <- sweep(x, 2, mind)
  x1 <- pmax(x1,0)
  return(x1)
}

nopool.lulu.ntc.blankC2 <- lapply(nopool.lulu.ntc.blankC1, blank.change)
pool.lulu.ntc.blankC2 <- lapply(pool.lulu.ntc.blankC1, blank.change)
pspool.lulu.ntc.blankC2 <- lapply(pspool.lulu.ntc.blankC1, blank.change)

## inf values because trying to subtract 0 from 0 .. converting back

nopool.lulu.ntc.blankC2 <- lapply(nopool.lulu.ntc.blankC2,function(x) replace(x, !is.finite(x), 0))
pool.lulu.ntc.blankC2 <- lapply(pool.lulu.ntc.blankC2, function(x) replace(x, !is.finite(x), 0))
pspool.lulu.ntc.blankC2 <- lapply(pspool.lulu.ntc.blankC2, function(x) replace(x, !is.finite(x), 0))

##########


# No field blanks in this study. Merge extraction-blank-subtracted groups back
# with do.call(rbind) — preserves sample rownames needed for downstream merge.
nopool.lulu.ntc.blankC2 <- lapply(nopool.lulu.ntc.blankC2, as.data.frame)
pool.lulu.ntc.blankC2   <- lapply(pool.lulu.ntc.blankC2,   as.data.frame)
pspool.lulu.ntc.blankC2 <- lapply(pspool.lulu.ntc.blankC2, as.data.frame)

# unname() prevents do.call(rbind) from prefixing rownames with the batch label
nopool.lulu.ntc.blank.fielddone <- do.call(rbind, unname(nopool.lulu.ntc.blankC2))
pool.lulu.ntc.blank.fielddone   <- do.call(rbind, unname(pool.lulu.ntc.blankC2))
pspool.lulu.ntc.blank.fielddone <- do.call(rbind, unname(pspool.lulu.ntc.blankC2))

# Remove extraction blank rows (B_ prefix, sampletype E) — they were needed
# inside each batch group for blank.change() but must not enter the biological
# analysis objects.
drop_blanks <- function(mat) mat[!grepl("^B", rownames(mat)), , drop = FALSE]
nopool.lulu.ntc.blank.fielddone <- drop_blanks(nopool.lulu.ntc.blank.fielddone)
pool.lulu.ntc.blank.fielddone   <- drop_blanks(pool.lulu.ntc.blank.fielddone)
pspool.lulu.ntc.blank.fielddone <- drop_blanks(pspool.lulu.ntc.blank.fielddone)

controlledblanks.to.samplelist <- function(x, sample_index) {
  minusNTCBLANK.sample <- x
  sample.info <- build_sample_info(sample_index)

  minusNTCBLANK.sample <- merge(minusNTCBLANK.sample, sample.info, by = 'row.names', all.x = TRUE)
  rownames(minusNTCBLANK.sample) <- minusNTCBLANK.sample$Row.names
  minusNTCBLANK.sample <- split(minusNTCBLANK.sample, minusNTCBLANK.sample$sample)

  # Drop columns with names shorter than 30 characters
  dropnames <- colnames(minusNTCBLANK.sample[[2]][, c(which(nchar(colnames(minusNTCBLANK.sample[[2]])) < 30))])
  minusNTCBLANK.sample <- lapply(minusNTCBLANK.sample, function(x) x[!(names(x) %in% dropnames)])

  return(minusNTCBLANK.sample)
}

nopool.lulu.controlled <- controlledblanks.to.samplelist(nopool.lulu.ntc.blank.fielddone, nopool.lulu.index)
pool.lulu.controlled   <- controlledblanks.to.samplelist(pool.lulu.ntc.blank.fielddone,   pool.lulu.index)
pspool.lulu.controlled <- controlledblanks.to.samplelist(pspool.lulu.ntc.blank.fielddone, pspool.lulu.index)

str(nopool.lulu.controlled)

cat("Controlled list sizes:\n")
cat("  nopool:", length(nopool.lulu.controlled), "\n")
cat("  pool  :", length(pool.lulu.controlled), "\n")
cat("  pspool:", length(pspool.lulu.controlled), "\n")

cat("OTU columns in first sample matrix:\n")
cat("  nopool:", ncol(nopool.lulu.controlled[[1]]), "\n")
cat("  pool  :", ncol(pool.lulu.controlled[[1]]), "\n")
cat("  pspool:", ncol(pspool.lulu.controlled[[1]]), "\n")

# =============================================================================
# SECTION 5 — PHYLOSEQ OBJECT CONSTRUCTION

# ---- PCR replicate filtering: rg2 (OTU must appear in both replicates) ------
rep.groups2 <- function(x){
  result <- apply(x, 2, function(c) replace(c, sum(c!=0)<2, 0))
  # apply() drops to a vector when x has only 1 row; restore matrix shape
  if (is.vector(result))
    matrix(result, nrow=1L, dimnames=list(rownames(x), names(result)))
  else result
}

rg2.nopool.lulu.controlled  <- lapply(nopool.lulu.controlled, rep.groups2)
rg2.pool.lulu.controlled    <- lapply(pool.lulu.controlled,   rep.groups2)
rg2.pspool.lulu.controlled  <- lapply(pspool.lulu.controlled, rep.groups2)

# ---- Build phyloseq objects -------------------------------------------------
library(dada2)
nopool.lulu.controlled  <- lapply(nopool.lulu.controlled,  as.matrix)
pool.lulu.controlled    <- lapply(pool.lulu.controlled,    as.matrix)
pspool.lulu.controlled  <- lapply(pspool.lulu.controlled,  as.matrix)

uniquesToFasta(as.matrix(nopool.lulu.controlled[[1]]),  fout="nopool.fasta",  ids=paste0("OTU", seq(length(getSequences(nopool.lulu.controlled[[1]])))))
uniquesToFasta(as.matrix(pool.lulu.controlled[[1]]),    fout="pool.fasta",    ids=paste0("OTU", seq(length(getSequences(pool.lulu.controlled[[1]])))))
uniquesToFasta(as.matrix(pspool.lulu.controlled[[1]]),  fout="pspool.fasta",  ids=paste0("OTU", seq(length(getSequences(pspool.lulu.controlled[[1]])))))

summary(nchar(getSequences("nopool.fasta")))

# Read FASTA files back in
nopoolseqs  <- read.fasta("nopool.fasta");  row.names(nopoolseqs)  <- nopoolseqs$id
poolseqs    <- read.fasta("pool.fasta");    row.names(poolseqs)    <- poolseqs$seq.name
pspoolseqs  <- read.fasta("pspool.fasta");  row.names(pspoolseqs)  <- pspoolseqs$seq.name

# Helper: collapse replicates into one matrix per sample
to.one.matrix <- function(x){
  lah <- do.call(rbind.data.frame, x)
  rownames(lah) <- names(x)
  colnames(lah) <- names(x[[1]])
  as.matrix(lah)
}

row.names(exdata) <- exdata$sample

# Helper: make a phyloseq object from OTU data, sample data, and sequences
make.phylo <- function(x, z, k){
  single  <- lapply(x, function(w) colSums(w))
  test    <- to.one.matrix(single)
  colnames(test) <- names(k)
  wanted  <- phyloseq(otu_table(test, taxa_are_rows=FALSE), sample_data(z))
  dna     <- Biostrings::DNAStringSet(colnames(x[[1]]))
  names(dna) <- names(k)
  merge_phyloseq(wanted, dna)
}

# ---- Raw data ----------------------------------------------------------------
nopoolps  <- make.phylo(nopool.lulu.controlled, exdata, nopoolseqs)
poolps    <- make.phylo(pool.lulu.controlled,   exdata, poolseqs)
pspoolps  <- make.phylo(pspool.lulu.controlled, exdata, pspoolseqs)

# ---- PCR-replicate-filtered (rg2: OTU in both replicates) -------------------
rg2.nopoolps <- make.phylo(rg2.nopool.lulu.controlled, exdata, nopoolseqs)
rg2.poolps   <- make.phylo(rg2.pool.lulu.controlled,   exdata, poolseqs)
rg2.pspoolps <- make.phylo(rg2.pspool.lulu.controlled, exdata, pspoolseqs)

# ---- Save phyloseq objects (without taxonomy - added in 3_taxonomic_assignment.R) ----
saveRDS(nopoolps,  "nopool_phyloseq.rds")
saveRDS(poolps,    "pool_phyloseq.rds")
saveRDS(pspoolps,  "pspool_phyloseq.rds")

saveRDS(rg2.nopoolps, "rg2.nopoolps.rds");  saveRDS(rg2.poolps, "rg2.poolps.rds");  saveRDS(rg2.pspoolps, "rg2.pspoolps.rds")

save.image("merged_ITS2.RData")
load("merged_ITS2.RData")

# =============================================================================
# SECTION 5b — FULL-COMPLEXITY PHYLOSEQ (pre-PCR-collapse)
#
# Purpose
# -------
# Build a phyloseq object where EVERY PCR replicate from EVERY root sample
# is stored as an independent sample row. This is the master input object for
# Monte Carlo resampling in Monte_Carlo.R — DO NOT collapse PCRs or roots here.
#
# 
# Naming convention for PCR replicate sample IDs (= otu_table row names)
# -----------------------------------------------------------------------
# Format: {sample_id}r{pcr_rep}
# Example: S_6_2_P1r1..r4
#
# The "plate" field (P1–P2) is the sequencing PLATE. There is no root/field
# replicate concept in this study (unlike the reference root study) — each
# individual x TimeBlock x Year has one dung sample, amplified via 2 PCR
# technical replicates (r1/r2, from Section 1's rep labeling).
#
# Metadata columns in the expanded PCR-level tables
# --------------------------------------------------
# All columns from exdata (ITS_metadata_ptarmigan_clean.xlsx, loaded in
# Section 4) are carried through plus:
#
# Auto-derived by expand_to_pcr_metadata():
#   pcr_sample_id  : exact PCR replicate sample ID (= otu_table rowname)
#                    format: paste0(sample, "r", pcr_rep_num)  e.g. "S_1_4_P1_1Br1"
#   root_id        : biological sample ID (= sample column; both PCRs share this)
#   pcr_rep_num    : PCR replicate number as integer (1/2)
#   pcr_rep_label  : PCR replicate as character ("r1"/"r2")
#
# Passed through from exdata:
#   sample         : dung sample ID (e.g. "S_1_4_P1_1B")
#   sampletype     : "S" = real dung sample (filtered to this below)
#   indivID, TimeBlock, Season, Year, sex, species, iNext_plate, ...
#
# Output CSV:  ITS_pcr_metadata.csv
# =============================================================================

# ---- Helper: stack all PCR replicate rows into one flat OTU matrix ----------
# pcr_list: named list (one element per root sample; rows = PCR reps, cols = OTU sequences)
# otu_ids : character vector of short OTU IDs (OTU1, OTU2, ...) matching the OTU seq columns
make_flat_otu_mat <- function(pcr_list, otu_ids) {
  flat <- do.call(rbind, lapply(pcr_list, as.matrix))
  storage.mode(flat) <- "numeric"
  flat[is.na(flat)] <- 0
  colnames(flat) <- otu_ids
  flat
}

# ---- Restrict to real dung samples (sampletype "S") before flattening -------
# PC_ (positive control) and T_ (tissue/other) samples were also sequenced and
# are still present as keys in *.lulu.controlled (only NTC/extraction blanks
# were dropped earlier), but the full-complexity object models real biological
# samples only (root_id/PCR replicate structure for GLLVM/HMSC).
dung_sample_ids <- exdata$sample[exdata$sampletype == "S"]

nopool.lulu.controlled.S <- nopool.lulu.controlled[names(nopool.lulu.controlled) %in% dung_sample_ids]
pool.lulu.controlled.S   <- pool.lulu.controlled[names(pool.lulu.controlled)     %in% dung_sample_ids]
pspool.lulu.controlled.S <- pspool.lulu.controlled[names(pspool.lulu.controlled) %in% dung_sample_ids]

# ---- OTU IDs (one per column-sequence, consistent with FASTA files above) ---
nopool_otu_ids  <- paste0("OTU", seq_len(ncol(nopool.lulu.controlled.S[[1]])))
pool_otu_ids    <- paste0("OTU", seq_len(ncol(pool.lulu.controlled.S[[1]])))
pspool_otu_ids  <- paste0("OTU", seq_len(ncol(pspool.lulu.controlled.S[[1]])))

# ---- Build flat OTU matrices ------------------------------------------------
fullmat_nopool  <- make_flat_otu_mat(nopool.lulu.controlled.S,  nopool_otu_ids)
fullmat_pool    <- make_flat_otu_mat(pool.lulu.controlled.S,    pool_otu_ids)
fullmat_pspool  <- make_flat_otu_mat(pspool.lulu.controlled.S,  pspool_otu_ids)

cat("Full-complexity matrix dimensions (PCR replicates × OTUs):\n")
cat("  nopool :", nrow(fullmat_nopool),  "PCR samples ×", ncol(fullmat_nopool),  "OTUs\n")
cat("  pool   :", nrow(fullmat_pool),    "PCR samples ×", ncol(fullmat_pool),    "OTUs\n")
cat("  pspool :", nrow(fullmat_pspool),  "PCR samples ×", ncol(fullmat_pspool),  "OTUs\n")

# ---- Reference sequences (DNAStringSet) — one per OTU -----------------------
dna_nopool  <- Biostrings::DNAStringSet(colnames(nopool.lulu.controlled.S[[1]]))
names(dna_nopool)  <- nopool_otu_ids

dna_pool    <- Biostrings::DNAStringSet(colnames(pool.lulu.controlled.S[[1]]))
names(dna_pool)    <- pool_otu_ids

dna_pspool  <- Biostrings::DNAStringSet(colnames(pspool.lulu.controlled.S[[1]]))
names(dna_pspool)  <- pspool_otu_ids

# ---- Build PCR-level metadata from exdata (Section 4) -----------------------
# exdata has sample IDs like: S_1_4_P1_1B
# We expand each biological sample row to 2 PCR rows:
#   S_1_4_P1_1Br1, S_1_4_P1_1Br2

expand_to_pcr_metadata <- function(meta_collapsed, reps = 1:2, sampletype_keep = "S") {
  required_cols <- c("sample", "sampletype")
  missing_cols <- setdiff(required_cols, colnames(meta_collapsed))
  if (length(missing_cols) > 0) {
    stop(sprintf("Collapsed metadata missing columns: %s", paste(missing_cols, collapse = ", ")))
  }

  base <- meta_collapsed %>%
    dplyr::filter(.data$sampletype %in% sampletype_keep) %>%
    dplyr::mutate(root_id = .data$sample)

  expanded <- base[rep(seq_len(nrow(base)), each = length(reps)), , drop = FALSE]
  expanded$pcr_rep_num <- rep(reps, times = nrow(base))
  expanded$pcr_rep_label <- paste0("r", expanded$pcr_rep_num)
  expanded$pcr_sample_id <- paste0(expanded$sample, expanded$pcr_rep_label)

  rownames(expanded) <- expanded$pcr_sample_id
  expanded
}

align_meta_to_otu <- function(meta_df, otu_mat) {
  missing <- setdiff(rownames(otu_mat), rownames(meta_df))
  if (length(missing) > 0) {
    stop(sprintf("Metadata missing %d PCR sample IDs from OTU matrix. Example: %s",
                 length(missing), paste(head(missing, 6), collapse = ", ")))
  }
  meta_df[rownames(otu_mat), , drop = FALSE]
}

meta_pcr <- expand_to_pcr_metadata(exdata, reps = 1:2, sampletype_keep = "S")

pcr_meta_nopool <- align_meta_to_otu(meta_pcr, fullmat_nopool)
pcr_meta_pool   <- align_meta_to_otu(meta_pcr, fullmat_pool)
pcr_meta_pspool <- align_meta_to_otu(meta_pcr, fullmat_pspool)

write.csv(pcr_meta_nopool, "ITS_pcr_metadata.csv", row.names = TRUE, na = "")
message("Wrote PCR-level metadata: ITS_pcr_metadata.csv")

# ---- Build full-complexity phyloseq (pre-PCR-collapse) ----------------------
make_fullcomplexity_ps <- function(flat_mat, meta_df, dna_set) {
  meta_ordered <- meta_df[rownames(flat_mat), , drop = FALSE]
  otu_ps  <- otu_table(flat_mat, taxa_are_rows = FALSE)
  samp_ps <- sample_data(meta_ordered)
  merge_phyloseq(phyloseq(otu_ps, samp_ps), dna_set)
}

fullps_nopool  <- make_fullcomplexity_ps(fullmat_nopool,  pcr_meta_nopool, dna_nopool)
fullps_pool    <- make_fullcomplexity_ps(fullmat_pool,    pcr_meta_pool,   dna_pool)
fullps_pspool  <- make_fullcomplexity_ps(fullmat_pspool,  pcr_meta_pspool, dna_pspool)

save.image("merged_ITS2.RData")

# Taxonomy is added in 3_taxonomic_assignment.R (same UNITE tax_table).
saveRDS(fullps_nopool,  "fullps_nopool.rds")
saveRDS(fullps_pool,    "fullps_pool.rds")
saveRDS(fullps_pspool,  "fullps_pspool.rds")

# Next step: 3_taxonomic_assignment.R
