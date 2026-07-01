# =============================================================================
# 3_taxonomic_assignment.R
# Fungal ITS taxonomy assignment using DADA2 + UNITE database, and
# integration of taxonomy tables into phyloseq objects produced in
# 2_DADA2_lulu.R
#
# Input:  phyloseq .rds and FASTA files from 2_DADA2_lulu.R
#         UNITE general FASTA release (sh_general_release_dynamic_s_19.02.2025)
# Output: taxonomy-annotated phyloseq objects (.rds)
#
# Reference: 3_taxonomic_assignment.R @ 65fbffa (Root_fungi_DADA2)
# =============================================================================

library(dada2)
library(Biostrings)
library(phyloseq)

setwd("/home/daniel/Ptarmigan/trimmed/mergedPlates/")

# =============================================================================
# SECTION 1 — ASSIGN TAXONOMY WITH UNITE DATABASE
# =============================================================================

# Path to UNITE general release FASTA
unite <- "/home/daniel/Ptarmigan/reference_db/sh_general_release_dynamic_s_19.02.2025.fasta"

# Read OTU sequences for each pooling strategy
nopool_fasta  <- readDNAStringSet("nopool.fasta")
pool_fasta    <- readDNAStringSet("pool.fasta")
pspool_fasta  <- readDNAStringSet("pspool.fasta")

summary(nchar(getSequences(nopool_fasta)))

# ---- Assign with minBoot = 80 (recommended) ---------------------------------
taxa_nopool_80  <- assignTaxonomy(nopool_fasta,  unite, multithread=50, tryRC=TRUE, minBoot=80)
taxa_pool_80    <- assignTaxonomy(pool_fasta,    unite, multithread=50, tryRC=TRUE, minBoot=80)
taxa_pspool_80  <- assignTaxonomy(pspool_fasta,  unite, multithread=50, tryRC=TRUE, minBoot=80)

# ---- Assign with minBoot = 60 (permissive alternative) ---------------------
taxa_nopool_60  <- assignTaxonomy(nopool_fasta,  unite, multithread=50, tryRC=TRUE, minBoot=60)
taxa_pool_60    <- assignTaxonomy(pool_fasta,    unite, multithread=50, tryRC=TRUE, minBoot=60)
taxa_pspool_60  <- assignTaxonomy(pspool_fasta,  unite, multithread=50, tryRC=TRUE, minBoot=60)

save.image(file="dada2_taxa.RData")
# load("dada2_taxa.RData")

# ---- Assignment success summary ---------------------------------------------
taxa_summary <- function(taxa) {
  sapply(1:ncol(taxa), function(i) sum(!is.na(taxa[, i])))
}

summary_df <- data.frame(
  Rank       = colnames(taxa_nopool_80),
  Nopool_80  = taxa_summary(taxa_nopool_80),
  Pool_80    = taxa_summary(taxa_pool_80),
  PsPool_80  = taxa_summary(taxa_pspool_80),
  Nopool_60  = taxa_summary(taxa_nopool_60),
  Pool_60    = taxa_summary(taxa_pool_60),
  PsPool_60  = taxa_summary(taxa_pspool_60)
)
print(summary_df)

# =============================================================================
# SECTION 2 — INTEGRATE TAXONOMY INTO PHYLOSEQ OBJECTS
# =============================================================================

# ---- Load phyloseq objects (from 2_DADA2_lulu.R) ----------------------------
nopoolps  <- readRDS("nopool_phyloseq.rds")
poolps    <- readRDS("pool_phyloseq.rds")
pspoolps  <- readRDS("pspool_phyloseq.rds")

rg2.nopoolps <- readRDS("rg2.nopoolps.rds");  rg2.poolps <- readRDS("rg2.poolps.rds");  rg2.pspoolps <- readRDS("rg2.pspoolps.rds")

# ---- Convert taxonomy to phyloseq format ------------------------------------
tax_table_nopoolps  <- tax_table(taxa_nopool_80)
tax_table_poolps    <- tax_table(taxa_pool_80)
tax_table_pspoolps  <- tax_table(taxa_pspool_80)

# Sanity: taxa IDs are sequences; need to map to OTU IDs
head(taxa_names(nopoolps), 2)
head(taxa_names(tax_table_nopoolps), 2)

# Create mapping between sequences and OTU identifiers
seqs_tax_nopool  <- rownames(tax_table_nopoolps)
seqs_tax_pool    <- rownames(tax_table_poolps)
seqs_tax_pspool  <- rownames(tax_table_pspoolps)

seq_to_otu_nopool  <- setNames(taxa_names(nopoolps),  seqs_tax_nopool)
seq_to_otu_pool    <- setNames(taxa_names(poolps),    seqs_tax_pool)
seq_to_otu_pspool  <- setNames(taxa_names(pspoolps),  seqs_tax_pspool)

rownames(tax_table_nopoolps)  <- seq_to_otu_nopool[rownames(tax_table_nopoolps)]
rownames(tax_table_poolps)    <- seq_to_otu_pool[rownames(tax_table_poolps)]
rownames(tax_table_pspoolps)  <- seq_to_otu_pspool[rownames(tax_table_pspoolps)]

# Confirm alignment
all(rownames(tax_table_nopoolps) == taxa_names(nopoolps))
all(rownames(tax_table_poolps)   == taxa_names(poolps))
all(rownames(tax_table_pspoolps) == taxa_names(pspoolps))

# ---- Merge taxonomy into raw phyloseq objects -------------------------------
nopoolps.dada2  <- merge_phyloseq(nopoolps,  tax_table_nopoolps)
poolps.dada2    <- merge_phyloseq(poolps,    tax_table_poolps)
pspoolps.dada2  <- merge_phyloseq(pspoolps,  tax_table_pspoolps)

saveRDS(nopoolps.dada2, "nopoolps.dada2.rds")
saveRDS(poolps.dada2,   "poolps.dada2.rds")
saveRDS(pspoolps.dada2, "pspoolps.dada2.rds")

# ---- Merge taxonomy into rg2-filtered phyloseq objects ----------------------
rg2.nopoolps.tax  <- rg2.nopoolps;  rg2.poolps.tax  <- rg2.poolps;  rg2.pspoolps.tax  <- rg2.pspoolps

tax_table(rg2.nopoolps.tax)  <- tax_table_nopoolps;  tax_table(rg2.poolps.tax)  <- tax_table_poolps;  tax_table(rg2.pspoolps.tax)  <- tax_table_pspoolps

# Quick check
head(tax_table(rg2.nopoolps.tax))
ntaxa(nopoolps) == ntaxa(nopoolps.dada2)

# Save taxonomy-annotated objects
saveRDS(rg2.nopoolps.tax, "rg2.nopoolps.tax.rds");  saveRDS(rg2.poolps.tax, "rg2.poolps.tax.rds");  saveRDS(rg2.pspoolps.tax, "rg2.pspoolps.tax.rds")

# =============================================================================
# SECTION 3b — ATTACH TAXONOMY TO FULL-COMPLEXITY PHYLOSEQ OBJECTS
#
# The full-complexity phyloseq objects (fullps_nopool/pool/pspool) were built
# in 2_DADA2_lulu.R Section 5b with rows = individual PCR replicates.
# Here we attach the SAME tax_table from UNITE to these objects to ensure
# full consistency with the main study.
# =============================================================================

fullps_nopool  <- readRDS("fullps_nopool.rds")
fullps_pool    <- readRDS("fullps_pool.rds")
fullps_pspool  <- readRDS("fullps_pspool.rds")

# Confirm OTU ID alignment before attaching
stopifnot(all(taxa_names(fullps_nopool)  %in% rownames(tax_table_nopoolps)))
stopifnot(all(taxa_names(fullps_pool)    %in% rownames(tax_table_poolps)))
stopifnot(all(taxa_names(fullps_pspool)  %in% rownames(tax_table_pspoolps)))

tax_table(fullps_nopool)  <- tax_table_nopoolps[taxa_names(fullps_nopool),  ]
tax_table(fullps_pool)    <- tax_table_poolps[taxa_names(fullps_pool),      ]
tax_table(fullps_pspool)  <- tax_table_pspoolps[taxa_names(fullps_pspool),  ]

# Quick check
cat("Full-complexity phyloseq (nopool): ",
    nsamples(fullps_nopool), "PCR samples, ",
    ntaxa(fullps_nopool),    "OTUs, taxonomy attached:", !is.null(tax_table(fullps_nopool, errorIfNULL=FALSE)), "\n")

saveRDS(fullps_nopool,  "fullps_nopool.rds")   # overwrites version without taxonomy
saveRDS(fullps_pool,    "fullps_pool.rds")
saveRDS(fullps_pspool,  "fullps_pspool.rds")

# Next step: 4_data_prep.R
