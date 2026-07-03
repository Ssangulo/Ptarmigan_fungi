# =============================================================================
# 5_community_composition.R
# Beta-diversity ordinations, community composition statistics, and
# taxonomic composition visualisations
#
# Requires objects from 4_data_prep.R:
#   alldat_full (list: nopool/pool/pspool) -- eco_analysis.RData
#
# Primary object: alldat_full[[1]] (== alldat_full$nopool), the full-
# complexity, per-PCR-replicate object (non-fungal OTUs removed, PlutoF
# taxonomy upgrade applied, taxonomy-cleaned; NOT depth-filtered/rarefied --
# see 4_data_prep.R header). Every analysis below needs one independent row
# per biological sample, so PCR replicates (2 per dung sample, sharing the
# same `sample` value) are summed into `ps_collapsed` ONCE in Section 1 and
# reused throughout. This is a fresh sum from alldat_full, not the
# pre-collapsed `alldat$nopool` from 4_data_prep.R -- note that unlike
# `alldat$nopool`, no min-depth filter has been applied here.
#
# This study has no site/habitat/elevation gradient (unlike the reference
# root study this pipeline derives from). The design is Season (winter/
# summer) x Year (2022-2024) -- see H1 in project notes. Reference-script
# sections built around habitat/site/elevation have been mapped onto
# Season/Year throughout; Section 5 (originally an elevation-gradient test)
# is adapted into a Season-centroid/Year-consistency check.
# =============================================================================

library(phyloseq)
library(microViz)
library(vegan)
library(ggplot2)
library(ggordiplots)
library(fantaxtic)
library(dplyr)
library(stringr)
library(patchwork)

setwd("/home/daniel/Ptarmigan/trimmed/mergedPlates/")
load("eco_analysis.RData")

out_dir  <- "models"
plot_dir <- "Plots2"
dir.create(out_dir,  showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE)

# =============================================================================
# SECTION 1 — COLLAPSE PCR REPLICATES (once) + QUICK METADATA CHECK
# =============================================================================

# ---- Collapse PCR replicates by summing --------------------------------------
# group_var = "sample": the dung-sample ID shared by both PCR replicates.
collapse_pcr_replicates <- function(ps, group_var = "sample") {
  otu <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) otu <- t(otu)           # -> samples (rows) x taxa (cols)

  md <- as(sample_data(ps), "data.frame")
  stopifnot(identical(rownames(otu), rownames(md)))
  grp <- factor(md[[group_var]])

  otu_collapsed <- rowsum(otu, group = grp)[levels(grp), , drop = FALSE]

  # PCR replicates of the same sample share identical metadata (Season,
  # TimeBlock, Year, indivID, ...) except the per-replicate pcr_* columns --
  # keep the first row per group and drop those replicate-specific columns.
  md_collapsed <- md[!duplicated(grp), , drop = FALSE]
  rownames(md_collapsed) <- md_collapsed[[group_var]]
  md_collapsed <- md_collapsed[levels(grp), , drop = FALSE]
  md_collapsed <- md_collapsed[, !names(md_collapsed) %in%
                                  c("pcr_sample_id", "pcr_rep_num", "pcr_rep_label"),
                                drop = FALSE]

  merge_phyloseq(
    otu_table(otu_collapsed, taxa_are_rows = FALSE),
    sample_data(md_collapsed),
    tax_table(ps),
    refseq(ps)
  )
}

ps_collapsed <- collapse_pcr_replicates(alldat_full[[1]], group_var = "sample")

cat("PCR-collapsed object (summed from alldat_full$nopool):\n")
cat("  ", nsamples(ps_collapsed), "samples,", ntaxa(ps_collapsed), "taxa\n")

k <- as.data.frame(sample_data(ps_collapsed))
head(k)
table(k$Season, k$Year)
table(k$Season, k$TimeBlock)

# =============================================================================
# SECTION 2 — PCA ON AITCHISON DISTANCES (robust CLR transform)
# =============================================================================

ps <- ps_collapsed

# decostand()/vegdist() need a plain matrix -- passing the phyloseq otu_table
# S4 object directly breaks internal subscripting (seen empirically: "logical
# subscript too long" inside vegan's rclr imputation / robust.aitchison).
otu_mat <- as(otu_table(ps), "matrix");  if (taxa_are_rows(ps)) otu_mat <- t(otu_mat)

rfyrg2 <- decostand(otu_mat, "rclr", MARGIN=1)
dimnames(rfyrg2) <- dimnames(otu_mat)  # vegan's rclr imputation path drops dimnames
my.rda <- rda(rfyrg2)
gdata  <- sample_data(ps)

p <- gg_ordiplot(my.rda, groups=gdata$Season, pt.size=3, spiders=TRUE, ellipse=FALSE)
p <- p$plot
ggsave(file.path(plot_dir, "PCA_Aitchison_Season.png"),
       plot=p, device="png", width=10, height=6, dpi=600)

# =============================================================================
# SECTION 3 — PERMANOVA + BETA-DISPERSION (Season x Year; H1)
# =============================================================================

ps <- ps_collapsed
sampledf <- data.frame(sample_data(ps))
otu_mat <- as(otu_table(ps), "matrix");  if (taxa_are_rows(ps)) otu_mat <- t(otu_mat)
dists <- vegdist(otu_mat, binary=FALSE, method="robust.aitchison")

# PERMANOVA: Season * Year (H1: seasonal turnover, reproducible across years)
adonis2(dists ~ Season * Year, by="terms", data=sampledf)

# Beta-dispersion (homogeneity of multivariate spread)
beta  <- betadisper(dists, factor(sampledf$Season))
beta1 <- betadisper(dists, factor(sampledf$Year))
permutest(beta)
permutest(beta1)

# =============================================================================
# SECTION 4 — CONSTRAINED ORDINATION (dbRDA; Season conditioned on Year)
# =============================================================================

ps      <- ps_collapsed
otu_mat <- as(otu_table(ps), "matrix");  if (taxa_are_rows(ps)) otu_mat <- t(otu_mat)
d_ait   <- vegdist(otu_mat, binary=FALSE, method="robust.aitchison")
sam     <- data.frame(sample_data(ps))
sam$Season <- factor(sam$Season)
sam$Year   <- factor(sam$Year)

## dbRDA model: Season, conditioned on Year
mod_ait <- capscale(d_ait ~ Season + Condition(Year), data=sam)

cat("\n== dbRDA overall ==\n");    print(anova.cca(mod_ait))
cat("\n== dbRDA by term ==\n");    print(anova.cca(mod_ait, by="terms"))
cat("\n== dbRDA by axis ==\n");    print(anova.cca(mod_ait, by="axis"))

eig_con  <- mod_ait$CCA$eig
cap1_pct <- round(100 * eig_con[1] / sum(eig_con), 1)
cap2_pct <- round(100 * eig_con[2] / sum(eig_con), 1)

# Publication-quality dbRDA plot
cols <- c("winter"="#1f77b4", "summer"="#d62728")
pchv <- c("winter"=16, "summer"=17)

png(file.path(plot_dir, "dbRDA_collapsed.png"),
    width=7.5, height=6.5, units="in", res=600)
par(mfrow=c(1,1), mar=c(4,4,1,1), xpd=NA, cex=1)

plot(mod_ait, display=c("sites","cn"), type="n", xaxt="n", yaxt="n", xlab="", ylab="")
with(sam, points(scores(mod_ait, display="sites"),
                 col=cols[Season], pch=pchv[Season], cex=0.6))
with(sam, ordiellipse(mod_ait, Season, kind="sd", level=0.95, draw="polygon",
                      col=adjustcolor(unname(cols[levels(Season)]), alpha.f=0.2),
                      border=unname(cols[levels(Season)]), lwd=1.2))

cen    <- with(sam, ordiellipse(mod_ait, Season, kind="sd", level=0.95, plot=FALSE))
cen_xy <- do.call(rbind, lapply(cen, `[[`, "centroid"))
points(cen_xy, pch=4, cex=1.2, lwd=1.4, col=unname(cols[rownames(cen_xy)]))

legend("topright", inset=0.01, bty="n",
       legend=c("Winter","Summer"),
       col=cols, pch=c(16,17), pt.cex=0.9, cex=0.9)
mtext(paste0("dbRDA1 (", cap1_pct, "%)"), side=1, line=2.5)
mtext(paste0("dbRDA2 (", cap2_pct, "%)"), side=2, line=2.5)
dev.off()

# =============================================================================
# SECTION 5 — DISTANCE-TO-SEASON-CENTROID TEST (consistency across years)
# Reference script tested whether beta-diversity escalated with elevation
# consistently across sites (parallel-slopes test). No elevation gradient
# exists here, so this is adapted to test whether a sample's dispersion from
# its Season centroid (in rCLR/Aitchison space) differs by Year -- a check
# on H1's "reproducible across years" claim.
# =============================================================================

ps   <- ps_collapsed
X    <- as(otu_table(ps), "matrix");  if (taxa_are_rows(ps)) X <- t(X)
rfyrg2 <- decostand(X, "rclr", MARGIN=1)
dimnames(rfyrg2) <- dimnames(X)  # vegan's rclr imputation path drops dimnames
meta <- data.frame(sample_data(ps), stringsAsFactors=FALSE, check.names=FALSE)
stopifnot(identical(rownames(rfyrg2), rownames(meta)))

# Season centroids in rCLR space
centroids <- rowsum(rfyrg2, group=meta$Season) / as.vector(table(meta$Season))

euclid <- function(x, y) sqrt(sum((x-y)^2))
dist_centroid <- vapply(
  seq_len(nrow(rfyrg2)),
  function(k) euclid(rfyrg2[k,], centroids[meta$Season[k],]),
  numeric(1)
)

meta$dist_centroid <- dist_centroid
sample_data(ps)$dist_centroid <- dist_centroid

df <- meta
df$Season <- factor(df$Season)
df$Year   <- factor(df$Year)

# Test: does distance-to-Season-centroid differ by Year?
mod <- lm(dist_centroid ~ Season * Year, data=df)
anova(mod); summary(mod)

# =============================================================================
# SECTION 6 — PCoA ON BRAY-CURTIS DISSIMILARITIES
# =============================================================================

ps <- ps_collapsed

mat     <- as(otu_table(ps), "matrix");  if (taxa_are_rows(ps)) mat <- t(mat)
mat_rel <- vegan::decostand(mat, method="total", MARGIN=1)
bc      <- vegan::vegdist(mat_rel, method="bray")
ord     <- vegan::capscale(bc ~ 1)

eig  <- pmax(ord$CA$eig, 0)
ax1p <- round(100 * eig[1] / sum(eig), 1)
ax2p <- round(100 * eig[2] / sum(eig), 1)

sc <- as.data.frame(vegan::scores(ord, display="sites", choices=1:2, scaling=1))
colnames(sc)[1:2] <- c("Axis1","Axis2")
sc$sample <- rownames(sc)

meta_p <- as.data.frame(sample_data(ps));  meta_p$sample <- rownames(meta_p)
sc$Season <- meta_p$Season[match(sc$sample, meta_p$sample)]
sc$Year   <- meta_p$Year[match(sc$sample, meta_p$sample)]
sc$Season <- droplevels(factor(sc$Season))
sc$Year   <- droplevels(factor(sc$Year))

cent <- aggregate(cbind(Axis1,Axis2) ~ Season, data=sc, FUN=mean)
names(cent) <- c("Season","cX","cY")
sc <- merge(sc, cent, by="Season", all.x=TRUE)

gp <- ggplot(sc, aes(Axis1, Axis2, colour=Season, shape=Year)) +
  geom_segment(aes(xend=cX, yend=cY, group=interaction(Season,sample)),
               alpha=0.5, linewidth=0.3, show.legend=FALSE) +
  geom_point(size=3) +
  geom_point(data=cent, aes(x=cX,y=cY,colour=Season), inherit.aes=FALSE,
             size=3, shape=4, stroke=1) +
  coord_equal() +
  xlab(paste0("PCoA1 (",ax1p,"%)")) +
  ylab(paste0("PCoA2 (",ax2p,"%)")) +
  guides(colour=guide_legend(title="Season",order=1),
         shape =guide_legend(title="Year",   order=2)) +
  theme_classic(base_size=12)

ggsave(file.path(plot_dir, "PCoA_BrayCurtis_Season.png"),
       plot=gp, width=12, height=8, dpi=800)

# =============================================================================
# SECTION 7 — TAXONOMIC COMPOSITION BAR PLOTS (microViz)
# =============================================================================

filter_taxa_by_rank <- function(physeq_obj, rank_prefix="o__", exclude_term="Incertae_sedis") {
  taxa_assigned <- grepl(rank_prefix, tax_table(physeq_obj)[,"Order"])
  exclude_inc   <- !grepl(exclude_term, tax_table(physeq_obj)[,"Order"])
  prune_taxa(taxa_assigned & exclude_inc, physeq_obj)
}

ps_filtered <- filter_taxa_by_rank(ps_collapsed)

# comp_barplot()/tax_agg() need every taxon to have a non-NA, unique label at
# the plotted rank -- tax_fix() fills NA ranks from the nearest classified
# parent for display purposes only, scoped to this local plotting copy. This
# does NOT touch ps_collapsed/alldat's tax_table (script 4 deliberately keeps
# real NAs there -- see 4_data_prep.R Section 7 comment).
ps_filtered <- ps_filtered %>% tax_fix()

# Clean rank prefixes (o__, g__, etc.)
tax_table(ps_filtered) <- apply(tax_table(ps_filtered), 2,
                                function(x) gsub("^[a-z]__","", x))

# Recode Season labels
sample_data(ps_filtered)$Season <- dplyr::recode(
  as.character(sample_data(ps_filtered)$Season),
  "winter"="Winter", "summer"="Summer",
  .default=NA_character_)
sample_data(ps_filtered)$Season <- factor(
  sample_data(ps_filtered)$Season, levels=c("Winter","Summer"), ordered=TRUE)

ps_filtered <- ps_filtered %>% ps_arrange(Year)

# Seriate samples within each Year by Order composition
year_levels <- sample_data(ps_filtered) %>% as.data.frame() %>% pull(Year) %>% unique()
samp_order  <- unlist(lapply(year_levels, function(y) {
  ps_filtered %>% ps_filter(Year==y) %>% ps_seriate(rank="Order") %>% sample_names()
}))

# ---- Phylum bar plot --------------------------------------------------------
phy <- ps_filtered %>%
  comp_barplot(tax_level="Phylum", n_taxa=5, label="Year",
               bar_outline_colour="grey5", facet_by="Season",
               sample_order=samp_order, merge_other=FALSE, other_name="Other phyla") +
  coord_flip() +
  theme(legend.text=element_text(size=6), legend.title=element_text(size=7),
        legend.key.size=unit(0.3,"cm"), legend.spacing.x=unit(0.2,"cm"))

ggsave(file.path(plot_dir, "phylum_barplot_ordered_season.png"),
       plot=phy, device="png", width=12, height=8, dpi=800)

# ---- Order bar plot ---------------------------------------------------------
p <- ps_filtered %>%
  comp_barplot(tax_level="Order", n_taxa=10, label="Year",
               bar_outline_colour="grey5", facet_by="Season",
               sample_order=samp_order, merge_other=FALSE, other_name="Other orders") +
  coord_flip() +
  theme(legend.text=element_text(size=6), legend.title=element_text(size=7),
        legend.key.size=unit(0.3,"cm"), legend.spacing.x=unit(0.2,"cm"))

ggsave(file.path(plot_dir, "Order_barplot_ordered_season.png"),
       plot=p, device="png", width=12, height=8, dpi=800)

# ---- Class bar plot ---------------------------------------------------------
class <- ps_filtered %>%
  comp_barplot(tax_level="Class", n_taxa=10, label="Year",
               bar_outline_colour="grey5", facet_by="Season",
               sample_order=samp_order, merge_other=FALSE, other_name="Other classes") +
  coord_flip() +
  theme(legend.text=element_text(size=6), legend.title=element_text(size=7),
        legend.key.size=unit(0.3,"cm"), legend.spacing.x=unit(0.2,"cm"))

ggsave(file.path(plot_dir, "class_barplot_ordered_season.png"),
       plot=class, device="png", width=12, height=8, dpi=800)

# ---- Genus bar plot ---------------------------------------------------------
genus <- ps_filtered %>%
  comp_barplot(tax_level="Genus", n_taxa=20, label="Year",
               bar_outline_colour="grey5", facet_by="Season",
               sample_order=samp_order, merge_other=FALSE, other_name="Other genera") +
  coord_flip() +
  theme(legend.text=element_text(size=6), legend.title=element_text(size=7),
        legend.key.size=unit(0.3,"cm"), legend.spacing.x=unit(0.2,"cm"))

ggsave(file.path(plot_dir, "genus_barplot_ordered_season.png"),
       plot=genus, device="png", width=12, height=8, dpi=800)

# ---- Panel plot (Phylum + Genus + Order) ------------------------------------
axis_tweak  <- theme(axis.text.x=element_text(size=6), axis.text.y=element_text(size=5),
                     axis.title=element_text(size=8))
legend_tweak <- theme(legend.position="right", legend.key.size=unit(0.5,"cm"),
                      legend.text=element_text(size=11), legend.title=element_text(size=11,face="bold"))

phy_lab   <- phy   + ggtitle("A) Phylum") + axis_tweak + legend_tweak +
  guides(fill=guide_legend(ncol=1)) + theme(plot.title=element_text(face="bold",size=12,hjust=0))
genus_lab <- genus + ggtitle("B) Genus")  + axis_tweak + legend_tweak +
  guides(fill=guide_legend(ncol=1)) + theme(plot.title=element_text(face="bold",size=12,hjust=0))
p_lab     <- p     + ggtitle("C) Order")  + legend_tweak +
  guides(fill=guide_legend(ncol=1)) +
  theme(plot.title=element_text(face="bold",size=12,hjust=0),
        axis.text.x=element_text(size=7), axis.text.y=element_text(size=6))

panel <- phy_lab + genus_lab + p_lab +
  plot_layout(design="AC\nBC", widths=c(1.3,1.65), heights=c(1,1), guides="keep")

ggsave(file.path(plot_dir, "Panel_PhylumGenus_Order_Season.png"),
       plot=panel, width=15, height=10, dpi=800)

# =============================================================================
# SECTION 8 — TAXONOMIC ASSIGNMENT SUMMARY (% reads/OTUs assigned per Season)
# =============================================================================

ps <- ps_collapsed

otu_mat <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)

tax <- as.data.frame(tax_table(ps));  tax$OTU_ID <- rownames(tax)
sam <- as.data.frame(sample_data(ps)); sam$SampleID <- rownames(sam)
sam$Season <- factor(sam$Season, levels=c("winter","summer"),
                      labels=c("Winter","Summer"))

tax <- tax %>% mutate(
  assigned_phylum = str_detect(Phylum,"^p__") & !str_detect(Phylum,regex("fungi kingdom",ignore_case=TRUE)),
  assigned_class  = str_detect(Class, "^c__") & !str_detect(Class, regex("fungi kingdom",ignore_case=TRUE)),
  assigned_order  = str_detect(Order, "^o__") & !str_detect(Order, regex("incertae",ignore_case=TRUE)),
  assigned_genus  = str_detect(Genus, "^g__") & !str_detect(Genus, regex("incertae",ignore_case=TRUE))
)

summarise_rank <- function(rank_flag) {
  # which() drops NA (unresolved ranks) instead of propagating them into the
  # subset, which would otherwise error on otu_mat[assigned_otus, ...] below.
  assigned_otus <- tax$OTU_ID[which(tax[[rank_flag]])]
  reads_df <- lapply(levels(sam$Season), function(s) {
    samp <- sam$SampleID[sam$Season == s]
    tibble(Season=s,
           total_reads   = sum(otu_mat[, samp, drop=FALSE]),
           reads_assigned = sum(otu_mat[assigned_otus, samp, drop=FALSE]))
  }) |> bind_rows()
  otus_df <- lapply(levels(sam$Season), function(s) {
    samp <- sam$SampleID[sam$Season == s]
    present          <- rowSums(otu_mat[, samp, drop=FALSE]) > 0
    present_assigned <- rowSums(otu_mat[assigned_otus, samp, drop=FALSE]) > 0
    tibble(Season=s, total_otus=sum(present), otus_assigned=sum(present_assigned))
  }) |> bind_rows()
  reads_df %>% left_join(otus_df, by="Season") %>%
    mutate(pct_reads_assigned = 100 * reads_assigned / total_reads,
           pct_otus_assigned  = 100 * otus_assigned  / total_otus)
}

tax_res_summary <- bind_rows(
  summarise_rank("assigned_phylum") %>% mutate(rank="Phylum"),
  summarise_rank("assigned_class")  %>% mutate(rank="Class"),
  summarise_rank("assigned_order")  %>% mutate(rank="Order"),
  summarise_rank("assigned_genus")  %>% mutate(rank="Genus")
) %>%
  mutate(rank=factor(rank, levels=c("Phylum","Class","Order","Genus")),
         Season=factor(Season, levels=c("Winter","Summer"))) %>%
  arrange(rank, Season) %>%
  select(rank, Season, total_reads, reads_assigned, pct_reads_assigned,
         total_otus, otus_assigned, pct_otus_assigned)

print(tax_res_summary, n=Inf, width=Inf)

writexl::write_xlsx(tax_res_summary, file.path(out_dir, "tax_assignment_summary_season.xlsx"))
