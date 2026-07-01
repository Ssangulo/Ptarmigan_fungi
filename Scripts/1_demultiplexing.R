# =============================================================================
# 1_demultiplexing.R
# Demultiplexing and primer trimming of ITS amplicon sequencing data
#
# Tools required (install via conda):
#   - sabre  (barcode-based demultiplexing)
#   - cutadapt (primer/adapter trimming)
#   - dos2unix (line-ending conversion)
# =============================================================================

# ---- CONDA ENVIRONMENT SETUP ------------------------------------------------
# Make sure you have access to a unix system
# install conda (or miniconda etc)
# create an environment in conda, call it a name and install:

# cutadapt
conda install -c bioconda cutadapt
# sabre
conda install -c bioconda sabre
# dos2unix
conda install -c conda-forge dos2unix
# R
conda install -c r r

# ---- RAW DATA SETUP ---------------------------------------------------------

#create data folder in my data directory
mkdir data
cd data 
mkdir seqfiles

# Open the terminal and navigate to your sequencing raw_data folder. For me it is here 
cd /data/bigexpansion/michadm/seqdata/2024-06-10_Novogene_NovaSeq_250PE_Pernettya_fungi_metabarcoding_/X204SC24050751-Z01-F001/01.RawData/
  cd /data/bigexpansion/michadm/seqdata/2024-06-20_Novogene_NovaSeq_250PE_Pernettya_fungi_metabarcoding__batch2/X204SC24050751-Z01-F002/01.RawData/
  
  
  # Copy raw data folder to my data directory 
  cp -r /data/bigexpansion//seqdata/2024-06-10_Novogene_NovaSeq_250PE_Pernettya_fungi_metabarcoding_/X204SC24050751-Z01-F001/01.RawData /data/lastexpansion//data/
  
  #Create folders for each plate and split samples by plate
  
  cd /home/daniel/Ptarmigan/raw_data/01.RawData/
  mkdir -p /home/daniel/Ptarmigan/data/newrawdata/P1 /home/daniel/Ptarmigan/data/newrawdata/P2
cp -r *_P1_*/ /home/daniel/Ptarmigan/data/newrawdata/P1/
  cp -r *_P2_*/ /home/daniel/Ptarmigan/data/newrawdata/P2/
  
  # ---- FILE CONCATENATION (two sequencing lanes per sample) -------------------
# My data was delivered as a single Novogene batch, but each sample was run
# across two lanes (L1 on flowcell HJ5VTDRX7, L2 on flowcell HKWGTDRX7).
# Need to concatenate L1+L2 forward reads together, and L1+L2 reverse reads
# together, per sample. This block is for P1 - same logic repeated for P2.

cd /home/daniel/Ptarmigan/data/newrawdata/P1/
  mkdir -p seqfiles
cp **/*.fq.gz /home/daniel/Ptarmigan/data/newrawdata/P1/seqfiles

cd /home/daniel/Ptarmigan/data/newrawdata/P1/seqfiles

# Get a list of unique sample prefixes (up to the MKDL part of naming)
# THIS WILL WORK ONLY FOR P1 -> CHANGE SCRIPT FOR OTHER PLATES

for sample in $(ls *_1.fq.gz | sed -E 's/(.+_P1_[^_]+)_.*_L[12]_[12]\.fq\.gz/\1/' | sort | uniq);
do

# Concatenate forward reads across lanes for each sample
cat ${sample}*_L1_1.fq.gz ${sample}*_L2_1.fq.gz > ${sample}_concat_1.fq.gz

# Concatenate reverse reads across lanes for each sample
cat ${sample}*_L1_2.fq.gz ${sample}*_L2_2.fq.gz > ${sample}_concat_2.fq.gz

#Verify the concatenation
echo "Forward reads for $sample combined into ${sample}_concat_1.fq.gz"
echo "Reverse reads for $sample combined into ${sample}_concat_2.fq.gz"
done

# Concatenation for P2:
cd /home/daniel/Ptarmigan/data/newrawdata/P2/
  mkdir -p seqfiles
cp **/*.fq.gz /home/daniel/Ptarmigan/data/newrawdata/P2/seqfiles
cd /home/daniel/Ptarmigan/data/newrawdata/P2/seqfiles

for sample in $(ls *_1.fq.gz | sed -E 's/(.+_P2_[^_]+)_.*_L[12]_[12]\.fq\.gz/\1/' | sort | uniq);
do

cat ${sample}*_L1_1.fq.gz ${sample}*_L2_1.fq.gz > ${sample}_concat_1.fq.gz
cat ${sample}*_L1_2.fq.gz ${sample}*_L2_2.fq.gz > ${sample}_concat_2.fq.gz

echo "Forward reads for $sample combined into ${sample}_concat_1.fq.gz"
echo "Reverse reads for $sample combined into ${sample}_concat_2.fq.gz"

done

# Move concatenated files to their destination (P1 and P2 combined: identical
# rep tags across plates in this run, so both plates' demux/cutadapt steps
# run together from one folder)
mkdir -p /home/daniel/Ptarmigan/data/newrawdata/concatAll/
  mv /home/daniel/Ptarmigan/data/newrawdata/P1/seqfiles/*_concat_1.fq.gz
/home/daniel/Ptarmigan/data/newrawdata/concatAll/
  mv /home/daniel/Ptarmigan/data/newrawdata/P1/seqfiles/*_concat_2.fq.gz
/home/daniel/Ptarmigan/data/newrawdata/concatAll/
  
  mv /home/daniel/Ptarmigan/data/newrawdata/P2/seqfiles/*_concat_1.fq.gz
/home/daniel/Ptarmigan/data/newrawdata/concatAll/
  mv /home/daniel/Ptarmigan/data/newrawdata/P2/seqfiles/*_concat_2.fq.gz
/home/daniel/Ptarmigan/data/newrawdata/concatAll/
  
  
  # ---- DEMULTIPLEXING WITH SABRE ----------------------------------------------
# Based on the forward barcode read, designated in the barcode file
# (euka01_repbarcodes.txt). Only 2 PCR replicates were used in this run
# (file trimmed to rep1/rep2 only). P1 and P2 share identical tags this run,
# so both plates were demultiplexed together from the merged concatAll/ folder.

cd /home/daniel/Ptarmigan/raw_data/data/newrawdata/concatAll/
  
  source ~/miniforge3/bin/activate
conda activate r_env   # r_env has sabre, cutadapt, dos2unix installed

dos2unix /home/daniel/Ptarmigan/euka01_repbarcodes.txt

for i in *1.fq.gz; do
bn=${i/1.fq.gz}
sabre pe -f ${bn}1.fq.gz -r ${bn}2.fq.gz -b /home/daniel/Ptarmigan/euka01_repbarcodes.txt -u
${bn}unassigned1.fq -w ${bn}unassigned.fq
mv rep1f.fq ${bn}_rep1f.fq
mv rep1r.fq ${bn}_rep1r.fq
mv rep2f.fq ${bn}_rep2f.fq
mv rep2r.fq ${bn}_rep2r.fq
done

# ---- PRIMER TRIMMING WITH CUTADAPT -----------------------------------------
# Doing replicate by replicate; reverse primer made an exact match to the
# corresponding reverse barcode (not removed by sabre).
# rc = reverse complement
#
# Base primers:
#   ITS86F Fwd: GTGAATCATCGAATCTTTGAA
#   ITS4   Rv:  TCCTCCGCTTATTGATATGC
#
# IMPORTANT — orientation differs from the original reference pipeline:
# in this run, read 1 (rep#f.fq) starts with ITS86F (no barcode, sabre-stripped),
# and read 2 (rep#r.fq) starts with the reverse barcode + ITS4. This was
# confirmed by direct inspection of raw reads and by testing trim rates
# (0% trimmed with the reference's original adapter orientation vs ~97%
# with the orientation below). Re-verify this against raw reads for any
# future sequencing run before reusing these adapters.

# --- Replicate 1 (forward barcode GGTAC, reverse barcode AGGAA) ---
for i in *rep1f.fq; do
bn=${i/rep1f.fq}
cutadapt -a GTGAATCATCGAATCTTTGAA...GCATATCAATAAGCGGAGGATTCCT -A
AGGAATCCTCCGCTTATTGATATGC...TTCAAAGATTCGATGATTCACGTACC --untrimmed-output ${bn}.rep1out1.fq.gz
--untrimmed-paired-output ${bn}.rep1out2.fq.gz -o ${bn}.rep1.trim1.fq.gz -p ${bn}.rep1.trim2.fq.gz
${bn}rep1f.fq ${bn}rep1r.fq
done

# --- Replicate 2 (forward barcode CAACAC, reverse barcode GAGTGG) ---
for i in *rep2f.fq; do
bn=${i/rep2f.fq}
cutadapt -a GTGAATCATCGAATCTTTGAA...GCATATCAATAAGCGGAGGACCACTC -A
GAGTGGTCCTCCGCTTATTGATATGC...TTCAAAGATTCGATGATTCACGTGTTG --untrimmed-output ${bn}.rep2out1.fq.gz
--untrimmed-paired-output ${bn}.rep2out2.fq.gz -o ${bn}.rep2.trim1.fq.gz -p ${bn}.rep2.trim2.fq.gz
${bn}rep2f.fq ${bn}rep2r.fq
done

# ---- MOVE TRIMMED FILES TO DADA2 INPUT DIRECTORY ---------------------------
mv *trim1.fq.gz /home/daniel/Ptarmigan/trimmed/
mv *trim2.fq.gz /home/daniel/Ptarmigan/trimmed/
  
# Then proceed to script: 2_DADA2_lulu.R
