##############
# May, 2020
# Zach Burcham, QIIME2 analysis on "A Pilot Study Characterizing Gravesoil Microbial Communities a Decade After Swine Decomposition"
# Samples taken from from soil cores

#####
# Sequences are uploaded to Qiita study 12465 data type 71318
# Sample could not be processed in Qiita so the individual FASTQ files 
# were imported to Qiime2 for analysis
#####


####
# QIIME2 2020.2 processing
####

### Data import from Qiita study 12465 data type 71318
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path metadata/DPSW_manifest.csv \
  --output-path seqs/paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33

# Visualize imported seqs
# 72 total samples
# Min seqs 34307, median 81148, mean 83535.4, max 135064, total 6014549
qiime demux summarize \
  --i-data seqs/paired-end-demux.qza \
  --o-visualization seqs/paired-end-demux.qzv
  
### Join paired end sequences
qiime vsearch join-pairs \
  --i-demultiplexed-seqs seqs/paired-end-demux.qza \
  --o-joined-sequences seqs/joined-seqs.qza \
  --verbose
  
# Visualize joined seqs
# 72 samples
# Min seqs 30637, median 73865, mean 74715.8, max 115792, total 5379539
qiime demux summarize \
  --i-data seqs/joined-seqs.qza \
  --o-visualization seqs/joined-seqs.qzv

### Deblur denoise
qiime deblur denoise-16S \
  --i-demultiplexed-seqs seqs/joined-seqs.qza \
  --p-trim-length 250 \
  --o-representative-sequences seqs/rep-seqs.qza \
  --o-table feature_tables/table.qza \
  --p-sample-stats \
  --o-stats deblur/deblur-stats.qza
  
# Visualize deblur stats
qiime deblur visualize-stats \
  --i-deblur-stats deblur/deblur-stats.qza \
  --o-visualization deblur/deblur-stats.qzv

### SEPP fragment insertion to obtain phylogenetic tree with gg13.8
wget \
  -O "sepp-refs-gg-13-8.qza" \
  "https://data.qiime2.org/2020.2/common/sepp-refs-gg-13-8.qza"

mv sepp-refs-gg-13-8.qza trees/

qiime fragment-insertion sepp \
  --i-representative-sequences seqs/rep-seqs.qza \
  --i-reference-database trees/sepp-refs-gg-13-8.qza \
  --p-threads 4 \
  --output-dir trees/fragment_insertion_out

### Feature table and data summaries
# Table
qiime feature-table summarize \
  --i-table feature_tables/table.qza \
  --o-visualization feature_tables/table.qzv \
  --m-sample-metadata-file metadata/DPSW_qiime2_metadata.txt
  
# Rep seqs
qiime feature-table tabulate-seqs \
  --i-data seqs/rep-seqs.qza \
  --o-visualization seqs/rep-seqs.qzv
   
### Get greengenes v4 classifier and move to appropriate directory
wget \
  -O "gg-13-8-99-515-806-nb-classifier.qza" \
  "https://data.qiime2.org/2020.2/common/gg-13-8-99-515-806-nb-classifier.qza"

mv gg-13-8-99-515-806-nb-classifier.qza taxonomy/

### GreenGenes 99% Classification (V4)
qiime feature-classifier classify-sklearn \
  --i-classifier taxonomy/gg-13-8-99-515-806-nb-classifier.qza \
  --i-reads seqs/rep-seqs.qza \
  --p-n-jobs -1 \
  --o-classification taxonomy/taxonomy-gg-99.qza

qiime metadata tabulate \
  --m-input-file taxonomy/taxonomy-gg-99.qza \
  --o-visualization taxonomy/taxonomy-gg-99.qzv

### Filter feature table of mitochondria and chloroplast and keep only bacteria
qiime taxa filter-table \
 --i-table feature_tables/table.qza \
 --i-taxonomy taxonomy/taxonomy-gg-99.qza \
 --p-include bacteria \
 --p-exclude mitochondria,chloroplast \
 --o-filtered-table feature_tables/table-gg-99-no-chlo-mito.qza

### Filter features based on frequency to remove features less than 20 times present
qiime feature-table filter-features \
 --i-table feature_tables/table-gg-99-no-chlo-mito.qza \
 --p-min-frequency 20 \
 --o-filtered-table feature_tables/table-gg-99-no-chlo-mito-minfreq20.qza

# Create taxa barplot
qiime taxa barplot \
  --i-table feature_tables/table-gg-99-no-chlo-mito-minfreq20.qza \
  --i-taxonomy taxonomy/taxonomy-gg-99.qza \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --o-visualization taxonomy/taxa-bar-plots-gg-99.qzv

## Removing sample that belongs to location B (0 meter north) since there is only one
qiime feature-table filter-samples \
  --i-table feature_tables/table-gg-99-no-chlo-mito-minfreq20.qza \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --p-where "sample_location='b'" \
  --p-exclude-ids \
  --o-filtered-table feature_tables/table-gg-99-no-chlo-mito-minfreq20-no-b.qza
  
## Removing sample that belongs to location D (0 meter south) since there is only one
qiime feature-table filter-samples \
  --i-table feature_tables/table-gg-99-no-chlo-mito-minfreq20-no-b.qza \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --p-where "sample_location='d'" \
  --p-exclude-ids \
  --o-filtered-table feature_tables/table-gg-99-no-chlo-mito-minfreq20-no-b-d.qza
  
## Removing DPSW 25 since it has no taxonomy assignment
qiime feature-table filter-samples \
  --i-table feature_tables/table-gg-99-no-chlo-mito-minfreq20-no-b-d-S25.qza \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --p-where "sample_number='25'" \
  --p-exclude-ids \
  --o-filtered-table feature_tables/table-gg-99-no-chlo-mito-minfreq20-no-b-d-S25.qza
  
### Summarize feature table
qiime feature-table summarize \
  --i-table feature_tables/table-gg-99-no-chlo-mito-minfreq20-no-b-d-S25.qza \
  --o-visualization feature_tables/table-gg-99-no-chlo-mito-minfreq20-no-b-d-S25.qzv \
  --m-sample-metadata-file metadata/DPSW_qiime2_metadata.txt  
  
# Create taxa barplot
qiime taxa barplot \
  --i-table feature_tables/table-gg-99-no-chlo-mito-minfreq20-no-b-d-S25.qza \
  --i-taxonomy taxonomy/taxonomy-gg-99.qza \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --o-visualization taxonomy/taxa-bar-plots-gg-99-final.qzv

############## STOP TO CHECK FEATURE TABLE AND PERFORM RAREFACTION   
  
### Perform alpha rarefaction to determine potential rarefaction depth
qiime diversity alpha-rarefaction \
  --i-table feature_tables/table.qza \
  --i-phylogeny trees/fragment_insertion_out/tree.qza \
  --p-max-depth 30000 \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --o-visualization rarefaction/alpha-rarefaction-30000.qzv  
  
############## STOP TO CHECK RAREFACTION AND PERFORM FILTER AT APPROPIATE DEPTH
  
#### Core metrics at 15161 rarefaction
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny trees/fragment_insertion_out/tree.qza \
  --i-table feature_tables/table-gg-99-no-chlo-mito-minfreq20-no-b-d-S25.qza \
  --p-sampling-depth 15161 \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --output-dir core-metrics-results-15161 

mkdir core-metrics-results-15161/alpha_diversity
mkdir core-metrics-results-15161/alpha_diversity/group_signif
mkdir core-metrics-results-15161/alpha_diversity/correlation
mkdir core-metrics-results-15161/beta_diversity
mkdir core-metrics-results-15161/beta_diversity/group_signif
mkdir core-metrics-results-15161/beta_diversity/mantel

mv core-metrics-results-15161/*_vector.qza core-metrics-results-15161/alpha_diversity
mv core-metrics-results-15161/*_emperor.qzv core-metrics-results-15161/beta_diversity
mv core-metrics-results-15161/*_pcoa*.qza core-metrics-results-15161/beta_diversity
mv core-metrics-results-15161/*_matrix.qza core-metrics-results-15161/beta_diversity

# Create taxa barplot
qiime taxa barplot \
  --i-table core-metrics-results-15161/rarefied_table.qza \
  --i-taxonomy taxonomy/taxonomy-gg-99.qza \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --o-visualization taxonomy/taxa-bar-plots-gg-99-final-rarefied.qzv

## Filter table at all taxonomic levels
# Species
qiime taxa collapse \
 --i-table feature_tables/table-gg-99-no-chlo-mito-minfreq20-no-b-d-S25.qza \
 --i-taxonomy taxonomy/taxonomy-gg-99.qza \
 --p-level 7 \
 --o-collapsed-table feature_tables/table-gg-99-no-chlo-mito-minfreq20-no-b-d-S25-L7.qza

# Genus
qiime taxa collapse \
 --i-table feature_tables/table-gg-99-no-chlo-mito-minfreq20-no-b-d-S25.qza \
 --i-taxonomy taxonomy/taxonomy-gg-99.qza \
 --p-level 6 \
 --o-collapsed-table feature_tables/table-gg-99-no-chlo-mito-minfreq20-no-b-d-S25-L6.qza

# Family
qiime taxa collapse \
 --i-table feature_tables/table-gg-99-no-chlo-mito-minfreq20-no-b-d-S25.qza \
 --i-taxonomy taxonomy/taxonomy-gg-99.qza \
 --p-level 5 \
 --o-collapsed-table feature_tables/table-gg-99-no-chlo-mito-minfreq20-no-b-d-S25-L5.qza

# Order
qiime taxa collapse \
 --i-table feature_tables/table-gg-99-no-chlo-mito-minfreq20-no-b-d-S25.qza \
 --i-taxonomy taxonomy/taxonomy-gg-99.qza \
 --p-level 4 \
 --o-collapsed-table feature_tables/table-gg-99-no-chlo-mito-minfreq20-no-b-d-S25-L4.qza

# Class
qiime taxa collapse \
 --i-table feature_tables/table-gg-99-no-chlo-mito-minfreq20-no-b-d-S25.qza \
 --i-taxonomy taxonomy/taxonomy-gg-99.qza \
 --p-level 3 \
 --o-collapsed-table feature_tables/table-gg-99-no-chlo-mito-minfreq20-no-b-d-S25-L3.qza

# Phylum
qiime taxa collapse \
 --i-table feature_tables/table-gg-99-no-chlo-mito-minfreq20-no-b-d-S25.qza \
 --i-taxonomy taxonomy/taxonomy-gg-99.qza \
 --p-level 2 \
 --o-collapsed-table feature_tables/table-gg-99-no-chlo-mito-minfreq20-no-b-d-S25-L2.qza
  
### Longitudinal Analysis
mkdir longitudinal longitudinal/all

# All locations
## Feature volatility
# Tells your important features 
## Add 0c
# ASV level
qiime longitudinal feature-volatility \
  --i-table feature_tables/table-gg-99-no-chlo-mito-minfreq20-no-b-d-S25.qza \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --p-state-column pmi_months \
  --p-individual-id-column pig_number \
  --p-cv 5 \
  --p-estimator 'RandomForestRegressor' \
  --p-parameter-tuning \
  --p-importance-threshold 'q1' \
  --p-feature-count 100 \
  --p-n-estimators 1000 \
  --p-random-state 999 \
  --p-n-jobs -1 \
  --output-dir longitudinal/all/pmi_months-feature-volatility 

# Genera level
qiime longitudinal feature-volatility \
  --i-table feature_tables/table-gg-99-no-chlo-mito-minfreq20-no-b-d-S25-L6.qza \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --p-state-column pmi_months \
  --p-individual-id-column pig_number \
  --p-cv 5 \
  --p-estimator 'RandomForestRegressor' \
  --p-parameter-tuning \
  --p-importance-threshold 'q1' \
  --p-feature-count 100 \
  --p-n-estimators 1000 \
  --p-random-state 999 \
  --p-n-jobs -1 \
  --output-dir longitudinal/all/pmi_months-feature-volatility-genera 
  
# Family level
qiime longitudinal feature-volatility \
  --i-table feature_tables/table-gg-99-no-chlo-mito-minfreq20-no-b-d-S25-L5.qza \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --p-state-column pmi_months \
  --p-individual-id-column pig_number \
  --p-cv 5 \
  --p-estimator 'RandomForestRegressor' \
  --p-parameter-tuning \
  --p-importance-threshold 'q1' \
  --p-feature-count 100 \
  --p-n-estimators 1000 \
  --p-random-state 999 \
  --p-n-jobs -1 \
  --output-dir longitudinal/all/pmi_months-feature-volatility-family 
  
# Order level
qiime longitudinal feature-volatility \
  --i-table feature_tables/table-gg-99-no-chlo-mito-minfreq20-no-b-d-S25-L4.qza \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --p-state-column pmi_months \
  --p-individual-id-column pig_number \
  --p-cv 5 \
  --p-estimator 'RandomForestRegressor' \
  --p-parameter-tuning \
  --p-importance-threshold 'q1' \
  --p-feature-count 100 \
  --p-n-estimators 1000 \
  --p-random-state 999 \
  --p-n-jobs -1 \
  --output-dir longitudinal/all/pmi_months-feature-volatility-order 
  
# Class level
qiime longitudinal feature-volatility \
  --i-table feature_tables/table-gg-99-no-chlo-mito-minfreq20-no-b-d-S25-L3.qza \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --p-state-column pmi_months \
  --p-individual-id-column pig_number \
  --p-cv 5 \
  --p-estimator 'RandomForestRegressor' \
  --p-parameter-tuning \
  --p-importance-threshold 'q1' \
  --p-feature-count 100 \
  --p-n-estimators 1000 \
  --p-random-state 999 \
  --p-n-jobs -1 \
  --output-dir longitudinal/all/pmi_months-feature-volatility-class 

# Phylum level
qiime longitudinal feature-volatility \
  --i-table feature_tables/table-gg-99-no-chlo-mito-minfreq20-no-b-d-S25-L2.qza \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --p-state-column pmi_months \
  --p-individual-id-column pig_number \
  --p-cv 5 \
  --p-estimator 'RandomForestRegressor' \
  --p-parameter-tuning \
  --p-importance-threshold 'q1' \
  --p-feature-count 100 \
  --p-n-estimators 1000 \
  --p-random-state 999 \
  --p-n-jobs -1 \
  --output-dir longitudinal/all/pmi_months-feature-volatility-phylum 

## Create time groups for linear models
# Filter time series to 0-1 months at genera level
qiime feature-table filter-samples \
  --i-table feature_tables/table-gg-99-no-chlo-mito-minfreq20-no-b-d-S25-L6.qza \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --p-where "pmi_months<=1" \
  --o-filtered-table feature_tables/genera-01-table.qza

# Filter time series to 0-12 months at genera level
qiime feature-table filter-samples \
  --i-table feature_tables/table-gg-99-no-chlo-mito-minfreq20-no-b-d-S25-L6.qza \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --p-where "pmi_months<=12" \
  --o-filtered-table feature_tables/genera-012-table.qza

# Filter time series to 0-24 months at genera level
qiime feature-table filter-samples \
  --i-table feature_tables/table-gg-99-no-chlo-mito-minfreq20-no-b-d-S25-L6.qza \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --p-where "pmi_months<=24" \
  --o-filtered-table feature_tables/genera-024-table.qza

# Filter time series to 6-24 months at genera level
qiime feature-table filter-samples \
  --i-table feature_tables/table-gg-99-no-chlo-mito-minfreq20-no-b-d-S25-L6.qza \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --p-where "pmi_months>=6 AND pmi_months<=24" \
  --o-filtered-table feature_tables/genera-624-table.qza

# Filter time series to 24-84 months at genera level
qiime feature-table filter-samples \
  --i-table feature_tables/table-gg-99-no-chlo-mito-minfreq20-no-b-d-S25-L6.qza \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --p-where "pmi_months>=24 AND pmi_months<=84" \
  --o-filtered-table feature_tables/genera-2484-table.qza

# Convert to relative abundances
qiime feature-table relative-frequency \
	--i-table feature_tables/genera-01-table.qza \
	--o-relative-frequency-table feature_tables/genera-01-relfreq-table.qza

qiime feature-table relative-frequency \
	--i-table feature_tables/genera-012-table.qza \
	--o-relative-frequency-table feature_tables/genera-012-relfreq-table.qza

qiime feature-table relative-frequency \
	--i-table feature_tables/genera-024-table.qza \
	--o-relative-frequency-table feature_tables/genera-024-relfreq-table.qza

qiime feature-table relative-frequency \
	--i-table feature_tables/genera-624-table.qza \
	--o-relative-frequency-table feature_tables/genera-624-relfreq-table.qza

qiime feature-table relative-frequency \
	--i-table feature_tables/genera-2484-table.qza \
	--o-relative-frequency-table feature_tables/genera-2484-relfreq-table.qza

### Create LME with rel freq table of time series and metric of important features 
# Genera chosen based on their relatively high importance, mean abundance, and cum avg change

mkdir longitudinal/all/genera_lme
## LME k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhodospirillales;f__Rhodospirillaceae;g__
# 0-12 months  
qiime longitudinal linear-mixed-effects \
  --m-metadata-file metadata/DPSW_qiime2_metadata_rename.txt \
  --i-table feature_tables/genera-012-relfreq-table.qza \
  --p-metric 'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhodospirillales;f__Rhodospirillaceae;g__' \
  --p-group-columns sample_detail \
  --p-state-column pmi_months \
  --p-individual-id-column pig_number \
  --o-visualization longitudinal/all/genera_lme/genera-rhodo-lme-012.qzv

# 24-84 months
qiime longitudinal linear-mixed-effects \
  --m-metadata-file metadata/DPSW_qiime2_metadata_rename.txt \
  --i-table feature_tables/genera-2484-relfreq-table.qza \
  --p-metric 'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhodospirillales;f__Rhodospirillaceae;g__' \
  --p-group-columns sample_detail \
  --p-state-column pmi_months \
  --p-individual-id-column pig_number \
  --o-visualization longitudinal/all/genera_lme/genera-rhodo-lme-2484.qzv

## LME k__Bacteria;p__Actinobacteria;c__Acidimicrobiia;o__Acidimicrobiales;f__EB1017;g__
# 0-24 month
qiime longitudinal linear-mixed-effects \
  --m-metadata-file metadata/DPSW_qiime2_metadata_rename.txt \
  --i-table feature_tables/genera-024-relfreq-table.qza \
  --p-metric 'k__Bacteria;p__Actinobacteria;c__Acidimicrobiia;o__Acidimicrobiales;f__EB1017;g__' \
  --p-group-columns sample_detail \
  --p-state-column pmi_months \
  --p-individual-id-column pig_number \
  --o-visualization longitudinal/all/genera_lme/genera-eb1017-lme-024.qzv

## LME k__Bacteria;p__Verrucomicrobia;c__[Spartobacteria];o__[Chthoniobacterales];f__[Chthoniobacteraceae];g__Chthoniobacter
# 0-1 month
qiime longitudinal linear-mixed-effects \
  --m-metadata-file metadata/DPSW_qiime2_metadata_rename.txt \
  --i-table feature_tables/genera-01-relfreq-table.qza \
  --p-metric 'k__Bacteria;p__Verrucomicrobia;c__[Spartobacteria];o__[Chthoniobacterales];f__[Chthoniobacteraceae];g__Chthoniobacter' \
  --p-group-columns sample_detail \
  --p-state-column pmi_months \
  --p-individual-id-column pig_number \
  --o-visualization longitudinal/all/genera_lme/genera-chthon-lme-01.qzv

# 0-24 month
qiime longitudinal linear-mixed-effects \
  --m-metadata-file metadata/DPSW_qiime2_metadata_rename.txt \
  --i-table feature_tables/genera-024-relfreq-table.qza \
  --p-metric 'k__Bacteria;p__Verrucomicrobia;c__[Spartobacteria];o__[Chthoniobacterales];f__[Chthoniobacteraceae];g__Chthoniobacter' \
  --p-group-columns sample_detail \
  --p-state-column pmi_months \
  --p-individual-id-column pig_number \
  --o-visualization longitudinal/all/genera_lme/genera-chthon-lme-024.qzv

# 24-84 month
qiime longitudinal linear-mixed-effects \
  --m-metadata-file metadata/DPSW_qiime2_metadata_rename.txt \
  --i-table feature_tables/genera-2484-relfreq-table.qza \
  --p-metric 'k__Bacteria;p__Verrucomicrobia;c__[Spartobacteria];o__[Chthoniobacterales];f__[Chthoniobacteraceae];g__Chthoniobacter' \
  --p-group-columns sample_detail \
  --p-state-column pmi_months \
  --p-individual-id-column pig_number \
  --o-visualization longitudinal/all/genera_lme/genera-chthon-lme-2484.qzv

## LME k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Hyphomicrobiaceae;g__Devosia
# 6-24 month
qiime longitudinal linear-mixed-effects \
  --m-metadata-file metadata/DPSW_qiime2_metadata_rename.txt \
  --i-table feature_tables/genera-624-relfreq-table.qza \
  --p-metric 'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Hyphomicrobiaceae;g__Devosia' \
  --p-group-columns sample_detail \
  --p-state-column pmi_months \
  --p-individual-id-column pig_number \
  --o-visualization longitudinal/all/genera_lme/genera-devo-lme-624.qzv

# 24-84 month 
qiime longitudinal linear-mixed-effects \
  --m-metadata-file metadata/DPSW_qiime2_metadata_rename.txt \
  --i-table feature_tables/genera-2484-relfreq-table.qza \
  --p-metric 'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Hyphomicrobiaceae;g__Devosia' \
  --p-group-columns sample_detail \
  --p-state-column pmi_months \
  --p-individual-id-column pig_number \
  --o-visualization longitudinal/all/genera_lme/genera-devo-lme-2484.qzv

## PCoA-based volatility
# Unweighted UniFrac
qiime longitudinal volatility \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --m-metadata-file core-metrics-results-15161/beta_diversity/unweighted_unifrac_pcoa_results.qza \
  --p-state-column pmi_months \
  --p-default-metric 'Axis 1' \
  --p-default-group-column sample_detail \
  --p-individual-id-column pig_number \
  --o-visualization longitudinal/all/pcoa-vol-unwunifrac.qzv

# Weighted UniFrac
qiime longitudinal volatility \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --m-metadata-file core-metrics-results-15161/beta_diversity/weighted_unifrac_pcoa_results.qza \
  --p-state-column pmi_months \
  --p-default-metric 'Axis 1' \
  --p-default-group-column sample_detail \
  --p-individual-id-column pig_number \
  --o-visualization longitudinal/all/pcoa-vol-wunifrac.qzv

## Alpha diversity volatility
# Faith's
qiime longitudinal volatility \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --m-metadata-file core-metrics-results-15161/alpha_diversity/faith_pd_vector.qza \
  --p-default-metric faith_pd \
  --p-default-group-column sample_detail \
  --p-state-column pmi_months \
  --p-individual-id-column pig_number \
  --o-visualization longitudinal/all/faith_pd-volatility.qzv
  
# Shannon
qiime longitudinal volatility \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --m-metadata-file core-metrics-results-15161/alpha_diversity/shannon_vector.qza \
  --p-default-metric shannon \
  --p-default-group-column sample_detail \
  --p-state-column pmi_months \
  --p-individual-id-column pig_number \
  --o-visualization longitudinal/all/shannon-volatility.qzv
  
# Richness
qiime longitudinal volatility \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --m-metadata-file core-metrics-results-15161/alpha_diversity/observed_otus_vector.qza \
  --p-default-metric observed_otus \
  --p-default-group-column sample_detail \
  --p-state-column pmi_months \
  --p-individual-id-column pig_number \
  --o-visualization longitudinal/all/observed_otus-volatility.qzv

# Evenness
qiime longitudinal volatility \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --m-metadata-file core-metrics-results-15161/alpha_diversity/evenness_vector.qza \
  --p-default-metric pielou_e \
  --p-default-group-column sample_detail \
  --p-state-column pmi_months \
  --p-individual-id-column pig_number \
  --o-visualization longitudinal/all/evenness-volatility.qzv

## Create time groups for alpha linear models
# Filter time series to 0-24 months at genera level
qiime feature-table filter-samples \
  --i-table core-metrics-results-15161/rarefied_table.qza \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --p-where "pmi_months<=24" \
  --o-filtered-table feature_tables/024-rarefied-table.qza

# Filter time series to 24-120 months at genera level
qiime feature-table filter-samples \
  --i-table core-metrics-results-15161/rarefied_table.qza \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --p-where "pmi_months>=24 AND pmi_months<=120" \
  --o-filtered-table feature_tables/24120-rarefied-table.qza

# Filter time series to 1-120 months at genera level
qiime feature-table filter-samples \
  --i-table core-metrics-results-15161/rarefied_table.qza \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --p-where "pmi_months>=1 AND pmi_months<=120" \
  --o-filtered-table feature_tables/1120-rarefied-table.qza

# Create time groups alpha metrics
mkdir longitudinal/all/alpha_time_groups
qiime diversity alpha \
	--i-table feature_tables/024-rarefied-table.qza \
	--p-metric observed_otus \
	--o-alpha-diversity longitudinal/all/alpha_time_groups/observed_otus_vector_024.qza

qiime diversity alpha \
	--i-table feature_tables/24120-rarefied-table.qza \
	--p-metric observed_otus \
	--o-alpha-diversity longitudinal/all/alpha_time_groups/observed_otus_vector_24120.qza

qiime diversity alpha \
	--i-table feature_tables/1120-rarefied-table.qza \
	--p-metric shannon \
	--o-alpha-diversity longitudinal/all/alpha_time_groups/shannon_vector_1120.qza

qiime diversity alpha \
	--i-table feature_tables/1120-rarefied-table.qza \
	--p-metric pielou_e \
	--o-alpha-diversity longitudinal/all/alpha_time_groups/evenness_vector_1120.qza
		
qiime diversity alpha-phylogenetic \
	--i-table feature_tables/024-rarefied-table.qza \
	--p-metric faith_pd \
	--i-phylogeny trees/fragment_insertion_out/tree.qza \
	--o-alpha-diversity longitudinal/all/alpha_time_groups/faith_pd_vector_024.qza

qiime diversity alpha-phylogenetic \
	--i-table feature_tables/24120-rarefied-table.qza \
	--p-metric faith_pd \
	--i-phylogeny trees/fragment_insertion_out/tree.qza \
	--o-alpha-diversity longitudinal/all/alpha_time_groups/faith_pd_vector_24120.qza

## Linear mixed effect model
# Faith's
# 0-24 month
qiime longitudinal linear-mixed-effects \
  --m-metadata-file metadata/DPSW_qiime2_metadata_rename.txt \
  --m-metadata-file longitudinal/all/alpha_time_groups/faith_pd_vector_024.qza \
  --p-metric faith_pd \
  --p-group-columns sample_detail \
  --p-state-column pmi_months \
  --p-individual-id-column pig_number \
  --o-visualization longitudinal/all/alpha_time_groups/faith_pd-024-lme.qzv

# 24-120 month
qiime longitudinal linear-mixed-effects \
  --m-metadata-file metadata/DPSW_qiime2_metadata_rename.txt \
  --m-metadata-file longitudinal/all/alpha_time_groups/faith_pd_vector_24120.qza \
  --p-metric faith_pd \
  --p-group-columns sample_detail \
  --p-state-column pmi_months \
  --p-individual-id-column pig_number \
  --o-visualization longitudinal/all/alpha_time_groups/faith_pd-24120-lme.qzv

# Richness
# 0-24 month
qiime longitudinal linear-mixed-effects \
  --m-metadata-file metadata/DPSW_qiime2_metadata_rename.txt \
  --m-metadata-file longitudinal/all/alpha_time_groups/observed_otus_vector_024.qza \
  --p-metric observed_otus \
  --p-group-columns sample_detail \
  --p-state-column pmi_months \
  --p-individual-id-column pig_number \
  --o-visualization longitudinal/all/alpha_time_groups/observed_otus-024-lme.qzv

# 24-120 month
qiime longitudinal linear-mixed-effects \
  --m-metadata-file metadata/DPSW_qiime2_metadata_rename.txt \
  --m-metadata-file longitudinal/all/alpha_time_groups/observed_otus_vector_24120.qza \
  --p-metric observed_otus \
  --p-group-columns sample_detail \
  --p-state-column pmi_months \
  --p-individual-id-column pig_number \
  --o-visualization longitudinal/all/alpha_time_groups/observed_otus-24120-lme.qzv

# Shannon
# 1-120 month
qiime longitudinal linear-mixed-effects \
  --m-metadata-file metadata/DPSW_qiime2_metadata_rename.txt \
  --m-metadata-file longitudinal/all/alpha_time_groups/shannon_vector_1120.qza \
  --p-metric shannon \
  --p-group-columns sample_detail \
  --p-state-column pmi_months \
  --p-individual-id-column pig_number \
  --o-visualization longitudinal/all/alpha_time_groups/shannon-1120-lme.qzv

# Evenness
# 1-120 month
qiime longitudinal linear-mixed-effects \
  --m-metadata-file metadata/DPSW_qiime2_metadata_rename.txt \
  --m-metadata-file longitudinal/all/alpha_time_groups/evenness_vector_1120.qza \
  --p-metric pielou_e \
  --p-group-columns sample_detail \
  --p-state-column pmi_months \
  --p-individual-id-column pig_number \
  --o-visualization longitudinal/all/alpha_time_groups/evenness-1120-lme.qzv

### LME Testing the liner trends detected in the PCoA volatility
mkdir longitudinal/all/beta_time_groups/

## Filter distance matrices
# Unweighted 0-12 months
qiime diversity filter-distance-matrix \
  --i-distance-matrix core-metrics-results-15161/beta_diversity/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --p-where "pmi_months<=12" \
  --o-filtered-distance-matrix longitudinal/all/beta_time_groups/unweighted_unifrac_distance_matrix_012.qza

# Unweighted 12-120 months
qiime diversity filter-distance-matrix \
  --i-distance-matrix core-metrics-results-15161/beta_diversity/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --p-where "pmi_months>=12 AND pmi_months<=120" \
  --o-filtered-distance-matrix longitudinal/all/beta_time_groups/unweighted_unifrac_distance_matrix_12120.qza

# Weighted 0-12 months
qiime diversity filter-distance-matrix \
  --i-distance-matrix core-metrics-results-15161/beta_diversity/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --p-where "pmi_months<=12" \
  --o-filtered-distance-matrix longitudinal/all/beta_time_groups/weighted_unifrac_distance_matrix_012.qza

# Weighted 12-120 months
qiime diversity filter-distance-matrix \
  --i-distance-matrix core-metrics-results-15161/beta_diversity/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --p-where "pmi_months>=12 AND pmi_months<=120" \
  --o-filtered-distance-matrix longitudinal/all/beta_time_groups/weighted_unifrac_distance_matrix_12120.qza

## Create PCOAs for each new time series
# Unweighted 0-12 months
qiime diversity pcoa \
	--i-distance-matrix longitudinal/all/beta_time_groups/unweighted_unifrac_distance_matrix_012.qza \
	--o-pcoa longitudinal/all/beta_time_groups/unweighted_unifrac_pcoa_012.qza

# Unweighted 12-120 months
qiime diversity pcoa \
	--i-distance-matrix longitudinal/all/beta_time_groups/unweighted_unifrac_distance_matrix_12120.qza \
	--o-pcoa longitudinal/all/beta_time_groups/unweighted_unifrac_pcoa_12120.qza

# Weighted 0-12 months
qiime diversity pcoa \
	--i-distance-matrix longitudinal/all/beta_time_groups/weighted_unifrac_distance_matrix_012.qza \
	--o-pcoa longitudinal/all/beta_time_groups/weighted_unifrac_pcoa_012.qza

# Weighted 12-120 months
qiime diversity pcoa \
	--i-distance-matrix longitudinal/all/beta_time_groups/weighted_unifrac_distance_matrix_12120.qza \
	--o-pcoa longitudinal/all/beta_time_groups/weighted_unifrac_pcoa_12120.qza

## PCOA Axis 1 LME
# Unweighted 0-12 months
qiime longitudinal linear-mixed-effects \
  --m-metadata-file metadata/DPSW_qiime2_metadata_rename.txt \
  --m-metadata-file longitudinal/all/beta_time_groups/unweighted_unifrac_pcoa_012.qza \
  --p-metric 'Axis 1' \
  --p-group-columns sample_detail \
  --p-state-column pmi_months \
  --p-individual-id-column pig_number \
  --o-visualization longitudinal/all/beta_time_groups/unwunifrac-pc1-distances-012-LME.qzv

# Unweighted 12-120 months
qiime longitudinal linear-mixed-effects \
  --m-metadata-file metadata/DPSW_qiime2_metadata_rename.txt \
  --m-metadata-file longitudinal/all/beta_time_groups/unweighted_unifrac_pcoa_12120.qza \
  --p-metric 'Axis 1' \
  --p-group-columns sample_detail \
  --p-state-column pmi_months \
  --p-individual-id-column pig_number \
  --o-visualization longitudinal/all/beta_time_groups/unwunifrac-pc1-distances-12120-LME.qzv

# Weighted 0-12 months
qiime longitudinal linear-mixed-effects \
  --m-metadata-file metadata/DPSW_qiime2_metadata_rename.txt \
  --m-metadata-file longitudinal/all/beta_time_groups/weighted_unifrac_pcoa_012.qza \
  --p-metric 'Axis 1' \
  --p-group-columns sample_detail \
  --p-state-column pmi_months \
  --p-individual-id-column pig_number \
  --o-visualization longitudinal/all/beta_time_groups/wunifrac-pc1-distances-012-LME.qzv

# Weighted 12-120 months
qiime longitudinal linear-mixed-effects \
  --m-metadata-file metadata/DPSW_qiime2_metadata_rename.txt \
  --m-metadata-file longitudinal/all/beta_time_groups/weighted_unifrac_pcoa_12120.qza \
  --p-metric 'Axis 1' \
  --p-group-columns sample_detail \
  --p-state-column pmi_months \
  --p-individual-id-column pig_number \
  --o-visualization longitudinal/all/beta_time_groups/wunifrac-pc1-distances-12120-LME.qzv

### Test if at 120 months PC1 or beta div and alpha divs are still significantly different
# Filter time series 
qiime feature-table filter-samples \
  --i-table core-metrics-results-15161/rarefied_table.qza \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --p-where "pmi_months=120" \
  --o-filtered-table feature_tables/120-rarefied-table.qza

# Core metrics at 15161 rarefaction
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny trees/fragment_insertion_out/tree.qza \
  --i-table feature_tables/120-rarefied-table.qza \
  --p-sampling-depth 15161 \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --output-dir core-metrics-results-15161-only12mo 

mkdir core-metrics-results-15161-only12mo/alpha_diversity
mkdir core-metrics-results-15161-only12mo/alpha_diversity/group_signif
mkdir core-metrics-results-15161-only12mo/alpha_diversity/correlation
mkdir core-metrics-results-15161-only12mo/beta_diversity
mkdir core-metrics-results-15161-only12mo/beta_diversity/group_signif
mkdir core-metrics-results-15161-only12mo/beta_diversity/mantel

mv core-metrics-results-15161-only12mo/*_vector.qza core-metrics-results-15161-only12mo/alpha_diversity
mv core-metrics-results-15161-only12mo/*_emperor.qzv core-metrics-results-15161-only12mo/beta_diversity
mv core-metrics-results-15161-only12mo/*_pcoa*.qza core-metrics-results-15161-only12mo/beta_diversity
mv core-metrics-results-15161-only12mo/*_matrix.qza core-metrics-results-15161-only12mo/beta_diversity


## Alpha group significance
# Faith's pd
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results-15161-only12mo/alpha_diversity/faith_pd_vector.qza \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --o-visualization core-metrics-results-15161-only12mo/alpha_diversity/group_signif/faith-pd-group-significance.qzv

# Peilou's evenness
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results-15161-only12mo/alpha_diversity/evenness_vector.qza \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --o-visualization core-metrics-results-15161-only12mo/alpha_diversity/group_signif/evenness-significance.qzv

# Shannon
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results-15161-only12mo/alpha_diversity/shannon_vector.qza \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --o-visualization core-metrics-results-15161-only12mo/alpha_diversity/group_signif/shannon-significance.qzv

# Observed OTUs (richness)
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results-15161-only12mo/alpha_diversity/observed_otus_vector.qza \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --o-visualization core-metrics-results-15161-only12mo/alpha_diversity/group_signif/observed_otus_significance.qzv
 
## Beta Diversity Group Significance
# Sample Location
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results-15161-only12mo/beta_diversity/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --m-metadata-column sample_detail \
  --o-visualization core-metrics-results-15161-only12mo/beta_diversity/group_signif/sample_location-unweighted-unifrac-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results-15161-only12mo/beta_diversity/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/DPSW_qiime2_metadata.txt \
  --m-metadata-column sample_detail \
  --o-visualization core-metrics-results-15161-only12mo/beta_diversity/group_signif/sample_location-weighted-unifrac-significance.qzv \
  --p-pairwise
  