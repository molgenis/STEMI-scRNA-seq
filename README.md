# STEMI-scRNA-seq
This repository contains the code that was used for the scRNA-seq study of PBMCs in a population of patients suffering from a STEMI. More information can be found in the original paper: https://doi.org/10.1101/2023.05.02.23289370

## data availability and description

Expression data is available here:
https://eqtlgen.org/sc/datasets/blokland2024-dataset.html

Expression data is available in three flavours:
- QC-ed and normalised (A)
- QC-ed without normalisation (B)
- pre-QC, STEMI only (C)

### A. normalized and QC-ed data

To use the normalized and QC-ed data, the following files are required:
- barcodes.tsv.gz
- features.tsv.gz
- matrix.mtx.gz

Given these three files are located in a given folder, with those exact filenames, they can be loaded into Seurat, using the following command:
```r
stemi_processed <- Read10X('/dir/to/three/files/', gene.column = 1, cell.column = 1)
```

or in Scanpy using the original filenames:

```python
# lead count data
stemi_counts = sc.read_mtx('/dir/to/three/files/matrix.mtx.gz')
# read barcodes
stemi_bc=pd.read_csv('/dir/to/three/files//barcodes.tsv.gz', header=None)
# read features
stemi_features=pd.read_csv('/dir/to/three/files/features.tsv.gz', header=None)
# transpose to scanpy format
stemi_counts = stemi_counts.T
# and make object
stemi_processed = sc.AnnData(stemi_counts)
# add barcodes and genes to obs and vars
stemi_processed.obs['cell_id']= stemi_bc[0].tolist()
stemi_processed.var['gene_name']= stemi_features[0].tolist()
# set indices for the obs and vars
stemi_processed.obs.index = stemi_processed.obs['cell_id']
stemi_processed.var.index = stemi_processed.var['gene_name']
```

### B. raw QC-ed data

To use the non-normalized counts, the following files are required:
- barcodes.tsv.gz
- features_raw.tsv.gz
- matrix_raw.mtx.gz

Given these three files are located in a given folder, and matrix_raw.mtx.gz is renamed to matrix.mtx.gz, and features_raw.tsv.gz is renamed to features.tsv.gz, they can be loaded into Seurat, using the following command:
```r
stemi_raw <- Read10X('/dir/to/three/files/', gene.column = 1, cell.column = 1)
```

or in Scanpy using the original filenames:

```python
# lead count data
stemi_counts = sc.read_mtx('/dir/to/three/files/matrix_raw.mtx.gz')
# read barcodes
stemi_bc=pd.read_csv('/dir/to/three/files//barcodes.tsv.gz', header=None)
# read features
stemi_raw=pd.read_csv('/dir/to/three/files/features_raw.tsv.gz', header=None)
# transpose to scanpy format
stemi_counts = stemi_counts.T
# and make object
stemi_raw = sc.AnnData(stemi_counts)
# add barcodes and genes to obs and vars
stemi_raw.obs['cell_id']= stemi_bc[0].tolist()
stemi_raw.var['gene_name']= stemi_features[0].tolist()
# set indices for the obs and vars
stemi_raw.obs.index = stemi_raw.obs['cell_id']
stemi_raw.var.index = stemi_raw.var['gene_name']

```

### C. pre-QC data

Data before QC is only available for the STEMI samples. The data before QC of the control samples is part of the 1m-BloodNL study: https://eqtlgen.org/sc/datasets/1m-scbloodnl.html

To use the pre-QC STEMI non-normalized counts, the following files are required:
- stemi_unfiltered_barcodes.tsv.gz
- stemi_unfiltered_features_raw.tsv.gz
- stemi_unfiltered_matrix_raw.mtx.gz

Given these three files are located in a given folder, and stemi_unfiltered_matrix.mtx.gz is renamed to matrix.mtx.gz, stemi_unfiltered_features.tsv.gz is renamed to features.tsv.gz and stemi_unfiltered_barcodes.tsv.gz is renamed to barcodes.tsv.gz, they can be loaded into Seurat, using the following command:

```r
stemi_unfiltered <- Read10X('/dir/to/three/files/', gene.column = 1, cell.column = 1)
```

or in Scanpy using the original filenames:

```python
# lead count data
stemi_counts = sc.read_mtx('/dir/to/three/files/stemi_unfiltered_matrix.mtx.gz')
# read barcodes
stemi_bc=pd.read_csv('/dir/to/three/files//stemi_unfiltered_barcodes.tsv.gz', header=None)
# read features
stemi_noqc=pd.read_csv('/dir/to/three/files/stemi_unfiltered_features.tsv.gz', header=None)
# transpose to scanpy format
stemi_counts = stemi_counts.T
# and make object
stemi_noqc = sc.AnnData(stemi_counts)
# add barcodes and genes to obs and vars
stemi_noqc.obs['cell_id']= stemi_bc[0].tolist()
stemi_noqc.var['gene_name']= stemi_features[0].tolist()
# set indices for the obs and vars
stemi_noqc.obs.index = stemi_noqc.obs['cell_id']
stemi_noqc.var.index = stemi_noqc.var['gene_name']
```

### metadata

metadata info on sample assignment, celltype assignment and timepoint/condition, is stored in the metadata.tsv.gz file. This data can be added in Seurat like this:

```r
stemi_metadata <- read.table('/dir/to/metadata.tsv.gz', header = T, row.names = 1, sep = '\t')
stemi_processed <- AddMetaData(stemi_processed, stemi_metadata[, setdiff(colnames(stemi_processed), 'orig.ident')])
```

or in scanpy like this:

```python
stemi_metadata = pandas.read_csv('/dir/to/metadata.tsv.gz', sep = '\t', header = 0, index_col = 0)
stemi_processed.obs = pandas.concat([stemi_processed.obs, stemi_metadata], axis=1).reindex(stemi_processed.obs.index)
```


The metadata contains several columns that are relevant in the contex of sample/condition/celltype assignment. The most important ones will be highlighted here:

- '*orig.ident*' the subset of the data the cell belongs to. stemi_v2 are cells from STEMI patients, processed with 10x v2 chemistry. stemi_v3 are cells from STEMI patients, processed with 10x v3 chemistry. hc_v2 are cells from controls, processed with 10x v2 chemistry. hc_v3 are cells from controls, processed with 10x v3 chemistry.

- '*chem*' the 10x version chemistry the cell was processed with

- '*lane*' and '*batch*' the 10x experimental batch the cells were processed in

- '*timepoint.final*' the condition the cell came from. UT are cells coming from controls. Baseline are cells coming from 0-4h post-STEMI. t24h are cells coming from 24h after STEMI. t8w are cells coming from t6-8 weeks after STEMI

- '*assignment.final*' the donors the cells came from. S1-S38 are from STEMI patients. C1-120 are from the controls. Note that the controls use the same sample numbers as in the 1m-BloodNL study (https://eqtlgen.org/sc/datasets/1m-scbloodnl.html), and as such do not run from 1-38, even though there are only 38 controls.

- '*cell_type*' the cell type assignment for the cell

- '*cell_type_lowerres*' the reclassification of the cell type to a lower granularity

## pipeline steps

If want to rerun any of the analysis steps in R, consider using the Singularity image used for most of the analyses: https://github.com/royoelen/single-cell-container-server

Relevant software versions:
- Seurat v4 (https://github.com/satijalab/seurat)
- Eagle v2.x (https://github.com/poruloh/Eagle)
- Souporcell v1.x (https://github.com/wheaton5/souporcell)
- in-house eQTL pipeline2 v1.4.0 (https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline)
- R v4.x (https://svn.r-project.org/R/)

Below we will outline the computational steps taken

### cell viability

'*experiment_qc/stemi_plot_viabilities.Rmd*' show the cell viabilities, to check if there is a difference between the controls and STEMI, as well as between timepoints for the STEMI samples

### alignment to reference genome

'*preprocessing/stemi_create_alignment_jobs.sh*' create SLURM job scripts to perform alignment of the sequence data to a human reference genome (note that the script uses the newer build 38, however for the original study, hg19 was used)

### doublet filtering and demultiplexing

'*demultiplexing/stemi_create_souporcell_jobs.sh*' create SLURM job scripts to do demultiplexing using the Souporcell tool

'*demultiplexing/stemi_create_v2_demux_jobs.sh*' create SLURM job scripts to do demultiplexing using the Demuxlet tool (note that this data is present in the final object, but the Souporcell tool was eventually chosen for the final assignments)

'*demultiplexing/stemi_create_v3_demux_jobs.sh*' create SLURM job scripts to do demultiplexing using the Demuxlet tool (note that this data is present in the final object, but the Souporcell tool was eventually chosen for the final assignments)

'*demultiplexing/stemi_scrublet.py*' run scrublet on each lane, to do doublet detection using the Scrublet tool (note that this data is present in the final object, but the Souporcell tool was eventually chosen for the final assignments)

'*demultiplexing/stemi_subset_genotypes*' subset the full genotype data to the genotype data of the participants present for a specific 10x lane, to reduce computational load

'*demultiplexing/stemi_correlate_soup_gts.R*' script to correlate the cluster genotypes that are output from Souporcell, to the genotypes of the participants present in the 10x lane, to do the final sample assignment

'*demultiplexing/stemi_correlate_soup_gts.sh*' shell script wrapper for the aforementioned R script, to perform the sample assignment for all the 10x lanes


### reading count data into Seurat

'*preprocessing/preprocess_stemi.R*' read the CellRanger output directories of the STEMI cohort into Seurat, do doublet filtering based on the Souporcell output mentioned above, and filter based mitochondrial percentage

'*preprocessing/preprocess_HC.R*' read the CellRanger output directories of the LifeLines cohort into Seurat, do doublet filtering based on the Souporcell output mentioned above, and filter based mitochondrial percentage, and filter on the participants we chose as our controls

### merge controls and STEMI cohort

'*dataset_integration/stemi_integrate_datasets.R*' merge the controls and the STEMI population, as well as the two different 10x chemistries

### cell type classification

'*cell_type_classification/stemi_predict_celtypes_azimuth.R*' use Azimuth from Seurat v4 to predict the celltypes based on a reference

'*cell_type_classification/stemi_celltype_classification.R*' add a lower resolution celltype annotation

### differential proportion analysis

'*cell_type_composition/stemi_cell_type_plots.R*' plot the cell type proportions present

'*cell_type_composition/stemi_cell_type_composition_differences.R*' use various methods to do differential proportion analysis, finally using the elife method

### differential gene expression analysis

'*differential_expression/stemi_mast.R*' perform differential gene expression analysis using the Seurat implementation of MAST

'*differential_expression/stemi_mast_output_analysis*' perform meta-analysis over the DE results of 10x v2 and 10x v3

'*differential_expression/stemi_limma_dream.R*' perform differential gene expression analysis limma dream

'*differential_expression/stemi_plot_de.Rmd*' plot the number of DE genes and their sharedness

'*differential_expression/stemi_plot_de_genes.Rmd*' plot the top DE genes for monocytes and NKs cells, and visualize these

'*differential_expression/stemi_gse_enrichr.Rmd*' perform gene set enrichment analysis on the differentially expressed gene

'*differential_expression/stemi_mast_to_excel.Rmd*' add DE results to Excel sheets

### cell-cell interactions

'*cell-cell_interactions/stemi_nichenet_analysis.Rmd*' perform cell-cell interaction analysis using Nichenet

'*cell-cell_interactions/stemi_nichenet_limma_input.R*' perform cell-cell interaction analysis using Nichenet, with the limma output

'*cell-cell_interactions/stemi_plot_nichenet_connections*' visualize and investigate the L-R and L-T interactions nominated in Nichenet

### differential protein analysis

'*stemi_differential_protein_expression.Rmd*' look at differential protein expression across timepoints

### qtl mapping

'*qtl_mapping/stemi_create_cis_configs_and_jobs_all_aragam2022confine.sh*' create SLURM jobs to do pQTL mapping for the STEMI samples, using only the GWAS significant SNPs of the Aragam 2022 study

'*qtl_mapping/stemi_create_cis_configs_and_jobs_lld_aragam2022confine.sh*' create SLURM jobs to do pQTL mapping for the LL samples, using only the GWAS significant SNPs of the Aragam 2022 study

'*qtl_mapping/stemi_map_interaction_qtls.Rmd*' do interaction analysis of the proteins with a pQTL effect

'*qtl_mapping/stemi_pqtl_to_excel.Rmd*' add pQTL results to Excel sheets

### clinical variable linkage

'*cell_type_composition/stemi_link_peakckmb_to_celltype_proportions.Rmd*' link cell type proportions to peak-ck-mb values

'*differential_expression/stemi_link_peakckmb_to_scrnaseq.Rmd*' link gene expression to peak-ck-mb values

'*protein_analysis/stemi_olink_peakckmb.Rmd*' link protein expression to peak-ck-mb values

### control cohort description

'*cohort_description/stemi_get_lifelines_variables.R*' get metadata table of controls described in table 1

### data upload

'*ega/stemi_create_md5_checksums.sh*' script to create md5 checksums of each sequence file

'*ega/stemi_create_ega_metadataa.Rmd*' create necessary metadata files for EGA upload

'*ega/stemi_create_pairing_file*' create pairing file necessary for EGA upload

'*ega/stemi_encrypt_ega.sh*' encrypt files for EGA upload

### anonymize data

'*preprocessing/stemi_anonimise_and_export_seurat.R*' anonymize sample names

### Heinig replication

'*factors/stemi_factor_replication_preprocess.Rmd*' preprocess scRNA-seq data for Factor analysis using MOFA

'*factors/stemi_factor_replication_mofa.Rmd*' perform factor analysis using MOFA


## License

The code availabe in this repository is available under GPLv2 License:
https://www.gnu.org/licenses/old-licenses/gpl-2.0.html
