# macs_perturbseq

Analysis code for Macs Perturb-seq.

## Analysis workflow

### 1. Ingestion from cellranger outputs
see
|-`basic_preprocessing_macs.py`
  |-`preprocess_macs.py`
  |-[src/preprocess.py](https://github.com/mangochiral/CRISPRa_Analysis_pipeline.git)
scripts 

Input files: 
- file1: sample_filter_feature_matrix.h5d
- file2: <description>

Output files:
- file1: `<sample_name>_gex_preprocessed.h5ad`
- file2: `<sample_name>_crispr_preprocessed.h5ad`

### 2. Guide assignment
see `script_name` (you can also add subtasks)
[add details on key substebs happening in the script]

Input files: 
- file1: <description> 
- file2: <description>

Output files:
- file1: <description>
- file2: <description>

### 3. Knockdown efficiency analysis
see `script_name`
[add details on key substebs happening in the script]

Input files: 
- file1: <description> 
- file2: <description>

Output files:
- file1: <description>
- file2: <description>

### 4. Comparison with FACS based screens

see `TNF_DE_vs_FACS.ipynb` - here we get measured expression of TNF in each cell and run a simple Wilcoxon test to assess which perturbations affect expression of TNF. The estimated logFC is compared to Mageck LFC in TNF protein expression from the FACS-based screen.

Input files: 
- `{exp_name}/Mac_{lane_id}_gex_multiguide.h5ad`: single cell counts data for each lane and experiment
- `metadata/Paired_hi_vs_low_2026-01-09.gene_summary.txt`: summary stats from FACS screen (mageck outputs)

Output files:
- `results/DE_t_on_TNF.t_vs_others.csv`: effects of perturbations on TNF expression (testing for cells with each perturbed gene vs rest of cells)
- `results/DE_t_on_TNF.t_vs_singleNTC.csv`: effects of perturbations on TNF expression (testing for cells with each perturbed gene vs cells with single NTC guide)

### 5. Selection of single-guide cells 

### 6. DESeq2 analysis

#### 6.1 Pseudobulking

#### 6.2 Filtering perturbed genes and measured genes

#### 6.3 DESeq2 test

