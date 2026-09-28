Dataset and Sample Selection

The experiment was conducted using a balanced subset of 162 BRCA1 SNV sequence samples, consisting of:

- 81 benign samples
- 81 pathogenic samples
- 162 samples in total

The samples were selected from the BRCA1 sequence dataset used in the study. Only single-nucleotide variants (SNVs) were considered for this experiment, and control sequences were excluded.

To avoid evaluating the representations on a class-imbalanced dataset, an equal number of benign and pathogenic samples was selected. The resulting balanced subset is documented in:

BRCA1_SNV_emparejadas_162.xlsx

This file contains the filename, class label, genomic coordinates, and sequence length associated with each selected sample.

The corresponding nucleotide sequences are provided in:

SeqBRCA1_Benigns/
SeqBRCA1_Pathogenics/

Each entry in `BRCA1_SNV_emparejadas_162.xlsx` is linked to its original sequence file through the `file_name` field. The comparison script uses this field to recover the sequence directly from the corresponding class directory and verifies that its length matches the metadata stored in the dataset.

Experimental Control

The same 162 sequences are used to generate all three representations evaluated in the experiment:

- WDK
- One-Hot Encoding (OHE)
- 1-mer nucleotide frequencies

Therefore, differences observed during the Gradient Boosting comparison cannot be attributed to different samples being used for each encoding.

In addition to class balancing, genomic windows are grouped using chromosome, start position, end position, and sequence length during cross-validation. Samples belonging to the same genomic window are prevented from appearing simultaneously in the training and test partitions.

Genomic coordinates and sequence length are used only for dataset organization and cross-validation grouping; they are not included as predictive features in the WDK, OHE, or 1-mer representations.

Data Source

The BRCA1 variants and sequence windows were obtained from publicly available genomic resources used in the study.

**[Exact source/database and genome assembly to be specified.]**

