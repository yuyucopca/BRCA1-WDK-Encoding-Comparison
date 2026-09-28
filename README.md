# BRCA1 WDK Encoding Comparison

This repository contains the sequence data and Python scripts used to reproduce the computational comparison of DNA sequence representations performed in our study of a position-weighted DNA encoding (WDK).

The experiment compares three representations of BRCA1 sequence windows:

- Position-weighted DNA encoding (WDK)
- One-Hot Encoding (OHE)
- 1-mer nucleotide frequencies

A Gradient Boosting classifier is used as a common downstream model to evaluate the three representations under identical cross-validation partitions.

The purpose of this experiment is to compare the computational behavior of the representations under controlled conditions. It is not intended to establish clinical pathogenicity prediction performance.

---

## Repository Structure

```text
BRCA1-WDK-Encoding-Comparison/
│
├── README.md
├── requirements.txt
│
├── code/
│   ├── EncodingNucleotidsV2.py
│   └── EncodingComparison_WDK_OHE_1mer_v2.py
│
└── data/
    ├── BRCA1_SNV_emparejadas_162.xlsx
    │
    ├── SeqBRCA1_Benigns/
    │   └── sequence files
    │
    └── SeqBRCA1_Pathogenics/
        └── sequence files
```

---

## Dataset Construction

The experimental dataset consists of **162 BRCA1 single-nucleotide variant (SNV) samples**, balanced between the two analyzed classes:

- **81 benign variants**
- **81 pathogenic variants**
- **162 samples in total**

Control sequences were excluded from the final comparison.

### Variant Identification

BRCA1 variants were first identified using the **NCBI Variation Viewer**. Benign and pathogenic SNVs were selected as the variants of interest for the experiment.

### Sequence Retrieval

After identifying the variants of interest, the corresponding genomic regions were located using the **Ensembl Genome Browser**.

The nucleotide sequences associated with these regions were then retrieved and downloaded as individual sequence files.

The resulting dataset contains sequence windows of **variable length**, reflecting the regions of interest retrieved for the selected variants.

The complete experimental workflow for data construction can therefore be summarized as:

```text
NCBI Variation Viewer
        ↓
Identification of BRCA1 SNVs
        ↓
Selection of benign and pathogenic variants
        ↓
Localization of the corresponding regions in Ensembl
        ↓
Retrieval of nucleotide sequence windows
        ↓
Balanced experimental subset
        ↓
81 benign + 81 pathogenic sequences
```

---

## Balanced Experimental Subset

An equal number of benign and pathogenic samples was used to avoid class imbalance in the representation comparison.

The final subset contains:

| Class | Number of samples |
|---|---:|
| Benign | 81 |
| Pathogenic | 81 |
| **Total** | **162** |

The samples included in the final experiment are documented in:

```text
data/BRCA1_SNV_emparejadas_162.xlsx
```

This file contains the metadata required to associate each selected sample with its corresponding nucleotide sequence, including:

- file name
- class
- chromosome
- start position
- end position
- sequence length

The corresponding nucleotide sequences are stored in:

```text
data/SeqBRCA1_Benigns/
data/SeqBRCA1_Pathogenics/
```

The sequence filenames are preserved because they are used to associate each entry in the selected-sample dataset with its corresponding sequence file.

---

## WDK Representation

The WDK representation uses a position-dependent weighting function defined as:

\[
w(i)=2^{-\frac{1}{2(1+i)}}
\]

where \(i\) represents the nucleotide position in the sequence.

The implementation uses positions starting at:

\[
i=1
\]

Each nucleotide contributes its positional weight to one of four nucleotide-specific components:

\[
E=(v_A,v_C,v_G,v_T)
\]

where each component contains the accumulated positional contributions associated with the corresponding nucleotide.

For the machine-learning comparison, the four WDK components are normalized by the actual sequence length.

The resulting representation therefore contains **four features per sequence**, independently of the original sequence length.

---

## Compared Representations

Exactly the same set of 162 nucleotide sequences is used to generate all three representations.

### WDK

The proposed position-weighted representation produces a four-dimensional vector:

```text
[vA, vC, vG, vT]
```

Positional information is incorporated through the weighting function \(w(i)\).

### One-Hot Encoding (OHE)

Each nucleotide is represented by four binary components:

```text
A = [1, 0, 0, 0]
C = [0, 1, 0, 0]
G = [0, 0, 1, 0]
T = [0, 0, 0, 1]
```

Because the analyzed sequence windows have variable lengths, zero-padding is applied to obtain a common dimensionality.

The resulting matrix is then flattened before being provided to the classifier.

### 1-mer

The 1-mer representation contains the normalized nucleotide frequencies:

```text
[freq_A, freq_C, freq_G, freq_T]
```

This representation also produces four features per sequence but does not explicitly encode nucleotide position.

---

## Gradient Boosting Comparison

The three representations are evaluated using the same machine-learning procedure.

The experiment uses:

- `GradientBoostingClassifier`
- 5-fold cross-validation
- `StratifiedGroupKFold`
- fixed random seed: `42`

The same cross-validation partitions are generated once and subsequently used for WDK, OHE, and 1-mer.

This ensures that differences between the representations are not caused by different train/test partitions.

### Genomic Grouping

To reduce information leakage between training and test partitions, genomic groups are constructed using:

```text
chromosome
start_position
end_position
sequence_length
```

Samples belonging to the same genomic window are kept within the same group during cross-validation.

The implementation verifies that genomic groups assigned to the training partition do not appear in the corresponding test partition.

Importantly, genomic coordinates and sequence length are used for **grouping and experimental control only**.

They are **not provided to Gradient Boosting as predictive features**.

The classifier receives only the features generated by the corresponding sequence representation.

---

## Evaluation Metrics

The comparison reports predictive and computational characteristics for each representation.

The machine-learning metrics are:

- Accuracy
- F1-score
- Cohen's Kappa

The computational comparison additionally records:

- Number of features
- Memory usage
- Memory per sample
- Training time
- Prediction time
- Total cross-validation time

These measurements are used to characterize the behavior of each representation under the same computational pipeline.

The Gradient Boosting experiment is included as a controlled comparison of sequence representations and should not be interpreted as evidence of clinical pathogenicity prediction capability.

---

## Scripts

### `EncodingNucleotidsV2.py`

This script processes the original sequence files and generates the WDK-encoded dataset.

Its main operations include:

- reading sequence files
- validating nucleotide content
- extracting genomic coordinates from sequence headers
- generating One-Hot Encoding
- applying the WDK positional weighting function
- generating the four WDK components
- storing genomic metadata
- excluding control sequences
- exporting the encoded dataset

### `EncodingComparison_WDK_OHE_1mer_v2.py`

This script performs the final representation comparison.

It:

1. Loads the selected balanced dataset.
2. Retrieves the corresponding nucleotide sequences.
3. Validates sequence lengths.
4. Recalculates the sequence representations.
5. Generates WDK, OHE, and 1-mer features.
6. Creates grouped and stratified cross-validation partitions.
7. Uses identical partitions for all representations.
8. Trains a Gradient Boosting classifier.
9. Calculates evaluation metrics.
10. Records computational characteristics.
11. Exports the comparison results.

---

## Requirements

The experiments were conducted using:

```text
Python 3.10.9
NumPy 1.26.4
pandas 1.5.3
scikit-learn 1.2.1
```

The scripts also require `openpyxl` for Excel file handling.

A minimal `requirements.txt` is:

```text
numpy==1.26.4
pandas==1.5.3
scikit-learn==1.2.1
openpyxl
```

Dependencies can be installed using:

```bash
pip install -r requirements.txt
```

---

## Reproducibility

The repository provides:

- the selected BRCA1 sequence samples,
- the balanced sample metadata,
- the WDK encoding implementation, and
- the script used for the WDK/OHE/1-mer Gradient Boosting comparison.

A fixed random seed (`42`) is used for both cross-validation generation and the Gradient Boosting classifier.

The same sequence samples and cross-validation partitions are used for all three representations.

Local file paths may need to be adjusted depending on the location where the repository is cloned.

---

## Scope

This repository accompanies a study investigating the mathematical and computational properties of a position-weighted DNA representation.

The broader study examines properties such as positional weighting, numerical behavior, dimensionality, computational cost, and the behavior of the representation relative to conventional DNA encoding approaches.

The machine-learning experiment contained in this repository serves specifically as a controlled comparison between WDK, OHE, and 1-mer representations.

It is **not intended as a clinical diagnostic tool or as evidence of clinical pathogenicity prediction performance**.

---

## Data Sources

The genomic resources used during dataset construction were:

- **NCBI Variation Viewer** — used to identify BRCA1 variants of interest.
- **Ensembl Genome Browser** — used to locate the corresponding genomic regions and retrieve the nucleotide sequences used in the experiment.

The processed subset and corresponding sequence files used in the final comparison are provided in the `data/` directory for reproducibility.

---

## Citation

If you use the code or data provided in this repository, please cite the associated publication.

Citation information will be added after publication.

---

## License

See the repository license for information regarding reuse of the source code.
