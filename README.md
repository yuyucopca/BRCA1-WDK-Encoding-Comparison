# BRCA1 WDK Encoding Comparison

This repository contains the data and Python scripts used to reproduce the computational comparison of DNA sequence representations performed in our study of a position-weighted encoding (WDK).

The experiment compares three representations of BRCA1 sequence windows:

- Position-weighted DNA encoding (WDK)
- One-Hot Encoding (OHE)
- 1-mer nucleotide frequencies

A Gradient Boosting classifier is used as a common downstream model to compare the representations under the same cross-validation partitions.

