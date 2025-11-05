# 🧬 Genomic String Matching

This repository implements two DNA sequence matching algorithms in Python:

1. **Naïve Matching with up to 2 Mismatches** — alignment algorithm that tolerates sequencing errors or mutations.  
2. **Boyer–Moore Matching with Comparison Counting** — string search algorithm that reports the number of character comparisons and alignments performed.

The above aim to demonstrate computational biology concepts used in genome analysis, such as **read alignment**, **approximate matching**, and **algorithmic performance evaluation**.

---

## 📖 Table of Contents

- [Overview](#-overview)
- [Project Structure](#-project-structure)
- [Algorithms](#-algorithms)
  - [Naïve 2-Mismatch Matcher](#1-naïve-2-mismatch-matcher)
  - [Boyer–Moore with Counts](#2-boyer–moore-with-counts)
- [Usage](#-usage)
- [Example Output](#-example-output)
- [Dependencies](#-dependencies)
- [License](#-license)
- [Author](#-author)

---

## 🧠 Overview

Modern sequencing technologies produce millions of short DNA fragments ("reads") that must be aligned to a **reference genome**.  
This project provides two simple but educational approaches for this problem:

- A **naïve matcher** that finds reads in a genome, allowing up to **2 mismatches** per read.
- A **Boyer–Moore matcher** that tracks **algorithmic efficiency** by counting alignments and character comparisons.

---

## 📁 Project Structure


> Large genome files are not included; small demo files can be added to the `data/` folder.

---

## 🚀 Usage

### Run Naïve 2-Mismatch Matcher
```bash
python naive_2mm.py data/myReads.fastq data/lambda_virus.fa


