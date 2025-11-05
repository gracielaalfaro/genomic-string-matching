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

> ⚠️ Large genome files are excluded from the repository — only small test files are included for demonstration purposes.

---

## 🔬 Algorithms

### 1️⃣ Naïve 2-Mismatch Matcher
**File:** `naive_2mm.py`

This algorithm scans the reference genome and compares each substring of equal length to the read, allowing up to **2 mismatches**.  
Only reads where **all base quality scores > 5** are considered.

**Output:**  
For each read that passes the quality filter, the program prints:
and reports how many reads were filtered out due to low quality.

---

### 2️⃣ Boyer–Moore with Counts
**File:** `bm_with_counts.py`

An extension of the classic **Boyer–Moore** string search algorithm.  
It returns:
- The list of **match positions**
- The number of **alignments** attempted
- The number of **character comparisons**

**Example:**
```python
p = 'word'
t = 'there would have been a time for such a word'

lowercase_alphabet = 'abcdefghijklmnopqrstuvwxyz '
p_bm = BoyerMoore(p, lowercase_alphabet)
occurrences, n_alignments, n_comparisons = boyer_moore_with_counts(p, p_bm, t)

print(occurrences, n_alignments, n_comparisons)
# Output: [40], 12, 15
GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG
python naive_2mm.py data/myReads.fastq data/lambda_virus.fa
python bm_with_counts.py
Read_1 found at positions [1023, 29810]
Read_2 found at positions [1045]
Read_3 did not pass quality filter
---
Total reads processed: 3
Reads filtered out (quality ≤ 5): 1
Occurrences: [40]
Alignments tried: 12
Character comparisons: 15

