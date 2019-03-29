# Assignment 1 "Project"


# Aim

This project aims to practice most of the technologies and software you've learned so far the course, starting from data download to annotation and differential expression.

# Important notes

- For each step, document your workflow and commands to reproduce the same expirement.
- **Everything** must be reported and documented on Github.
- Contribution

---

# Project details

## 1- Data Download

Download ~25M Human RNA-Seq reads from the SRA.

## 2- Prepare the data

- Construct 5 Samples, each sample contains two parts
    - **Part 1:** 1M unique reads from the main reads file.
    - **Part 2:** Shuffle the main reads file, and take random 1M Reads.
    - **Example:** Sample 1 will be divided into S1_1 & S1_2, Unique and Random selected after shuffeling respectively.

## 3- FASTQ Quality Control

- For **Sample 1** only, use FASTQC to report the difference between between S1_1 and S1_2

## 4- Trimming

- For **Sample 1** only, Apply:
    - Mild Trimming for S1_1.
    - Aggressive Trimming for S1_2

## 5- Alignment

- Align all the samples (1:5) using **BWA** and **Hisat** against the human reference file.
- Export some useful statistics report for **each sample indvidually**.

## 6- Assembly

- Apply reference-based trasncriptome assembly using stringTie.
    -  **Step1.** For the 5 samples (without shuffeling) parts.
    -  **Step2.** For the 5 samples (with shuffeling applied) parts.

## 7- Apply Differential Expression.

## 8- Using GTF-Compare Compare the reference and the generated annotation files.

