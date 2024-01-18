
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Analysis of lncRNAs after acute exercise

## Background

This repository serves to analyze the changes in long non-coding RNAs
(lncRNAs) after an acute bout of aerobic or resistance exercise.

Exercise induces a multitude of physiological changes, and many of these
are downstream of changes in cellular gene expression. Though there has
been a great deal of study into protein-coding genes, there has been
little research into the response of non-coding RNAs in response to
exercise. Here, we seek to address this gap. We examine this problem
along two different dimensions - exercise modality and exercise time.

We first show similarities and differences in lncRNA expression after
different durations and types of exercise.

<img src="images/Exercise_Modality.png" width="100%" />

For a more granular view, we can examine the pairwise venn diagrams
between different durations of exercise, and between resistance and
aerobic exercise.

<img src="images/Aerobic_Exercise.png" width="50%" height="20%" /><img src="images/Resistance_Exercise.png" width="50%" height="20%" /><img src="images/1_Hour_Exercise.png" width="50%" height="20%" /><img src="images/4_Hour_Exercise.png" width="50%" height="20%" />

------------------------------------------------------------------------

We can expand upon the differences in lncRNA expression after
aerobic exercise compared to baseline using a volcano plot.

<img src="images/lncRNA_aerobic_volcano.png" width="100%" />

------------------------------------------------------------------------

We can do the same for resistance exercise.

<img src="images/lncRNA_resistance_volcano.png" width="100%" />

------------------------------------------------------------------------

## Directories

The ‘proc’ folder contains post-processed results and tables. The
‘images’ folder contains generated figures from the manuscript and
repository

## Usage

To use this package, simply run the `main.R` script in R version 4.3.0
or higher. To generate venn diagrams, run the `generateVennDiagrams.R`
script after.

## Attribution

This project is a collaboration between the [Muscle Physiology
Lab](https://skhs.queensu.ca/qmpl/) at Queen’s University and the
[Clarke Laboratory for Quantitative Exercise
Biology](https://www.sfu.ca/clarkelab-bpk.html) at Simon Fraser
University.
