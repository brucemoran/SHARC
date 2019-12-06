# SHARC NextFlow Pipeline Rewrite

## Intro
SHARC is a pipeline for somatic SV calling and filtering from tumor-only Nanopore sequencing data. It performs mapping, SV calling, SV filtering, random forest classification, blacklist filtering and SV prioritization, followed by automated primer design for PCR amplicons of 80-120 bp that are useful to track cancer ctDNA molecules in liquid biopsies.

##Issues
```
1: Is splitting VCFs necessary?
   Is this just to get mean coverages for each Chr?
```
