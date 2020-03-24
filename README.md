# Tertiary Analysis on tumour

This pipeline involves:
1. normalisation of vcf file (split multiallelic sites; left-align indels)
2. remove non-pass variants
3. annotating rest of variants (including population AF, prediction values from various software, CIVic and etc...)
4. remove benign & likely benign variants

## Dependencies

Please make sure following software were installed and added to your environment. 

- BCFtools
- Annovar
- Ensembl-VEP

Also, you need to download following database in order to annotate

- Reference Fasta file
- VEP database 
- refFlat 
- CIVic (ClinicalEvidenceSummaries.tsv)
- gnomad mnv database
- clinic_transcript.tsv (generated from parse_lrg.py)

## Usage

When you first time use this pipeline, you need to open the `tertiary_analysis_pipeline.sh` and type in paths for required database or scripts. Then you are ready to run:

```bash
bash tertiary_analysis_pipeline.sh <vcf>
```



