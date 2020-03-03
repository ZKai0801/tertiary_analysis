# Tertiary Analysis on tumour

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
- clinic_transcript.tsv (generated from parse_lrg.py)

## Usage

When you first time use this pipeline, you need to open the `tertiary_analysis_pipeline.sh` and type in paths for required database or scripts. Then you are ready to run:

```bash
bash tertiary_analysis_pipeline.sh <vcf>
```



