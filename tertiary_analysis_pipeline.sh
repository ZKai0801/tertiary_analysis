#!/usr/bin/bash

#  
#
#
# Input file: raw VCF file generated from variant-calling software
# Output file: detailed annotated Excel file 

input_vcf=$1

# database
reference="/public/database/GATK_Resource_Bundle/hg19/ucsc.hg19.fasta"
vep_dir=""
civic=""
clinic_transcripts=""

# scripts
merge_mnv=""
anno_hgvs=""
anno_civic=""


# step1: filter failed variant 
echo "** Step1: filtering starts at `date`"
grep "PASS\|#" $input_vcf > ${input_vcf%%.vcf*}.step1_filtered.vcf

# step2: normalisation of vcf file
echo "** Step2: normalisation starts at `date`"
bcftools norm -m -both -f ${reference} ${input_vcf%%.vcf*}.step1_filtered.vcf -o ${input_vcf%%.vcf*}.step2_normalised.vcf

# step3: merge multi-nucleotide polymorphism
echo "** Step3: merging starts at `date`"
python3 $merge_mnv ${input_vcf%%.vcf*}.step2_normalised.vcf ${reference} -o ${input_vcf%%.vcf*}.step3_merged.vcf

# step4: annotate population freq, software prediction results etc.
echo "** Step3: annotation starts at `date`"
table_annovar.pl ${input_vcf%%.vcf*}.step3_merged.vcf ${annovar_dir} \
--outfile ${input_vcf%%.vcf*}.step4_temp --buildver hg19  \
--protocol 1000g2015aug_all,esp6500siv2_all,exac03nonpsych,cosmic70,clinvar_20191013,intervar_20180118,dbnsfp35a,avsnp150,dbscsnv11  \
--operation f,f,f,f,f,f,f,f,f  --vcfinput --remove

vep --dir ${vep_dir} --cache --offline --cache_version 98 \
--use_given_ref --refseq --assembly GRCh37  \
--fa ${reference} --force_overwrite --vcf --variant_class --gene_phenotype --vcf_info_field ANN \
--hgvs --hgvsg --transcript_version \
-i ${input_vcf%%.vcf*}.step4_temp.hg19_multianno.vcf -o ${input_vcf%%.vcf*}.step4_anno.vcf

rm ${input_vcf%%.vcf*}.step4_temp.hg19_multianno.vcf ${input_vcf%%.vcf*}.step4_temp.hg19_multianno.txt

# step5: select for correct transcript
python3 $anno_hgvs ${input_vcf%%.vcf*}.step4_anno.vcf ${clinic_transcripts} ${refFlat} -o ${input_vcf%%.vcf*}.step5_hgvs.vcf

# step6: annotate clinical significant associated with each variant
python3 $anno_civic ${input_vcf%%.vcf*}.step5_hgvs.vcf ${civic} -o ${input_vcf%%.vcf*}.step6_civic.vcf



