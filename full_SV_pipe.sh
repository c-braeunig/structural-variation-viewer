#!/bin/bash

#assign in/out and params
infasta=../raw_data/OEVEW1_p1p2_csa.fasta
windowsize=10000
stepsize=500
outbed=../OEVEW1_p1p2_csa.fasta.bed
outfasta=../OEVEW1_p1p2_csa.LR
reference=../raw_data/OEVE_p1p2_csa.fasta
outtemp=../OEVEW1vOEVE


#index infasta file
samtools faidx $infasta

#generate bedfile for simulated long reads
python makebed.py --fi $infasta.fai --sz $windowsize --st $stepsize --bed $outbed

#use bedtools to generate long read fasta file w/ bed file and infasta
bedtools getfasta -fi $infasta -bed $outbed -fo $outfasta

#align reads against reference w/ NGMLR
ngmlr -t 4 -r $reference -q $outfasta -o $outtemp.sam

#delete intermediates
rm $infasta.fai
rm $outbed
rm $outfasta

#convert sam to bam, sort
samtools view -S -b $outtemp.sam > $outtemp.bam
samtools sort $outtemp.bam -o $outtemp.sorted.bam

#variant calling w/ sniffles
sniffles -m $outtemp.sorted.bam -v $outtemp.vcf

#delete intermediate files
rm $outtemp.sam
rm $outtemp.bam
rm $outtemp.sorted.bam
