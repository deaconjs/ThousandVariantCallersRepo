### Swi9 
[2017 Study](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1417-7) 

compares 9 SNV: deepSNV, GATK HaplotypeCaller, GATK UnifiedGenotyper, JointSNVMix2, MuTect, SAMtools, SiNVICT, SomaticSniper, and VarScan2 and does parameter selection on some of the callers

 * Virmid - top 4. similar to strelka, best in deep seq
 * Strelka - top 4. similar to virmid
 * EBCall - top 4. best in deep sequencing
 * Mutect - top 4
 * Seurat - largest # calls. better for high sequencing depths
 * Sniper - largest # calls, not as good in deep seq
 * Shimmer, Varscan, DeepSnV - mediocre to poor performance

### Gor4 
[2016 Study](https://www.nature.com/articles/srep36540) 

compares 4 SNV: Varscan, SomaticSniper, Strelka and MuTect2 

 * Strelka - top 2
 * Shearwater - top 2

### Den9 
[2016 Study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4803342/)

compares 9 SNV/indel callers: EBCall, Mutect, Seurat, Shimmer, Indelocator, Somatic Sniper, Strelka, VarScan 2 and Virmid 

 * JointSNVMix - 2 most sens to LFV. 2 most unique FPs
 * DeepSNV - 2 most sens to LFV. 2 most unique FPs. always does better with increased depth. germline filter improves LFV
 * Mutect - germline filter improves LFV
 * Varscan2 - use minvarfreq 0.02 for LVF


### Bcb8 
[2015](http://bcb.io/2015/03/05/cancerval/) 

compares 8 variant callers: MuTect, VarDict, FreeBayes, Varscan, Scalpel, LUMPY, DELLY, WHAM

 * Mutect performs best

[evaluated](https://bcbio-nextgen.readthedocs.io/en/latest/contents/testing.html#cancer-tumor-normal) cancer tumor/normal variant calling with [synthetic dataset 3](https://www.synapse.org/#!Synapse:syn312572/wiki/62018) from the DREAM challenge, using multiple approaches to detect SNPs, indels and structural variants 


### Bcb3 
2014
compares 3 SV variant callers: delly, lumpy, cm.mops


### Bro5 
[2014 Study](https://academic.oup.com/bioinformatics/article/30/20/2843/2422145) 
compares 5 haplotyping variant callers with different mappers: GATK HC, GATK UG, Platypus, FreeBayes, SAMTools


### Qia5 
[2014 Study](https://www.ncbi.nlm.nih.gov/pubmed/24678773) 
compares 5 somatic mutation callers: GATK UG with subtraction, MuTect, Strelka, SomaticSniper and VarScan2 

 * Strelka - top 2 sens
 * Mutect - top 2 sens


### Aus4 
[2013 Study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3753564/) 
compares 4 variant callers on cancer/normal exome somatic tumor data: VarScan, SomaticSniper, JointSNVMix and Strelka 

 * JointSNVMix - inconvenient, little benefit. can have germline FPs at high depths
 * Strelka - highest TP pass rate through custom filters, lowest germline SNP FP (dbSNP presence)
 * Sniper - credible with low and medium MAF
 * Varscan2 - returns lots of germline FPs. MAF 20-75%. doesn't add anything useful


### Van6 
[2013 Study](https://genomemedicine.biomedcentral.com/articles/10.1186/gm495) 
compares 6 SNV/indel callers: JointSNVMix, SomaticSniper, Strelka, and VarScan 2, MuTect and EBCall

 * Mutect - more sensitive for LAF variants
 * Varscan2 - more sensitive for high allele frequency variants
 * EBCall - unmatched error distribution between normal references and target samples


### CSH5 
[2013 Study](https://genomemedicine.biomedcentral.com/articles/10.1186/gm432) 
compares 5 pipelines: SOAP, BWA-GATK, BWA-SNVer, GNUMAP, and BWA-SAMtools 

### Bcb3 
variant caller [comparison](https://bcbio.wordpress.com/2013/10/21/updated-comparison-of-variant-detection-methods-ensemble-freebayes-and-minimal-bam-preparation-pipelines/) (2013): 
compares 3 pipelines: GATK UG, GATK HC both with and without Picard and GATK adjustments, and FreeBayes
