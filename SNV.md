## SNP Variant Callers

|caller|pubyear|from|study|source|algorithm|tictac|
|------|-------|----|-----|------|---------|------|
|[graphtyper](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#graphtyper)|2017|deCODE genetics|[study](https://www.nature.com/articles/ng.3964)|[source](https://github.com/DecodeGenetics/graphtyper)|Population-scale genotyping using pangenome graphs|
|[muse](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#muse)|2016|MD Anderson Cancer Center|[study](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1029-6)|[source](http://bioinformatics.mdanderson.org/main/MuSE)|F81 Markov Substitution Model|
|[sinvict](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#sinvict)|2016|Simon Frasiser University, Canada|[study](https://academic.oup.com/bioinformatics/article-abstract/33/1/26/2755714/SiNVICT-ultra-sensitive-detection-of-single)|[source](https://sfu-compbio.github.io/sinvict/)|
|[multigems](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#multigems)|2016|University of California, Riverside|[study](https://academic.oup.com/bioinformatics/article/32/10/1486/1742567/MultiGeMS-detection-of-SNVs-from-multiple-samples)|[source](https://github.com/cui-lab/multigems)|Multinomial Bayesian, base and alignment quality priors|
|[somaticseq](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#somaticseq)|2015|Roche Bina|[study](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0758-2)|[source](http://bioinform.github.io/somaticseq/)|meta-caller, decision tree|
|[discosnp](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#discosnp)|2015|Genscale France|[study](https://www.ncbi.nlm.nih.gov/pubmed/25404127)|[source](http://colibread.inria.fr/software/discosnp/)|reference-free, de bruijn graph|
|[2kplus2](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#2kplus2)|2015|Norwich Research Park, UK, Sainsbury lab|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4341063/)|[source](https://github.com/danmaclean/2kplus2)|reference-free, de bruijn graph|
|[exscalibur](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#exscalibur)|2015|University of Chicago|[study](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0135800)|[source](https://github.com/cribioinfo)||
|[multisnv](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#multisnv)|2015|Cambridge Tavare|[study](https://academic.oup.com/nar/article/43/9/e61/1113203/multiSNV-a-probabilistic-approach-for-improving)|[source](https://bitbucket.org/joseph07/multisnv/wiki/Home)|joint paired, timepoint pooling|
|[rarevator](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#rarevator)|2015|University of Florence|[study](http://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-1481-9)|[source](https://sourceforge.net/projects/rarevator/)|Fisher's exact test, conserved loci only|
|[snv-ppilp](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#snv-ppilp)|2015|University of Helsinki, Finland|[study](https://academic.oup.com/bioinformatics/article/31/7/1133/180437/SNV-PPILP-refined-SNV-calling-for-tumor-data-using)|[source](https://www.cs.helsinki.fi/en/gsa/snv-ppilp/)|perfect phylogeny/integer linear programming|
|[platypus](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#platypus)|2014|U Oxford|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4753679/)|[source](https://github.com/andyrimmer/Platypus)|Haplotype, bayesian, multi-sample, local realignment|
|[baysic](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#baysic)|2014|Baylor/Genformatic LLC|[study](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-104)|[source](http://genformatic.com/baysic/)|Meta-caller, Bayesian, unsupervised|
|[hapmuc](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#hapmuc)|2014|Kyoto University, Japan|[study](https://www.ncbi.nlm.nih.gov/pubmed/25123903)|[source](https://github.com/usuyama/hapmuc)|Haplotype, Bayesian HMM|
|[snpest](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#snpest)|2014|U Copenhagen|[study](http://bmcresnotes.biomedcentral.com/articles/10.1186/1756-0500-7-698)|[source](https://github.com/slindgreen/SNPest)|reference-free, generative probabilistic|
|[variantmaster](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#variantmaster)|2014|Geneva Medical School, Switzerland|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3912425/)|[source](https://sourceforge.net/projects/variantmaster/)|reference-free, pedigree inference|
|[mutect](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#mutect)|2013|Broad Getz|[study](https://www.ncbi.nlm.nih.gov/pubmed/23396013)|[source](https://github.com/broadinstitute/mutect)|Beta-binomial, Variable Allele Fraction, filter population SNPs|
|[niks](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#niks)|2013|Max Planck Institute for Plant Breeding Research, Germany|[study](https://www.nature.com/articles/nbt.2515)|[source](https://sourceforge.net/projects/niks/)|
|[ebcall](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#ebcall)|2013|Vanderbilt Zhao|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3627598/)|[source](https://github.com/friend1ws/EBCall)|Heuristic, multiple feature|
|[shearwater](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#shearwater)|2013|U Cambridge/Welcome Trust|[study](http://bioinformatics.oxfordjournals.org/content/early/2014/01/31/bioinformatics.btt750.full)|[source](https://bioconductor.org/packages/release/bioc/html/deepSNV.html)|Beta-binomial, DeepSNV with aggregate control counts|
|[shimmer](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#shimmer)|2013|NHGRI Larsen|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3673219/)|[source](https://github.com/nhansen/Shimmer)|Fisher's exact test, variant read count > N|
|[bubbleparse](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#bubbleparse)|2013|Norwich Research Park Sainsbury Lab, UK|[study](https://www.ncbi.nlm.nih.gov/pubmed/23536903/)|[source](https://github.com/richardmleggett/bubbleparse)|Reference-free, de Bruijn graph|
|[cake](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#cake)|2013|Welcome Trust Adams|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3740632/)|[source](http://cakesomatic.sourceforge.net/)|Meta-caller, simple 2x consensus, post-filter|
|[denovogear](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#denovogear)|2013|WashU St Louis Conrad|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4003501/)|[source](https://github.com/denovogear/denovogear)|Beta-binomial, pedigree|
|[qsnp](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#qsnp)|2013|U Queensland|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3826759/)|[source](https://sourceforge.net/p/adamajava/wiki/qSNP%201.0/)|Heuristic, min 3 reads, post-filter|
|[rvd](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#rvd)|2013|Stanford University School of Medicine|[study](https://bmcresnotes.biomedcentral.com/articles/10.1186/1756-0500-6-206)|[source](http://dna-discovery.stanford.edu/software/rvd/)|Beta-binomial|
|[seurat](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#seurat)|2013|Translational Genomics Research Institute|[study](http://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-14-302)|[source](https://sites.google.com/site/seuratsomatic/)|Joint-paired, beta-binomial|
|[snptools](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#snptools)|2013|Baylor College of Medicine|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3638139/)|[source](https://sourceforge.net/projects/snptools/)|Haplotype, Bayesian HMM|
|[vcmm](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#vcmm)|2013|RIKEN Japan|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3703611/)|[source](http://www.mybiosoftware.com/vcmm-variant-caller-with-multinomial-probabilistic-model.html)|Multinomial Bayesian, priors corrected Illumina q-score|
|[vip](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#vip)|2013|Case Western, Li lab|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530907/)|[source](http://cbc.case.edu/VIP/)|Overlapping Pools|
|[virmid](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#virmid)|2013|UCSD Bafna|[study](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-8-r90)|[source](https://sourceforge.net/p/virmid/wiki/Home/)|Joint-paired, Beta-binomial, purity estimation|
|[varscan2](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#varscan2)|2012|WashU St Louis Wilson|[study](https://www.ncbi.nlm.nih.gov/pubmed/22300766)|[source](http://varscan.sourceforge.net/)|Heuristic, min 3 reads, filter|
|[jointsnvmix](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#jointsnvmix)|2012|U British Columbia Vancouver|[study](https://www.ncbi.nlm.nih.gov/pubmed/22285562)|[source](https://code.google.com/archive/p/joint-snv-mix/)|Joint-paired, Beta-binomial|
|[lofreq](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#lofreq)|2012|Genome Institute of Singapore|[study](http://nar.oxfordjournals.org/content/40/22/11189.long)|[source](http://csb5.github.io/lofreq/)|Joint-paired, Poisson-binomial|
|[strelka](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#strelka)|2012|Illumina|[study](http://bioinformatics.oxfordjournals.org/content/28/14/1811)|[source](https://sites.google.com/site/strelkasomaticvariantcaller/home/)|Joint-paired, multinomial Bayesian|
|[atlas2](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#atlas2)|2012|Baylor Yu|[study](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-8)|[source](https://www.hgsc.bcm.edu/software/atlas-2)|Heuristic, reads ratio plus filter|
|[conan-snv](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#conan-snv)|2012|British Columbia Cancer Agency, Canada|[study](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0041551)|[source](http://compbio.bccrc.ca/software/conan-snv/)|CNV-informed SNV calls|
|[cortex](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#cortex)|2012|Welcome Trust, University of Oxford, UK|[study](https://www.ncbi.nlm.nih.gov/pubmed/22231483)|[source](http://cortexassembler.sourceforge.net/index_cortex_var.html)|Reference-free, de Bruijn|
|[deepsnv](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#deepsnv)|2012|ETH Zurich|[study](http://www.nature.com/articles/ncomms1814)|[source](https://bioconductor.org/packages/release/bioc/html/deepSNV.html)|Probabilistic beta-binomial|
|[gems](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#gems)|2012|UCal Riverside|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3338331/)|[source](https://github.com/cui-lab/multigems)|Multinomial Bayesian, base- and alignment-quality priors|
|[impute2](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#impute2)|2012|University of Chicago|[study](https://www.ncbi.nlm.nih.gov/pubmed/22820512)|[source](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html)|Haplotype|
|[somatic_sniper](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#somatic_sniper)|2011|Wash U St Louis Ding|[study](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btr665)|[source](https://github.com/genome/somatic-sniper)|Joint-paired, multinomial Bayesian|
|[bambino](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#bambino)|2011|"NCI L Population Genetics Buetow"|[study](https://www.ncbi.nlm.nih.gov/pubmed/21278191)|[source](https://github.com/NCIP/cgr-bambino)|Heuristic, multiple features|
|[freebayes](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#freebayes)|2011|Erik Garrison|[study](https://arxiv.org/abs/1207.3907)|[source](https://github.com/ekg/freebayes)|Haplotype, multi-allelic, non-uniform copy number|
|[mutationseq](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#mutationseq)|2011|"U British Columbia Vancouver Shah"|[study](http://bioinformatics.oxfordjournals.org/content/28/2/167.abstract?keytype=ref&ijkey=oj0Wpkhils4hmyC)|[source](http://compbio.bccrc.ca/software/mutationseq/)|Machine learning|
|[snver](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#snver)|2011|New Jersey I of T|[study](https://www.ncbi.nlm.nih.gov/pubmed/21813454/)|[source](http://snver.sourceforge.net/)|Overlapping pools|
|[syzygy](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#syzygy)|2011|broad|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3378381/)|[source](http://software.broadinstitute.org/software/syzygy/home)|Probabilistic, strand, sequence context, neighborhood quality score priors|
|[vipr](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#vipr)|2011|Max Plank Institute|[study](http://bioinformatics.oxfordjournals.org/content/27/13/i77.full)|[source](https://sourceforge.net/projects/htsvipr/files/vipR/)|Overlapping pools|
|[crisp](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#crisp)|2010|Scripp's Translational Science Institute, US|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2881398/)|[source](https://sites.google.com/site/vibansal/software/crisp)|Overlapping pools|
|[indelocator](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#indelocator)||CGA/Broad||[source](http://archive.broadinstitute.org/cancer/cga/indelocator)||

### graphtyper

Validated vs: GATK (UG, UGLite, HC, HC joint), Samtools, Platypus, FreeBayes

Used by: deCODE genetics on WGS of > 28,000 Icelanders

Notes: fast, highly scalable, includes HLA typing

Algorithm: Iterative creation of local pangenome graphs for accurate re-alignment of reads to all possible haplotypes.

### muse

Validated vs: mutect, sniper, strelka

Used by: GDC, SomaticSeq

Algorithm: co-local realignment of paired normal/tumor reads, pre-filter, estimate allele equilibrium frequences and evolutionary disance with F81 Markov substitution model, weighs frequencies against sample-specific error model, requires higher stringency at dbSNP locations

Notes: produces somatic calls. should give competitive performance on impure samples

Description: Markov Substitution model for Evolution (MuSE), which models the evolution of the reference allele to the allelic composition of the tumor and normal tissue at each genomic locus. We further adopt a sample-specific error model to identify cutoffs, reflecting the variation in tumor heterogeneity among samples.


### sinvict 

Validated vs: MuTect, VarScan2, Freebayes

Notes: Captures very low allele frequences, to detect mutations in free floating tumor dna. Can do time-series analysis. No confidence score assigned.


### multigems

Validated vs: freebayes, gatk, samtools, varscan

Notes: assumes diploid. Multiple-sample version of GeMS. 

Description: estimates sample genotypes and genotype probabilities for possible SNV sites. Unlike other popular multiple sample SNV callers, the MultiGeMS statistical model accounts for enzymatic substitution sequencing errors. Also, in consideration of the multiple testing problem associated with SNV calling, SNVs are called using a local false discovery rate (lFDR) estimator. Further, MultiGeMS utilizes high performance computing (HPC) techniques for computational efficiency and is robust to low-quality sequencing and alignment data.


### somaticseq

Validated vs: mutect/indellocator, varscan2, somaticsniper, jointsnvmix2, vardict

Notes: calls SNV/indel. produces somatic calls. meta caller, AI consensus

Description: Collects up to 72 features per mutation by SAMtools, HaplotypeCaller, and five orthogonal variant callers. The Adaptive Boosting model constructs a decision tree classifier that yields P for each variant.


### discosnp

Validated vs: niks, bubbleparse, cortex

Notes: Calls SNV/indel. FastQ input. de novo calls. Uses Cortex. Ranks predictions. Compute efficient.

Description:  designed to call isolated SNPs directly from sequenced reads, without a reference genome. ... DISCOSNP finds and ranks high quality isolated heterozygous or homozygous SNPs from any number of read sets, from 1 to n. It introduces new features to distinguish SNPs from sequencing errors or false positives due to approximate repeats. DISCOSNP can be used for finding high-quality isolated SNPs, either heterozygous, e.g. to build databases of high-quality markers within and across populations, or homozygous between individuals/strains, e.g. to create discriminant markers.


### 2kplus2

Notes: detects SNV/SV. Cortex input. Reference free de novo de Bruijn graph. 

Description: begins by producing a tree of kmers
for an input read set picking a seed k-mer and assuming that it
lies on one path through a SNP and then looks for an opposite kmer,
one substitution different, which would lie on another path
through the bubble. If this can be found in the k-mer tree, then a recursive
algorithm builds paths left and right of each k-mer until they
join or no k-mer can be found. Further to graph structure, the attributes
of the sample and sampled sequence reads can be used.


### exscalibur

Notes: Reports the union of multiple pipelines

Description: WES analysis pipelines for the detection of germline and somatic mutations,
with the implementation of three aligners, six germline callers, and six somatic callers. It
automates the full analysis workflow from raw sequencing reads to annotated variants and provides
an interactive visualization of the results


### multisnv

Notes: somatic calls

Validated vs: SomaticSniper, MuTect, UnifiedGenotyper and Platypus

Algorithm: probabilistic, timepoint-pooling, multiple samples from same patient

Description: a somatic variant caller that extends pairwise analysis of tumour-normal pairs to joint analysis of multiple samples from the same patient ... multiSNV calls somatic SNVs across all available same-patient samples without pooling reads. It is based on a Bayesian framework that captures the relatedness between samples by modelling the probability of a mutation in a given sample, conditioned on the somatic status of all other samples. 


### rarevator

Validated vs: mutect, varscan2

Notes: outputs SNV/indel. filter only, GATK UG input. somatic calls. validation only mentions how many new variants were called by rarevator, not how many were missed

Algorithm: Fisher exact test on conserved loci from hg19


### snv-ppilp

Notes: Filter only, gatk ug vcf input.

Algorithm: perfect phylogeny/integer linear programming

Description: a tool for refining GATK’s Unified Genotyper SNV calls for multiple samples. We assume these samples form a character-based phylogeny, the characters being the SNVs reported by GATK. As in Salari et al. (2013), we work with the perfect phylogeny model; however, we have a new problem formulation for fitting GATK’s calls to such a phylogeny, which we solve exactly using integer linear programming (ILP).


[top](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#snp-variant-callers)


### platypus

Validated vs: gatk ug/hc, samtools

Compared in: [Bro5](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#bro5)

Used by: bcbio, bioconda

Notes: [docs](http://www.well.ox.ac.uk/platypus-doc). calls SNV/SV/indel. somatic calls. haplotype-based. No dependencies, fast. Also see [Somatypus](https://github.com/baezortega/somatypus).


### baysic

Notes: input is paired vcfs only. 

Algorithm: meta, unsupervised bayesian consensus

Description: uses a Bayesian statistical method based on latent class analysis to combine variant sets produced by different bioinformatic packages (e.g., GATK, FreeBayes, Samtools) into a high-confidence set of genome variants.


### hapmuc

Validated vs: VarScan 2, SomaticSniper, Strelka and MuTect

Notes: SNV/indel. pileups input. somatic calls

Algorithm: bayesian model on haplotype inference

Description: two generative models under a Bayesian statistical framework: one represents true somatic mutations and the other regards candidate somatic mutations as errors. In our generative models, we prepared four candidate haplotypes by combining a candidate mutation and a heterozygous germ line variant, if available. The alignment probabilities of the observed reads given each candidate haplotypewere then computed by using profile hidden Markov models. Next, we inferred the haplotype frequencies and calculated the marginal likelihoods by using a variational Bayesian algorithm. Finally, we derived a Bayes factor, which is the ratio of the marginal likelihoods of these two models, to evaluate the possibility of the presence of somatic mutations.


### snpest

Validated vs: GeMS, freebayes, GATK HC, samtools

Notes: SNV/indel. pileups input. confidence score ranks predictions, does not model aneuploidy

Algorithm: reference-free probablistic model, generative probabilistic graphical model

Description: models the genotyping and SNP calling from the raw read sequences in a fully probabilistic framework. The problem is described using a generative probabilistic graphical model


### variantmaster

Notes: SNV/indel. somatic, do novo calls

Algorithm: reference-free probiblistic model, inference through inheritance

Description: uses raw sequence data information available in BAM files (binary sequence alignment/map format) in addition to the variants reported in the VCF files. BAM and VCF files can be generated using standard tools such as BWA, SAMtools, or GATK with default parameters. More specifically, for each variant in each affected individual, the algorithm estimates the strand bias and the probability that each family member is a carrier accounting for the respective fraction of supporting reads and the corresponding base call error rate


### mutect (1 & 2)

Validated vs: somatic sniper, jointSNVmix, strelka

Used by: GDC, SomaticSeq, bcbio, rave

Compared in: [Den9](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#den9), [Wash7](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#wash7), [Bcb8](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#bcb8), [Barc2](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#barc2), [Van6](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#van6), [Gor4](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#gor4), [Swi9](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#swi9)

Algorithm: bayesian with variable allele fraction, filter variants appearing in normal pool unless they are known variants. paired but not joint calling. No confidence score.

Notes: population-based calls. sensitive for low allelic frequency

Description: Bayesian classifier designed to detect somatic mutations with very low allele-fractions, requiring only a few supporting reads, followed by a set of carefully tuned filters


### niks

Notes: Identifies mutagen-induced mutations in paired samples.

Algorithm: reference free. map all reads to k-mers and count frequency changes between paired samples, reassemble

Description: reference-free genome comparison based solely on the frequencies of short subsequences within whole-genome sequencing data. It is geared toward identifying mutagen-induced, small-scale, homozygous differences between two highly related genomes, independent of their inbred or outbred background, and provides a route to identification of mutations without requiring any prior information about reference sequences or genetic maps


### ebcall

Notes: outputs SNV/indel. exome-only? population based calls. doesn't output vcfs. sensitive for low allelic frequency

Validated vs: varscan 2, somatic sniper

Compared by: [Den9](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#den9), [Van6](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#den9)

Algorithm: heuristic with beta-binomial error model from pooled normal bams, simply subtract germ line variants for somatic

Description: empirically estimating the distribution of sequencing errors by using a set of non-paired normal samples. Using this approach, we can directly evaluate the discrepancy between the observed allele frequencies and the expected scope of sequencing errors


### shearwater

Validated vs: caveman, mutect, deepsnv

Notes: population-based calls. for targeted sequencing.

Algorithm: beta-binomial model for variant calling with multiple samples

Description: exploits the power of a large sample set for precisely defining the local error rates and which uses prior information to call variants with high specificity and sensitivity.


### shimmer

Compared by: [Den9](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#den9), [Wash7](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#wash7)

Validated vs: varscan 2, somatic sniper, deepsnv, jointSNVmix 2

Algorithm: Fisher's exact test with multiple testing correction

Notes: somatic calls. employs a statistical model quite similar to that of Varscan 2, but in addition to this it performs a correction for multiple testing

Description: If the total number of reads displaying a non-reference allele in the two samples is greater than a minimum threshold nvar, a Fisher’s exact test is performed to test the null hypothesis that variant alleles are distributed randomly between the two samples 


### bubbleparse

Validated vs: cortex, samtools

Notes: Outputs SNV/SV. Cortex input

Description: minimal error-cleaning routine, followed by a depth-first search in the graph to find bubbles. 

[top](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#snp-variant-callers)


### cake

Notes: somatic calls

Validated vs: bambino, caveman, mpileup, varscan2

Algorithm: meta-caller - merge, consensus, filter. 

Description: integrates four publicly available somatic variant-calling algorithms to identify single nucleotide variants, Bambino, CaVEMan, SAMtools mpileup, and VarScan 2 with extra filtering


### denovogear

Notes: outputs SNV/indel. works with trios. do novo calls.

Validated vs: gatk, polymutt, samtools

Used in: biocondor

Algorithm: joint statistical analysis over multiple samples

Description: model consists of individual genotype likelihoods, transmission probabilities, and priors on the probability of observing a polymorphism or a de novo mutation at any given site in the genome


### qsnp

Validated vs: GATK, strelka

Algorithm: heuristic; minimum of 3 reads, compare to in house database of variants

Notes: somatic calls. fast, easy to run on a cluster

Description: Classification into germline and somatic calls follows a number of simple rules that were designed to accommodate for the expected low mutant allele ratio in low purity tumors


### rvd

Notes: MatLab. Detects very low frequency alleles

Description: [See here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3245950/). We use a multi-reference, indexed experimental design to minimize experimental variance and characterize a position-specific error distribution. We employ a rigorous statistical model to estimate the position-specific error rate distribution for reference sequences and thus the probability of a true mutation at each position in the sample. The statistical model provides a rigorous framework for hypothesis testing and estimation that minimizes false positives in variant calling.


### seurat

Compared by: [Den9](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#den9), [Wash7](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#wash7)

Notes: outputs SNV/indel/LOH/SV. somatic calls

Validated vs: varscan 2, strelka, somatic sniper

Description: calculates the joint posterior probability that a variant exists in the tumor sample and not in the normal sample. The resulting VCF file contains both SNVs and indels.


### snptools

Algorithm: haplotype imputation, effective base depth, binomial mixture modeling

Notes: includes genotype liklihood estimation

Description: SNPTools is organized by functionality into four modules ... EBD calculation: It summarizes mapping and base quality information to improve computational performance and reduce storage space.  SNP site discovery: The variance ratio statistic utilizes EBD information to provide high-quality SNP variant calls. ... GL estimation: BAM-specific parameter estimation allows this algorithm to overcome data heterogeneity due to platforms reference bias (from mapping or capture), and low-quality data. Genotype/haplotype imputation: A constrained Li-Stephens population haplotype sampling schema 


### vcmm

Notes: Detects SNV/indel/SV. pileups input. 

Validated vs: gatk, samtools

Algorithm: multinomial bayesian from paper in notes & strand bias filter

Description: The SNV calls were distinguished by the ratio of the probabilities that the minor allele at a nucleotide site is an error Perror and a major allele Pallele as described previously


### vip

Validated vs: dna sudoku, overlap log

Algorithm: overlapping pools

Description: A complete data analysis framework for overlapping pool designs, with novelties in all three major steps: variant pool and variant locus identification, variant allele frequency estimation and variant sample decoding. VIP is very flexible and can be combined with any pool design approaches and sequence mapping/alignment tools.


### virmid

Compared by: [Den9](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#den9)

Notes: measures purity. Exome only. somatic calls

Validated vs: jointSNVmix 2, strelka, varscan 2, 

Algorithm: Estimate purity, bayesian inference with estimated joint genotype probability matrix as the prior distribution

Description: estimate α, the level of impurity, i.e. the admixture of stromal cells in the cancer sample. A maximum likelihood estimation method is used. Next, the most probable genotype is estimated in the somatic variant caller step, using a Bayesian algorithm.


### varscan2

Compared by: [Den9](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#den9), [Wash7](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#wash7), [Bcb8](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#bcb8), [Aus4](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#aus4), [Van6](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#van6), [Gor4](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#gor4), [Swi9](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#swi9)

Notes: Calls SNV/indel/CNV. Mpileups input. somatic calls.

Validated vs: Somatic Sniper

Used by: GDC, SomaticSeq, bcbio, rave, bioconda

Algorithm: fisher's exact test, CBS alg for cnv, filters snps by heuristic criteria

Description: heuristic pairwise comparisons of base calls and normalized sequence depths at each position. Variants are classified into germline, somatic, LOH and unknown


[top](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#snp-variant-callers)


### jointsnvmix

Compared by: [Aus4](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#aus4), [Van6](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#van6), [Swi9](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#swi9)

Validated vs: compared to identical but non-joint and joint with fisher's exact test

Used in: SomaticSeq, rave

Algorithm: bayesian joint genotype of the samples

Notes: somatic calls. ranks the mutations. true joint calling

Description: probabilistic graphical model to analyse sequence data from tumour/normal pairs. allows statistical strength to be borrowed across the samples and therefore amplifies the statistical power to identify and distinguish both germline and somatic events in a unified probabilistic framework.


### lofreq

Validated vs: snver, breseq, samtools, some custom methods

Used in: somaticseq, biocondor

Algorithm: poisson binomial with bernoulli trials

Notes: somatic calls. uses bonferoni correction, tries for deep sequence <0.05 low allele frequency

Description: models sequencing run-specific error rates to accurately call variants occurring in <0.05% of a population


### strelka

Compared by: [Den9](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#den9), [Wash7](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#wash7), [Barc2](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#barc2), [Aus4](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#aus4), [Van6](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#van6), [Gor4](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#gor4)

Validated vs: varscan, samtools

Used in: SomaticSeq

Algorithm: bayesian joint probability of normal and somatic, indel realign

Notes: Outputs SNV/indel. somatic calls. works in presence of impurities, joint calling

Description: Bayesian approach wherein the tumor and normal allele frequencies are treated as continuous values. Search for candidate indels, realign, produce somatic variant probabilities. Strelka uses allele frequencies rather than diploid genotypes


### atlas2

Validated vs: gatk ug, dindel, samtools mpileup

Algorithm: logistic regression model includes reference/variant reads ratio for calling, and variety of features for filtering

Notes: outputs SNV/indel. exome only? part of Genboree. fast. may be a windows app. works for SOLiD, Illumina, and Roche 454

Description: Est. error as 11bp window rolling average. Filter variants on uni-directional reads.


### conan-snv

Notes: CNV-informed SNV. binomial mixture model, one per copy

Description: integrates information about copy number state of different genomic segments into the inference of single nucleotide variants. CoNAn-SNV requires as input a pileup file (either Maq or Samtools format) and model parameters, as well as a file demarcating segmentation boundaries of copy number amplifications


### cortex

Notes: Detects SNV/SV.

Algorithm: reference-free de bruijn graph

Description: extend classical de Bruijn graphs37,38 by colouring the nodes and edges in the graph by the samples in which they are observed. This approach accommodates information from multiple samples, including one or more reference sequences and known variants.


### deepsnv

Compared by: [Swi9](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#swi9)

Notes: SNV/indel. for targeted sequencing. population-based calls

Validated vs: varscan 2, crisp, vipr

Used by: Biocondor

Algorithm: beta-binomial model, error model from population data

Notes: uses population data, fast due to C implementation

Description: Model for error distribution is based on the observation that sequencing artifacts are recurrent on specific loci. In a large cohort this allows to define a background error distribution on each locus, above which true variants can be called.


### gems

Notes: pileups input. 

Validated vs: varscan2, snvmix2, freebayes, maq, samtools, gatk, atlas, soapsnp

Algorithm: bayesian multinomial, base- and alignment-quality priors, Dixon's Q-test

Notes: max of 2 alleles

Description: statistical model accounts for enzymatic substitution sequencing errors, addresses the multiple testing problem 


### impute2

Used by: biocondor

Algorithm: haplotype imputation

Description: statistically estimate the haplotypes underlying the GWAS genotypes (“pre-phasing”), then impute into these haplotypes as if they were correct


### somatic_sniper

Compared by: [Den9](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#den9), [Wash7](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#wash7), [Aus4](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#aus4), [Van6](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#van6), [Gor4](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#gor4), [Swi9](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#swi9)

Notes: somatic calls.

Validated vs: snvmix 2

Used by: GDC, SomaticSeq, rave, bioconda

Algorithm: basic joint probability bayesian genotyping

Description: is like Mutect based on a Bayesian posterior possibility. Somatic Sniper reports a somatic score (SSC), a Phred-scaled probability between 0 and 255, that the tumor and normal genotypes are different


[top](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#snp-variant-callers)


### bambino

Notes: outputs SNV/indel. somatic calls

Algorithm: basic. some filters

Description: Bambino's variant detector and assembly viewer are capable of pooling and analyzing data from multiple BAM files simultaneously.


### freebayes

Compared by: [Bcb8](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#bcb8), [Bro5](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#bro5)

Notes: outputs SNV/indel/MNPs

Used by: bcbio, biocondor

Algorithm: haplotype-aware bayesian inference with  multiallelic loci and non-uniform copy number across the samples

Description: generalize the Bayesian statistical method described by Marth to allow multiallelic loci and non-uniform copy number across the samples under consideration.


### mutationseq

Notes: Somatic calls. Available at command line for JointSNVMix

Validated vs: samtools, gatk ug

Algorithm: classic machine learning for somatic calling

Description: Comparison of four classic machine learning algorithms toward SNV calling 


### snver

Validated vs: CRISP, samtools, gatk

Algorithm: model minor alleles from pooled cancer/normal samples, using binomial dist

Notes: somatic calls. fast, early paired model, reports p-val

Description: statistical tool SNVer for calling SNPs in analysis of pooled or individual NGS data. Different from the previous models employed by CRISP, it analyzes common and rare variants in one integrated model, which considers and models all relevant factors including variant distribution and sequencing errors simultaneously.


### syzygy

Notes: mpileup input. population-based calls.

Algorithm: multinomial bayesian with filters

Description: empirical modeling of the sequencing error processes and filters to remove sites with strand inconsistency or clusters of variants suggestive of read misalignment


### vipr

Notes: mpileup input. outputs SNV/deletions. population-based calls. uses pooled samples

Validated vs: crisp, poisson, varscan

Description: vipR identifies sequence positions that exhibit significantly different minor allele frequencies in at least two DNA pools using the Skellam distribution.


### crisp

Notes: outputs SNV/indel

Algorithm: Fisher's exact test

Description: identify rare variants by comparing the distribution of allele counts across multiple DNA pools using contingency tables. To detect common variants, we utilize individual base-quality values to compute the probability of observing multiple non-reference base calls due to sequencing errors alone. Additionally, we incorporate information about the distribution of reads on the forward and reverse strands and the size of the pools to filter out false variants.


### indelocator
Compared by: [Den9](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#den9)

Used by: SomaticSeq

Notes: not published


[top](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SNV#snp-variant-callers)
