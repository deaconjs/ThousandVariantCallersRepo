## Structural Variant Callers


|Caller|Year|From|Study|Source|Algorithm|
|------|----|----|-----|------|---------|
|[popins](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#popins)|2017|deCODE genetics|[study](https://www.nature.com/articles/ng.3801)|[source](https://github.com/bkehr/popins)|Assembly of unmapped reads across multiple samples, placing of contigs into reference genome, genotyping|
|[svaba](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#svaba)|2017|Broad Institute|[study](http://biorxiv.org/content/early/2017/02/01/105080)|[source](https://github.com/walaj/svaba)|Discordant reads, classify, follow split reads to pool reads, assemble contigs|
|[valor](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#valor)|2017|Bilkent University|[study](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-3444-1)|[source](https://github.com/BilkentCompGen/Valor)|Long-range sequencing|
|[novobreak](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#novobreak)|2017|University of Texas Maryland Anderson Cancer Center|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5199621/)|[source](https://github.com/abjonnes/novoBreak)|de bruijn kmer hash to filter and assemble mutant contigs|
|[vardict](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#vardict)|2016|Astra Zeneca|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4914105/)|[source](https://github.com/AstraZeneca-NGS/VarDict)|split-read start, paired end, soft clipping, explicit alignment|
|[gridss](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#gridss)|2016|Walter and Eliza Hall Institute of Medical Research|[study](http://genome.cshlp.org/content/early/2017/11/02/gr.222109.117.short)|[source](https://github.com/PapenfussLab/gridss/)|split-read, paired end, soft clipping, genome-wide positional de bruijn graph break-end assembly|
|[sv-bay](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#sv-bay)|2016|Curie Institute|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4896370/)|[source](https://github.com/InstitutCurie/SV-Bay)|discordant reads, classify, read coverage, bayesian caller|
|[cosmos](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#cosmos)|2016|Advanced Industrial Science and Technology (AIST), Tokyo|[study](http://nar.oxfordjournals.org/content/44/8/e78)|[source](http://wall-lab.stanford.edu/projects/cosmos/)|discordant reads, classify, read coverage, score by binom dist on flanking read depths|
|[svstat](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#svstat)|2016|Baylor College of Medicine|[study](https://scfbm.biomedcentral.com/articles/10.1186/s13029-016-0051-0)|[source](https://github.com/abjonnes/novoBreak)|
|[svelter](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#svelter)|2016|University of Michigan|[study](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0993-1)|[source](https://github.com/mills-lab/svelter)|
|[skald](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#skald)|2016|University of Kansas Medical Center|[study](http://www.nature.com/articles/npjgenmed201626?utm_source=feedburner&utm_medium=feed&utm_campaign=Feed%3A+npjgenmed%2Frss%2Fcurrent+(npj+Genomic+Medicine))|N/A|
|[devro](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#devro)|2016|Uppsala University|[study](http://biorxiv.org/content/early/2016/12/15/094474)|N/A|
|[seq2c](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#seq2c)|2015|AstraZeneca|[study](http://cancerres.aacrjournals.org/content/76/14_Supplement/5268.short)|[source](https://libraries.io/github/AstraZeneca-NGS/Seq2CJava)|
|[manta](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#manta)|2015|Illumina|[study](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btv710)|[source](https://github.com/Illumina/manta)|split-read start, paired end, soft clipping, classification, assembly|
|[metasv](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#metasv)|2015|Bina/Roche|[study](https://www.ncbi.nlm.nih.gov/pubmed/25861968)|[source](http://bioinform.github.io/metasv/)|meta-caller|
|[wham](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#wham)|2015|U of Utah|[study](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004572)|[source](https://github.com/zeeev/wham)|soft-clipping, cluster alternative alignments from bwa output|
|[breakmer](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#breakmer)|2015|Broad, MacConaill|[study](http://nar.oxfordjournals.org/content/43/3/e19)|[source](https://github.com/ccgd-profile/BreaKmer)|split-read start, paired end, soft clipping, assemble kmer collections of sv reads|
|[indelminer](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#indelminer)|2015|Penn State U|[study](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0483-6)|[source](https://github.com/aakrosh/indelMINER)|split-read start, paired end, soft clipping, explicit alignment (indel)|
|[breakseek](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#breakseek)|2015|Chinese Academy of Sciences|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4538813/)|[source](https://sourceforge.net/projects/breakseek/)|soft-clipping, paired-end classification, sophisticated probabilistic scoring for indels|
|[raptr-sv](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#raptr-sv)|2015|USDA|[study](https://academic.oup.com/bioinformatics/article/31/13/2084/195773/RAPTR-SV-a-hybrid-method-for-the-detection-of)|[source](https://github.com/njdbickhart/RAPTR-SV)|split-read start, paired end, soft clipping|
|[scanindel](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#scanindel)|2015|University of Minnesota|[study](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-015-0251-2)|[source](https://github.com/cauyrd/ScanIndel)|soft-clipping, binomial distribution for breakpoints, assembly/mapping, freeBayes|
|[speedseq](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#speedseq)|2015|WashU St Louis, Hall|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4589466/)|[source](https://github.com/hall-lab/speedseq)|
|[lumpy](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#lumpy)|2014|U Virginia|[study](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-6-r84)|[source](https://github.com/arq5x/lumpy-sv)|breakpoint probability map, read-pair, split read, read-depth|
|[scalpel](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#scalpel)|2014|Cold Spring Harbor, Simons Center for Quantitative Biology|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4180789/)|[source](http://scalpel.sourceforge.net/)|de bruijn graph, iterative k-mer adjustment|
|[gindel](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#gindel)|2014|U Conneticut|[study](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0113324)|[source](https://sourceforge.net/projects/gindel/)|Support Vector Machine (SVM)|
|[gustaf](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#gustaf)|2014|Freie U, Berlin|[study](http://bioinformatics.oxfordjournals.org/content/30/24/3484)|[source](http://www.seqan.de/apps/gustaf/)|split-read start, align unmapped reads for breakpoints|
|[smufin](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#smufin)|2014|Barcelona Supercomputing Center|[study](http://www.nature.com/nbt/journal/v32/n11/full/nbt.3027.html)|[source](http://cg.bsc.es/smufin/)|somatic reference free, quaternary sequence tree|
|[socrates](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#socrates)|2014|The Walter and Eliza Hall Institute of Medical Research|[study](http://bioinformatics.oxfordjournals.org/content/30/8/1064)|[source](https://github.com/PapenfussLab/socrates)|split-read start, merge clusters|
|[ulysses](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#ulysses)|2014|Lab of Computational and Quantitative Biology, Paris|[study](http://bioinformatics.oxfordjournals.org/content/31/6/801.full?keytype=ref&%2520ijkey=fiLPQEO731TMaNt)|[source](http://www.lcqb.upmc.fr/ulysses/#citeUlysses)|
|[vivar](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#vivar)|2014|Ghent U, Belgium|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4264741/)|[source](https://www.cmgg.be/vivar/)|
|[bellerophon](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#bellerophon)|2013|Case Western|[study](https://www.ncbi.nlm.nih.gov/pubmed/23734783)|[source](http://cbc.case.edu/Bellerophon/)|discordant reads, interchromosomal, soft-clipping breakpoints|
|[sv-m](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#sv-m)|2013|Max Planck Institute|[study](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-14-132)|[source](https://github.com/BorgwardtLab/Indel-Prediction-SV-M)|
|[pesv-fisher](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#pesv-fisher)|2013|Center for Genomic Regulation, Spain|[study](https://www.ncbi.nlm.nih.gov/pubmed/23704902/)|[source](http://gd.crg.eu/tools/)|
|[isvp](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#isvp)|2013|Tohoku University, Japan|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4029547/)|[source](http://nagasakilab.csml.org/en/isvp)|meta-caller|
|[meerkat](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#meerkat)|2013|Harvard, Park|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3704973/)|[source](http://compbio.med.harvard.edu/Meerkat/)|discordant reads, classify, recognize specific complex events e.g. repair pathways|
|[soapindel](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#soapindel)|2013|BGI Shenzhen|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530679/)|[source](http://soap.genomics.org.cn/soapindel.html)|de bruijn graph, identifies breakpoints from discordant reads|
|[softsearch](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#softsearch)|2013|Mayo Clinic, Kocher|[study](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0083356)|[source](http://bioinformaticstools.mayo.edu/research/softsearch/)|soft-clipping, heuristic, number of soft-clipped reads|
|[tigra](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#tigra)|2013|WashU MD Anderson Cancer Center, Weinstock|[study](http://genome.cshlp.org/content/early/2013/12/04/gr.162883.113)|[source](http://bioinformatics.mdanderson.org/main/TIGRA)|de bruijn graph, requires input break points|
|[delly](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#delly)|2012|EMBL|[study](http://bioinformatics.oxfordjournals.org/content/28/18/i333.abstract)|[source](https://github.com/dellytools/delly)|discordant reads, classify, adds long-range mate pairs, split reads to get breakpoint|
|[cn.mops](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#cn.mops)|2012|Johannes Kepler U, Australia|[study](http://nar.oxfordjournals.org/content/40/9/e69)|[source](http://www.bioconductor.org/packages/release/bioc/html/cn.mops.html)|
|[battenberg](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#battenberg)|2012|Welcome Trust Sanger Institute, UK|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3428864/)|N/A|
|[breakpointer](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#breakpointer)|2012|Max Plank Institute, Haas|[study](http://bioinformatics.oxfordjournals.org/content/28/7/1024)|[source](https://github.com/ruping/Breakpointer)|
|[clever](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#clever)|2012|Life Sciences Group, Amsterdam|[study](http://bioinformatics.oxfordjournals.org/content/28/22/2875.long)|[source](http://www.mybiosoftware.com/clever-2-0rc1-clique-enumerating-variant-finder.html)|discordant reads, cluster on concordant pairs|
|[forestsv](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#forestsv)|2012|University of California San Diego, Sebat lab|[study](https://www.ncbi.nlm.nih.gov/pubmed/22751202)|[source](http://sebatlab.ucsd.edu/index.php/software-data)|
|[gasvpro](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#gasvpro)|2012|Brown U|[study](http://genomebiology.biomedcentral.com/articles/10.1186/gb-2012-13-3-r22)|[source](https://code.google.com/archive/p/gasv/)|
|[hugeseq](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#hugeseq)|2012|Stanford University|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4720384/)|[source](https://github.com/StanfordBioinformatics/HugeSeq)|plural caller with simple aggregation and voting|
|[prism](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#prism)|2012|U Toronto,|[study](http://bioinformatics.oxfordjournals.org/content/28/20/2576)|[source](http://compbio.cs.toronto.edu/prism/)|
|[splitread](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#splitread)|2012|Howard Hughs Medical Institute|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3269549/)|[source](http://splitread.sourceforge.net/)|split-read start, hamming distance, de novo via read depth|
|[svm2](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#svm2)|2012|Univerisity of Milan|[study](https://academic.oup.com/nar/article/40/18/e145/2411257/SVM-2-an-improved-paired-end-based-tool-for-the)|[source](http://elixir-italy.org/milano/en/archives/servizi/svm2-1-2)|
|[clipcrop](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#clipcrop)|2011|U Tokyo|[study](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-S14-S7)|[source](https://github.com/shinout/clipcrop)|
|[crest](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#crest)|2011|St Jude, Zhang|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3527068/)|[source](http://www.stjuderesearch.org/site/lab/zhang)|soft-clipping, binomial distribution, assembly/mapping|
|[genomestrip](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#genomestrip)|2011|Broad, McCarroll|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5094049/)|[source](http://software.broadinstitute.org/software/genomestrip/)|discordant reads, reassemble by allele, read-depth, breakpoint database|
|[ingap-sv](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#ingap-sv)|2011|Chinese Academy of Sci, Zhao|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3125812/)|[source](http://ingap.sourceforge.net/)|depth of coverage, paired end|
|[hydra](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#hydra)|2010|U Va|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2860164/)|[source](https://code.google.com/archive/p/hydra-sv/)|split-read start, paired end|
|[age](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#age)|2010|Yale, Gerstein|[study](http://bioinformatics.oxfordjournals.org/content/27/5/595.abstract)|[source](http://sv.gersteinlab.org/age/)|
|[slope](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#slope)|2010|WashU St Louis, Pfiefer|[study](https://academic.oup.com/bioinformatics/article/26/21/2684/214667/SLOPE-a-quick-and-accurate-method-for-locating-non)|[source](http://www-genepi.med.utah.edu/suppl/SLOPE/slope_guide.txt)|
|[svdetect](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#svdetect)|2010|Curie Institute|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2905550/)|[source](http://svdetect.sourceforge.net/Site/Manual.html)|discordant reads, classify, sliding window clustering, mate-pair|
|[svmerge](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#svmerge)|2010|Sanger Institute|[study](http://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-12-r128)|[source](http://svmerge.sourceforge.net//)|meta-caller|
|[pindel](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#pindel)|2009|EMBL, Ning|[study](https://www.ncbi.nlm.nih.gov/pubmed/19561018)|[source](http://gmt.genome.wustl.edu/packages/pindel/index.html)|split-read start|
|[breakdancer](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#breakdancer)|2009|WashU St Louis, Mardis|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3661775/)|[source](http://gmt.genome.wustl.edu/packages/breakdancer/)|discordant reads, classify, maq-based|
|[breakseq](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#breakseq)|2009,15|Yale, Gerstein|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2951730/)|[source](http://sv.gersteinlab.org/breakseq/)|map reads to breakpoints|
|[pemer](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#pemer)|2009|Yale, Gerstein|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2688268/)|[source](http://sv.gersteinlab.org/pemer/)|split-read start, merge clusters|

### popins

Notes: Non-reference sequence insertions from short-reads, population-scale

Validated vs: [MindTheGap](http://gatb.inria.fr/software/mind-the-gap/), [Pamir](https://academic.oup.com/bioinformatics/article/33/14/i161/3953969)

Used on WGS data of > 15,000 Icelanders and included in [Graph Genome Pipeline](https://www.biorxiv.org/content/early/2017/09/29/194530) and [Illumina Polaris](https://github.com/Illumina/Polaris)

Algorithm: joint assembly of unmapped reads across samples

Description: Collects reads without good alignment to the reference genome, filters these reads for contamination, and assembles the remaining ones into contigs. Next it merges the contigs across samples, which improves the assembly of non-reference sequence insertions shared by several individuals. The merged contigs are anchored into the reference genome using paired-end information and exact breakpoints positions are determined using split alignment. Popins finishes by computing genotype likelihoods for all anchored contig ends in all samples.

### svaba 

Notes: low memory, fast

Algorithm: identifies clipped, discordant, and unmapped reads, split pairs, and reads with
deletions or insertions in the CIGAR string. Discordant reads are re-aligned, filtered, and clustered. Split read partners are identified and reads are pooled. Then realign grouped reads into contigs with specialized SGA assembler, align these contigs to reference with bwa-mem. Re-align all constituent reads to the contig or to reference keeping reads that match contig better. 

Description:  perform local assembly to create consensus contigs from sequence reads with
divergence from the reference, and to apply this procedure to every region of the genome. The
contigs are then compared to the reference to annotate the variants. By uniting the different
classes of variant-supporting reads into a single framework, we further expect that this assemblyfirst
approach would be effective for variants of all sizes and require few parameters

### valor

Notes: long range sequencing e.g. 10X Genomics linked-read sequencing, pooled clone sequencing

Description:  (variation using long range information). Briefly, valor searches for both read pair and split clone sequence signatures using the mapping locations of long range sequencing reads, and requires split clones from different pools to cluster at the same putative inversion breakpoints. Ambiguity due to multiple possible pairings of split clones are resolved using an approximation algorithm for the maximal quasi clique problem.

### novobreak

Description: novoBreak algorithm divides tumor reads into shorter sections that are k nucleotides long (k-mers). By hashing and filtering out k-mers that match the reference and normal genomes, the algorithm identifies k-mers that are unique to the tumor and indicate breakpoints. Essentially a de Bruijn graph. Reads containing the unique k-mers are assembled into contigs local to the breakpoints, which are then aligned to the reference genome in order to infer exact breakpoints and their associated structural variants.

Description: obtains genome-wide local assembly of breakpoints from clusters of reads sharing a set of k-mers uniquely present in a subject genome but not in the reference genome or any control data... constructs a hash table from the tumor reads, containing all the k-mers, their host reads and frequencies in the set. Next, it filters out k-mers representing reference alleles or sequencing errors, and retains those representing variants... classifies the k-mers into 1) germline k-mers, those present in both the tumor and the normal genome, and 2) somatic k-mers, those present in the tumor but not the normal genome. Then, novoBreak identifies clusters of read pairs spanning each somatic breakpoint, and assembles each cluster of reads into contigs. By comparing the resulting high-quality contigs with the reference, novoBreak identifies breakpoints and associated SVs. Finally, novoBreak quantifies the amount of the supporting evidences at each breakpoint and outputs a final report.

### vardict

Notes: outputs SV/indel/SNV/LOH. somatic caller. Efficient with ultra deep seq. calls complex variants (same read, multiple vars). filters PCR artifacts. Estimates SV allele frequency. 

Validated vs: GATK UG/HC, Freebayes, varscan, pindel, scalpel, manta, lumpy

Used by: SomaticSeq, RAVE, bcbio

Compared by: [bcb8](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#bcb8)

Algorithm: consensus on realigned soft-clipped reads used as search query 

Description: Calls SNV, MNV, InDels, complex and structural variants, performs local realignments on the fly. Performance scales linearly to sequencing depth. Performs amplicon aware variant calling for polymerase chain reaction (PCR)-based targeted sequencing often used in diagnostic settings. Is able to detect PCR artifacts. Detects differences in somatic and loss of heterozygosity variants between paired samples.


### gridss

Used by: bcbio

Description: performs alignment-constrained whole genome breakend assembly using a novel positional de Bruijn graph algorithm and a probabilistic structural variant caller that combines assembly, split read, and read pair evidence in a unified variant scoring model. 


### sv-bay

Validated vs: gasvpro, breakdancer, lumpy, delly

Notes: works with paired end or mate pair; analyzes tumor/normal pair concurrently

Algorithm: PEM & read coverage with bayesian testing for adjacency

Descripton: we combine both PEM signatures and information about changes in DOC in regions flanking each candidate rearrangement. Our method takes into account GC-content and mappability. The use of a Bayesian framework based on both PEM and DOC information allows us to significantly decrease the level of false positive predictions while retaining high sensitivity. Additionally, SV-Bay infers 15 different types of structural variant from the detected novel genomic adjacencies


### cosmos

Validated vs: Breakdancer, GasVPro, Delly, Lumpy

Notes: somatic caller. compares the statistics of the mapped read pairs in tumor samples with isogenic normal control samples in a distinct asymmetric manner. fast. mouse model for validation plus synthetic

Algorithm: discordant pair reads, classify, DOC binomial

Description: compares the mapping read status of paired-end short reads in a tumor sample with a normal sample in an asymmetric manner: groups of discordant read pairs, which are indicative of SVs, are generated from the tumor sample, following which the groups are filtered against individual discordant read pairs, instead of the group equivalents, in the normal sample to eliminate false positives. Next, we introduce the concept of strand-specific read depth, which allows prioritization of candidate SVs more efficiently than the conventional strand-independent read depth. 


### svstat

Description: we explored methods for quantifying support for nucleotide-resolved breakpoints of SVs without PE, SR, or DN. all reads are aligned to the reference genome... Recurrent alignment stop or start coordinates indicate candidate breakpoints... Candidate breakpoint regions are paired with each other to form a sequence “library” of candidate junctions. Stack reads are then aligned to the library of candidate SVs, and evidence for each candidate junction (C) is calculated based on 1) the number of bases in the tails aligned to the partner region, and 2) the quality scores of the alignments.

### svelter

Notes: Handles complex rearrangements

Description: accurately resolve complex structural genomic rearrangements in whole genomes. Unlike previous “bottom up” strategies that search for deviant signals to infer structural changes, our “top down” approach works by virtually rearranging segments of the genomes in a randomized fashion and attempting to minimize such aberrations relative to the observed characteristics of the sequence data. In this manner, SVelter is able to interrogate many different types of rearrangements, including multi-deletion and duplication-inversion-deletion events as well as distinct overlapping variants on homologous chromosomes

### skald

Notes: detects >50nt deletion structural variants

Description: combines calls from two tools (Breakdancer and GenomeStrip) with calibrated filters and clinical interpretation rules.


[top](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#structural-variant-callers)


### devro

Notes: population: paired end/depth of coverage


### seq2c

Notes: From AstraZeneca, inactive


### manta

Notes: produces SV and indel calls. runs on trios. Parallelized for clusters. handles degraded FFPE samples, uses pedigree-consistency and cosmic for validation. faster than delly

Validated vs: pindel, delly

Used by: bcbio, bioconda

Algorithm: breakend graph/assembly

Description: provides scoring models for germline analysis of diploid individuals and somatic analysis of tumor-normal sample pairs, with additional applications under development for RNA-Seq, de novo variants, and unmatched tumors.  less than a tenth of the time that comparable methods require


### metasv

Validated vs: pindel, BreakSeq2, LUMPY, BreakDancer, Delly, CNVNator, MindTheGap

Used by: bioconda, [bcbio](http://blog.bina.com/read/metasv-integration-in-bcbio)

Algorithm: merge then assemble: pindel, BreakSeq2, LUMPY, BreakDancer, Delly, CNVNator, MindTheGap

Notes: consensus caller

Description: merging SVs from multiple tools for all types of SVs. It also analyzes soft-clipped reads from alignment to detect insertions accurately since existing tools underestimate insertion SVs. Local assembly in combination with dynamic programming is used to improve breakpoint resolution. Paired-end and coverage information is used to predict SV genotypes.


### wham

Notes: de novo calls. requires bwa. does association testing

Validated vs: lumpy, delly, softsearch

Used by: bcbio, biocondor

Compared by: [Bcb8](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#bcb8)

Algorithm: mate-pair & split read mapping, soft-clipping, alternative alignment, consensus sequence based evidence

Description: pinpoint SVs in pooled and genotypic data associated with phenotypic variation. uses split-read, mate-pair, and alternative alignments to find the other SV breakpoint. Positions in the pileup where three or more primary reads share the same breakpoint are interrogated as a putative SV. Use SA and XA cigar tags as alternative alignment locations, cluster those. SW align clipped consensus to alternative locations. intra-chromosomal require min 2 reads. 


### breakmer

Notes: For targeted sequencing. somatic calls. 

Validated vs: crest, meerkat, breakdancer, pindel

Algorithm: soft clip, identify kmers in reads but not in reference

Description: uses a ‘kmer’ strategy to assemble misaligned sequence reads for predicting insertions, deletions, inversions, tandem duplications and translocations at base-pair resolution in targeted resequencing data. Variants are predicted by realigning an assembled consensus sequence created from sequence reads that were abnormally aligned to the reference genome


### indelminer

Notes: outputs indels only. simple de novo, somatic calls. validatin against synthetic variants introduced to chr22 and the na18507 data set. recommended to align  with gatk indelRealigner

Validated vs: samtools, pindel, prism

Algorithm: split-read, paired-end, soft-clipped. align unmapped reads at both ends to look for indels

Description: uses a split-read approach to identify the precise breakpoints for indels of size less than a user specified threshold, and supplements that with a paired-end approach to identify larger variants that are frequently missed with the split-read approach


### breakseek

Notes: describes parameters for competitors, works reasonable well for all size indels, estimates level of heterozygosity

Validated vs: pindel, lumpy, crest, soapindel, breakdancer, prism, delly

Algorithm: soft-clipping/breakread break points, paired-end span to validate indel, sophisticated probabilistic scoring model

Description: unbiasedly and efficiently detect both homozygous and heterozygous INDELs, ranging from several base pairs to over thousands of base pairs, with accurate breakpoint and heterozygosity rate estimations


### raptr-sv

Notes: sensitivity for tandem duplications

Algorithm: discordant read-pair, split read, soft-clip, filter

Description: combining their predictions to generate highly confident SV calls, which can be filtered at runtime for improved accuracy.


### scanindel

Algorithm: Identify quality soft-clipped reads, cluster and use binomial distribution to identify breakpoints, remap these breakpoint contigs and unmapped reads to reference with BLAT, reassemble breakpoint contigs de novo with inchworm/Trinity and map with BLAT, both remapped and reassembled paths produce BAM files freeBayes calls indels on.

Description:  integrates multiple signals from all three sources (gapped alignment, split reads and de novo assembly) allows for more sensitive indel discovery than methods examining merely one or two signals. Our framework scans the initial mapping file from a gapped NGS aligner and refines the alignment of the soft-clipped reads meeting tiered criteria. Next, de novo assembly is performed for the selected soft-clipped reads and unmapped reads. Subsequent to the re-alignment and assembly, we have applied a Bayesian haplotype-based variant caller to detect indels. 


### speedseq

Notes: assumes diploid. uses pedigree+ to call validation variants. discusses parameters. reports confidence score.

Algorithm: meta. freebayes, lumpy, cnvnator, and custom caller svtyper. 


[top](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#structural-variant-callers)

### lumpy

Notes: sensitive to low MAF. not great for small dels

Validated vs: gasvpro, delly, pindel

Used by: bcbio, metasv, biocondor

Compared by: [bcb8](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#bcb8), [bcb3](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#bcb3)

Algorithm: merge read-pair, split read, read-depth in a breakpoint probability map. classify and cluster 

Description: LUMPY integrates disparate signals by converting them to a common format in which the two predicted breakpoint intervals in the reference genome are represented as paired probability distributions.


### scalpel

Notes: indel only. exome-capture data. de novo calls. slow, not for wgs. does indel normalization

Validated vs: gatk hc, SOAPindel

Used by: somaticseq, bcbio, biocondor

Compared by: [Bcb8](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#bcb8)

Algorithm: de bruijn graph traversal, local assembly with iterative k-mer k value reassessment to eliminate repeats. 

Description: localized micro-assembly of specific regions of interest with the goal of detecting mutations with high accuracy and increased power. It is based on the de Bruijn graph assembly paradigm and implements an on-the-fly repeat composition analysis coupled with a self-tuning k-mer strategy


### gindel

Notes: indels > 50bp only. efficient

Validated vs: pindel, cleversv

Algorithm: SVM on 7 features: discordant pair, split-read, read depth, concordant encompassing pair, single-end-mapped pair, partially mapped reads, fully-mapped spanning reads, 

Description: An approach for calling genotypes of both insertions and deletions from sequence reads. GINDEL uses a machine learning approach which combines multiple features extracted from next generation sequencing data. It performs well for insertion genotyping on both simulated and real data. GINDEL can not only call genotypes of insertions and deletions (both short and long) for high and low coverage population sequence data, but also is more accurate and efficient than other approaches.


### gustaf

Notes: small validation set. claims best for small SVs 30-100bp. small validation set. might work with FFPE

Validated vs: delly, pindel

Used by: biocondor

Algorithm: local alignment between unmapped reads shows breakpoints

Description: based on a generic multi-split alignment strategy that can identify SV breakpoints with base pair resolution.


### smufin

Notes: produces SV and SNV calls. paired exome fastq input. somatic calls. parallelized

Validated vs: mutect, breakdancer, pindel, delly, crest

Algorithm: directly compares reads "quaternary sequence tree"

Description: directly compares sequence reads from normal and tumor genomes to accurately identify and characterize a range of somatic sequence variation, from single-nucleotide variants (SNV) to large structural variants at base pair resolution.


### socrates

Notes: somatic calls. needs parameterized. fast. split-read only algorithms are better for short read FFPE data. "On real tumour data without additional information, we find it impractical to run at its most sensitive settings, but it is easily tuned."

Algorithm: re-aligns and clusters soft-clipped reads

Description: uses split reads to find breakpoints. It is optimized to be fast and extremely sensitive.


### ulysses

Notes: mate-pair only

Description: assessing, in a principled manner, the statistical significance of each possible variant (duplications, deletions, translocations, insertions and inversions) against an explicit model for the generation of experimental noise.


### vivar

Notes: needs a reference set for sequencing error model

Description: facilitates the processing, analysis and visualization, of structural variation based on massive parallel sequencing data


### bellerophon

Notes: detects interchromosomal translocations. really basic caller

Validated vs: GASV, Breakdancer, SVDetect, and CREST

Algorithm: discordant reads, soft-clipped

Description: uses discordant read pairs and "soft-clipped" reads to predict the location of the precise breakpoints.  for each chimeric breakpoint, attempts to classify it as a participant in an unbalanced translocation, balanced translocation, or interchromosomal insertion.


### sv-m

[top](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#structural-variant-callers)

### pesv-fisher

### isvp

Notes: meta-caller. deletions only

Validated vs: breakdancer, delly, pindell, haplotypecaller

Algorithm: limit each method's calls to their optimal size ranges



### meerkat

Notes: somatic calls. specific complex rearrangement events like dna repair pathways. 

Algorithm: discordant read pair clustering with refinement

Description: considers local clusters of discordant read pairs to recognize specific complex events. uses split, clipped, and multiple-aligned reads


### soapindel

Notes: calls indels. similar sensitivity and specificity for small indels, higher sensitivity for large indels. might be slow. should call SNPs too. assigns confidence q-scores. weird validation looking at hg19 vs venter genome and chimpanzee vs hg19

Validated vs: dindel, pindel, gatk

Algorithm: identifies breakpoints from discordant reads, multi-path de bruijn graph assembly

Description: assign all unmapped reads with a mapped partner to their expected genomic positions and then perform extensive de novo assembly on the regions with many unmapped reads to resolve homozygous, heterozygous, and complex indels by exhaustive traversal of the de Bruijn graph


### softsearch

Notes: slow but high TP rate. works at low depth but needs parameters adjusted. Levenstein Distance "confidence scores"

Validated vs: breakdancer, delly, crest, svseq

Algorithm: soft-clipping heuristic, number of soft-clipped reads per position

Description: Assuming soft clipping delineates the exact breakpoint position and direction, DRPs overlapping such soft-clipped areas should already contain the information about the type and size of SV, obviating the need for secondary alignments.


### tigra

Notes: calls breakpoints. population calls locate common alleles. de novo? uses population data, low FDR

Algorithm: iterative breakpoint collection, de bruijn graph, assembles the breakpoints


### delly

Notes: calls SV and CNV. high sensitivity and specificity, lower sensitivity to small deletions

Validated vs: pindel, breakdancer, gasv, hydra

Used by: bcbio, metasv, biocondor

Compared by: [bcb8](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#bcb8), [bcb3](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/Compared#bcb3)

Algorithm: integrates short insert paired-ends and long-range mate-pairs to identify discordant pairs, then uses split-read alignments to identify breakpoints

Description: integrates short insert paired-ends, long-range matepairs and split-read alignments to accurately delineate genomic rearrangements


### cn.mops

Notes: calls CNV (move to CNV). population-based calls. 

Used by: biocondor

Compared by: bcb3

Description: decomposes variations in the depth of coverage across samples into integer copy numbers and noise by means of its mixture components and Poisson distributions


### cpgbattenberg

### breakpointer

Notes: single-end read breakpoint locator

Validated vs: pindel

Description: By taking advantage of local non-uniform read distribution and misalignments created by SVs, Breakpointer scans the alignment of single-end reads to identify regions containing potential breakpoints.


[top](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#structural-variant-callers)


### clever

Notes: better at 20-100bp size range

Validated vs: gasv, variationhunter, breakdancer, hydra

Used by: biocondor

Algorithm: clustering on concordant pairs

Description: enumerates all max-cliques and statistically evaluates them for their potential to reflect insertions or deletions.


### forestsv

Algorithm: random forests


### gasvpro

Notes: validates vs HuRef, NA18705, and NA12878, has ROC curves. models uncertainty in call/reference overlap for truth calls to be more precise 

Validated vs: hydra, breakdancer, CNVer

Algorithm: joint P on paired read and read depth

Description: Combines read depth information along with discordant paired-read mappings into a single probabilistic model two common signals of structural variation.


### hugeseq

Notes: SNP/SV/CNV caller. "whole pipeline" includes PCR artifact removal, GATK realignment, variant calling, functional annotation with Annovar.

Algorithm: SNPs with UnifiedGenotyper and SAMtools, indels with Dindel, CNV/SV with BreakDancer, Pindel, CNVnator, BreakSeq. basic vote 2+ gives "high confidence"


### prism

Description: uses a split-alignment approach informed by the mapping of paired-end reads, hence enabling breakpoint identification of multiple SV types, including arbitrary-sized inversions, deletions and tandem duplications


### splitread

Notes: exome only. de novo calls. compares read depth between parents and child to identify de novo mutations

Algorithm: discordant pairs clustering, map with mrsFast and hamming distance, call anomolous mappings, split unmapped reads, search.

Description: searches for clusters of mate pairs where one end maps to the reference genome but the other end does not because it traverses a breakpoint creating a mapping inconsistency with respect to the reference sequence


### clipcrop

Notes: not for somatic. doesn't recognize useful mutation types according to socrates authors

Validated vs: breakdancer, cnvnator, pindel

Description: A soft-clipped sequence is an unmatched fragment in a partially mapped read


### crest

Notes: somatic calls. made for somatic comparisons, lower performance for small deletions

Algorithm: single read soft clipping, classification

Description: uses the soft-clipping reads to directly map the breakpoints of structural variations


### genomestrip

Notes: calls SV/CNV. population-based calls. validated by 1kGP, "most sensitive and accurate". less sensitive to small SVs. 2.0 adds CNV detection

Validated vs: spanner, pindel, breakdancer, pemer, cnvnator

Algorithm: discordant read pair clustering, reassemble breakpoint-spanning reads by allele, read-depth for copy number estimate, align unmapped reads to breakpoint database

Description: designed to find shared variation using data from multiple individuals. Genome STRiP looks both across and within a set of sequenced genomes to detect variation.


### ingap-sv

Notes: good validation against 12878, differentiates homo- and hetero-zygous variants

Validated vs: Breakdancer, variationhunter, spanner, PEMer, cortex, pindel

Algorithm: paired-end mapping & depth of coverage


### hydra

Notes: sanger split-read & illumina paired-end input. 

Validated vs: None. 

Algorithm: split-read + paired end


[top](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#structural-variant-callers)

### age

Notes: improved alignment algorithm, but does not call variants. narrow scope validation, proof of theory. has been modified for metasv inclusion

Used by: biocondor

Algorithm: read-depth

Description: AGE for Alignment with Gap Excision, finds the optimal solution by simultaneously aligning the 5′ and 3′ ends of two given sequences and introducing a ‘large-gap jump’ between the local end alignments to maximize the total alignment score. We also describe extensions allowing the application of AGE to tandem duplications, inversions and complex events involving two large gaps.


### slope

Notes: basic simulated data

Validated vs: pindel, breakdancer

Description: detect sequence breakpoints from only one side of a split read, and therefore does not rely on the insert size for detection.


### svdetect

Notes: paired-end/mate pairs input. 

Validated vs: GasV

Algorithm: discordant read pairs, adds mate-pairs, sliding window for clustering

Description: anomalously mapped read pairs provided by current short read aligners to localize genomic rearrangements and classify them according to their type, e.g. large insertions– deletions, inversions, duplications and balanced or unbalanced interchromosomal translocations.


### svmerge

Algorithm: meta-caller - BDMax, Pindel, SECluster, RetroSeq, RDXplorer


### pindel

Notes: slow, high FP rate

Used by: metasv, biocondor

Algorithm: split-read clustering, pattern growth algorithm to search local space for unmapped (split) read

Description: detect breakpoints of large deletions (1bp-10kbp) and medium sized insertions (1-20bp) from paired-end short reads


### breakdancer

Notes: somatic calls. confidence scores, use Q>80

Used by: metasv

Validated vs: MoDIL, VariationHunter

Algorithm: paired-end MAQ calls; classify, cluster, multi-nomial Poisson-based confidence score

Descriptions:  predicts large and small 10-100bp indels, inversions and translocations


### breakseq

Used by: metasv

Algorithm: map reads to known breakpoints from a database

Description: scanning the reads from short-read sequenced genomes against our breakpoint library to accurately identify previously overlooked SVs


### pemer

Validated vs: PEM

Algorithm: split-read clustering, merges clusters


[top](https://github.com/deaconjs/ThousandVariantCallersRepo/wiki/SV#structural-variant-callers)
