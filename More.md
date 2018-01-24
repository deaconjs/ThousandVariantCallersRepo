## MSI CALLERS

|Caller|PubYear|Specialty|From|Study|Source|
|------|-------|---------|----|-----|------|
|msisensor|2014|Microsatellite Instability(MSI)|Memorial Sloan Kettering Cancer Center|[study](http://ascopubs.org/doi/pdf/10.1200/PO.17.00084)|[source](https://github.com/ding-lab/msisensor)|
|msings|2014|MSI|University of Washington, Seattle|[study](https://www.ncbi.nlm.nih.gov/pubmed/24987110)|[source](https://bitbucket.org/uwlabmed/msings)|[study](http://www.sciencedirect.com/science/article/pii/S1525157815001531)|
|msiseq|2015|MSI|National Cancer Center Singapore|[study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4549793/)|[source](https://cran.r-project.org/web/packages/MSIseq/index.html)|
|mirmmr|2017|MSI with methylation and mutation|Siteman Cancer Center, Wash U in St. Louis|[study](https://www.ncbi.nlm.nih.gov/pubmed/28961932)|[source](https://github.com/ding-lab/MIRMMR)|
|mantis|2017|MSI|The Ohio State University|[study](https://www.ncbi.nlm.nih.gov/pubmed/27980218)|[source](https://github.com/OSU-SRLab/MANTIS)|
|no name|2017|MSI classification|Harvard Medical School|[study](https://www.nature.com/articles/ncomms15180)|contact authors|


## Other Devices & de novo Filters

|caller|orig pub|caller class|input type|somatic/denovo|from|validated vs.|cited|used by|compared by|algorithm|features|description|installation|study|source|
|------|--------|------------|----------|--------------|----|-------------|-----|-------|-----------|---------|--------|-----------|-----------|-----|------|
|msisensor|2014|MSI instability||somatic|||15|||||||https://academic.oup.com/bioinformatics/article/30/7/1015/236553/MSIsensor-microsatellite-instability-detection|https://github.com/ding-lab/msisensor|
|hipSTR|2016|STR|||polymut||0|||haplotype|STR|||http://www.biorxiv.org/content/early/2016/09/27/077727|https://hipstr-tool.github.io/HipSTR/|
|otg-snpcaller (available?)|2014||torrent||||8|||||||http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0097507|available?|
R453Plus1toolbox|2011||454||||18|||||||https://academic.oup.com/bioinformatics/article/27/8/1162/228803/R453Plus1Toolbox-an-R-Bioconductor-package-for|http://www.bioconductor.org/packages/2.10/bioc/html/R453Plus1Toolbox.html|
|eSNV-Detect|2014|eSNV|rna-only||||9|||||||https://academic.oup.com/nar/article/42/22/e172/2410988/The-eSNV-detect-a-computational-system-to-identify|http://bioinformaticstools.mayo.edu/research/esnv-detect/|
|graphmap|2016|SNV|nanopore||||8|||||||http://www.nature.com/articles/ncomms11307|https://github.com/benedictpaten/marginAlign/blob/master/src/margin/marginCaller.py|
marginAlign|2015|SNV|nanopore||||149|||||||https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4907500/|https://github.com/benedictpaten/marginAlign|
|radia|2014|SNV/indel|rna & dna||||12|||||||http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0111516|https://github.com/aradenbaugh/radia/|
|UNCeqR|2014|SNV/indel|rna & dna||||20|||||||https://academic.oup.com/nar/article/42/13/e107/1277201/Integrated-RNA-and-DNA-sequencing-improves|http://lbg.med.unc.edu/~mwilkers/unceqr_dist/|
|pyroHMMvar|2013|SNV/indel|torrent/454||Tsinghua U, Beijing||6|||||||https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3888126/|https://github.com/homopolymer/PyroTools/|
|SMRT-SV|2016|SV|pacbio||||||||||||http://genome.cshlp.org/content/early/2017/03/31/gr.214007.116|https://github.com/EichlerLab/pacbio_variant_caller|
|HySA|2017|SV/indel|pacbio + illumina||||0||||||http://genome.cshlp.org/content/early/2017/01/19/gr.214767.116.abstract#corresp-1|https://bitbucket.org/xianfan/hybridassemblysv|
|Parliament|2015|SV|bio nano||Baylor, Gibbs||37|||BioNano assembly|Integrates Illumina, PacBio, BioNano reads|adds bio nano||https://www.ncbi.nlm.nih.gov/pubmed/25886820|https://sourceforge.net/projects/parliamentsv/|
|Multibreak-SV|2014|SV|PacBio||Brown U, Raphael||15|||||https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4253835/|https://github.com/raphael-group/multibreak-sv|
|PBHoney|2014|SV|PacBio||Baylor, Reid||16|||||http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-180|https://sourceforge.net/projects/pb-jelly/|
|BreakTrans|2013|SV|wgs/rna-seq|somatic|Wash U St Louis, Ding||10|||limits breakpoint graph to transcribed regions|adds RNA-seq data|systematically maps predicted gene fusions to structural rearrangements||https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-8-r87|http://bioinformatics.mdanderson.org/main/BreakTrans|
|MindTheGap|2014|insertions|c elegans||||19|||||||https://academic.oup.com/bioinformatics/article/30/24/3451/2422179/MindTheGap-integrated-detection-and-assembly-of|http://gatb.inria.fr/software/mind-the-gap/|
|Factera|2015|SV|||Stanford, Alizadeh||11|||fusion detection. soft-clip, kmer|trained on real tumor data|Fusion And Chromosomal Translocation Enumeration and Recovery Algorithm||https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4296148/|https://factera.stanford.edu/|
|denovolyzer|2015|SNV/indel|vcf|de novo|||12||||de novo filter|||https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4606471/||
|Triodenovo (former polymutt)|2015|SNV|vcfs|germ|Vanderbilt, Li|denovogear|8|||bayesian with flexible post hoc priors|works well for low coverage|||https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4410659/|http://genome.sph.umich.edu/wiki/Triodenovo#Download|
|MendelScan|2014|SNV/indel|vcfs exome|germ|U Texas|Inginuity Variant Analysis (Qiagen)|21||||exome|||https://www.ncbi.nlm.nih.gov/pubmed/24560519|https://github.com/genome/mendelscan|
|trioCaller|2013|SNV|vcfs|denovo|||34|||HMM||||https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530674/|http://genome.sph.umich.edu/wiki/TrioCaller|
|FamSeq|2012|SNV/indel|vcf |denovo|U Texas MD Anderson Cancer Center|gatk, samtools mpileup both unfiltered|24|||Bayesian network/markov chain-monte carlo||||http://www.pnas.org/content/110/10/3985.long|http://odin.mdacc.tmc.edu/~wwang7/FamSeqIndex.html|
|ForestDNM|2012|SNV/indel|vcf|denovo|||236||||de novo filter|||https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3712641/|http://sebatlab.ucsd.edu/index.php/software-data|
|Var-MD (available?)|2012|SNV/indel|vcfs exome|denovo|NHGRI, Boerkoel|n/a|22|||variety of mendelian/quality filters|sorts variants by disease likelihood|||https://www.ncbi.nlm.nih.gov/pubmed/22290570|https://research.nhgri.nih.gov/software/Var-MD/|
|kggseq|2012|SNV|mendelian filter/annotation engine|https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3326332/||

## Unsorted Callers

|caller|orig pub|caller class|input type|study|source|
|------|--------|------------|----------|-----|------|
|dindel|2011|indels|Sanger, slow|http://europepmc.org/abstract/MED/20980555|N/A|
|piCALL|2011|indel|scripps, exon|https://www.ncbi.nlm.nih.gov/pubmed/21653520||
|svseq 1 & 2|2011|SV|split-read|https://www.ncbi.nlm.nih.gov/pubmed/21994222||
|Spanner|2011|SV||https://www.ncbi.nlm.nih.gov/pubmed/21293372||
|SPLINTER|2010|SNV|SNPSeeker + indels|https://www.ncbi.nlm.nih.gov/pubmed/21041413/||
|bam2mpg|2010|SNV||https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2945191/||
|svmerge|2010|
|snvmix|2010|SNV|early cancer-specific, with low purity expectations|https://academic.oup.com/bioinformatics/article/26/6/730/245170/SNVMix-predicting-single-nucleotide-variants-from||
|cnD|2010|cnv||https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2820678/||
|dna sudoku|2009|SNV|overlapping pooling design|https://www.ncbi.nlm.nih.gov/pubmed/19447965/||
|breseq|2009|SNV|binomial model|https://www.researchgate.net/publication/232776704_Genome_evolution_and_adaptation_in_a_long-term_experiment_with_Escherichia_coli||
|maq|2009|SNV|bayesian model includes sequencing errors|https://www.ncbi.nlm.nih.gov/pubmed/18714091||
|SNPSeeker|2009|SNV|pooling. compares observed allele frequencies against the distribution of sequencing errors as measured by the Kullback Leibler (KL) distance|https://www.ncbi.nlm.nih.gov/pubmed/19252504/||
|modil|2009|indels||http://www.nature.com/nmeth/journal/v6/n7/full/nmeth.f.256.html||
|VariationHunter|2009|SV||http://genome.cshlp.org/content/19/7/1270.short||
|mrcanavar|2009|cnv||https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2875196/||
|maq|2008|SNV||https://www.ncbi.nlm.nih.gov/pubmed/18714091||
|PEM|2007|SV |paired-end|http://science.sciencemag.org/content/318/5849/420||
|PolyPhred|2006|SNV||https://www.ncbi.nlm.nih.gov/pubmed/16493422/||
|SNPDetector|2005|SNV||https://www.ncbi.nlm.nih.gov/pubmed/16261194/||
|novoSNP|2005|SNV||https://www.ncbi.nlm.nih.gov/pubmed/15741513/||
|ssahaSNP|2001|SNV||https://www.ncbi.nlm.nih.gov/pubmed/11591649/||
|polybayes|1999|SNV||https://www.ncbi.nlm.nih.gov/pubmed/10581034/||

## UNPUBLISHED CALLERS

|caller|orig pub|caller class|input type|study|source|
|------|--------|------------|----------|-----|------|
|caveman||||||
|Sentieon (LLC)||||https://peerj.com/preprints/1672.pdf|http://www.sentieon.com/products.html|
|gatk hc||||||
|gatk uv||||||
|samtools mpileup||||||
|BREPA||SV|||https://bitbucket.org/xianfan/brepa|
|snippy||||n/a|http://www.vicbioinformatics.com/software.snippy.shtml|
|cna-seq||||||
|copycat||||?|https://github.com/chrisamiller/copyCat|
|bassovac|2011|||n/a|http://tvap.genome.wustl.edu/tools/bassovac/|
|SoapSNV|2013|SNV/SV/indel|||http://soap.genomics.org.cn/SOAPsnv.html|
|glfTools|2010|SNV||n/a|http://genome.sph.umich.edu/wiki/GlfSingle|
|CNValidator|2012|CNV|||https://code.google.com/archive/p/cnvalidator/|
|SeqCNVCBS|2012|CNV||?|http://www.mybiosoftware.com/seqcnvcbs-1-0-scan-statistics-cnv-detection-cbs.html|
|BreakDown|2015|SV & VAF||https://scholarship.rice.edu/handle/1911/87870||
|kissnp||||https://hal.inria.fr/inria-00514887/||
|QuadGT|2013|SNV|exome bam trio (quad)|http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-S5-S3|http://www.iro.umontreal.ca/~csuros/quadgt/|
|excaliburSMD|2014||||https://github.com/cribioinfo/ExScaliburSMD|
|UPS-indel|2016|indel||http://ieeexplore.ieee.org/document/7802793/|https://sourceforge.net/projects/ups-indel/|
|splazers available?|2012|indel||https://academic.oup.com/bioinformatics/article/28/5/619/248213/Detecting-genomic-indel-variants-with-exact|?|
|sv-m structural variant machine|2012|indel|||https://www.bsse.ethz.ch/mlcb/research/bioinformatics-and-computational-biology/structural-variant-machine--sv-m-.html|
|takeabreak|2014|inversion break points||http://link.springer.com/chapter/10.1007%2F978-3-319-07953-0_10|https://colibread.inria.fr/software/takeabreak/|
|mogul|2012|||https://link.springer.com/chapter/10.1007/978-3-642-12683-3_23||
|concod|2016|sv||http://ieeexplore.ieee.org/abstract/document/7822495/||
|cnndel|2016|sv||http://ieeexplore.ieee.org/abstract/document/7822793/||
|sprites|2016|sv||https://academic.oup.com/bioinformatics/article/32/12/1788/1743630/Sprites-detection-of-deletions-from-sequencing||
|svmod|2016|sv||http://link.springer.com/article/10.1007/s00180-016-0674-2||