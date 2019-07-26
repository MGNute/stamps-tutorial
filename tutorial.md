SEPP and TIPP Tutorial
======================
Let's analyze metagenomic datasets using SEPP (SATe-Enabled Phylogenetic Placement) and TIPP (Taxonomic Identification and Phylogenetic Profiling)! Recall that both SEPP and TIPP require the following inputs:
+ Set of query sequences, i.e., fragments/reads of unknown origin
+ Reference alignment and tree or taxonomy

**Query sequences**   
In this tutorial, we will analyze (subsets of) metagenomic datasets acquired through (16S) amplicon and shotgun sequencing. The 16S sample ([SRR1219742](https://www.ncbi.nlm.nih.gov/biosample/SAMN02725485) -- Lemur Vaginal Sample -- 454 GS FLX Titanium) comes from Yildirim et al., 2014. The shotgun sample ([SRR059421](https://www.ncbi.nlm.nih.gov/sra/SRX022983[accn]) -- Human Stool Sample -- Illumina Genome Analyzer II) comes from the [Human Microbiome Project](http://www.hmpdacc.org). All runs were downloaded from the NCBI database using this [fastq-dump command](tools/fastq_dump.sh).

**Reference alignments and trees**   
Phylogenetic placement (with SEPP) and taxonomic identification (with TIPP) of 16S samples can be performed using a [reference alignment and tree](refpkgs/RDP_2016_Bacteria.refpkg) built on the 11,988 bacterial reference sequences from the [Ribosomal Database Project (RDP) database](https://rdp.cme.msu.edu/) by selecting >1200 site-length, type isolates with quality "Good", and the NCBI taxonomy. 

In this tutorial, we will use only a subset of the sequences from the RDP Bacteria reference as the size of the reference as well as the number of query sequences affect the running time of SEPP and TIPP. For example, the number of sequences in the reference alignment and the alignment subset size determine how many profile HMMs must be built over the reference alignment. Then each query sequence in the sample must be aligned (and scored) to each of these profile HMMs. To reduce the size of the reference, TIPP was run on the first 2,500 sequences from the 16S sample (SRR1219742). The majority of classified reads (929 reads) were identified as Clostridia. This subset of reads can now be classified at the family, genus, and species levels using TIPP with the RDP Bacteria reference constrained to the 707 sequences in the Clostridia class.

Taxonomic identification and abundance profiling (with TIPP) of whole shotgun samples can be performed using a collection of reference alignment and trees built from the 30 marker genes (i.e., "use universal housekeeping genes that are unlikely to undergo duplication or horizontal gene transfer") in MetaPhyler.

Part I: Taxonomic Identification using TIPP
-------------------------------------------
To begin, log in to your the cloud instance using the [Jupyter link provided by the STAMPS organizers](https://hackmd.io/@astrobiomike/stamps2019). Once you are into the Jupyter main screen, click the "terminal" link to start a command-line session. Once that has started, clone the repository for this tutorial:
```
git clone https://github.com/MGNute/stamps-tutorial.git
```
Change into the tipp directory,
```
cd stamps-tutorial/tipp
```
and run TIPP using the following command
```
python /opt/sepp/run_tipp.py \
    -a ../refpkgs/RDP_2016_Clostridia.refpkg/pasta.fasta \
    -t ../refpkgs/RDP_2016_Clostridia.refpkg/pasta.taxonomy \
    -r ../refpkgs/RDP_2016_Clostridia.refpkg/RAxML_info.taxonomy \
    -tx ../refpkgs/RDP_2016_Clostridia.refpkg/taxonomy.table \
    -txm ../refpkgs/RDP_2016_Clostridia.refpkg/species.mapping \
    -A 100 \
    -P 1000 \
    -at 0.95 \
    -pt 0.95 \
    -f ../samples/16S/SRR1219742_RDP_2016_Clostridia.fasta \
    -o TIPP-RDP-CLOSTRIDIA-95-SRR1219742 \
    --tempdir tmp \
    --cpu 2
```
This will take 5-6 minutes to finish. In the meantime, let's break down the command. The first five options specify files included for the reference
+ `-a [`[`reference multiple sequence alignment -- fasta format`](refpkgs/RDP_2016_Clostridia.refpkg/pasta.fasta)`]`
+ `-t [`[`reference taxonomy -- newick format`](refpkgs/RDP_2016_Clostridia.refpkg/pasta.taxonomy)`]`
+ `-r [`[`reference tree model parameters -- RAxML info file`](refpkgs/RDP_2016_Clostridia.refpkg/RAxML_info.taxonomy)`]`
+ `-tx [`[`mapping taxonomic id to taxonomy information -- csv`](refpkgs/RDP_2016_Clostridia.refpkg/taxonomy.table)`]`
+ `-txm [`[`mapping sequence names to taxonomic IDs -- csv`](refpkgs/RDP_2016_Clostridia.refpkg/species.mapping)`]`

The next two options specify the decomposition of the reference alignment and tree into subsets.
+ `-A [alignment subset size]`
+ `-P [placement subset size]`

TIPP was run with an alignment subset size of 100 (slightly less than 10% of the Clostridia reference package) and a placement subset size of 1000 (greater than the entire Clostridia reference package). Recall that running SEPP/TIPP with larger placement subset sizes can increase accuracy but is more computationally intensive. The default alignment subset size follows the 10% rule (i.e., 10% of the placement subset). 

The next two options specify the support thresholds used by TIPP.
+ `-at [alignment support threshold]`
+ `-pt [placement support threshold]`

TIPP was run with support thresholds of 0.95, which is the default. (More recent experiments have suggested this might be conservative, so relaxing this is an option to get broader classification at the risk of a slight increase in false-positives.)

The next two options specify the input and output.
+ `-f [`[`fragment file -- fasta`](samples/16S/SRR1219742_RDP_2016_Clostridia.fasta)`]`
+ `-o [prefix of output files]`

The final options are set specifically for STAMPS tutorial to prevent temporary files from being written all over the MBL servers and limit the number of CPUs per user.

To see all of the [TIPP options](tipp-help.md), run the following command. There are a few practical pointers on these options in the [reference slides](https://github.com/MGNute/stamps-tutorial/blob/master/sepp_tipp_tutorial_slides.pdf) for this tutorial.
```
python /opt/sepp/run_tipp.py -h
```
By now TIPP may have finished and written the following files
+ [classification information -- csv](tipp/out/TIPP-RDP-CLOSTRIDIA-95-SRR1219742_classification.txt)
+ [phylogenetic placement information -- json](tipp/out/TIPP-RDP-CLOSTRIDIA-95-SRR1219742_placement.json)
+ [alignment on both the reference and query sequences -- fasta](tipp/out/TIPP-RDP-CLOSTRIDIA-95-SRR1219742_alignment.fasta.gz)

The classification file shows the support of classifying sequences at each taxonomic rank. Check out the support for each read classified at the species level
```
grep ",species," TIPP-RDP-CLOSTRIDIA-95-SRR1219742_classification.txt
```
Computing the number of reads classified at each taxonomic rank
```
python ../tools/restructure_tipp_classification.py \
    -i TIPP-RDP-CLOSTRIDIA-95-SRR1219742_classification.txt \
    -o FINAL-TIPP-RDP-CLOSTRIDIA-95-SRR1219742
```
and examining the read count for species-level classification
```
cat FINAL-TIPP-RDP-CLOSTRIDIA-95-SRR1219742_species.csv
```
shows the majority of reads are unclassified (545 reads). Classified reads are largely Fastidiosipila sanguinis (213 reads) and Anaerovorax odorimutans (144 reads). What do read counts look like at the genus and family level? 

*Before moving on, repeat this portion of the tutorial running TIPP with a lower alignment/placement support threshold (e.g., 0.50). What do the support values look like for reads classified at the species level? How does the number of reads unclassified at the species level compare to TIPP run with an alignment/placement support threshold of 0.95?*

**NOTE: In general, SEPP/TIPP should be run on reads and their reverse complement; however, this has already been handled for this tutorial.**

Part I.a: Visualizing the TIPP Placements
------------------------------------------
Remember that part of TIPP involves using SEPP to place the reads into a *taxonomy* (i.e., a phylogenetic tree that is forced to conform to a pre-determined taxonomy). Sometimes it is helpful to visualize the placement of the full (or near-full) sample at once. This can sometimes make the composition of each sample clearer in a way that a table of percentages can't always.

Here are two pictures of the reference taxonomy that we used as the backbone for TIPP. The first labels each leaf by its [order](images/RDP_2016_Clostridia_taxonomy_by_order.pdf), and the second goes further and labels them by [family](images/RDP_2016_Clostridia_taxonomy_by_family.pdf). This gives us an idea of what the backbone taxonomy looks like, for reference when we look at the placement. 

Now here is a picture of where the reads that *were* classified were [placed in the taxonomy](images/TIPP-RDP-CLOSTRIDIA-95-SRR1219742_phylo_heat_map.pdf). In this picture, a more vivid (or "saturated") color in one spot means more reads were placed there, whereas a faint color means fewer reads. The hue (e.g. red/orange/yellow/...) indicates the branch length. Recall that these reads were taken from a lemur sample.

Compare this heat map to the reference image by [family](images/RDP_2016_Clostridia_taxonomy_by_family.pdf), and take a look at the TIPP classifications by family:
```
cat FINAL-TIPP-RDP-CLOSTRIDIA-95-SRR1219742_family.csv
```
- *What is the most abundant family in this sample?*
- *Is that family well represented or is there one subgroup that is a heavy majority? What about the second most abundant?*

Now let's go down to the genus level. [Here is a partially labelled picture](https://github.com/MGNute/stamps-tutorial/blob/master/images/RDP_2016_Clostridia_taxonomy_by_genus_partial.pdf) of the reference tree showing several key genera. Compare that to TIPP classifications at the genus level:
```
cat FINAL-TIPP-RDP-CLOSTRIDIA-95-SRR1219742_genus.csv
```
- *What is the most abundant genus in this sample?* 
- *Where are the reads in the heat map located in relation to the leaves? What does that tell you about these microbes?*
- *What about the second most abundant?*


Part II: Phylogenetic Placement using SEPP
------------------------------------------
Now let's take a closer look at some select reads using phylogenetic placement. Based on the TIPP classification (using an alignment/placement support threshold of 0.50), many sequences were identified as the Ruminococcaceae family including
+ GEQJ1S112HF5CU
+ GEQJ1S110GHR11
+ GEQJ1S110GEAZV 
+ GEQJ1S112HNELY  
+ GBEHU2E07D5RLY

You can examine the support at which these reads are classified at the family level by using grep, e.g.,
```
grep "GEQJ1S112HF5CU" TIPP-RDP-CLOSTRIDIA-95-SRR1219742_classification.txt 
```

SEPP can be used to place these five query sequences into the RDP Bacteria reference package; however, we will use SEPP with an RDP Bacteria reference package constrained to the Ruminococcaceae family (55 sequences) for visualization purposes. Change into the sepp directory
```
cd ../sepp
```
and run SEPP using the following command
```
python /opt/sepp/run_sepp.py \
    -a ../refpkgs/RDP_2016_Ruminococcaceae.refpkg/pasta_labeled.fasta \
    -t ../refpkgs/RDP_2016_Ruminococcaceae.refpkg/pasta_labeled.taxonomy \
    -r ../refpkgs/RDP_2016_Ruminococcaceae.refpkg/RAxML_info.taxonomy \
    -A 25 \
    -P 100 \
    -f ../samples/16S/READS-Ruminococcaceae.fasta \
    -o SEPP-RDP-RUMINO-READS \
    --tempdir tmp \
    --cpu 2
```
Both the command and the output are nearly identical to TIPP; however, alignment and placement support thresholds are not specified, and the classification file is not written. Use the [placement file from SEPP](sepp/out/SEPP-RDP-RUMINO-READS_placement.json) to rank the five reads by the branch length connecting the read to the Ruminococcaceae tree. You may need read more about the json file format [here](https://matsen.github.io/pplacer/generated_rst/pplacer.html) (search for JSON format specification) or [here](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0031009). 

For this example, we will try out some other tools for visualizing phylogenetic placement. We will first be converting the placement file into a tree format (e.g., newick or xml) with [guppy](https://matsen.github.io/pplacer/generated_rst/guppy.html)
```
/opt/sepp/.sepp/bundled-v4.3.8/guppy tog \
    --xml \
    SEPP-RDP-RUMINO-READS_placement.json
```
Ordinarily you would download the xml file onto your personal computer. Since we are logged in Jupyter this year, copying this to your local machine is more trouble than it's worth, so download [this pre-baked version](https://www.dropbox.com/s/sbsi4xthmg2rgp7/SEPP-RDP-RUMINO-READS_placement.tog.xml?dl=0) instead. You can get it right from a web browser, or you can use wget from the terminal if you are using Linux:
```
wget https://www.dropbox.com/s/sbsi4xthmg2rgp7/SEPP-RDP-RUMINO-READS_placement.tog.xml?dl=0
```
Now open the file with your favorite tree viewer. [EvolView](http://www.evolgenius.info/evolview) can be used to visualize the placement of query sequences in the reference tree with [colored branches](http://evolview.codeplex.com/wikipage?title=DatasetBranchColor) and [colored leaves](https://evolview.codeplex.com/wikipage?title=DatasetLeafColor). After the uploading the file, hover over "Annatotation upload", click on the buttons that are second (branch color) and third (leaf color) from the left, and add the text
```
GEQJ1S112HF5CU red ad   
GEQJ1S110GHR11 red ad   
GEQJ1S110GEAZV red ad   
GEQJ1S112HNELY red ad   
GBEHU2E07D5RLY red ad   
```
**NOTE: that the above spaces need to be tabs!**

Examining the [cladogram](images/rumino-cladogram.pdf), we can notice that GEQJ1S112HNELY was placed closer to the root than the GEQJ1S112HF5CU, which was placed sister to Saccharofermentans acetigenes -- however based on the [branch length](images/rumino-phylogeny.pdf) GEQJ1S112HF5CU may not necessarily have very high sequence identity to Saccharofermentans accetigenes. Branch lengths can also be very short, see placement of the query sequence GEQJ1S112HN8VO on the Heliobacterium reference package [here](https://github.com/MGNute/stamps-tutorial/blob/master/images/helio-phylogeny.pdf).

Use the [alignment file from SEPP](sepp/out/SEPP-RDP-RUMINO-READS_alignment.fasta) to compare the read GEQJ1S112HF5CU to the reference sequence Saccharofermentans_acetigenes_1. Extract the two sequences into a new fasta file
```
grep -A1 "GEQJ1S112HF5CU" SEPP-RDP-RUMINO-READS_alignment.fasta > SEPP-RDP-RUMINO-READS_subset.fasta
grep -A1 "Saccharofermentans_acetigenes_1" SEPP-RDP-RUMINO-READS_alignment.fasta >> SEPP-RDP-RUMINO-READS_subset.fasta
```
and download it onto your personal computer. Again, copying it is more trouble that it's worth, so download it from [here](https://www.dropbox.com/s/nailtpu67rt8dkv/SEPP-RDP-RUMINO-READS_subset.fasta?dl=0), or use the terminal if you are using Linux:
```
wget https://www.dropbox.com/s/nailtpu67rt8dkv/SEPP-RDP-RUMINO-READS_subset.fasta?dl=0
```
[MSAViewer](http://msa.biojs.net) can be used to visualize the multiple sequence alignment. Scroll down and click the little arrow icon under "Use It". Then click "Import" followed by "From file".

*Before moving on, let's consider the relationship between alignment, placement, and classification. Use the [cladogram](images/helio-cladogram.pdf) to identify reference sequences near GEQJ1S112HN8VO (e.g., Heliobacterium_modesticaldum_11 -- taxon ID is in the TIPP classification file 498761 -- and the sequence ID is S000469502 -- grep for S000469502). Go back to the TIPP directory, and extract these sequences from the alignment file from TIPP. Visualize the alignment, and compare it to the visualization of GEQJ1S112HF5CU. Now examine the placement file from TIPP. What are branch lengths and maximum likelihood scores for placements of GEQJ1S112HN8VO onto the Clostridia reference tree? Based on this alignment and placement information, discuss TIPP (with 0.50 support thresholds) classifying GEQJ1S112HN8VO as Heliobacterium modesticaldum Ice1 (below species level) versus TIPP (with 0.95 support thresholds) classifying GEQJ1S112HN8VO as Clostridiales order and Unclassified at the family, genus, and species levels.*

**JUST A REMINDER: Small reference alignments and trees are used in this tutorial to save time and make visualization easier; however, the benefits of using SEPP/TIPP are greatest when trees have a large evolutionary diameters and when query sequences are potentially novel or divergent from the any reference sequences. New tools for visualizing phylogenetic placements in the graphical format used in Part I are on the way. A *very* preliminary version of the code to generate them is available on [Github](https://github.com/MGNute/pican_pi) and an accompanying manuscript has been submitted.**

Part III: Phylogenetic (Abundance) Profiling with TIPP
------------------------------------------------------
All prior analyses are on 16S, which is not single copy. TIPP can be used for phylogenetic (abundance) profiling by using a collection of marker genes (i.e., single copy universal genes) as reference alignments and trees. First, BLAST is used to identify whether a read is a match for a specific marker gene. If so, TIPP is used to classify the read. To run this analysis (in the future), create an output directory
```
mkdir TIPP-COGS-95-SRR059420
mkdir TIPP-COGS-95-SRR059420/markers
```
and run TIPP using the following command
```
python /opt/sepp/run_abundance.py \
    -G cogs \
    -c /opt/sepp/.sepp/tipp.config \
    -at 0.95 \
    -pt 0.95 \
    -f ../samples/shotgun/SRR059420_pass_1_1-25000.fasta \
    -d TIPP-COGS-95-SRR059420 \
    --tempdir tmp \
    --cpu 2
```
The above command includes the following new options,
+ `-G [whether marker genes or COGs should be used]`
+ `-c [configuration file for TIPP]`
+ `-d [name of output directory]`

and the [output](tipp/out/TIPP-95-COGS-SRR059420) will include abundance profile for each taxonomic rank, for example, the [genus-level abundance profile](tipp/out/TIPP-95-COGS-SRR059420/abundance.genus.csv) shows that the 95% of reads (that match to a COG) are classified as Bacteroides. Finally, the [markers folder](https://github.com/MGNute/stamps-tutorial/tree/master/tipp/out/TIPP-95-COGS-SRR059420/markers) contains the output from running TIPP on each of the markers.

*Before leaving Woods Hole, consider using the commands in this portion of the tutorial to analyze a shotgun dataset of interest to you. Let us know how it goes!* **Thank you for taking the time to do this tutorial!**

Part IV (Bonus Round): Multiple Sequence Alignment (MSA) with PASTA
------------------------------------------------------
Since TIPP and SEPP rely on a reference package, essentially a set of full-length genes that have been aligned along with a phylogenetic tree, it is helpful to run through an example of producing an alignment, which we will do with PASTA. PASTA includes a useful GUI, although it is sometimes easier to run through the command line.

The data for this is a set of 96 amino acid sequences pulled from NCBI. (Note: since SEPP and TIPP are used for placing sequencing reads, a more applicable data set would use DNA sequences, but it is still a helpful example.) The sequences are DNA Gyrase A (gyrA, for short) and can be found [here](https://www.ncbi.nlm.nih.gov/Structure/cdd/COG0188).

Lets first move to the subfolder for the PASTA example:
```
cd ../pasta
```
Note the files containing the unaligned sequences `gyrA_raw.fasta` and the MSA provided by NCBI `gyrA_ncbi_aln.fasta`. We are going to give the former to PASTA as input. To run PASTA on these sequences, run the following:
```
python /opt/pasta-code/pasta/run_pasta.py 
    -i gyrA_raw.fasta \
    --temporaries "./tmp/" \
    -d Protein \
    -j gyrA_pasta \
    --num-cpus $(nproc) \
    -o "./out"
```

This should take approximately 2 minutes on the Cloud instances. Here, we are specifying the folder `out/` as the destination for the various files that PASTA creates. Once it has completed, the file ending with `.aln` will be non-empty, and that is the final alignment.

When PASTA is done, let's do a quick comparison of the alignment provided by NCBI and the alignment we've just created. Run these two commands to get a few summary statistics about each alignment:
```
python ../tools/alignment_stats.py out/gyrA_pasta.marker001.gyrA_raw.aln
python ../tools/alignment_stats.py gyrA_ncbi_aln.fasta
```
Here are some questions to consider:
+ Which alignment has more gaps?
+ Which alignment finds more matches between sites?
+ Which alignment is *better*? Why?

Since this is biological sequence data, that last one is kind of a trick question. Presumably the alignment provided by NCBI has been curated to ensure that no obvious homologous sites are missed, but humans make mistakes too. If there were 10,000 sequences instead of 96, how accurate would you expect a curated alignment to be? 

The only way to say whether one alignment is "better" than another is to compare them on simulated sequences, where the true alignment is known from the simulation. Here is a link to [some simulated alignment data](https://sites.google.com/eng.ucsd.edu/datasets/alignment/sate-i) that we have used in the past. The disadvantage of simulated data is that it may not reflect realistic biological sequence data. In one recent paper, we found a notable [divergence in aligner performance](https://academic.oup.com/sysbio/article/68/3/396/5133545) between biological and simulated data. 

Citations
---------
Cole, J. R., Q. Wang, J. A. Fish, B. Chai, D. M. McGarrell, Y. Sun, C. T. Brown, A. Porras-Alfaro, C. R. Kuske, and J. M. Tiedje. 2014. Ribosomal Database Project: data and tools for high throughput rRNA analysis. *Nucl. Acids Res.* 42(Database issue):D633-D642. doi:[10.1093/nar/gkt1244](https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gkt1244)

Liu, B., T. Gibbons, M. Ghodsi, T. Treangen, and M. Pop. Accurate and fast estimation of taxonomic profiles from metagenomic shotgun sequences. 2011. *BMC Genomics* 12(Suppl 2):S4. doi:[10.1186/1471-2164-12-S2-S4](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-12-S2-S4)

Mirarab, S., N. Nguyen, and T. Warnow. (2012). SEPP: SATe-Enabled Phylogenetic Placement. *Proceedings of the 2012 Pacific Symposium on Biocomputing (PSB 2012)* 17:247-258. doi:[10.1142/9789814366496_0024](http://www.worldscientific.com/doi/abs/10.1142/9789814366496_0024)

Nguyen, N., S. Mirarab, B. Liu, M. Pop, and T. Warnow (2014). TIPP:Taxonomic Identification and Phylogenetic Profiling. *Bioinformatics* 30(24):3548-3555. doi:[10.1093/bioinformatics/btu721](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btu721)

Yildirim, S., C. J. Yeoman, S. C. Janga, S. M. Thomas, M. Ho, and S. R. Leigh; Primate Microbiome Consortium, White BA4, Wilson BA2, Stumpf RM3. 2014. Primate vaginal microbiomes exhibit species specificity without universal Lactobacillus dominance. *The ISME Journal* 8(12):2431-44.
doi:[10.1038/ismej.2014.90](http://www.nature.com/ismej/journal/v8/n12/full/ismej201490a.html)

Nute M., Saleh, E., and Warnow, T. (2018). Evaluating statistical multiple sequence alignment in comparison to other alignment methods on protein data sets. *Systematic Biology*, Volume 68, Issue 3, May 2019, Pages 396-411, [https://doi.org/10.1093/sysbio/syy068](https://academic.oup.com/sysbio/article/68/3/396/5133545)

Mirarab, S., Nguyen N., Guo, S., Wang, L.-S., Kim, J., and Warnow, T. (2014) PASTA: Ultra-Large Multiple Sequence Alignment for Nucleotide and Amino-Acid Sequences. *Journal of Computational Biology* 22(5):377-386.
