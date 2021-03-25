ExplorATE - *Explore Active Transposable Elements* -
====
### Description

**ExplorATE** (*Explore Active Transposable Elements*) is an R package for the exploration and identification of active transposons in RNA-seq data. The package offers functions to manipulate the RepeatMasker output files, and allows to discriminate TEs-coding transcripts from those repeats that are co-transcribed with genes coding non-transposon proteins. Through a simple pipeline you can solve overlaps of the repetitions that RepeatMasker cannot solve based on either higher score (`HS`), longer length (`LE`) or lower Kimura’s distances (`LD`). The transposons are finally annotated in the reference file using a selection criterion based on the percentage of identity, the percentage of the length of the repeat in accordance with the transcript, and a minimum length in base pairs (Wicker's rule or user's defined). In addition, a decoy file and a transcriptome salmon-formated are created to be used for indexing and quantification with Salmon. Finally, a function is incorporated to import the quantification estimates of the transcripts into the R environment for their subsequent differential expression analysis.


### Requirements

ExplorATE requires some previously installed packages, select those packages that are not available in your environment:

```{r eval=FALSE}

install.packages(c("stringr","seqinr","foreach","doParallel"))
install.packages(c("BiocManager","devtools")) 

BiocManager::install(c("readr","GenomicRanges", "IRanges","csaw", "edgeR","SummarizedExperiment","DESeq2", "tximport"))

```

### Install ExplorATE

Install ExplorATE from GitHub

```{r eval=FALSE}

devtools::install_github("FemeniasM/ExplorATEproject")

```

### Make input files

The pipeline requires a transposable elements annotations file and a transcripts annotations file. The TEs annotations file can be created with [RepeatMasker](http://www.repeatmasker.org) by comparing the transcriptome against curated libraries such as [Dfam](https://www.dfam.org/) (HMM library profile derived from Repbase sequences) and [Repbase](https://www.girinst.org/). Also, a *de novo* library can be created with [RepeatModeler2](http://www.repeatmasker.org/RepeatModeler/) that incorporate a module to the identification of LTR elements. Since TEs are highly variable between species, *de novo* annotations combined with consensus libraries reduce the age-related TEs bias. The transcripts annotations file, for example generated by [BLAST](https://blast.ncbi.nlm.nih.gov/), must be in `.outfmt6` format and all TEs-related proteins must be eliminated (*e. g.* reverse transcriptases, endonucleases, transposases, tyrosine recombinases, etc.). ExplorATE will only take the first two columns of the `.outfmt6` file. Also, a GFF3 file can be used to summarize the information for each function (optional). If working with a *de novo* transcriptome, the `.gff` file can be generated with [TransDecoder](https://github.com/TransDecoder/TransDecoder/wiki) keeping only the best ORF. Finally, one can also use a RepeatMasker alignment file (optional) `.align` to resolve overlaps. A shell [script](https://github.com/FemeniasM/ExplorATEproject/quick_start_guide/make_input_file.sh) is provided with the package to help run these programs.

### User guide

Check the [vignettes](https://github.com/FemeniasM/ExplorATEproject/tree/master/vignettes) and the [quick_start_guide](https://github.com/FemeniasM/ExplorATEproject/tree/master/quick_start_guide) with simulated data
