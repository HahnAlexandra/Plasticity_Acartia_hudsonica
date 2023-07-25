# pipeline to analyze mitotype data

This pipeline holds an automated method to analyze mitotype data for a tonsa roughly based on the method from of Figueroa et al., 2020: Phylogeography ofAcartia tonsaDana, 1849 (Calanoida: Copepoda)and phylogenetic reconstruction of the genusAcartiaDana, 1846.

`all_samples.fasta` contains the sequence data for all samples used in this analysis. See below for details of how to reproduce the analysis. 


---

View the help file with: `./mitotype_pipeline.sh --help` or `./mitotype_pipeline.sh -h`

To reproduce the pipeline as in the manuscript, run the following:

```bash
./mitotype_pipeline.sh input.fasta tree_plot

```

or generally for other data:


```bash
./mitotype_pipeline.sh your.fasta plot_name

```

where `your.fasta` contains your sequence data of unknown samples and `plot_name` is the name you want for your output plot.

All output will be directed to `output/`, which will be created if it doesn't already exist. 

I assume some R packages are installed: `treeio`, `tidyverse`, `ggtree`, `ape`, `stringr`.

This script requires a number of fasta files containing the reference sequences in addition to the scripts in the main directory. Do not modify these files. These files are:

- `reference_haplotypes.fasta`: the reference sequences of species to compare with
- `make_mr_bayes.R`: R script to make Mr Bayes input
- `plot_tree.R`: R script to generate final tree
 
You will also need to ensure the appropriate R packages are installed as well as MrBayes and muscle and their correct locations referenced in the scripts. These locations are currently hardcoded and need to be manually edited.
