# Plasticity_Acartia_hudsonica
Contains the code used to generate graphs and statistics for "Phenotypic Plasticity drives Seasonal Thermal Tolerance in a Baltic Copepod"

`Plasticity_Acartia_main_figures.R` is used to produce the main figures, it requires `data.csv` and `temperature.csv`
`temperature.csv` uses data from Hiebenthal et al. 2023 (https://doi.pangaea.de/10.1594/PANGAEA.963281)

`Plasticity_Acartia_stats.R` is used to produce the main statistics and statistic for individual plots, it requires `data.csv` and `temperature.csv`

`Plasticity_Acartia_supplemental_figures.R` is used to produce the supplemental figures and some corresponding stats, it requires `data.csv`, `temp_running.csv` and `CTD_data.csv`

The directory `genotyping` contains all files necessary to generate a phylogenetic tree of experimental animals

NCBI Accession Numbers for all animals genotyped in this study can be found in `NCBI_Accessions.txt`
