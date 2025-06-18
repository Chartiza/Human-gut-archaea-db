# PROJECT: Human gut archaea collection from Estonian population

✨ The aim of this project is to explore the archaeal diversity within the Estonian population. The final EstMB MAGdb Archaea-273 collection comprises 144 strains representing 21 archaeal species.
![alt text](db.png)
**Fig.1.** Visual representation of EstMB MAGdb Archaea-273 collection.

## Data availible for download here

- Paper draft: [Download link](https://www.biorxiv.org/content/10.1101/2024.07.06.602324v1.full)
- EstMB MAGdb Archaea-273 collection: [Download link](https://www.ebi.ac.uk/ena/browser/view/PRJEB81541)
- Functional profiling results: [Download link](https://figshare.com/articles/dataset/Prokka_annotation/29329166)

## Scriprs

Raw data processing and MAG reconstruction on cluster:
| Script Name                         | Description                                                                                                            |
| ----------------------------------- | ---------------------------------------------------------------------------------------------------------------------- |
| `Reads_preprocessing_RemoveLQ.sh`   | Removes low-quality reads from sequencing data using **fastp**.                                                        |
| `Reads_preprocessing_RemoveHost.sh` | Removes host-derived reads from sequencing data using **Bowtie2**.                                                     |
| `MAG_assembly.sh`                   | Assembles contigs from filtered reads using **MEGAHIT**.                                                       |
| `MAG_binning.sh`                    | Bins contigs into MAGs using **MetaBAT2**, **MaxBin2**, and **VAMB**, followed by refinement with **DAS Tool**.        |
| `MAG_clustering.sh`                 | Dereplicates MAGs at the species level using **dRep** (ANI threshold > 95%).                                           |
| `Taxonomic_annotation.sh`           | Assigns taxonomic classifications to MAGs using **GTDB-Tk v2** with **GTDB release 226**.                              |
| `MAGs_quality.sh`                   | Estimates MAG quality and completeness using **CheckM2**.                                                              |
| `Prevalence_and_abundance.sh`       | Estimates the prevalence and relative abundance of representative species in the Estonian population using **CoverM**. |
| `parse_CoverM_results.py`           | Aggregates CoverM output into summary tables across all samples.                                                       |
| `Functional_annotation.sh`          | Performs functional annotation of MAGs using **Prokka**.                                                               |

Follow-up Analysis (Jupyter Notebooks)

| Script Name                         | Description                                                                                                            |
| ----------------------------------- | ---------------------------------------------------------------------------------------------------------------------- |
| `Paper_figures_and_stats.ipynb`   | Generates all figures and summary statistics reported in the manuscript.                                                        |
