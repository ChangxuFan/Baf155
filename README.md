# The core code for the publication: Jun Wu and Changxu Fan, et al. Cell Reports. 2024
## TODO:
* add release numbers to referenced R packages

Sequence of analyses:
* Data processing starts at `sth/count_exonOnly/count.sh`, where cellranger was used to align the reads.
* Scripts under `fast_check` were run first to perform clustering on a single-sample basis
* Then, scripts under `sync_all` were used to generate R objects containing all samples with batch correction
* `day9_only` contains scripts that generated all figures in the paper. Although it's named "day9_only", it also contains code that generated plots for day11 post-5FU data.

Stochasticity of the pipeline:
* We noticed that the doublet detection process has some stochasticity: if you run it twice, the doublet scores will not be exactly the same, although they are quite similar. Therefore, after doublet filtering, you are likely to end up with a slightly different set of cells. Unfortunately, such slight differences might cause ArchR to assign different cluster names. We offer the integrated Seurat/ArchR objects (referred to as `soi` and `aoi` in the code, respectively) at GSE240585, which contain the cells we obtained after our run of doublet filtering, and their associated cluster assignments.

About dead symlinks:
* R_scFanc and R_liteRnaSeqFanc links to the /R directory you can download from the corresponding repositories for these R packages, as mentioned above.
* Other dead symlinks point to large files/directories that github cannot handle. But you can generate them by running the code.
  
