# The core code for the publication: Jun Wu and Changxu Fan, et al. Cell Reports. 2024
## TODO:
* add release numbers to referenced R packages

Sequence of analyses:
* Data processing starts at `sth/count_exonOnly/count.sh`, where cellranger was used to align the reads.
* Scripts under `fast_check` were run first to perform clustering on a single-sample basis
* Then, scripts under `sync_all` were used to generate R objects containing all samples with batch correction
* `day9_only` contains scripts that generated all figures in the paper. Although it's named "day9_only", it also contains code that generated plots for day11 post-5FU data.

About dead symlinks:
* R_scFanc and R_liteRnaSeqFanc links to the /R directory you can download from the corresponding repositories for these R packages, as mentioned above.
* Other dead symlinks point to large files/directories that github cannot handle. But you can generate them by running the code.
  
