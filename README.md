# The core code for the publication: Jun Wu and Changxu Fan, et al. Cell Reports. 2024
## Dependencies:
The analysis workflow heavily depends on the following packages I wrote:
* [scFanc](https://github.com/ChangxuFan/scFanc) Release v1.0.1
* [liteRnaSeqFanc](https://github.com/ChangxuFan/liteRnaSeqFanc) Release v1.0.1

To begin analysis, you need to copy/softlink the .R files from these packages to the corresponding directories (R_scFanc, R_liteRnaSeqFanc).

Although not critical, some functions in this repo might call functions from other packages I wrote, such as `utilsFanc`, `abaFanc2`, etc. These can also be found in the corresponding github repos under the account `ChangxuFan`. 
In addition, some stand-alone scripts (`*.R`, `*.sh`) are used, and can be found in the `R_for_bash` or `scripts` repos under the account `ChangxuFan`.

## Sequence of analyses:
* Data processing starts at `sth/count_exonOnly/count.sh`, where cellranger was used to align the reads.
* Scripts under `fast_check` were run first to perform clustering on a single-sample basis
* Then, scripts under `sync_all` were used to generate R objects containing all samples with batch correction
* `day9_only` contains scripts that generated all figures in the paper. Although it's named "day9_only", it also contains code that generated plots for day11 post-5FU data.
* In general, scripts should be run sequentially, with step1.1 before step1.2. However, this is not true for scripts labeled "revisionx.x", because they were numbered based on reviewers' questions.

## Stochasticity of the pipeline:
* We noticed that the doublet detection process has some stochasticity: if you run it twice, the doublet scores will not be exactly the same, although they are quite similar. Therefore, after doublet filtering, you are likely to end up with a slightly different set of cells. Unfortunately, such slight differences might cause ArchR to assign different cluster names. We offer the integrated Seurat/ArchR objects (referred to as `soi` and `aoi` in the code, respectively) at GSE240585, which contain the cells we obtained after our run of doublet filtering, and their associated cluster assignments.

## About dead symlinks:
* R_scFanc and R_liteRnaSeqFanc link to the /R directories under scFanc and liteRnaSeqFanc, respectively. You can download them from the corresponding repositories, as mentioned above.
* Other dead symlinks point to large files/directories that github cannot handle. But you can generate them by running the code.
