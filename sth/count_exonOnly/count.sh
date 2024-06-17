#!/bin/bash

#SBATCH --array=1-4%4 --cpus-per-task=20 --mem=100G --time=10-00:00:00

sampleArray=("WT_rep2" "WT_rep3" "KO_rep2" "KO_rep3")
sample=${sampleArray[${SLURM_ARRAY_TASK_ID}-1]}
fastqDir=/scratch/twlab/fanc/jun/2022-01-25/fastq/

cellranger-arc count --id ${sample}_ExonOnly --reference /scratch/twlab/fanc/software/cellrangerArc_2.0/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
--localcores 18 --localmem 90 \
--gex-exclude-introns \
--libraries ${fastqDir}/${sample}/${sample}.csv 1>count.log.${sample} 2>&1

