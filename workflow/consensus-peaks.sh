#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=30G
#SBATCH --time=8:00:00

#SBATCH --job-name=consensusPeak
#SBATCH --output=../logs/consensus.out

# set directories
cd scripts

cell_summit="../../data/rawdata/ATACSeq/cell_lines/summits"
#pdx_summit="../../data/rawdata/ATACSeq/PDXs/summits"
#cf_summit="/cluster/projects/bhklab/rawdata/BCaATAC/ATAC-cfDNA/summits"
#cf_all_summit="/cluster/projects/bhklab/projects/BCaATAC/BCa_ARCHE_Scoring/data/rawdata/ATAC-cfDNA/summits-allSamples"
tcga_summit="/cluster/projects/bhklab/rawdata/BCaATAC/TCGA/summits"

#cell_outdir="../../data/results/cell_lines/consensus"
#pdx_outdir="../../data/results/PDXs/consensus"
#cf_outdir="../../data/results/cfDNA/consensus"
#cf_all_outdir="../../data/results/cfDNA/consensus-allSamples"
tcga_cell_outdir="../../data/results/TCGA/consensus_cell"
#tcga_pdx_outdir="../../data/results/TCGA/consensus_pdx"

module load R
#Rscript consensus_corces.R $cell_summit $cell_outdir
#Rscript consensus_corces.R $pdx_summit $pdx_outdir
#Rscript consensus_corces.R $cf_all_summit $cf_all_outdir

Rscript consensus_corces.R $tcga_summit $cell_summit $tcga_cell_outdir
#Rscript consensus_corces.R $tcga_summit $pdx_summit $tcga_pdx_outdir