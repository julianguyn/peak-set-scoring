#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=30G
#SBATCH --time=8:00:00

#SBATCH --mail-user=jlm.nguyen@mail.utoronto.ca
#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --job-name=consensusPeak
#SBATCH --output=../logs/consensus.out

# set directories
cd scripts

cell_summit="../../data/rawdata/ATACSeq/cell_lines/summits"
pdx_summit="../../data/rawdata/ATACSeq/PDXs/summits"
pdo_summit="../../data/rawdata/ATACSeq/PDO/summits"
#cf_summit="/cluster/projects/bhklab/rawdata/BCaATAC/ATAC-cfDNA/summits"
#cf_all_summit="/cluster/projects/bhklab/projects/BCaATAC/BCa_ARCHE_Scoring/data/rawdata/ATAC-cfDNA/summits-allSamples"
tcga_summit="/cluster/projects/bhklab/rawdata/BCaATAC/TCGA/summits"
#pancancer_summit="/cluster/projects/bhklab/projects/BCaATAC/peak-set-scoring/data/rawdata/TCGA-pancancer/summits"

#cell_outdir="../../data/results/cell_lines/consensus"
#pdx_outdir="../../data/results/PDXs/consensus"
#cf_outdir="../../data/results/cfDNA/consensus"
#cf_all_outdir="../../data/results/cfDNA/consensus-allSamples"
#tcga_cell_outdir="../../data/results/TCGA/consensus_cell"
#tcga_pdx_outdir="../../data/results/TCGA/consensus_pdx"
#tcga_cell_pdx_outdir="../../data/results/TCGA_cell_PDX/consensus/"
tcga_cell_pdx_pdo_outdir="../../data/results/TCGA_cell_PDX_PDO/consensus/"
#tcga_pancancer_outdir="../../data/results/TCGA-pancancer/consensus"

#Rscript consensus_corces.R $cell_summit $cell_outdir
#Rscript consensus_corces.R $pdx_summit $pdx_outdir
#Rscript consensus_corces.R $cf_all_summit $cf_all_outdir

conda activate ARCHEscoring

Rscript consensus_corces.R $tcga_summit $cell_summit $pdx_summit $pdo_summit $tcga_cell_pdx_pdo_outdir
#Rscript consensus_corces.R $pancancer_summit $tcga_pancancer_outdir