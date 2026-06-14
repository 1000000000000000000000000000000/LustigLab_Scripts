#!/bin/bash

#SBATCH -n 1
#SBATCH --cpus-per-task=10
#SBATCH --time=100:00:00
#SBATCH -J Q96EB6
#SBATCH --mem=0
#SBATCH -D /home/rpearson/AlphaFoldv2/alphafold-2.2.0/research/Q96EB6
#SBATCH --partition=gpus
#SBATCH --gres=gpu:2
#SBATCH --mem-per-gpu=256G

DATEVAR=$(date +"%Y-%m-%d")
echo --------------------------------------------------------------
echo Today is: $DATEVAR
echo --------------------------------------------------------------

source activate alphafold

echo The python being used is the following:
which python
echo Please make sure this python is in the env you want to use!
echo --------------------------------------------------------------

module load cuda-11.4
module load gnu/6.3.0
nvcc --version
nvidia-smi

echo --------------------------------------------------------------
echo XLA_PYTHON_CLIENT_PREALLOCATE is set to: $XLA_PYTHON_CLIENT_PREALLOCATE
export XLA_PYTHON_CLIENT_PREALLOCATE=false
echo XLA_PYTHON_CLIENT_PREALLOCATE is set to: $XLA_PYTHON_CLIENT_PREALLOCATE

echo --------------------------------------------------------------
echo The TF_FORCE_GPU_ALLOW_GROWTH is set to: $TF_FORCE_GPU_ALLOW_GROWTH
export TF_FORCE_GPU_ALLOW_GROWTH=true
echo TF_FORCE_GPU_ALLOW_GROWTH is set to: $TF_FORCE_GPU_ALLOW_GROWTH
echo --------------------------------------------------------------

# The Following is required for larger proteins >1200aa
export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_MEM_FRACTION=0.5
export XLA_PYTHON_CLIENT_ALLOCATOR=platform

#export CUDA_DEVICE_ORDER="PCI_BUS_ID"

echo $CUDA_DEVICE_ORDER
echo $CUDA_VISIBLE_DEVICES

cd /home/rpearson/AlphaFoldv2/alphafold-2.2.0/

srun bash ./run_alphafold.sh \
  -d /home/rpearson/AlphaFoldv2/alphafold_databases \
  -o /home/rpearson/AlphaFoldv2/alphafold-2.2.0/research/Q96EB6 \
  -f /home/rpearson/AlphaFoldv2/alphafold-2.2.0/research/Q96EB6/Q96EB6.fasta \
  -t $DATEVAR
