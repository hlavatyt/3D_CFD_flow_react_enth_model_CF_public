#!/bin/bash
#PBS -S /bin/bash
#PBS -M hlavatyo@vscht.cz
#PBS -m bae
#PBS -N reactionCaseTest_uIn_2.040
#PBS -q batch
#PBS -l nodes=1:ppn=8
#PBS -l walltime=300:00:00

# specify the case folder
caseDir='reactionCaseTest_uIn_2.040'

# copy the necessary files to scratch
rsync -a ~/$caseDir/ /scratch/hlavatyo/$caseDir/

# set the environment for OpenFOAM
export PATH=/usr/lib64/mpi/gcc/openmpi/bin:$PATH
export OMPI_MCA_btl="sm,self"

export FOAM_INST_DIR=/opt/OpenFOAM/
source /opt/OpenFOAM/OpenFOAM-4.x/etc/bashrc WM_LABEL_SIZE=64 WM_COMPILER_TYPE=ThirdParty FOAMY_HEX_MESH=yes

# run the calculation
cd /scratch/hlavatyo/$caseDir
bash Allrun-parallel

# copy all the results back to the home folder (no need for the 'processor*' folders)
cd ~
rsync -a /scratch/hlavatyo/$caseDir/ ~/$caseDir/ --exclude "processor*"
