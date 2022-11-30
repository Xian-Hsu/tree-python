#!/bin/bash

#SBATCH --job-name=tree_sim
#SBATCH -o slurm.out
#BATCH -J Python
#SBATCH -p cpu
#SBATCH --ntasks=1
#SBATCH --nodes=1
##SBATCH --nodelist c[13]
##SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=48:00:00
##SBATCH --export=NONE
#******************************************

# Please change the following accordingly 

#export EXEC=main.py
export EXEC="main.py"
#export FILEDIR="magp dens_pres_temp magfield magfield_cylindrical magpres pres3 velocity radial_velocity animation"

export INPUTS="classes.py constants.py energy.py evolve.py initial.py IO.py search_tree.py tree_structure.py visualize.py"
export SCR=/data/xianhsu/others/computational_astrophysics/project/plots/solar_01  # working directory
#export SCR=/data/xianhsu/CRBUB/mhd_injectiontest041_new_model_ref10_10field/velocity  # working directory

# The following directories are where the above files are located

export EXECDIR=/data/xianhsu/others/computational_astrophysics/project/plots
#export PARMDIR=/data/yhlin
export INPUTDIR=/data/xianhsu/others/computational_astrophysics/project/plots

#++++++++++++++++++++++++++++++++++++++++++

# These are the modules loaded in order to compile the FLASH code

#source ~/.bash_profile
module purge
module load python
source activate astro

# Create the scratch directory for the job and cd to it; make it private
if [ "$SCR" != "" ]; then
    echo Scratch directory is ${SCR}
    mkdir -p $SCR
    cd $SCR
    chmod o-rwx .
#    echo Creating output directory $FILEDIR ...
#    mkdir $FILEDIR
else
    echo Could not get scratch directory!
    exit
fi

cp $EXECDIR/$EXEC $SCR

# Get input files from home directory.
if [ "$INPUTS" != "" ]; then
    echo Getting files from home directory...
    for file in $INPUTS; 
    do
        echo Copying input file $file ...
        /bin/cp -r $INPUTDIR/$file .
    done
fi

cd $SCR

echo Running at `date`
echo about to run $EXEC

#for file in $EXEC;
#    do
#        echo about to run $file...
#        srun -c 1 python $SCR/$file  
#    done
srun -c 1 python $SCR/$EXEC
#srun -c 1 python test.py
#python $EXEC

echo Run complete at `date`
