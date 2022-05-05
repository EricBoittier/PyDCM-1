#!/bin/bash

#SBATCH --job-name=g_scn_2
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1600MB
#SBATCH --partition=vshort
#SBATCH --exclude=node[109-124]
#$ -S /bin/bash

# Define Parameters
hostname
export WORKDIR=/home/toepfer/PyDCM/input
export DATADIR=/home/toepfer/PyDCM/data
export INPFILE=gaussian_2_scn.com
export OUTFILE=gaussian_2_scn.out
export CHKFILE=gaussian_2_scn.chk
export FCKFILE=gaussian_2_scn.fchk
export DNSFILE=gaussian_2_scn_dens.cube
export ESPFILE=gaussian_2_scn_esp.cube

# Make scratch temporary directory
export STMPDIR=/scratch/$USER/jobs/tmp.$SLURM_JOBID
if [ -d $STMPDIR ]; then
  echo "$STMPDIR exists; double job start; exit"
  exit 1
fi
mkdir -p $STMPDIR

cd $STMPDIR

# Copy input file
cp $WORKDIR/$INPFILE $STMPDIR/$INPFILE

# Prepare Gaussian parameters
export g16root=/opt/cluster/programs/g16-c.01
source $g16root/g16/bsd/g16.profile
export GAUSS_SCRDIR=$STMPDIR

# Execute Gaussian jobs
$g16root/g16/g16 $STMPDIR/$INPFILE $STMPDIR/$OUTFILE
    
# Copy result file to output directory
cp $STMPDIR/$OUTFILE $DATADIR/$OUTFILE

# Delete output file if successfully copied
status=$?
if [ $status -eq 1 ]; then
    echo "$STMPDIR/$OUTFILE could not be transferred to $OUTPDIR/$OUTFILE"
fi
        
# Write formchk file
$g16root/g16/formchk $CHKFILE $FCKFILE

# Extract density cube file
$g16root/g16/cubegen 0 Density $FCKFILE $DNSFILE -2 h

# Extract esp cube file
$g16root/g16/cubegen 0 Potential $FCKFILE $ESPFILE -2 h

# Copy result to data directory
cp $CHKFILE $DNSFILE $ESPFILE $DATADIR/.

if [ $status -eq 0 ] 
then 
    # Delete temporary directory
    rm -rf $STMPDIR
else
    echo "$STMPDIR results could not be transferred to $DATADIR"
fi
