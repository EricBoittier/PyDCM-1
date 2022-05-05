#!/usr/bin/python

# Basics
import os
import numpy as np

# Load FDCM modules
from pydcm import Scan, DCM, mdcm

# Write dictionary for scan class
config_scan = {
    "system_label"              : "scn",
    "system_file_coord"         : "sources/scn.xyz",
    "system_file_type"          : "xyz",
    "system_total_charge"       : -1,
    "system_spin_multiplicity"  : 1,
    "scan_dofs"                 : [
        [0, 1, 2]
        ],
    "scan_steps"                : [
        [179, 175, 171, 167, 163, 159, 155, 151, 147, 143, 139, 135, 131]
        ],
    "scan_qm_program"           : "Gaussian",
    "scan_qm_method"            : "MP2",
    "scan_qm_basis_set"         : "aug-cc-pVTZ",
    "scan_constrained_opt"      : True,
    "scan_parallel_tasks"       : 4,
    "scan_cpus_per_task"        : 4,
    "scan_memory_per_task"      : 1600,
    "scan_overwrite"            : False,
    "scan_time_check_tasks"     : 10
    }

file_config_scan = "config_scan_scn.txt"

# Initialize scan class either by:
scan = Scan()
scan.initialize_scan(config_scan)
# or direct:
#scan = Scan(config=config_scan)

# After initialization, a config file is written that can be used to read config
#scan.initialize_scan(file_config_scan)

# Print documentation of config dictionary
scan.print_doc()

# Optional: Prepare all input files for scan to check for correct creation
scan.prepare_scan()

# Execute scan (prepare_scan, submit jobs, evaluate_scan) according to 
# config_scan definition
scan.execute_scan()

# Optional: Evaluate results, already done at the end of execute_scan
scan.evaluate_scan()

# Get final list of ESP and density cube files
dnslist = scan.get_files_cube_dens()
esplist = scan.get_files_cube_esp()

# Identify energetically lowest scan step
scan_smin = np.nanargmin(scan.get_potential())



#------------------


# MDCM:
#-------

# Just test initialization of still empty DCM class
fdcm = DCM()


mdcm_cxyz = "sources/8-charges_mike.xyz"
mdcm_clcl = "sources/8-charges_mike.dcm"

# Prepare some cube file list
scan_fesp = [
    "data/gaussian_0_scn_esp.cube", 
    "data/gaussian_1_scn_esp.cube",
    "data/gaussian_2_scn_esp.cube",
    "data/gaussian_3_scn_esp.cube"]
scan_fdns = [
    "data/gaussian_0_scn_dens.cube", 
    "data/gaussian_1_scn_dens.cube",
    "data/gaussian_2_scn_dens.cube",
    "data/gaussian_3_scn_dens.cube"]

Nfiles = len(scan_fesp)
Nchars = int(np.max([
    len(filename) for filelist in [scan_fesp, scan_fdns] 
    for filename in filelist]))

esplist = np.empty([Nfiles, Nchars], dtype='c')
dnslist = np.empty([Nfiles, Nchars], dtype='c')

for ifle in range(Nfiles):
    esplist[ifle] = "{0:{1}s}".format(scan_fesp[ifle], Nchars)
    dnslist[ifle] = "{0:{1}s}".format(scan_fdns[ifle], Nchars)

# Load cube files, read MDCM global and local files
mdcm.load_cube_files(Nfiles, Nchars, esplist.T, dnslist.T)
mdcm.load_clcl_file(mdcm_clcl)
mdcm.load_cxyz_file(mdcm_cxyz)

# Write MDCM global from local and Fitted ESP cube files
mdcm.write_cxyz_files()
mdcm.write_mdcm_cube_files()

# Get and set local MDCM array (to check if manipulation is possible)
clcl = mdcm.mdcm_clcl
mdcm.set_clcl(clcl)
clcl = mdcm.mdcm_clcl

# Get and set global MDCM array (to check if manipulation is possible)
cxyz = mdcm.mdcm_cxyz
mdcm.set_cxyz(cxyz)
cxyz = mdcm.mdcm_cxyz

# Get RMSE, averaged or weighted over ESP files, or per ESP file each
rmse = mdcm.get_rmse()
print(rmse)
wrmse = mdcm.get_rmse_weighted(Nfiles, [1.]*Nfiles)
print(wrmse)
srmse = mdcm.get_rmse_each(Nfiles)
print(srmse)


# Define simple charge constrained function returning RMSE for given local MDCM
# configuration
def mdcm_rmse(clcl, qtot=-1.0):
    
    # Constraint total charge
    qsum = 0.0
    for i in range(3, len(clcl), 4):
        qsum += clcl[i]
    diff = qtot - qsum
    clcl[i] += diff
    Nchg = float(len(clcl)//4)
    qsum = 0.0
    for i in range(3, len(clcl), 4):
        qsum += clcl[i]
        
    mdcm.set_clcl(clcl)
    rmse = mdcm.get_rmse()
    
    print(qsum, rmse)
    
    return rmse
    
# Apply simple minimization without any feasibility check (!)
# Leads to high amplitudes of MDCM charges and local positions
from scipy.optimize import minimize
res = minimize(mdcm_rmse, clcl, tol=1e-5)
print(res)

# Recompute final RMSE each
srmse = mdcm.get_rmse_each(Nfiles)
print(srmse)

# Not necessary but who knows when it become important to deallocate all 
# global arrays
mdcm.dealloc_all()

