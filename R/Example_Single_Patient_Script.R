## Example single patient script

########### Set up cfg structure ###########
cfg=list()
files=list.files('/Users/jordandworkin/LQT/R')
for(i in files[-c(1:2)]){
  source(paste0('/Users/jordandworkin/LQT/R/',i))
}
library(fslr);library(R.matlab)
library(neurobase);library(igraph)

#### Assign relevant paths ####
# Path to DSI_Studio program
cfg$dsi_path = '/Applications/dsi_studio.app/Contents/MacOS/dsi_studio'
# Path to HCP842 .fib template and tractography atlas
cfg$source_path = '/Users/jordandworkin/LQT/inst/extdata/Tractography_Atlas'
# Path to output directory (patient and atlas result directories will be created within the output directory)
cfg$out_path = '/Users/jordandworkin/Desktop/TestLes'
# Path to lesion (pre-registered to MNI template)
cfg$lesion_path = '/Users/jordandworkin/Desktop/MSles/gold_stand_les.nii.gz'
# Path to parcellation (should have identical dimensions to lesion and be in MNI template space)
cfg$parcel_path = '/Users/jordandworkin/LQT/inst/extdata/Schaefer_Yeo_Plus_Subcort/100Parcels7Networks.nii.gz'

#### Output Filename Options ####
# Patient ID (used as prefix for output files)
cfg$pat_id = 'Subject1'
# File suffix -- used as suffix for output files. Atlas name is recommended (e.g. AAL, Power, Gordon, etc.).
cfg$file_suffix = 'Yeo7100'

#### Connectivity Options ####
# Connectivity type ('end' or 'pass'): if 'end', connections are defined based on streamline endpoints. If 'pass', connections are defined based on streamline pass-throughs. 'End' is recommended (most conservative).
cfg$con_type = 'end'
# Percent spared threshold for computing SSPLs (e.g. 100 means that only fully spared regions will be included in SSPL calculation; 1 means that all regions with at least 1% spared will be included. Default is 50)
cfg$sspl_spared_thresh = 50

#### Connectivity Output Options ####
t=read.csv('/Users/jordandworkin/LQT/inst/extdata/Schaefer_Yeo_Plus_Subcort/100Parcels7Networks.csv',header=T)
cfg$node_label = t$RegionName # n_regions vector of strings corresponding to node labels (i.e. parcel names)
cfg$node_color = t$NetworkID # n_regions vector of integer values corresponding to e.g. network assignments or partitions (used to color nodes in external viewers)

#### Summary Figure Options####
# show summary figures for damage and disconnection measures (single-subject runs only).
cfg$show_summary_figs = 1
# subsequent options are only relevant if cfg.show_summary_figs == 1.
# Parcel coordinates. Used for plotting ball and stick brain graph. If not supplied, they will be estimated from the parcel file
cfg$parcel_coords = cbind(t$X, t$Y, t$Z)
# Percentile threshold for displaying SSPL increases. Will only display SSPL increases above percentile threshold (Default = 90).
cfg$delta_sspl_thresh = 99
# Parcel damage threshold for displaying SSPL increases. Will not display SSPL increases for parcels with % damage greater than or equal to the threshold (default = 100)
cfg$parcel_dmg_thresh = 100
# Tract disconnection threshold for displaying tracts in summary figure. Will not display data for tracts with % disconnection below threshold (Default = 5)
cfg$tract_sdc_thresh = 5
cfg$smooth = 2

########### Create Damage and Disconnection Measures ###########
# Get parcel damage for patient
get_parcel_damage(cfg)
# Get tract SDC for patient
get_tract_discon(cfg)
# Get parcel SDC and SSPL measures for patient
get_parcel_cons(cfg)
