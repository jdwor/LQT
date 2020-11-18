## Example single patient script

########### Set up cfg structure ###########
cfg=list()
files=list.files('/Users/jordandworkin/LQT/R')
for(i in files[-c(1:2)]){
  source(paste0('/Users/jordandworkin/LQT/R/',i))
}
library(fslr);library(R.matlab)
library(neurobase);library(igraph)

#### Create cfg object ####
pat_id = "Subject1"
out_path = '~/Desktop/TestLes'
lesion_path = '~/LQT/inst/extdata/Example_Lesions/ExampleLesion1.nii.gz'

cfg = create_cfg_object(pat_id=pat_id,out_path=out_path,
                        lesion_path=lesion_path)

########### Create Damage and Disconnection Measures ###########
# Get parcel damage for patient
get_parcel_damage(cfg)
# Get tract SDC for patient
get_tract_discon(cfg)
# Get parcel SDC and SSPL measures for patient
get_parcel_cons(cfg)
