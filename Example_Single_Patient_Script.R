## Example single patient script
library(LQT)

########### Set up cfg structure ###########
pat_id = "Subject1"
out_path = '/Users/jordandworkin/Desktop/TestLes'
lesion_path = '/Users/jordandworkin/LQT/inst/extdata/Example_Lesions/ExampleLesion1.nii.gz'

cfg = create_cfg_object(pat_id=pat_id,out_path=out_path,
                        lesion_path=lesion_path)

########### Create Damage and Disconnection Measures ###########
# Get parcel damage for patient
get_parcel_damage(cfg)
# Get tract SDC for patient
get_tract_discon(cfg)
# Get parcel SDC and SSPL measures for patient
get_parcel_cons(cfg)
