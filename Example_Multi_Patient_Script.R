## Example single patient script

library(LQT)

########### Set up cfg structure ###########
pat_ids = c("Subject1","Subject2","Subject3","Subject4")
lesion_paths = rep('/Users/jordandworkin/Desktop/s01/lesion.nii.gz',4)
out_path = '/Users/jordandworkin/Desktop/s01'

cfg = create_cfg_object(pat_ids=pat_ids,
                        lesion_paths=lesion_paths,
                        out_path=out_path)

########### Create Damage and Disconnection Measures ###########
# Get parcel damage for patient
get_parcel_damage(cfg, cores=2)
# Get tract SDC for patient
get_tract_discon(cfg, cores=2)
# Get parcel SDC and SSPL measures for patient
get_parcel_cons(cfg, cores=2)

########### Build and View Summary Plots ###########
plots = plot_subject_summaries(cfg, "Subject1")
plots$parcel.damage / plots$tract.discon
plots$parcel.discon
plots$parcel.sspl




