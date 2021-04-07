## Example mutli-patient script

library(LQT)

########### Set up config structure ###########
pat_ids = paste0("Subject", 1:45)
lesion_paths = list.files('/Users/JaneGoodall/Study/LesionMasks',
                          full.names = TRUE)
parcel_path = system.file("extdata","Schaefer_Yeo_Plus_Subcort",
                          "100Parcels7Networks.nii.gz",package="LQT")
out_path = '/Users/JaneGoodall/Study/Results'

cfg = create_cfg_object(pat_ids=pat_ids,
                        lesion_paths=lesion_paths,
                        parcel_path=parcel_path,
                        out_path=out_path)

########### Create Damage and Disconnection Measures ###########
# Get parcel damage for patients
get_parcel_damage(cfg, cores=2)
# Get tract SDC for patients
get_tract_discon(cfg, cores=2)
# Get parcel SDC and SSPL measures for patients
get_parcel_cons(cfg, cores=2)

########### Build and View Summary Plots ###########
plot_lqt_subject(cfg, "Subject1", "parcel.damage")
plot_lqt_subject(cfg, "Subject1", "tract.discon")
plot_lqt_subject(cfg, "Subject1", "parcel.discon")
plot_lqt_subject(cfg, "Subject1", "parcel.sspl")

########### Compile Datasets for Analysis ###########
data = compile_data(cfg, cores = 2)
list2env(data, .GlobalEnv); rm(data)

