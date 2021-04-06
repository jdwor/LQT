## Example single patient script

library(LQT)

########### Set up config structure ###########
pat_id = "Subject1"
lesion_path = "/Users/JaneGoodall/Study/Images/Subject1/lesion_mask.nii.gz"
parcel_path = system.file("extdata","Schaefer_Yeo_Plus_Subcort",
                          "100Parcels7Networks.nii.gz",package="LQT")
out_path = "/Users/JaneGoodall/Study/Results"

cfg = create_cfg_object(pat_ids=pat_id,
                        lesion_paths=lesion_path,
                        parcel_path=parcel_path,
                        out_path=out_path)

########### Create Damage and Disconnection Measures ###########
# Get parcel damage for patient
get_parcel_damage(cfg)
# Get tract SDC for patient
get_tract_discon(cfg)
# Get parcel SDC and SSPL measures for patient
get_parcel_cons(cfg)

########### Build and View Summary Plots ###########
plot_lqt_subject(cfg, "parcel.damage")
plot_lqt_subject(cfg, "tract.discon")
plot_lqt_subject(cfg, "parcel.discon")
plot_lqt_subject(cfg, "parcel.sspl")

