% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_parcel_discon.R
\name{get_parcel_discon}
\alias{get_parcel_discon}
\title{Get parcel-based disconnection}
\usage{
get_parcel_discon(cfg, cores = 1)
}
\arguments{
\item{cfg}{a pre-made cfg structure (as list object).}

\item{cores}{an integer value that indicates how many parallel cores the function should be run on.}
}
\value{
An .RData file with the suffix .disconnectivity.RData. This contains the structural disconnection matrix (disconnectivity);
an .RData file with the suffix .network_measures.RData, which contains various graph measures for the SC matrix;
an .RData file with the suffix _percent_parcel_mats.RData. This file contains a disconnection adjacency matrix (pct_sdc_matrix) and a spared connection adjacency matrix (pct_spared_sc_matrix);
a .txt file with the suffix .connectogram.txt. This file contains a connectogram that can be viewed on http://mkweb.bcgsc.ca/tableviewer/visualize/ by checking the two size options in step 2A (col with row size, row with col size);
a .node file with the suffix _percent_parcel_SDC.node. This file contains the node information for external connectome viewers (e.g. MRIcroGL). Node sizes are proportional to the number of affected connections. Node colors can be pre-assigned in the .cfg file (cfg.node_color), but if not, they will be proportional to the amount of disconnection sustained analogous to node size;
a .edge file with the suffix _percent_parcel_SDC.edge. This file contains the percent SDC matrix in a format that can be loaded into external viewers (e.g. MRICroGL);
a .node file with the suffix _percent_parcel_spared_SC.node. This is analogous to (7), but for the spared SC matrix;
a .edge file with the suffix _percent_parcel_spared_SC.edge. This is analogous to (8), but for the spared SC matrix;
a .tdi.nii.gz file named the same way as the .trk.gz file. This contains a nifti image volume with track density imaging (TDI) values from the .trk.gz file at each voxel. It is essentially a way of converting the .trk.gz file into voxel space. Higher values indicate higher streamline densities at each grid element (voxel);
a .nii file with the suffix _percent_tdi.nii.gz. For each voxel, values correspond the % reduction in streamline density relative to the atlas when accounting for the effects of the lesion.
}
\description{
This function computes parcel-based direct disconnection measures using an
MNI-registered lesion file and an MNI-registered brain parcellation.
}
