# RSIGWAS

This is the script repo for performing RSI GWAS, as documented in the preprint Fan et al., Jun 2021,  DOI:10.1101/2021.06.24.449723

# Codes for performing multivariate combined principal components GWAS
The main code is mCPC.m while it depends on the following function:
nancorr.m
PlinkRead_binary2.m
svdsecon.m

# Codes for performing replication analyses
mCPC_diffusion_rep_ABCD.m includes the PVS calculation on ABCD data, while evoking R script ABCD_ukb_rep.R to do the linear mixed effects analyses to account for complex design of ABCD

# Codes for regional enrichment analyses
mCPC_regional_enrichment.m

# Codes for calculating the multi-tissue compartment metrics (RSI model)
The process involves two steps. First, FOD is calculated given diffusion signals from multi-shell diffusion images:
rsi_fit_MFOD.m

Second, the multiple compartment FOD is then used to generate the corresponding metrics, ND, N0, and NF
rsi_calc_MFmeas.m

  
