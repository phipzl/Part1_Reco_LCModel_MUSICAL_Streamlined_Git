# Concept_HAM_LUKAS \
/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part1_Reco_LCModel_copy/Part1_ProcessMRSI.sh \
-c /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/ComparisonRollerVsSpiral/Concept_Spiral_Vol002/CONCEPT/meas_MID128_bs_RollerMUSI5_Hamming_64x64_FID16766.dat \
-b "/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Basis_Sets/Basis_Set_0.033333ms_WithLacAc_PChConc1_NoNAAGGluMoiety_CorrectedGABA_CorrectLip_DwnfldTests/LCModelOutput_GovindGABA_D55NAA_D100GSH_D190Cr_D130Gln/fid_1.300000ms.basis" \
-o /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/Concept_HAM_LUKAS \
-t /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/ComparisonRollerVsSpiral/Concept_Spiral_Vol002/MP2RAGE/UNI \
-a /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/ComparisonRollerVsSpiral/Concept_Spiral_Vol002/MP2RAGE/INV2 \
-j /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part1_Reco_LCModel/ControlFiles/LCModel_Control_Template.m \
-m "BET" \
-g "" \

/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part2_MRSI_Evaluation/Part2_EvaluateMRSI.sh \
-o /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/Concept_HAM_LUKAS \
