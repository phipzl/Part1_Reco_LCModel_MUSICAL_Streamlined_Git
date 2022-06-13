#!/bin/bash
#SBATCH -p mrce
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --partition mrce
#SBATCH --mem 38192
## Above lines are for slurm

## MS_PATIENT_004_SSL (A single slice undersampled dataset) \
#/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/bstrasser/Sourcecode/MRSI_Processing/Part1_Reco_LCModel/ActualCode/Part1_ProcessMRSI.sh \
#-c /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/Patients/MS_PATIENT_004/DAT/meas_MID758_bs_gh_ghsi_30b_R6_CAIPI_SSL_FID10087.dat \
#-b "/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Basis_Sets/Basis_Set_0.033333ms_WithLacAc_PChConc1_NoNAAGGluMoiety_CorrectedGABA_CorrectLip/LCModelOutput_GovindGABA/fid_1.300000ms.basis" \
#-o /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Process_Results/MS_Patients/MS_Patient_004_4SL_2014-10-03 \
#-g "" \
#-j /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Sourcecode/MRSI_Processing_ReleaseVersions/Part1_Reco_LCModel_v0.9/ControlFiles/LCModel_Control_Volunteers.m \
#-m "threshold" \
#-h 100 \
#-I "[0.5 0.5 1]" \





## CSI_PROBAND_22_Downfield \
#./Part1_ProcessMRSI.sh \
#-c /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/Multichannel_combination/CSI_PROBAND_22/DAT/meas_85908.dat \
#-b "/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Basis_Sets/Basis_Set_0.033333ms_WithLacAc_PChConc1_NoNAAGGluMoiety_CorrectedGABA_CorrectLip_DwnfldTests/LCModelOutput_GovindGABA_D55NAA_D100GSH_D190Cr_D130Gln/fid_1.300000ms.basis" \
#-o /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/bstrasser/Process_Results/DownfieldTests/LowRes/CSI_PROBAND_22_DownfieldD55NAA_D100GSH_D190Cr_D110Gln_Res8x8 \
#-i "/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/Multichannel_combination/CSI_PROBAND_22/DAT/meas_85903.dat" \
#-t /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/Multichannel_combination/CSI_PROBAND_22/IMA/T1 \
#-g "" \
#-j /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Sourcecode/MRSI_Processing_ReleaseVersions/Part1_Reco_LCModel_v0.9.5/ControlFiles/LCModel_Control_Volunteers_Downfield.m \
#-m "BET" \
#-h 100 \
#-I "[8 8 1]" \

#/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Sourcecode/MRSI_Processing_ReleaseVersions/Part2_Registration_Evaluation_v1.3.2/evaluation_split_real_3D.sh \
#-o /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/bstrasser/Process_Results/DownfieldTests/LowRes/CSI_PROBAND_22_DownfieldD55NAA_D100GSH_D190Cr_D110Gln_Res1x1 \


# Concept_philipp_64 \
/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part1_Reco_LCModel_copy/Part1_ProcessMRSI.sh \
-c /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/Rollercoaster/CONCEPT_TESTS/ConceptTest_Vol04/meas_MID490_bs_RollerMUSI5_Natural_VC_FID16422.dat \
-b "/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Basis_Sets/Basis_Set_0.033333ms_WithLacAc_PChConc1_NoNAAGGluMoiety_CorrectedGABA_CorrectLip_DwnfldTests/LCModelOutput_GovindGABA_D55NAA_D100GSH_D190Cr_D130Gln/fid_1.300000ms.basis" \
-o /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/lukas_VC_memory \
-t /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/Rollercoaster/CONCEPT_TESTS/ConceptTest_Vol04/MP2RAGE/UNI \
-a /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/Rollercoaster/CONCEPT_TESTS/ConceptTest_Vol04/MP2RAGE/INV2 \
-j /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part1_Reco_LCModel/ControlFiles/LCModel_Control_Template.m \
-m "BET" \
-g "" \
-h 100 \

/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part2_MRSI_Evaluation/Part2_EvaluateMRSI.sh \
-o /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/lukas_VC_memory \



## Blue Phantom PHaseRollTest \
#/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part1_Reco_LCModel/Part1_ProcessMRSI.sh \
#-c /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/Rollercoaster/CONCEPT_TESTS/test58_BluePhantomForPhaseRollTests/meas_MID568_bs_RollerMUSI5_Natural_64x64_FID17202.dat \
#-b "/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Basis_Sets/Basis_Set_0.033333ms_WithLacAc_PChConc1_NoNAAGGluMoiety_CorrectedGABA_CorrectLip_DwnfldTests/LCModelOutput_GovindGABA_D55NAA_D100GSH_D190Cr_D130Gln/fid_1.300000ms.basis" \
#-o /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/lukas_VC_memory \
#-t /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/Rollercoaster/CONCEPT_TESTS/ConceptTest_Vol04/MP2RAGE/UNI \
#-a /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/Rollercoaster/CONCEPT_TESTS/ConceptTest_Vol04/MP2RAGE/INV2 \
#-j /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part1_Reco_LCModel/ControlFiles/LCModel_Control_Template.m \
#-m "BET" \
#-g "" \

#/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part2_MRSI_Evaluation/Part2_EvaluateMRSI.sh \
#-o /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/lukas_VC_memory \

## Concept_philipp_64 \
#/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part1_Reco_LCModel/Part1_ProcessMRSI.sh \
#-c /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/Rollercoaster/CONCEPT_TESTS/Concept_Vol04/meas_MID574_bs_Rlold_Hamming_64x64_TR550_FID12033.dat \
#-b "/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Basis_Sets/Basis_Set_0.033333ms_WithLacAc_PChConc1_NoNAAGGluMoiety_CorrectedGABA_CorrectLip_DwnfldTests/LCModelOutput_GovindGABA_D55NAA_D100GSH_D190Cr_D130Gln/fid_1.300000ms.basis" \
#-o /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/lukas_64_ham_999_NAAtrue_FIX \
#-t /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/Rollercoaster/CONCEPT_TESTS/Concept_Vol04/MP2RAGE/UNI \
#-a /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/Rollercoaster/CONCEPT_TESTS/Concept_Vol04/MP2RAGE/INV2 \
#-j /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part1_Reco_LCModel/ControlFiles/LCModel_Control_Template.m \
#-m "BET" \

#/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part2_MRSI_Evaluation/Part2_EvaluateMRSI.sh \
#-o /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/lukas_64_ham_999_NAAtrue_FIX \


## Volunteer_42_FullRef \
#./Part1_ProcessMRSI.sh \
#-c /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/Parallel_Imaging/CSI_PROBAND_42/DAT/meas_57150.dat \
#-b "/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Basis_Sets/Basis_Set_0.033333ms_WithLacAc_PChConc1_NoNAAGGluMoiety_CorrectedGABA_CorrectLip/LCModelOutput_GovindGABA/fid_1.300000ms.basis /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Basis_Sets/Basis_Set_0.033333ms_WithLacAc_PChConc1_NoNAAGGluMoiety_CorrectedGABA_CorrectLip/LCModelOutput_GovindGABA/fid_2.300000ms.basis" \
#-o /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/bstrasser/Process_Results/MRSI_ParallelImaging/Step1_Simulation/2+1DCaipi/gFactorComp_new/Volunteer_42_FullRef_2015-02-26_DELETEME \
#-t /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/Parallel_Imaging/CSI_PROBAND_42/T1 \
#-a /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/Parallel_Imaging/CSI_PROBAND_42/INV2 \
#-g "" \
#-j /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Sourcecode/MRSI_Processing_ReleaseVersions/Part1_Reco_LCModel_v0.9/ControlFiles/LCModel_Control_Volunteers.m \
#-m "BET" \
#-h 100 \
#-r "InPlaneCaipPattern = [1; 0]; VD_Radius = 0;" \
#-R "FoVShifts_x = [0 0.3]; FoVShifts_y = [0 0.3];" \
#-I "[0.5 0.5 1]" \

#/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Sourcecode/MRSI_Processing_ReleaseVersions/Part2_Registration_Evaluation_v1.3.2/evaluation_split_real_3D.sh \
#-o /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/bstrasser/Process_Results/MRSI_ParallelImaging/Step1_Simulation/2+1DCaipi/gFactorComp_new/Volunteer_42_FullRef_2015-02-26 \
#-k "['fullrange']" \




## 3D_spiral_dataset \
#./Part1_ProcessMRSI.sh \
#-c /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/udydak/PD-CSI-anonymous/AD40/Scan1/CSI_corr/RETRO_RECON.MR.PHYSIKER_HNILICOVA.0004.0003.2015.05.19.08.57.25.906250.48801725.IMA \
#-b "/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Basis_Sets/BASIS_MEGA_68ms_newSchema/basis_MEGA_LASER_diff_TE=68_131101.BASIS" \
#-o /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/bstrasser/Process_Results/TheOthers/zf_tests/PD-CSI_zf8 \
#-t /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/udydak/PD-CSI-anonymous/AD40/Scan1/T1long \
#-j /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Sourcecode/MRSI_Processing_ReleaseVersions/Part1_Reco_LCModel_v1.9.x18/ControlFiles/LCModel_Control_Volunteers_4_2_to_1_4.m \
#-m "VOI" \

