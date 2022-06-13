#!/bin/bash
#SBATCH -p mrce
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --partition mrce
#SBATCH --mem 38192
## Above lines are for slurm

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Concept_NAT_BERNI \
/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part1_Reco_LCModel_copy/Part1_ProcessMRSI.sh \
-c /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/ComparisonRollerVsSpiral/Concept_Spiral_Vol001/CONCEPT/meas_MID35_bs_RollerMUSI5_Natural_64x64_FID16673.dat \
-b "/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Basis_Sets/Basis_Set_0.033333ms_WithLacAc_PChConc1_NoNAAGGluMoiety_CorrectedGABA_CorrectLip_DwnfldTests/LCModelOutput_GovindGABA_D55NAA_D100GSH_D190Cr_D130Gln/fid_1.300000ms.basis" \
-o /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/Concept_NAT_BERNI_WITHMIRROR \
-t /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/ComparisonRollerVsSpiral/Concept_Spiral_Vol001/MP2RAGE/UNI \
-a /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/ComparisonRollerVsSpiral/Concept_Spiral_Vol001/MP2RAGE/INV2 \
-j /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part1_Reco_LCModel_copy/ControlFiles/LCModel_Control_Template.m \
-m "BET" \
-g "" \
-h 100 \

/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part2_MRSI_Evaluation/Part2_EvaluateMRSI.sh \
-o /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/Concept_NAT_BERNI_WITHMIRROR \


# Concept_HAM_BERNI \
/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part1_Reco_LCModel_copy/Part1_ProcessMRSI.sh \
-c /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/ComparisonRollerVsSpiral/Concept_Spiral_Vol001/CONCEPT/meas_MID36_bs_RollerMUSI5_Hamming_64x64_FID16674.dat \
-b "/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Basis_Sets/Basis_Set_0.033333ms_WithLacAc_PChConc1_NoNAAGGluMoiety_CorrectedGABA_CorrectLip_DwnfldTests/LCModelOutput_GovindGABA_D55NAA_D100GSH_D190Cr_D130Gln/fid_1.300000ms.basis" \
-o /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/Concept_HAM_BERNI_WITHMIRROR \
-t /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/ComparisonRollerVsSpiral/Concept_Spiral_Vol001/MP2RAGE/UNI \
-a /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/ComparisonRollerVsSpiral/Concept_Spiral_Vol001/MP2RAGE/INV2 \
-j /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part1_Reco_LCModel_copy/ControlFiles/LCModel_Control_Template.m \
-m "BET" \
-g "" \

/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part2_MRSI_Evaluation/Part2_EvaluateMRSI.sh \
-o /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/Concept_HAM_BERNI_WITHMIRROR \

##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Concept_NAT_LUKAS \
/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part1_Reco_LCModel_copy/Part1_ProcessMRSI.sh \
-c /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/ComparisonRollerVsSpiral/Concept_Spiral_Vol002/CONCEPT/meas_MID127_bs_RollerMUSI5_Natural_64x64_FID16765.dat \
-b "/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Basis_Sets/Basis_Set_0.033333ms_WithLacAc_PChConc1_NoNAAGGluMoiety_CorrectedGABA_CorrectLip_DwnfldTests/LCModelOutput_GovindGABA_D55NAA_D100GSH_D190Cr_D130Gln/fid_1.300000ms.basis" \
-o /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/Concept_NAT_LUKAS_WITHMIRROR \
-t /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/ComparisonRollerVsSpiral/Concept_Spiral_Vol002/MP2RAGE/UNI \
-a /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/ComparisonRollerVsSpiral/Concept_Spiral_Vol002/MP2RAGE/INV2 \
-j /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part1_Reco_LCModel_copy/ControlFiles/LCModel_Control_Template.m \
-m "BET" \
-g "" \
-h 100 \

/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part2_MRSI_Evaluation/Part2_EvaluateMRSI.sh \
-o /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/Concept_NAT_LUKAS_WITHMIRROR \


## Concept_HAM_LUKAS \
#/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part1_Reco_LCModel_copy/Part1_ProcessMRSI.sh \
#-c /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/ComparisonRollerVsSpiral/Concept_Spiral_Vol002/CONCEPT/meas_MID128_bs_RollerMUSI5_Hamming_64x64_FID16766.dat \
#-b "/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Basis_Sets/Basis_Set_0.033333ms_WithLacAc_PChConc1_NoNAAGGluMoiety_CorrectedGABA_CorrectLip_DwnfldTests/LCModelOutput_GovindGABA_D55NAA_D100GSH_D190Cr_D130Gln/fid_1.300000ms.basis" \
#-o /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/Concept_HAM_LUKAS_WITHMIRROR \
#-t /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/ComparisonRollerVsSpiral/Concept_Spiral_Vol002/MP2RAGE/UNI \
#-a /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/ComparisonRollerVsSpiral/Concept_Spiral_Vol002/MP2RAGE/INV2 \
#-j /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part1_Reco_LCModel_copy/ControlFiles/LCModel_Control_Template.m \
#-m "BET" \
#-g "" \

#/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part2_MRSI_Evaluation/Part2_EvaluateMRSI.sh \
#-o /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/Concept_HAM_LUKAS_WITHMIRROR \

##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Concept_NAT_SCHOKO \
/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part1_Reco_LCModel_copy/Part1_ProcessMRSI.sh \
-c /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/ComparisonRollerVsSpiral/Concept_Spiral_Vol003/CONCEPT/meas_MID226_bs_RollerMUSI5_Natural_64x64_FID16864.dat \
-b "/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Basis_Sets/Basis_Set_0.033333ms_WithLacAc_PChConc1_NoNAAGGluMoiety_CorrectedGABA_CorrectLip_DwnfldTests/LCModelOutput_GovindGABA_D55NAA_D100GSH_D190Cr_D130Gln/fid_1.300000ms.basis" \
-o /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/Concept_NAT_SCHOKO_WITHMIRROR \
-t /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/ComparisonRollerVsSpiral/Concept_Spiral_Vol003/MP2RAGE/UNI \
-a /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/ComparisonRollerVsSpiral/Concept_Spiral_Vol003/MP2RAGE/INV2 \
-j /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part1_Reco_LCModel_copy/ControlFiles/LCModel_Control_Template.m \
-m "BET" \
-g "" \
-h 100 \

/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part2_MRSI_Evaluation/Part2_EvaluateMRSI.sh \
-o /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/Concept_NAT_SCHOKO_WITHMIRROR \


# Concept_HAM_SCHOKO \
/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part1_Reco_LCModel_copy/Part1_ProcessMRSI.sh \
-c /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/ComparisonRollerVsSpiral/Concept_Spiral_Vol003/CONCEPT/meas_MID227_bs_RollerMUSI5_Hamming_64x64_FID16865.dat \
-b "/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Basis_Sets/Basis_Set_0.033333ms_WithLacAc_PChConc1_NoNAAGGluMoiety_CorrectedGABA_CorrectLip_DwnfldTests/LCModelOutput_GovindGABA_D55NAA_D100GSH_D190Cr_D130Gln/fid_1.300000ms.basis" \
-o /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/Concept_HAM_SCHOKO_WITHMIRROR \
-t /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/ComparisonRollerVsSpiral/Concept_Spiral_Vol003/MP2RAGE/UNI \
-a /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/ComparisonRollerVsSpiral/Concept_Spiral_Vol003/MP2RAGE/INV2 \
-j /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part1_Reco_LCModel_copy/ControlFiles/LCModel_Control_Template.m \
-m "BET" \
-g "" \

/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part2_MRSI_Evaluation/Part2_EvaluateMRSI.sh \
-o /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/Concept_HAM_SCHOKO_WITHMIRROR \

##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Concept_NAT_WOLF-DIETRICH \
/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part1_Reco_LCModel_copy/Part1_ProcessMRSI.sh \
-c /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/ComparisonRollerVsSpiral/Concept_Spiral_Vol004/CONCEPT/meas_MID254_bs_RollerMUSI5_Natural_64x64_FID16892.dat \
-b "/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Basis_Sets/Basis_Set_0.033333ms_WithLacAc_PChConc1_NoNAAGGluMoiety_CorrectedGABA_CorrectLip_DwnfldTests/LCModelOutput_GovindGABA_D55NAA_D100GSH_D190Cr_D130Gln/fid_1.300000ms.basis" \
-o /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/Concept_NAT_WOLFDIETRICH_WITHMIRROR \
-t /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/ComparisonRollerVsSpiral/Concept_Spiral_Vol004/MP2RAGE/UNI \
-a /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/ComparisonRollerVsSpiral/Concept_Spiral_Vol004/MP2RAGE/INV2 \
-j /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part1_Reco_LCModel_copy/ControlFiles/LCModel_Control_Template.m \
-m "BET" \
-g "" \
-h 100 \

/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part2_MRSI_Evaluation/Part2_EvaluateMRSI.sh \
-o /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/Concept_NAT_WOLFDIETRICH_WITHMIRROR \


# Concept_HAM_WOLF-DIETRICH \
/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part1_Reco_LCModel_copy/Part1_ProcessMRSI.sh \
-c /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/ComparisonRollerVsSpiral/Concept_Spiral_Vol004/CONCEPT/meas_MID255_bs_RollerMUSI5_Hamming_64x64_FID16893.dat \
-b "/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Basis_Sets/Basis_Set_0.033333ms_WithLacAc_PChConc1_NoNAAGGluMoiety_CorrectedGABA_CorrectLip_DwnfldTests/LCModelOutput_GovindGABA_D55NAA_D100GSH_D190Cr_D130Gln/fid_1.300000ms.basis" \
-o /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/Concept_HAM_WOLFDIETRICH_WITHMIRROR \
-t /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/ComparisonRollerVsSpiral/Concept_Spiral_Vol004/MP2RAGE/UNI \
-a /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/ComparisonRollerVsSpiral/Concept_Spiral_Vol004/MP2RAGE/INV2 \
-j /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part1_Reco_LCModel_copy/ControlFiles/LCModel_Control_Template.m \
-m "BET" \
-g "" \

/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part2_MRSI_Evaluation/Part2_EvaluateMRSI.sh \
-o /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/Concept_HAM_WOLFDIETRICH_WITHMIRROR \


