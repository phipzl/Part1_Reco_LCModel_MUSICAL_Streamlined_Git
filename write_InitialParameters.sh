#!/bin/bash
#####################################################
#### WRITE INFO TO PARAMETER FILE FOR MATLAB USE ####
#####################################################

# Open Parameter file
Par="${tmp_dir}/InitialParameters.m"
touch $Par			
chmod 755 $Par


# Flags
echo "image_normal_flag = ${image_normal_flag};" > $Par
echo "image_flip_flag = ${image_flip_flag};" >> $Par
echo "image_VC_flag = ${image_VC_flag};" >> $Par
echo "T1w_flag = ${T1w_flag};" >> $Par
echo "WaterReference_flag = ${WaterReference_flag};" >> $Par
echo "mask_flag = ${mask_flag};" >> $Par
echo "hamming_flag = ${hamming_flag};" >> $Par
echo "LipidDecon_flag = ${LipidDecon_flag};" >> $Par
echo "NuisRem_flag = ${NuisRem_flag};" >> $Par
echo "exponential_filter_Hz_flag = ${exponential_filter_Hz_flag};" >> $Par
echo "TwoDCaipParallelImaging_flag = ${TwoDCaipParallelImaging_flag};" >> $Par
echo "SliceParallelImaging_flag = ${SliceParallelImaging_flag};" >> $Par
echo "noisedecorrelation_flag = ${noisedecorrelation_flag};" >> $Par
echo "FirstOrderPhaseCorr_flag = ${FirstOrderPhaseCorr_flag};" >> $Par
echo "FirstOrderPhaseModulation_flag = ${FirstOrderPhaseModulation_flag};" >> $Par
echo "use_phantom_flag = ${use_phantom_flag};" >> $Par
echo "ZeroFillMetMaps_flag = ${ZeroFillMetMaps_flag};" >> $Par
echo "InterpolateCSIResolution_flag = ${InterpolateCSIResolution_flag};" >> $Par
echo "TimeInterpolation_flag = ${TimeInterpolation_flag};" >> $Par
echo "AlignFreq_flag = ${AlignFreq_flag};" >> $Par
echo "dont_compute_LCM_flag = ${dont_compute_LCM_flag};" >> $Par
echo "LCM_ControlPath_flag = ${LCM_ControlPath_flag};" >> $Par
echo "LCM_ControlPath_Water_flag = ${LCM_ControlPath_Water_flag};" >> $Par
echo "phase_encoding_direction_is_RL_flag = ${phase_encoding_direction_is_RL_flag};" >> $Par
echo "basis_echo_flag = ${basis_echo_flag};" >> $Par
echo "control_echo_flag = ${control_echo_flag};" >> $Par
echo "XPACE_motion_correction_flag = ${XPACE_motion_correction_flag};" >> $Par
#BOW - for phase & frequency prior knowledge
echo "priors_flag = ${priors_flag};" >> $Par
echo "NonCartTraj_flag = ${NonCartTraj_flag};" >> $Par
echo "DebugAdditionalInput_flag = ${NonCartTraj_flag};" >> $Par
echo "GradientDelay_flag = ${GradientDelay_flag};" >> $Par


# Mandatory Input files, Output directory
Loop_index=0
for csi_path_dummy in ${csi_path}; do
	((Loop_index = loop_index + 1))
	echo "csi_path{${Loop_index}} = '${csi_path_dummy}';" >> $Par
done

NumberOfCSIFiles=$(echo ${csi_path} | grep -oi  "\.dat\|\.IMA" | wc -l)
echo ${csi_path}
if [[ $NumberOfCSIFiles -eq 0 ]]; then
	NumberOfCSIFiles=$(ls ${csi_path}/* | grep -oi -m 1  "\.dat\|\.IMA" | wc -l)
fi

Loop_index=0
for basis_path_dummy in ${basis_path}; do
	((Loop_index = loop_index + 1))
	echo "basis_path{${Loop_index}} = '${basis_path_dummy}';" >> $Par
done
echo "out_path = '${out_path}';" >> $Par


# Optional Input files
Loop_index=0
for image_normal_path_dummy in ${image_normal_path}; do
	((Loop_index = loop_index + 1))
	echo "image_normal_path{${Loop_index}} = '${image_normal_path_dummy}';" >> $Par
done

Loop_index=0
for image_flip_path_dummy in ${image_flip_path}; do
	((Loop_index = loop_index + 1))
	echo "image_flip_path{${Loop_index}} = '${image_flip_path_dummy}';" >> $Par
done

Loop_index=0
for image_VC_path_dummy in ${image_VC_path}; do
	((Loop_index = loop_index + 1))
	echo "image_VC_path{${Loop_index}} = '${image_VC_path_dummy}';" >> $Par
done

if [[ $noisedecorrelation_flag -eq 1 && -n $noisedecorrelation_path ]]; then				# -n tests for non-emptiness
	echo "noisedecorrelation_path = '${noisedecorrelation_path}';" >> $Par
fi
Loop_index=0
if [[ $AlignFreq_flag -eq 1 && -n $AlignFreq_MethodAndPath ]]; then				# -n tests for non-emptiness
	export AlignFreq_method=$(echo ${AlignFreq_MethodAndPath} | cut -d, -f1)	# The first field is the Method
	export AlignFreq_path=$(echo ${AlignFreq_MethodAndPath} | cut -d, -f2)		# The second field is the path
	if [[ $AlignFreq_path == $AlignFreq_method ]]; then					# If the user only inputs "Alignment" without comma, the cutting doesnt work
		export AlignFreq_path=""
	fi
	echo "AlignFreq_method = '${AlignFreq_method}';" >> $Par
	for AlignFreq_path_dummy in ${AlignFreq_path}; do
		((Loop_index = loop_index + 1))
		echo "AlignFreq_path{${Loop_index}} = '${AlignFreq_path_dummy}';" >> $Par
	done
fi
if [[ $T1w_flag -eq 1 ]]; then
	echo "T1w_path = '${T1w_path}';" >> $Par
fi
if [[ $WaterReference_flag -eq 1 ]]; then
	echo "WaterReference_MethodAndFile = '${WaterReference_MethodAndFile}';" >> $Par
fi
if [[ $LCM_ControlPath_Water_flag -eq 1 ]]; then
	echo "LCM_Control_Water_path= '${LCM_Control_Water_path}';" >> $Par
fi
if [[ $LCM_ControlPath_flag -eq 1 ]]; then
	echo "LCM_ControlPath= '${LCM_ControlPath}';" >> $Par
fi
echo "LCM_Path= '${LCM_Path}';" >> $Par
if [[ ${priors_flag} -eq 1 ]]; then
	echo "priors_path = '${priors_path}';" >> $Par
fi

# Additional User Input
if [[ $hamming_flag -eq 1 ]]; then
	echo "hamming_factor = ${hamming_factor};" >> $Par
fi
if [[ $exponential_filter_Hz_flag -eq 1 ]]; then
	echo "exponential_filter_Hz = ${exponential_filter_Hz};" >> $Par
fi
if [[ $TwoDCaipParallelImaging_flag -eq 1 ]]; then
	if [[ $(echo $InPlaneCaipPattern_And_VD_Radius | grep -c "Skip_Matrix") -eq 1 ]]; then
		InPlaneCaipPattern_And_VD_Radius=$(echo $InPlaneCaipPattern_And_VD_Radius | sed s/Skip_Matrix/InPlaneCaipPattern/)
	fi
	echo ${InPlaneCaipPattern_And_VD_Radius} >> $Par
fi
if [[ $SliceParallelImaging_flag -eq 1 ]]; then
	echo ${SliceAliasingPattern} >> $Par
fi
if [[ $image_normal_flag -eq 1 || $image_VC_flag -eq 1 ]]; then
	echo "phase_encod_dir_is = '${phase_encoding_direction_is}';" >> $Par
fi
if [[ $LipidDecon_flag -eq 1 ]]; then
	echo "LipidDecon_MethodAndNoOfLoops = '${LipidDecon_MethodAndNoOfLoops}';" >> $Par
fi
if [[ $NuisRem_flag -eq 1 ]]; then
	echo "NuisRem_ControlPath = '${NuisRem_ControlPath}';" >> $Par
fi
if [[ $mask_flag -eq 1 ]]; then
	echo "mask_method = '${mask_method}';" >> $Par
fi
if [[ $ZeroFillMetMaps_flag -eq 1 ]]; then
	echo "ZeroFillMetMaps = ${ZeroFillMetMaps};" >> $Par
fi
if [[ $InterpolateCSIResolution_flag -eq 1 ]]; then
	echo "InterpolateCSIResolution = '${InterpolateCSIResolution}';" >> $Par
fi

if [[ $TimeInterpolation_flag -eq 1 ]]; then
	echo "TimeInterpolationFactor = ${TimeInterpolationFactor};" >> $Par
fi
Loop_index=0
if [[ $NonCartTraj_flag -eq 1 ]]; then
	for dummy in ${NonCartTrajFile_path}; do
		((Loop_index = loop_index + 1))
		echo "NonCartTrajFile_path{${Loop_index}} = '${dummy}';" >> $Par
	done
fi
if [[ $DebugAdditionalInput_flag -eq 1 ]]; then
	echo "DebugAdditionalInput = '${DebugAdditionalInput}';" >> $Par
fi
if [[ $GradientDelay_flag -eq 1 ]]; then
	echo "GradientDelay = '${GradientDelay}';" >> $Par
else
	echo "GradientDelay = 0;" >> $Par
fi

if [[ $basis_echo_flag -eq 1 ]]; then
	echo "basis_echo_path = '${basis_echo_path}';" >> $Par
fi
if [[ $control_echo_flag -eq 1 ]]; then
	echo "control_echo_path = '${control_echo_path}';" >> $Par
fi

if [[ $XPACE_motion_correction_flag -eq 1 ]]; then
	echo "XPACE_motion_correction_path = '${XPACE_motion_correction_path}';" >> $Par
fi

#Gh added 2019
if [[ $FirstOrderPhaseModulation_flag -eq 1 ]]; then
	echo "FID_Truncation_in_ms = ${FID_Truncation_in_ms};" >> $Par
fi

if ! [[ ${RunLCModelOn} = "" ]]; then
	echo "RunLCModelOn = '${RunLCModelOn}';" >> $Par
	echo "RunLCModelAs = '${RunLCModelAs}';" >> $Par

fi
if ! [[ ${RunLCModel_CPUCores} = "" ]]; then
	echo "RunLCModel_CPUCores = ${RunLCModel_CPUCores};" >> $Par
fi
