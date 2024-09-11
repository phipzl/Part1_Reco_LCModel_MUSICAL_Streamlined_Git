#

FoundMinc=$(basename ${AlignFreq_path} | grep -ci ".mnc")
FoundMagRemember=0
if [[ $FoundMinc -eq 0 ]]; then
	# Convert from dicom to minc
	for AlignFreq_dummy in ${AlignFreq_path}; do
		FoundMag=$(basename $AlignFreq_dummy | grep -ci "mag")
		FoundPha=$(basename $AlignFreq_dummy | grep -ci "pha")
		if [[ $FoundMag -eq 1 ]]; then
			FoundMagRemember=1
			dcm2mnc -clobber $AlignFreq_dummy -dname ${tmp_dir} -fname AlignFreq_mag_dum .
		else
			dcm2mnc $AlignFreq_dummy -dname ${tmp_dir} -fname AlignFreq_pha .
		fi
	done

	# Need that mask later
	mincresample -clobber -nearest_neighbour -like ${tmp_dir}/AlignFreq_template.mnc ${tmp_dir}/mask_brain_hires.mnc ${tmp_dir}/mask_brain_B0res.mnc	
	# Run the prepared script for assigning the DeltaTE, and the LarmorFrequency
	 . ${tmp_dir}/B0Map_Dummy.sh
	# PHASE UNWRAPPING
	# Only if magnitude image is provided --> Do phase unwrapping
	if [[ ${FoundMagRemember} -eq 1 ]]; then
		# Make the hi-res brain-mask smaller, to be sure to only include brain voxels
		mincmorph -clobber -successive EE ${tmp_dir}/mask_brain_hires.mnc ${tmp_dir}/mask_brain_hires_dum.mnc
		#cp ${tmp_dir}/mask_brain_hires.mnc ${tmp_dir}/mask_brain_hires_dum.mnc

		# Mask the magnitude map
		mincresample -clobber -nearest_neighbour -like ${tmp_dir}/AlignFreq_mag_dum.mnc ${tmp_dir}/mask_brain_hires_dum.mnc ${tmp_dir}/mask_brain_hires.mnc
		mincmath -mult ${tmp_dir}/AlignFreq_mag_dum.mnc ${tmp_dir}/mask_brain_hires.mnc ${tmp_dir}/AlignFreq_mag.mnc

		# Get phase map from range [-4096, 4095] to [0 4096]
		mincmath ${tmp_dir}/AlignFreq_pha.mnc -add -const 4096 ${tmp_dir}/AlignFreq_pha_dum.mnc
		mincmath -clobber ${tmp_dir}/AlignFreq_pha_dum.mnc -div -const 2 ${tmp_dir}/AlignFreq_pha.mnc

		# Convert to Nifti
		mnc2nii ${tmp_dir}/AlignFreq_pha.mnc ${tmp_dir}/AlignFreq_pha_nii.nii
		mnc2nii ${tmp_dir}/AlignFreq_mag.mnc ${tmp_dir}/AlignFreq_mag_nii.nii

		# Phase-unwrap data with fsl_prepare_fieldmap (output in rad/s)
		fsl_prepare_fieldmap SIEMENS ${tmp_dir}/AlignFreq_pha_nii.nii ${tmp_dir}/AlignFreq_mag_nii.nii ${tmp_dir}/AlignFreq_B0FieldMap_nii.nii $AlignFreq_DeltaTE_ms
		gunzip ${tmp_dir}/AlignFreq_B0FieldMap_nii.nii.gz 
		nii2mnc ${tmp_dir}/AlignFreq_B0FieldMap_nii.nii ${tmp_dir}/AlignFreq_B0FieldMap.mnc
		cp ${tmp_dir}/AlignFreq_B0FieldMap.mnc ${tmp_dir}/AlignFreq_pha.mnc
	fi

	# Resample fieldmap to csi-data & rescale field map
	# Rescale factor is 10^6/(2*pi*B0Map_Par.LarmorFreq) in case of phase unwrapping, and 10^12/(8192*dTE_us*B0Map_Par.LarmorFreq) otherwise
	# Just to be sure, multiply with brain mask
	mincresample -tricubic -clobber -like ${tmp_dir}/AlignFreq_template.mnc ${tmp_dir}/AlignFreq_pha.mnc ${tmp_dir}/AlignFreq_B0FieldMap_dum.mnc 
	mincmath -clobber -mult -const -${AlignFreq_RescaleFactor} ${tmp_dir}/AlignFreq_B0FieldMap_dum.mnc ${tmp_dir}/AlignFreq_B0FieldMap_dum2.mnc
	mincmath -nocheck_dimensions -clobber -mult ${tmp_dir}/AlignFreq_B0FieldMap_dum2.mnc ${tmp_dir}/mask_brain_B0res.mnc ${tmp_dir}/AlignFreq_B0FieldMap.mnc
else
	cp ${AlignFreq_path} ${tmp_dir}/AlignFreq_B0FieldMap.mnc
fi


# Convert to raw and copy to output folder
mkdir -p ${out_path}/AlignFreq
cp ${tmp_dir}/AlignFreq_B0FieldMap.mnc ${out_path}/AlignFreq/AlignFreq_B0FieldMap.mnc
minctoraw ./${tmp_dir}/AlignFreq_B0FieldMap.mnc -nonormalize -float > ./${tmp_dir}/AlignFreq_B0FieldMap.raw 
# Remove temporary files
rm -f ${tmp_dir}/*_dum* ${tmp_dir}/*.gz ${tmp_dir}/*.nii ${tmp_dir}/AlignFreq_mag.mnc ${tmp_dir}/AlignFreq_pha.mnc
rm -f ${tmp_dir}/mask_brain_hires.mnc ${tmp_dir}/mask_brain_B0res.mnc


