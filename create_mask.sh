# 2.
############ USE MASK OR CREATE MASK OUT OF (IN PRIORITY ORDER): MASK, T1_MAP, VOLUME COIL, SUM OVER CHANNELS OF IMAGING-MAGNITUDES, SUM OVER CHANNELS OF CSI-MAGNITUDES ############


# SCHEME OF FOLLOWING IF'S (OUTDATED):


#if [[ $mask_flag -eq 1 ]]; then									
#	# DO MASK STUFF
#else
#	if [[ $T1w_flag -eq 1 ]]; then								
#		# DO T1-STUFF
#	else
#		if [[ $image_VC_flag -eq 1 ]]; then	
#			# DO VC STUFF
#		elif [[ $image_flag -eq 1 ]]; then						
#			# DO AC STUFF
#		fi
#		# DO AC+VC STUFF
#	fi
#	# DO AC+VC+T1-STUFF:	
#	if [[ $use_phantom_flag -eq 1 ]]; then
#		PHANTOM
#	else
#		BRAIN
#	fi
#	if [[ $T1w_flag -eq 0 ]]; then											
#		DO AC+VC STUFF
#	fi
#	DO AC+VC+T1-STUFF
#fi
# DO VC+AC+T1+MASK - STUFF


#read -p "Stop before creating mask1."

if [[ $T1w_flag -eq 1 ]]; then				   											# if T1_map is inputted, create magnitude minc file									

	if [[ "$T1w_path" == */*.mnc ]]; then												# Copy magnitude.mnc if T1w_path is minc-file.
		cp $T1w_path ${tmp_dir}/magnitude.mnc
	else																				# Only perform the dcm2mnc stuff if T1w_path is NOT a minc-file
		mkdir ${tmp_dir}/mnc
		dcm2mnc $T1w_path ${tmp_dir}/mnc 
		mv ${tmp_dir}/mnc/*/* ${tmp_dir}/magnitude.mnc
		rm -r ${tmp_dir}/mnc
#		echo "File magnitude.mnc was created in $tmp_dir."

		if [[ $T1w_AntiNoise_flag -eq 1 ]]; then										# If there was inputted another 3D-measurement for removing the noise of the T1w-image around the head.
			mkdir ${tmp_dir}/mnc
			dcm2mnc $T1w_AntiNoise_path ${tmp_dir}/mnc 
			mv ${tmp_dir}/mnc/*/* ${tmp_dir}/magnitude_AntiNoise.mnc
			rm -r ${tmp_dir}/mnc
#			echo "File magnitude_AntiNoise.mnc was created in $tmp_dir."

			max_magnitude=`mincstats -quiet -max ${tmp_dir}/magnitude_AntiNoise.mnc`
			lower_threshold=$(echo "scale=6 ; ${max_magnitude}/25" | bc)
			mincmath -clobber -segment -const2 $lower_threshold $max_magnitude ${tmp_dir}/magnitude_AntiNoise.mnc ${tmp_dir}/mask_AntiNoise.mnc
			mincmath -clobber -mult ${tmp_dir}/magnitude.mnc ${tmp_dir}/mask_AntiNoise.mnc ${tmp_dir}/magnitude2.mnc

			##GH: added for new masking
			mnc2nii -quiet ${tmp_dir}/magnitude_AntiNoise.mnc ${tmp_dir}/magnitude_inversion_2.nii

			rm ${tmp_dir}/magnitude.mnc
			mv ${tmp_dir}/magnitude2.mnc ${tmp_dir}/magnitude.mnc
			rm ${tmp_dir}/*AntiNoise*
		fi
	fi
	
else 																  	# In this case create mask out of GRE image or even CSI data (done in GetPar_CreateTempl_MaskPart1)
	rawtominc -float -like ${tmp_dir}/mag_template.mnc -input ${tmp_dir}/magnitude.raw ${tmp_dir}/magnitude.mnc	# Done there: READ IN DAT OF ALL CHANNELS, SoS of CHANNELS, WRITE DAT AS .RAW
fi
cp ${tmp_dir}/magnitude.mnc ${out_path}/maps/magnitude.mnc


##### Create FLAIR minc
if [[ $FLAIR_flag -eq 1 ]]; then
echo -e "\n\n Do it with FLAIR! \n\n"	
mkdir ${tmp_dir}/mnc
dcm2mnc $FLAIR_path ${tmp_dir}/mnc 
mv ${tmp_dir}/mnc/*/* ${tmp_dir}/flair.mnc
rm -r ${tmp_dir}/mnc
echo "File flair.mnc was created in $tmp_dir."
cp ${tmp_dir}/flair.mnc ${out_path}/maps/flair.mnc
#else
#echo -e "\n\n You should not see this message. I am error. \n\n"
fi

mkdir -p ${out_path}/maps/Extra/


##### Create B1 minc
if [[ $B1corr_flag -eq 1 ]]; then
./B1_preparations.sh
fi

#read -p "Stop before creating mask2."

if [[ $mask_flag -eq 1 ]]; then

	# Create mask out of magnitude image, VoI-Info, threshold, or copy user-given mask.
	voi_found=$(echo $mask_method | grep -c -i "voi")
	bet_found=$(echo $mask_method | grep -c -i "bet")
	ThreeD_found=$(echo $mask_method | grep -c -i "dreid")
	thresh_found=$(echo $mask_method | grep -c -i "thresh")
	
	 if [[ $ThreeD_found > 0 ]]; then



		#VOI pfusch
		rawtominc -float -like ${tmp_dir}/csi_template.mnc -input ${tmp_dir}/mask_brain_VOI.raw ${tmp_dir}/mask_brain_VOI.mnc >/dev/null    # The VOI mask was created by the MATLAB script

		if [[ $InterpolateCSIResolution_flag -eq 1 ]]; then
		    rawtominc -float -like ${tmp_dir}/csi_template_BefInterpol.mnc -input ${tmp_dir}/mask_brain_VOI_BefInterpol.raw ${tmp_dir}/mask_brain_BefInterpol_VOI.mnc >/dev/null
		fi
		#end VOI pfusch


		#BET pfusch
		mnc2nii -quiet ${tmp_dir}/magnitude.mnc ${tmp_dir}/nii_magnitude.nii >/dev/null					# If brain was inputted, use BET2
		

		
		#if [[ $T1w_flag -eq 1 ]]; then
		#	${betp} ${tmp_dir}/nii_magnitude ${tmp_dir}/brain -f 0.33 -g 0 -n -m     #-f 0.5 -g 0 -n -m
		#else
		#	${betp} ${tmp_dir}/nii_magnitude ${tmp_dir}/brain -f 0.7 -g 0 -n -m      #-f 0.5 -g 0 -n -m ##################### -Z # Improve results if FoV is very small in z-direction
		#fi

		if [[ $T1w_AntiNoise_flag -eq 1 ]]; then
			${betp} ${tmp_dir}/magnitude_inversion_2.nii ${tmp_dir}/brain -f 0.1 -g 0.0 -n -m  #-B 
			${betp} ${tmp_dir}/magnitude_inversion_2.nii ${tmp_dir}/lipid $BetOptions -n -A
		elif [[ $T1w_flag -eq 1 ]]; then
			${betp} ${tmp_dir}/nii_magnitude ${tmp_dir}/brain $BetOptions -B -f 0.1 -g 0.0 -n -m #-n -m     #-f 0.5 -g 0 -n -m
			${betp} ${tmp_dir}/nii_magnitude ${tmp_dir}/lipid $BetOptions -n -A
		else
			${betp} ${tmp_dir}/nii_magnitude ${tmp_dir}/brain $BetOptions -n -m      #-f 0.5 -g 0 -n -m ##################### -Z # Improve results if FoV is very small in z-direction
		fi
		
		gunzip ${tmp_dir}/brain_mask.nii.gz		
		gunzip ${tmp_dir}/lipid_skull_mask.nii.gz

		gunzip ${tmp_dir}/lipid_outskin_mask.nii.gz
		gunzip ${tmp_dir}/lipid_outskull_mask.nii.gz
		gunzip ${tmp_dir}/lipid_inskull_mask.nii.gz
		
		nii2mnc -quiet ${tmp_dir}/lipid_outskin_mask.nii ${tmp_dir}/OUTSKIN.mnc >/dev/null
		nii2mnc -quiet ${tmp_dir}/lipid_outskull_mask.nii ${tmp_dir}/OUTSKULL.mnc >/dev/null
		nii2mnc -quiet ${tmp_dir}/lipid_inskull_mask.nii ${tmp_dir}/INSKULL.mnc >/dev/null



		rm ${tmp_dir}/mask_brain_unres.mnc
		rm ${tmp_dir}/mask_lipid_unres.mnc

		#generate brain mask in .mnc format
		nii2mnc -quiet ${tmp_dir}/brain_mask.nii ${tmp_dir}/mask_brain_unres.mnc >/dev/null
				
		#generate lipid mask by subtracting the brain mask from the OUTSKIN mask.
		#for definition of OUTSKIN, OUTSKULL and INSKULL see BETsurf (http://poc.vl-e.nl/distribution/manual/fsl-3.2/bet2/)
		mincmath -sub ${tmp_dir}/OUTSKIN.mnc ${tmp_dir}/mask_brain_unres.mnc ${tmp_dir}/mask_lipid_unres.mnc >/dev/null

		rm ${tmp_dir}/*.nii








		#### Solving some problem with the dircos (???)
		if [[ $T1w_flag -eq 0 ]]; then
			minctoraw ${tmp_dir}/mask_brain_unres.mnc -nonormalize -float > ${tmp_dir}/mask_brain_unres.raw
			rawtominc -float -clobber -like ${tmp_dir}/mag_template.mnc -input ${tmp_dir}/mask_brain_unres.raw ${tmp_dir}/mask_brain_unres.mnc >/dev/null
		
			minctoraw ${tmp_dir}/mask_lipid_unres.mnc -nonormalize -float > ${tmp_dir}/mask_lipid_unres.raw
			rawtominc -float -clobber -like ${tmp_dir}/mag_template.mnc -input ${tmp_dir}/mask_lipid_unres.raw ${tmp_dir}/mask_lipid_unres.mnc >/dev/null

		fi

		## Resample to CSI
		mincresample -clobber -nearest_neighbour -float -like ${tmp_dir}/csi_template.mnc ${tmp_dir}/mask_brain_unres.mnc ${tmp_dir}/mask_brain_BET.mnc >/dev/null
		mincresample -clobber -nearest_neighbour -like ${tmp_dir}/csi_template.mnc ${tmp_dir}/mask_lipid_unres.mnc ${tmp_dir}/mask_lipid.mnc >/dev/null
		
		if [[ $InterpolateCSIResolution_flag -eq 1 ]]; then
			mincresample -clobber -nearest_neighbour -like ${tmp_dir}/csi_template_BefInterpol.mnc ${tmp_dir}/mask_brain_unres.mnc ${tmp_dir}/mask_brain_BefInterpol_BET.mnc >/dev/null
			mincresample -clobber -nearest_neighbour -like ${tmp_dir}/csi_template_BefInterpol.mnc ${tmp_dir}/mask_lipid_unres.mnc ${tmp_dir}/mask_lipid_BefInterpol.mnc >/dev/null
		fi
		
		#read -p "Stop before creating mask2."
		minctoraw ${tmp_dir}/mask_lipid.mnc -nonormalize -float > ${tmp_dir}/mask_lipid.raw
		#read -p "Stop before creating mask22."
		
		# create common mask
		#mincmath -clobber -mult -nocheck_dimensions ${tmp_dir}/mask_brain_BET.mnc ${tmp_dir}/mask_brain_BET.mnc ${tmp_dir}/mask_brain.mnc

		mincmath -clobber -mult -nocheck_dimensions ${tmp_dir}/mask_brain_BET.mnc ${tmp_dir}/mask_brain_VOI.mnc ${tmp_dir}/mask_brain.mnc >/dev/null
		mincmath -clobber -mult -nocheck_dimensions ${tmp_dir}/mask_brain_BefInterpol_BET.mnc ${tmp_dir}/mask_brain_BefInterpol_VOI.mnc ${tmp_dir}/mask_brain_BefInterpol.mnc >/dev/null

		# For some reason the mask is flipped (because of the nii-stuff?). Undo this flip
		if [[ $T1w_flag -eq 0 ]]; then
			minctoraw -clobber ${tmp_dir}/mask_brain.mnc -nonormalize -float > ${tmp_dir}/mask_brain.raw
			minctoraw ${tmp_dir}/mask_lipid.mnc -nonormalize -float > ${tmp_dir}/mask_lipid.raw
			${matlabp} -r "tmp_dir = '${tmp_dir}'; ${MatlabStartupCommand}" -nodisplay -nojvm < ./flip_mask.m
			rawtominc -float -clobber -like ${tmp_dir}/csi_template.mnc -input ${tmp_dir}/mask_brain.raw ${tmp_dir}/mask_brain.mnc >/dev/null
			rawtominc -float -clobber -like ${tmp_dir}/csi_template.mnc -input ${tmp_dir}/mask_lipid.raw ${tmp_dir}/mask_lipid.mnc >/dev/null
			

			if [[ $InterpolateCSIResolution_flag -eq 1 ]]; then
				minctoraw ${tmp_dir}/mask_brain_BefInterpol.mnc -nonormalize -float > ${tmp_dir}/mask_brain_BefInterpol.raw
				minctoraw ${tmp_dir}/mask_lipid_BefInterpol.mnc -nonormalize -float > ${tmp_dir}/mask_lipid_BefInterpol.raw
				${matlabp} -r "tmp_dir = '${tmp_dir}'; ${MatlabStartupCommand}" -nodisplay -nojvm < ./flip_mask.m
				rawtominc -float -clobber -like ${tmp_dir}/csi_template_BefInterpol.mnc -input ${tmp_dir}/mask_brain_BefInterpol.raw ${tmp_dir}/mask_brain_BefInterpol.mnc >/dev/null
				rawtominc -float -clobber -like ${tmp_dir}/csi_template_BefInterpol.mnc -input ${tmp_dir}/mask_lipid_BefInterpol.raw ${tmp_dir}/mask_lipid_BefInterpol.mnc >/dev/null
			fi
		fi
		#end BET pfusch

	

	    fi


	#############################
	########   V  o  I   ########
	#############################
	if [[ $voi_found > 0 ]]; then

		rawtominc -float -like ${tmp_dir}/csi_template.mnc -input ${tmp_dir}/mask_brain_VOI.raw ${tmp_dir}/mask_brain.mnc >/dev/null	# The VOI mask was created by the MATLAB script

	if [[ $InterpolateCSIResolution_flag -eq 1 ]]; then
		rawtominc -float -like ${tmp_dir}/csi_template_BefInterpol.mnc -input ${tmp_dir}/mask_brain_VOI_BefInterpol.raw ${tmp_dir}/mask_brain_BefInterpol.mnc >/dev/null	
	fi






	#############################
	########   B  E  T   ########
	#############################
	elif [[ $bet_found > 0 ]]; then		# if mask is created with brain extraction tool


		mnc2nii -quiet ${tmp_dir}/magnitude.mnc ${tmp_dir}/nii_magnitude.nii >/dev/null					# If brain was inputted, use BET2

		if [[ $T1w_flag -eq 1 ]]; then
			${betp} ${tmp_dir}/nii_magnitude ${tmp_dir}/brain -f 0.33 -g 0 -n -m     #-f 0.5 -g 0 -n -m
		else
			${betp} ${tmp_dir}/nii_magnitude ${tmp_dir}/brain -f 0.7 -g 0 -n -m      #-f 0.5 -g 0 -n -m ##################### -Z # Improve results if FoV is very small in z-direction
		fi
		gunzip ${tmp_dir}/brain_mask.nii.gz
		rm ${tmp_dir}/mask_brain_unres.mnc
		nii2mnc -quiet ${tmp_dir}/brain_mask.nii ${tmp_dir}/mask_brain_unres.mnc
		rm ${tmp_dir}/*.nii


		#### Solving some problem with the dircos (???)
		if [[ $T1w_flag -eq 0 ]]; then
			minctoraw ${tmp_dir}/mask_brain_unres.mnc -nonormalize -float > ${tmp_dir}/mask_brain_unres.raw
			rawtominc -float -clobber -like ${tmp_dir}/mag_template.mnc -input ${tmp_dir}/mask_brain_unres.raw ${tmp_dir}/mask_brain_unres.mnc >/dev/null
		fi

		## Resample to CSI
		mincresample -clobber -nearest_neighbour -like ${tmp_dir}/csi_template.mnc ${tmp_dir}/mask_brain_unres.mnc ${tmp_dir}/mask_brain.mnc >/dev/null
		if [[ $InterpolateCSIResolution_flag -eq 1 ]]; then
			mincresample -clobber -nearest_neighbour -like ${tmp_dir}/csi_template_BefInterpol.mnc ${tmp_dir}/mask_brain_unres.mnc ${tmp_dir}/mask_brain_BefInterpol.mnc >/dev/null
		fi


		# For some reason the mask is flipped (because of the nii-stuff?). Undo this flip
		if [[ $T1w_flag -eq 0 ]]; then
			minctoraw ${tmp_dir}/mask_brain.mnc -nonormalize -float > ${tmp_dir}/mask_brain.raw
			${matlabp} -r "tmp_dir = '${tmp_dir}'; ${MatlabStartupCommand}" -nodisplay -nojvm < ./flip_mask.m
			rawtominc -float -clobber -like ${tmp_dir}/csi_template.mnc -input ${tmp_dir}/mask_brain.raw ${tmp_dir}/mask_brain.mnc >/dev/null
			
			if [[ $InterpolateCSIResolution_flag -eq 1 ]]; then
				minctoraw ${tmp_dir}/mask_brain_BefInterpol.mnc -nonormalize -float > ${tmp_dir}/mask_brain_BefInterpol.raw
				${matlabp} -r "tmp_dir = '${tmp_dir}'; ${MatlabStartupCommand}" -nodisplay -nojvm < ./flip_mask.m
				rawtominc -float -clobber -like ${tmp_dir}/csi_template_BefInterpol.mnc -input ${tmp_dir}/mask_brain_BefInterpol.raw ${tmp_dir}/mask_brain_BefInterpol.mnc >/dev/null
			fi
		fi





	#############################
	########  Threshold  ########
	#############################
	elif [[ $thresh_found > 0 ]]; then		# if mask is thresholded

		max_magnitude=`mincstats -quiet -max ${tmp_dir}/magnitude.mnc`
		lower_threshold=$(echo "scale=6 ; ${max_magnitude}/7" | bc)
		mincmath -clobber -segment -const2 $lower_threshold $max_magnitude ${tmp_dir}/magnitude.mnc ${tmp_dir}/mask_brain_unres.mnc >/dev/null

		## Resample to CSI
		mincresample -clobber -nearest_neighbour -like ${tmp_dir}/csi_template.mnc ${tmp_dir}/mask_brain_unres.mnc ${tmp_dir}/mask_brain.mnc >/dev/null
		if [[ $InterpolateCSIResolution_flag -eq 1 ]]; then
			mincresample -clobber -nearest_neighbour -like ${tmp_dir}/csi_template_BefInterpol.mnc ${tmp_dir}/mask_brain_unres.mnc ${tmp_dir}/mask_brain_BefInterpol.mnc >/dev/null
			minctoraw ${tmp_dir}/mask_brain_BefInterpol.mnc -nonormalize -float > ${tmp_dir}/mask_brain_BefInterpol.raw
		fi
		minctoraw ${tmp_dir}/mask_brain.mnc -nonormalize -float > ${tmp_dir}/mask_brain.raw
		${matlabp} -r "tmp_dir = '${tmp_dir}'; ${MatlabStartupCommand}" -nodisplay -nojvm < ./ExtractBrain_mask.m
		rawtominc -float -clobber -like ${tmp_dir}/csi_template.mnc -input ${tmp_dir}/mask_brain.raw ${tmp_dir}/mask_brain.mnc >/dev/null
		if [[ $InterpolateCSIResolution_flag -eq 1 ]]; then
			rawtominc -float -clobber -like ${tmp_dir}/csi_template_BefInterpol.mnc -input ${tmp_dir}/mask_brain_BefInterpol.raw ${tmp_dir}/mask_brain_BefInterpol.mnc >/dev/null	
		fi



	#############################
	########  ext. mask  ########
	#############################
	elif [[ "$mask_method" == */*.mnc ]]; then						# if mask is inputted, COPY MASK-FILE TO tmp-FOLDER WITH CORRECT NAME

		cp $mask_method ${tmp_dir}/mask_brain.mnc

	fi



	# Create .raw file and copy that to $out_path/maps
	minctoraw ${tmp_dir}/mask_brain.mnc -nonormalize -float > ${tmp_dir}/mask_brain.raw
	cp ${tmp_dir}/mask_brain.raw ${out_path}/maps/mask.raw
	if [[ $InterpolateCSIResolution_flag -eq 1 ]]; then
		minctoraw ${tmp_dir}/mask_brain_BefInterpol.mnc -nonormalize -float > ${tmp_dir}/mask_brain_BefInterpol.raw
		cp ${tmp_dir}/mask_brain_BefInterpol.raw ${out_path}/maps/mask_BefInterpol.raw
	fi

	# REMOVE ALL THE UNNECCESSARY STUFF
	rm ${tmp_dir}/*unres*





else

	rawtominc -float -like ${tmp_dir}/csi_template.mnc -input ${tmp_dir}/mask_brain.raw ${tmp_dir}/mask_brain.mnc
	cp ${tmp_dir}/mask_brain.raw ${out_path}/maps/mask.raw
	
	if [[ $InterpolateCSIResolution_flag -eq 1 ]]; then
		rawtominc -float -like ${tmp_dir}/csi_template_BefInterpol.mnc -input ${tmp_dir}/mask_brain_BefInterpol.raw ${tmp_dir}/mask_brain_BefInterpol.mnc >/dev/null
		cp ${tmp_dir}/mask_brain_BefInterpol.raw ${out_path}/maps/mask_BefInterpol.raw
	fi


fi	# if mask_flag = 1

















#if [[ $mask_flag -eq 1 ]]; then									# if mask is inputted, COPY MASK-FILE TO tmp-FOLDER WITH CORRECT NAME
#	cp $mask_method ${tmp_dir}/mask_brain.mnc
#	minctoraw ${tmp_dir}/mask_brain.mnc -nonormalize -float > ${tmp_dir}/mask_brain.raw
#	cp ${tmp_dir}/mask_brain.raw ${out_path}/maps/mask.raw
#else


#	if [[ $T1w_flag -eq 1 ]]; then					   # if T1_map is inputted, create magnitude minc file									
#		dcm2mnc $T1w_path -dname ${tmp_dir} -fname magnitude .
#	else 																   # In this case create mask out of GRE image or even CSI data (done in GetPar_CreateTempl_MaskPart1)
#		rawtominc -float -like ${tmp_dir}/mag_template.mnc -input ${tmp_dir}/magnitude.raw ${tmp_dir}/magnitude.mnc # Done there: READ IN DAT OF ALL CHANNELS, SoS of CHANNELS, WRITE DAT AS .RAW
#	fi

#	cp ${tmp_dir}/magnitude.mnc ${out_path}/maps/magnitude.mnc
#	rawtominc -float -like ${tmp_dir}/csi_template.mnc -input ${tmp_dir}/mask_brain_VOI.raw ${tmp_dir}/mask_brain_VOI.mnc	# The VOI mask was created by the MATLAB script


#	# Create mask
#	if [[ $use_phantom_flag -eq 1 ]]; then								# If data is phantom, don't use brain extraction tool BET2, but simple thresholding.
#		max_magnitude=`mincstats -quiet -max ${tmp_dir}/magnitude.mnc`
#		lower_threshold=$(echo "scale=6 ; ${max_magnitude}/7" | bc)
#		mincmath -clobber -segment -const2 $lower_threshold $max_magnitude ${tmp_dir}/magnitude.mnc ${tmp_dir}/mask_brain_unres.mnc

#	else

#	
#		mnc2nii -quiet ${tmp_dir}/magnitude.mnc ${tmp_dir}/nii_magnitude.nii 				# If brain was inputted, use BET2

#		if [[ $T1w_flag -eq 1 ]]; then
#			/usr/local/fsl/bin/bet ${tmp_dir}/nii_magnitude ${tmp_dir}/brain -f 0.33 -g 0 -n -m     #-f 0.5 -g 0 -n -m
#		else
#			/usr/local/fsl/bin/bet ${tmp_dir}/nii_magnitude ${tmp_dir}/brain -f 0.7 -g 0 -n -m      #-f 0.5 -g 0 -n -m
#		fi
#		gunzip ${tmp_dir}/brain_mask.nii.gz
#		rm ${tmp_dir}/mask_brain_unres.mnc
#		nii2mnc -quiet ${tmp_dir}/brain_mask.nii ${tmp_dir}/mask_brain_unres.mnc
#		rm ${tmp_dir}/*.nii
#	fi


#	# CUT SLICES OUT OF 3d-VOLUME OR JUST RESAMPLING
#	if [[ $T1w_flag -eq 0 && $use_phantom_flag -eq 0 ]]; then
#		minctoraw ${tmp_dir}/mask_brain_unres.mnc -nonormalize -float > ${tmp_dir}/mask_brain_unres.raw
#		rawtominc -float -clobber -like ${tmp_dir}/mag_template.mnc -input ${tmp_dir}/mask_brain_unres.raw ${tmp_dir}/mask_brain_unres.mnc
#	fi
#	
#	mincresample -clobber -nearest_neighbour -like ${tmp_dir}/csi_template.mnc ${tmp_dir}/mask_brain_unres.mnc ${tmp_dir}/mask_brain_noVOI.mnc

#	# For some reason the mask is flipped (because of the nii-stuff?). Undo this flip
#	if [[ $T1w_flag -eq 0 && $use_phantom_flag -eq 0 ]]; then
#		minctoraw ${tmp_dir}/mask_brain_noVOI.mnc -nonormalize -float > ${tmp_dir}/mask_brain.raw
#		${matlabp} -r "tmp_dir = '${tmp_dir}'; ${MatlabStartupCommand}" -nodisplay -nojvm < ./flip_mask.m
#		rawtominc -float -clobber -like ${tmp_dir}/csi_template.mnc -input ${tmp_dir}/mask_brain.raw ${tmp_dir}/mask_brain_noVOI.mnc
#	fi
#	
#	# Minc is so stupid to create different dircos values if I do "ratominc -like csi_template.mnc [...]" or "mincresample -like csi_template.mnc" !!! Therefor -nocheck_dimensions necessary.  
#	mincmath -nocheck_dimensions -clobber -mult ${tmp_dir}/mask_brain_noVOI.mnc ${tmp_dir}/mask_brain_VOI.mnc ${tmp_dir}/mask_brain.mnc

#	# Create .raw file and copy that to $out_path/maps
#	minctoraw ${tmp_dir}/mask_brain.mnc -nonormalize -float > ${tmp_dir}/mask_brain.raw
#	cp ${tmp_dir}/mask_brain.raw ${out_path}/maps/mask.raw


#	# REMOVE ALL THE UNNECCESSARY STUFF
#	rm ${tmp_dir}/*unres*
#	



#fi













#	##### erode outer part of mask if desired #####
#	mask_ok="n"
#	while [[ $mask_ok = "n" ]]
#	do
#		register ${tmp_dir}/mask_brain2.mnc ${tmp_dir}/magmap.mnc
#		echo -n -e "Is this mask ok? type y for continueing with processing, n for eroding mask.\n"
#		read mask_ok
#		echo -e "\nmask_ok = $mask_ok"
#		if [[ $mask_ok = "n" ]]; then
#			echo -e "\nERODING\n"
#			mincmorph -successive 'E' ${tmp_dir}/mask_brain2.mnc ${tmp_dir}/mask_brain3.mnc -clobber -2D08 -clobber	#erode outer part of brain mask;
#			#'E' means the outermost voxels get eroded. In some cases this is not useful.
#			rm ${tmp_dir}/mask_brain2.mnc
#			cp ${tmp_dir}/mask_brain3.mnc ${tmp_dir}/mask_brain2.mnc
#			rm ${tmp_dir}/mask_brain3.mnc 
#		fi
#	done
#	minctoraw ${tmp_dir}/mask_brain2.mnc -nonormalize -double > ${tmp_dir}/mask_brain.raw
#	#rm ${tmp_dir}/magmap.mnc ${tmp_dir}/res_magmap.mnc ${tmp_dir}/nii_magmap.nii ${tmp_dir}/nii_magmap.nii.gz ${tmp_dir}/res_magmap_brain_mask.nii ${tmp_dir}/mask_brain2.mnc




