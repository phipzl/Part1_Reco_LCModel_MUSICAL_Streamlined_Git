#! /bin/bash
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

# Required here because dcm2mnc only takes relative paths
part1_path=$(pwd)
cd "${tmp_dir}" || exit
tmp_dir="."

source "$part1_path/run_matlab.sh"

#read -p "Stop before creating mask1."

if [[ $T1w_flag -eq 1 ]]; then # if T1_map is inputted, create magnitude minc file

    if [[ "$T1w_path" == */*.mnc ]]; then # Copy magnitude.mnc if T1w_path is minc-file.
        cp $T1w_path ./${tmp_dir}/magnitude.mnc
    else # Only perform the dcm2mnc stuff if T1w_path is NOT a minc-file
        dcm2mnc $T1w_path -dname ./${tmp_dir} -fname magnitude .
        if [[ $T1w_AntiNoise_flag -eq 1 ]]; then # If there was inputted another 3D-measurement for removing the noise of the T1w-image around the head.
            dcm2mnc $T1w_AntiNoise_path -dname ./${tmp_dir} -fname magnitude_AntiNoise .
            max_magnitude=$(mincstats -quiet -max ./${tmp_dir}/magnitude_AntiNoise.mnc)
            lower_threshold=$(echo "scale=6 ; ${max_magnitude}/25" | bc)
            mincmath -clobber -segment -const2 $lower_threshold $max_magnitude ./${tmp_dir}/magnitude_AntiNoise.mnc ./${tmp_dir}/mask_AntiNoise.mnc
            mincmath -clobber -mult ./${tmp_dir}/magnitude.mnc ./${tmp_dir}/mask_AntiNoise.mnc ./${tmp_dir}/magnitude2.mnc

            ##GH: added for new masking. bstr: Not necessary anymore, I use magnitude_AntiNoise.mnc later for the mask creation if it exists
            #mnc2nii -quiet ./${tmp_dir}/magnitude_AntiNoise.mnc ./${tmp_dir}/magnitude_inversion_2.nii

            rm ./${tmp_dir}/magnitude.mnc
            mv ./${tmp_dir}/magnitude2.mnc ./${tmp_dir}/magnitude.mnc
            # rm ./${tmp_dir}/*AntiNoise*
        fi
    fi

else                                                                                                                    # In this case create mask out of GRE image or even CSI data (done in GetPar_CreateTempl_MaskPart1)
    $rawtomincp -float -like ./${tmp_dir}/mag_template.mnc -input ./${tmp_dir}/magnitude.raw ./${tmp_dir}/magnitude.mnc # Done there: READ IN DAT OF ALL CHANNELS, SoS of CHANNELS, WRITE DAT AS .RAW
fi
cp ./${tmp_dir}/magnitude.mnc ${out_path}/maps/magnitude.mnc

##### Create FLAIR minc
if [[ $FLAIR_flag -eq 1 ]]; then
    echo -e "\n\n Do it with FLAIR! \n\n"
    "$part1_path/Split_DICOM_Folder.sh" -f $FLAIR_path
    dcm2mnc $FLAIR_path ./${tmp_dir}/flair.mnc
    echo "File flair.mnc was created in $tmp_dir."
    cp ./${tmp_dir}/flair.mnc ${out_path}/maps/flair.mnc
fi

mkdir -p ${out_path}/maps/Extra/

##### Create B1 minc
if [[ $B1corr_flag -eq 1 ]]; then
    "$part1_path/B1_preparations.sh"
fi

#read -p "Stop before creating mask2."
echo "" && echo Mask flag: $mask_flag, method: $mask_method
if [[ $mask_flag -eq 1 ]]; then

    # Create mask out of magnitude image, VoI-Info, threshold, or copy user-given mask.
    voi_found=$(echo $mask_method | grep -c -i "voi")
    bet_found=$(echo $mask_method | grep -c -i "bet")
    thresh_found=$(echo $mask_method | grep -c -i "thresh")
    ThreeD_found=$(echo $mask_method | grep -c -i "dreid")
	if [[ $ThreeD_found > 0 ]]; then
		echo -e "\n\nWARNING: YOU ARE USING AN OUTDATED MASKING OPTION \"dreid\" WHICH NO LONGER EXISTS.\nSWITCHING TO \"bet\" INSTEAD."
		mask_method="bet"
		ThreeD_found=0; bet_found=1;
	fi


	#############################
	########   V  o  I   ########
	#############################
	if [[ -f "./${tmp_dir}/mask_brain_VOI.raw" ]]; then


		$rawtomincp -float -like ./${tmp_dir}/csi_template.mnc -input ./${tmp_dir}/mask_brain_VOI.raw ./${tmp_dir}/mask_brain_VOI.mnc >/dev/null    # The VOI mask was created by the MATLAB script

		mincresample -like ./${tmp_dir}/csi_template_zf.mnc -nearest_neighbour ./${tmp_dir}/mask_brain_VOI.mnc ./${tmp_dir}/mask_brain_VOI_zf.mnc

		if [[ $InterpolateCSIResolution_flag -eq 1 ]]; then
		    $rawtomincp -float -like ./${tmp_dir}/csi_template_BefInterpol.mnc -input ./${tmp_dir}/mask_brain_VOI_BefInterpol.raw ./${tmp_dir}/mask_brain_BefInterpol_VOI.mnc >/dev/null
		fi
	fi
	
	
	
	#############################
	########   B  E  T   ########
	#############################
	if [[ $bet_found > 0 ]]; then		# if mask is created with brain extraction tool
		
	# See if the user provided a -f and -g option for bet
		BetOptionsFound=$(echo $mask_method | grep -ci "bet,\s*")

		BetOptions=${mask_method#bet,}
		#BetOptions=$(echo $mask_method | grep -oi "\-f +\{0,1\}-\{0,1\}[0-9]\{1,\}\.*[0-9]* \-g +\{0,1\}-\{0,1\}[0-9]\{1,\}\.*[0-9]*") # Search for sth like "-f +-0.5 -g +-0.1" sd
		if [[ "$BetOptions" == "" ]]; then
			if [[ $BetOptionsFound -eq 1 ]]; then
				echo -e "\n\n\n\nWARNING: IT SEEMS YOU INPUTTED PARAMETERS FOR BET, BUT I COULD NOT RECOGNIZE THEM. DID YOU USE A WRONG FORMAT?\n\n\n"
			fi
			BetOptions="-f 0.33 -g 0"				# Default Bet option when T1w used
			if [[ $T1w_flag -eq 0 ]]; then
				BetOptions="-f 0.7 -g 0"			# Default Bet opotion otherwise
			fi
		fi
		#mnc2nii -quiet ./${tmp_dir}/magnitude.mnc ./${tmp_dir}/nii_magnitude.nii 					# If brain was inputted, use BET2
		if [[ $T1w_AntiNoise_flag -eq 1 ]]; then
        	echo -e "\nRunning bet on data provided with flag -A:\n$T1w_AntiNoise_path"		
			mnc2nii -quiet ./${tmp_dir}/magnitude_AntiNoise.mnc ./${tmp_dir}/nii_magnitude.nii
		else
			mnc2nii -quiet ./${tmp_dir}/magnitude.mnc ./${tmp_dir}/nii_magnitude.nii			
		fi
		
        echo "Running BET now. This will take some time. $(date)"
		#read -p "BeforeMask"

        echo "${betp} ./${tmp_dir}/nii_magnitude.nii ./${tmp_dir}/brain $BetOptions -n -m &"
        ${betp} ./${tmp_dir}/nii_magnitude.nii ./${tmp_dir}/brain $BetOptions -n -m & #-f 0.5 -g 0 -n -m
		if [[ $T1w_flag -eq 1 ]]; then
            echo "${betp} ./${tmp_dir}/nii_magnitude ./${tmp_dir}/lipid $BetOptions -n -A &"
            ${betp} ./${tmp_dir}/nii_magnitude ./${tmp_dir}/lipid $BetOptions -n -A &
		fi
        wait && echo "BET completed."
		
		#read -p "AfterMask"		
        for fileniigz in brain_mask.nii.gz lipid_skull_mask.nii.gz lipid_outskin_mask.nii.gz lipid_outskull_mask.nii.gz lipid_inskull_mask.nii.gz; do
            if [ -f "./${tmp_dir}/${fileniigz}" ]; then
                gunzip "./${tmp_dir}/${fileniigz}"
            fi
        done
		
		nii2mnc -quiet ./${tmp_dir}/lipid_outskin_mask.nii ./${tmp_dir}/OUTSKIN.mnc >/dev/null
		nii2mnc -quiet ./${tmp_dir}/lipid_outskull_mask.nii ./${tmp_dir}/OUTSKULL.mnc >/dev/null
		nii2mnc -quiet ./${tmp_dir}/lipid_inskull_mask.nii ./${tmp_dir}/INSKULL.mnc >/dev/null

		rm ./${tmp_dir}/mask_brain_unres.mnc
		rm ./${tmp_dir}/mask_lipid_unres.mnc

		#generate brain mask in .mnc format
		nii2mnc -quiet ./${tmp_dir}/brain_mask.nii ./${tmp_dir}/mask_brain_unres.mnc >/dev/null
				
		#generate lipid mask by subtracting the brain mask from the OUTSKIN mask.
		#for definition of OUTSKIN, OUTSKULL and INSKULL see BETsurf (http://poc.vl-e.nl/distribution/manual/fsl-3.2/bet2/)
		mincmath -sub ./${tmp_dir}/OUTSKIN.mnc ./${tmp_dir}/mask_brain_unres.mnc ./${tmp_dir}/mask_lipid_unres.mnc >/dev/null

		rm ./${tmp_dir}/*.nii

		#### Solving some problem with the dircos (???)
		if [[ $T1w_flag -eq 0 ]]; then
			minctoraw ./${tmp_dir}/mask_brain_unres.mnc -nonormalize -float > ./${tmp_dir}/mask_brain_unres.raw
			$rawtomincp -float -clobber -like ./${tmp_dir}/mag_template.mnc -input ./${tmp_dir}/mask_brain_unres.raw ./${tmp_dir}/mask_brain_unres.mnc >/dev/null
		
			minctoraw ./${tmp_dir}/mask_lipid_unres.mnc -nonormalize -float > ./${tmp_dir}/mask_lipid_unres.raw
			$rawtomincp -float -clobber -like ./${tmp_dir}/mag_template.mnc -input ./${tmp_dir}/mask_lipid_unres.raw ./${tmp_dir}/mask_lipid_unres.mnc >/dev/null
		fi

		## Resample to CSI
		mincresample -clobber -nearest_neighbour -float -like ./${tmp_dir}/csi_template.mnc ./${tmp_dir}/mask_brain_unres.mnc ./${tmp_dir}/mask_brain.mnc >/dev/null
		mincresample -clobber -nearest_neighbour -like ./${tmp_dir}/csi_template.mnc ./${tmp_dir}/mask_lipid_unres.mnc ./${tmp_dir}/mask_lipid.mnc >/dev/null
		mincresample -clobber -nearest_neighbour -like ./${tmp_dir}/csi_template_zf.mnc ./${tmp_dir}/mask_brain_unres.mnc ./${tmp_dir}/mask_brain_zf.mnc		
		if [[ $InterpolateCSIResolution_flag -eq 1 ]]; then
			mincresample -clobber -nearest_neighbour -like ./${tmp_dir}/csi_template_BefInterpol.mnc ./${tmp_dir}/mask_brain_unres.mnc ./${tmp_dir}/mask_brain_BefInterpol.mnc >/dev/null
			mincresample -clobber -nearest_neighbour -like ./${tmp_dir}/csi_template_BefInterpol.mnc ./${tmp_dir}/mask_lipid_unres.mnc ./${tmp_dir}/mask_lipid_BefInterpol.mnc >/dev/null
		fi

		# For some reason the mask is flipped (because of the nii-stuff?). Undo this flip
		if [[ $T1w_flag -eq 0 ]]; then
			minctoraw -clobber ./${tmp_dir}/mask_brain.mnc -nonormalize -float > ./${tmp_dir}/mask_brain.raw
			minctoraw ./${tmp_dir}/mask_lipid.mnc -nonormalize -float > ./${tmp_dir}/mask_lipid.raw
            run_matlab flip_mask
			$rawtomincp -float -clobber -like ./${tmp_dir}/csi_template.mnc -input ./${tmp_dir}/mask_brain.raw ./${tmp_dir}/mask_brain.mnc >/dev/null
			$rawtomincp -float -clobber -like ./${tmp_dir}/csi_template.mnc -input ./${tmp_dir}/mask_lipid.raw ./${tmp_dir}/mask_lipid.mnc >/dev/null
			
			if [[ $InterpolateCSIResolution_flag -eq 1 ]]; then
				minctoraw ./${tmp_dir}/mask_brain_BefInterpol.mnc -nonormalize -float > ./${tmp_dir}/mask_brain_BefInterpol.raw
				minctoraw ./${tmp_dir}/mask_lipid_BefInterpol.mnc -nonormalize -float > ./${tmp_dir}/mask_lipid_BefInterpol.raw
                run_matlab flip_mask
				$rawtomincp -float -clobber -like ./${tmp_dir}/csi_template_BefInterpol.mnc -input ./${tmp_dir}/mask_brain_BefInterpol.raw ./${tmp_dir}/mask_brain_BefInterpol.mnc >/dev/null
				$rawtomincp -float -clobber -like ./${tmp_dir}/csi_template_BefInterpol.mnc -input ./${tmp_dir}/mask_lipid_BefInterpol.raw ./${tmp_dir}/mask_lipid_BefInterpol.mnc >/dev/null
			fi
		fi

	#############################
	########  Threshold  ########
	#############################
	elif [[ $thresh_found > 0 ]]; then		# if mask is thresholded
		lower_threshold=$(echo $mask_method | grep -oi ",.*[0-9]\{1,\}\.*[0-9]*" | grep -oi "[0-9]\{1,\}\.*[0-9]*")

		max_magnitude=$(mincstats -quiet -max ./${tmp_dir}/magnitude.mnc)
		if [[ "$lower_threshold" == "" ]]; then
			lower_threshold=$(echo "scale=6 ; ${max_magnitude}/7" | bc)
		fi
		mincmath -clobber -segment -const2 $lower_threshold $max_magnitude ./${tmp_dir}/magnitude.mnc ./${tmp_dir}/mask_brain_unres.mnc

		## Resample to CSI
		mincresample -clobber -nearest_neighbour -like ./${tmp_dir}/csi_template.mnc ./${tmp_dir}/mask_brain_unres.mnc ./${tmp_dir}/mask_brain.mnc
		mincresample -clobber -nearest_neighbour -like ./${tmp_dir}/csi_template_zf.mnc ./${tmp_dir}/mask_brain_unres.mnc ./${tmp_dir}/mask_brain_zf.mnc

		if [[ $InterpolateCSIResolution_flag -eq 1 ]]; then
			mincresample -clobber -nearest_neighbour -like ./${tmp_dir}/csi_template_BefInterpol.mnc ./${tmp_dir}/mask_brain_unres.mnc ./${tmp_dir}/mask_brain_BefInterpol.mnc
			minctoraw ./${tmp_dir}/mask_brain_BefInterpol.mnc -nonormalize -float > ./${tmp_dir}/mask_brain_BefInterpol.raw
		fi
		minctoraw ./${tmp_dir}/mask_brain.mnc -nonormalize -float > ./${tmp_dir}/mask_brain.raw
		minctoraw ./${tmp_dir}/mask_brain_zf.mnc -nonormalize -float > ./${tmp_dir}/mask_brain_zf.raw
        run_matlab ExtractBrain_mask
		$rawtomincp -float -clobber -like ./${tmp_dir}/csi_template.mnc -input ./${tmp_dir}/mask_brain.raw ./${tmp_dir}/mask_brain.mnc
		$rawtomincp -float -clobber -like ./${tmp_dir}/csi_template_zf.mnc -input ./${tmp_dir}/mask_brain_zf.raw ./${tmp_dir}/mask_brain_zf.mnc
		if [[ $InterpolateCSIResolution_flag -eq 1 ]]; then
			$rawtomincp -float -clobber -like ./${tmp_dir}/csi_template_BefInterpol.mnc -input ./${tmp_dir}/mask_brain_BefInterpol.raw ./${tmp_dir}/mask_brain_BefInterpol.mnc
			cp ./${tmp_dir}/mask_brain_BefInterpol.mnc ./${tmp_dir}/mask_lipid_BefInterpol.mnc	
		fi		
		# For thresh method, the lipid is useless. Create lipid masks only to avoid if conditions later
		cp ./${tmp_dir}/mask_brain.mnc ./${tmp_dir}/mask_lipid.mnc				

	#############################
	########  ext. mask  ########
	#############################
	elif [[ "$mask_method" == */*.mnc ]]; then						# if mask is inputted, COPY MASK-FILE TO tmp-FOLDER WITH CORRECT NAME

		mincresample -nearest_neighbour -like ./${tmp_dir}/csi_template_zf.mnc $mask_method ./${tmp_dir}/mask_brain_zf.mnc
		cp $mask_method ./${tmp_dir}/mask_brain.mnc

	fi	
	
    #########################################
    ########   V  o  I (Once again)  ########
    #########################################
	if [[ -f "./${tmp_dir}/mask_brain_VOI.mnc" ]]; then

        if [[ -f "./${tmp_dir}/mask_brain.mnc" ]]; then
            mincmath -nocheck_dimensions -mult ./${tmp_dir}/mask_brain_VOI.mnc ./${tmp_dir}/mask_brain.mnc ./${tmp_dir}/mask_brain2.mnc
            mincmath -nocheck_dimensions -mult ./${tmp_dir}/mask_brain_VOI_zf.mnc ./${tmp_dir}/mask_brain_zf.mnc ./${tmp_dir}/mask_brain2_zf.mnc
            mincmath -nocheck_dimensions -mult ./${tmp_dir}/mask_brain_VOI.mnc ./${tmp_dir}/mask_lipid.mnc ./${tmp_dir}/mask_lipid2.mnc
            rm ./${tmp_dir}/mask_brain.mnc
            rm ./${tmp_dir}/mask_brain_zf.mnc
            rm ./${tmp_dir}/mask_lipid.mnc
            mv ./${tmp_dir}/mask_brain2.mnc ./${tmp_dir}/mask_brain.mnc
            mv ./${tmp_dir}/mask_brain2_zf.mnc ./${tmp_dir}/mask_brain_zf.mnc
            mv ./${tmp_dir}/mask_lipid2.mnc ./${tmp_dir}/mask_lipid.mnc
        else
            mv ./${tmp_dir}/mask_brain_VOI.mnc ./${tmp_dir}/mask_brain.mnc
            mv ./${tmp_dir}/mask_brain_VOI_zf.mnc ./${tmp_dir}/mask_brain_zf.mnc
        fi
    fi

    #########################################
    ########   postprocessing        ########
    #########################################
    if [[ -f "./${tmp_dir}/mask_brain_BefInterpol_VOI.mnc" ]]; then
        if [[ -f "./${tmp_dir}/mask_brain_BefInterpol.mnc" ]]; then
            mincmath -nocheck_dimensions -mult ./${tmp_dir}/mask_brain_BefInterpol_VOI.mnc ./${tmp_dir}/mask_brain_BefInterpol.mnc ./${tmp_dir}/mask_brain2.mnc
            mincmath -nocheck_dimensions -mult ./${tmp_dir}/mask_brain_BefInterpol_VOI.mnc ./${tmp_dir}/mask_lipid_BefInterpol.mnc ./${tmp_dir}/mask_lipid2.mnc
            rm ./${tmp_dir}/mask_brain_BefInterpol.mnc
            rm ./${tmp_dir}/mask_lipid_BefInterpol.mnc
            mv ./${tmp_dir}/mask_brain2.mnc ./${tmp_dir}/mask_brain_BefInterpol.mnc
            mv ./${tmp_dir}/mask_lipid2.mnc ./${tmp_dir}/mask_lipid_BefInterpol.mnc
        else
            mv ./${tmp_dir}/mask_brain_BefInterpol_VOI.mnc ./${tmp_dir}/mask_brain_BefInterpol.mnc
        fi
    fi

    # Create .raw file and copy that to $out_path/maps
    minctoraw ./${tmp_dir}/mask_brain.mnc -nonormalize -float >./${tmp_dir}/mask_brain.raw
    cp ./${tmp_dir}/mask_brain.raw "${out_path}/maps/mask.raw"

    if [[ -f ./${tmp_dir}/mask_brain_zf.mnc ]]; then
        minctoraw ./${tmp_dir}/mask_brain_zf.mnc -nonormalize -float >./${tmp_dir}/mask_brain_zf.raw
        cp ./${tmp_dir}/mask_brain_zf.raw "${out_path}/maps/mask_zf.raw"
    fi

    if [[ -f ./${tmp_dir}/mask_lipid.mnc ]]; then
        minctoraw ./${tmp_dir}/mask_lipid.mnc -nonormalize -float >./${tmp_dir}/mask_lipid.raw
        cp ./${tmp_dir}/mask_lipid.raw "${out_path}/maps/mask_lipid.raw"
    fi

    if [[ $InterpolateCSIResolution_flag -eq 1 ]]; then
        minctoraw ./${tmp_dir}/mask_brain_BefInterpol.mnc -nonormalize -float >./${tmp_dir}/mask_brain_BefInterpol.raw
        cp ./${tmp_dir}/mask_brain_BefInterpol.raw "${out_path}/maps/mask_BefInterpol.raw"
        if [[ -f ./${tmp_dir}/mask_lipid_BefInterpol.mnc ]]; then
            minctoraw ./${tmp_dir}/mask_lipid_BefInterpol.mnc -nonormalize -float >./${tmp_dir}/mask_lipid_BefInterpol.raw
            cp ./${tmp_dir}/mask_lipid_BefInterpol.raw "${out_path}/maps/mask_lipid_BefInterpol.raw"
        fi
    fi

    # REMOVE ALL THE UNNECCESSARY STUFF
    if [[ -f ./${tmp_dir}/mask_brain_unres.mnc ]]; then
        mv ./${tmp_dir}/mask_brain_unres.mnc ./${tmp_dir}/mask_brain_hires.mnc # Need this mask for unwrapping the B0-map. Will be deleted later.
        rm ./${tmp_dir}/*unres*
    else
        if [[ $AlignFreq_flag -eq 1 ]] && ! [[ ${AlignFreq_path} == "" ]]; then
            mincresample -clobber -nearest_neighbour -like ./${tmp_dir}/AlignFreq_template.mnc ./${tmp_dir}/mask_brain_zf.mnc ./${tmp_dir}/mask_brain_hires.mnc
        fi
    fi

else # If no mask-flag was used. There will be still a mask full of ones, created previously in Matlab

    $rawtomincp -float -like ./${tmp_dir}/csi_template.mnc -input ./${tmp_dir}/mask_brain.raw ./${tmp_dir}/mask_brain.mnc
    $rawtomincp -float -like ./${tmp_dir}/csi_template_zf.mnc -input ./${tmp_dir}/mask_brain_zf.raw ./${tmp_dir}/mask_brain_zf.mnc
    cp ./${tmp_dir}/mask_brain.raw ${out_path}/maps/mask.raw
    cp ./${tmp_dir}/mask_brain_zf.raw ${out_path}/maps/mask_zf.raw

    if [[ $InterpolateCSIResolution_flag -eq 1 ]]; then
        $rawtomincp -float -like ./${tmp_dir}/csi_template_BefInterpol.mnc -input ./${tmp_dir}/mask_brain_BefInterpol.raw ./${tmp_dir}/mask_brain_BefInterpol.mnc
        cp ./${tmp_dir}/mask_brain_BefInterpol.raw ${out_path}/maps/mask_BefInterpol.raw
    fi

fi # if mask_flag = 1
