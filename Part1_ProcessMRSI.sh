#!/bin/bash
###################################################################################################
###	Process 2D/3D/Multislice CSI data, prepare them for LCModel Processing and start LCModel	###
###	   Data can be miltichannel data in which case GRE/reference imaging data is used to		###
###	                             combine the individual channels.	      		  				###
###	Data can also be undersampled, in which case GRAPPA-based 2D-Caipirinha and 1D-Caipirinha	###
###								algorithms are used for Reco.									###
###################################################################################################

# -1. Preparations
# In case you hit ctrl-c, kill all processes, also background processes. Trap all exit signals and the SIGUSR1 signal.
trap "{ Trapped=1; echo Try to kill: PID: $$, PPID: $PPID; TerminateProgram; echo 'Kill all processes due to user request.'; sleep 3s; kill 0;}" SIGINT SIGUSR1
TerminateProgram() {

    echo "Stop tee."
    # close and restore backup; both stdout and stderr

    exec 1>&6 # duplicate 6 to 1 again
    exec 6>&- # close 6
    exec 2>&7 # duplicate 7 to 2 again
    exec 7>&- # close 7
    sleep 1

    echo -e "\n\nTerminate Program & Backup: "

    cd "$calldir" || exit

    # Copy the logfile
    echo "Copy logfile to $out_path/UsedSourcecode/logfile.log."
    if [ -d "$out_path/UsedSourcecode" ]; then
        cp "$logfile" "$out_path/UsedSourcecode"
        cp "$logfile" "$out_path"
        cp "$tmp_dir/ErrorFile.sh" "$out_path/UsedSourcecode"
    fi

    if [ -f "$out_path/TempServerDir/KillLCMProcesses.sh" ]; then
        if [[ $RunLCModelOn == "" ]] || [[ $CurrentComputer == "$RunLCModelOn" ]]; then
            "$out_path/TempServerDir/KillLCMProcesses.sh" # kill locally
        else
            if [[ $RunLCModelAs == "" ]]; then
                ssh -o PasswordAuthentication=no "$RunLCModelOn" "$out_path/TempServerDir/KillLCMProcesses.sh"
            else
                ssh -o PasswordAuthentication=no -l "$RunLCModelAs" "$RunLCModelOn" "$out_path/TempServerDir/KillLCMProcesses.sh"
            fi
        fi
    fi
    rm -fR "$out_path/TempServerDir"

    if [[ "$1" == "1" ]]; then # DebugFlag = 1
        rmtmpdir="n"
    else
        rmtmpdir="y"
    fi
    echo -e "\n$rmtmpdir"
    if [[ "$rmtmpdir" == "y" ]]; then

        if [ -n "${pid_list[1]}" ]; then
            echo "Wait for all LCModel Processes to terminate gracefully."
            sleep 40 # So that all lcmodel processes can close
        fi

        rm -R -f "$tmp_dir"
        if [ -d "$tmp_dir" ]; then
            sleep 10
            rm -R -f "$tmp_dir" # Try again if it didnt work
        fi
    fi
    echo "Stop now."
    echo -e "\n\n\n\t\tE N D\n\n\n"

    if [[ "$Trapped" == "0" ]]; then
        exit 0
    fi
}


# -1.0 Debug flag
DebugFlag=0
Trapped=0

# -1.1 Install Paths etc.
echo -e "Install Program\n\n"
source "$(dirname "${BASH_SOURCE[0]}")/InstallProgramPaths.sh"
source "$(dirname "${BASH_SOURCE[0]}")/run_matlab.sh" # convenience function to run matlab scripts, defined externally to be used in other scripts as well

# -1.2 Change to directory of script
calldir=$(pwd)
cd "$(dirname "${BASH_SOURCE[0]}")" || exit

# -1.3 Create directories
tmp_trunk="tmp"
tmp_num=1
tmp_dir="${tmp_folder}/${tmp_trunk}${tmp_num}"
while [ -d "$tmp_dir" ]; do
    ((tmp_num = tmp_num + 1))
    tmp_dir="${tmp_folder}/${tmp_trunk}${tmp_num}"
done
export tmp_dir
mkdir "$tmp_dir"
if [ ! -d "$tmp_dir" ]; then
    echo "Could not create temporary directory. Exiting."
    exit 1
fi
chmod 775 "$tmp_dir"
abs_tmp_dir=$(readlink -f "$tmp_dir")

# -1.4 Write the script output to a logfile
logfile=${tmp_dir}/logfile.log
hostyy=$(hostname)
echo -e "Run Script on $hostyy with parameters:\n$0 $*\n\n"
echo -e "Run Script on $hostyy with parameters:\n$0 $*" >"$logfile"

# backup the original filedescriptors, first
# stdout (1) into fd6; stderr (2) into fd7
exec 6<&1                        # Copy 1 to 6
exec 7<&2                        # Copy 2 to 7
exec > >(tee -a "$logfile")      # Copy 1 to tee which writes to logfile
exec 2> >(tee -a "$logfile" >&2) # Copy 2 to tee which writes to logfile and redirect to 1 and 2 [???]

# INTRODUCING COLORS
NC='\033[0m'    # No color, not bold
BOLDF='\033[1m' # Bold font
RED='\033[0;31m'
ORA='\033[0;33m'
GRE='\033[0;32m'
BLU='\033[0;34m'
PUR='\033[0;35m'

printf "\n\n0.\t\t${BOLDF}S T A R T ${NC}(PID $$, PPID $PPID)\n\n\n"
date
sleep 3s
# 0.
############# DEFINE ARGUMENTS/PARAMETER OPTIONS ####################
# export so that every child process (like readcsi.sh), grandchild-process etc. can use these variables

# FLAGS
# mandatory
export csi_flag=0
export basis_flag=0
export out_flag=0

# optional
export image_normal_flag=0
export image_flip_flag=0
export image_VC_flag=0
export T1w_flag=0
export T1w_AntiNoise_flag=0
export FLAIR_flag=0
export WaterReference_flag=0
export mask_flag=0
export hamming_flag=0
export LipidDecon_flag=0
export NuisRem_flag=0
export TwoDCaipParallelImaging_flag=0
export SliceParallelImaging_flag=0
export noisedecorrelation_flag=0
export FirstOrderPhaseCorr_flag=0
export FirstOrderPhaseModulation_flag=0
export use_phantom_flag=0
export ZeroFillMetMaps_flag=0
export InterpolateCSIResolution_flag=0
export TimeInterpolation_flag=0
export AlignFreq_flag=0
export dont_compute_LCM_flag=0
export LCM_ControlPath_flag=0
export LCM_ControlPath_Water_flag=0
export exponential_filter_Hz_flag=0
export control_echo_flag=0
export basis_echo_flag=0
export XPACE_motion_correction_flag=0
export julia_reconstruction=0
export B1corr_flag=0
export NonCartTraj_flag=0
export compiled_matlab_flag=0
#BOW
export priors_flag=0
export DebugAdditionalInput_flag=0
export GradientDelay_flag=0

# INITIALIZING
export phase_encoding_direction_is_RL_flag=0
export phase_encoding_direction_is="AP"
LipidDecon_MethodAndNoOfLoops="L1,10"
export julia_n_threads="auto"
export julia_mmap="false"

while getopts 'c:b:o:a:A:B:D:e:E:f:g:G:h:i:I:j:J:k:L:m:n:p:P:r:R:s:S:t:T:v:w:W:X:z:dFKlu?' OPTION; do
    case $OPTION in

    #mandatory
    c)
        export csi_flag=1
        export csi_path="$OPTARG"
        ;;
    b)
        export basis_flag=1
        export basis_path="$OPTARG"
        ;;
    o)
        export out_flag=1
        export out_path="$OPTARG"
        ;;

    #optional
    a)
        export T1w_AntiNoise_flag=1
        export T1w_AntiNoise_path="$OPTARG"
        ;;
    A)
        export AlignFreq_flag=1
        export AlignFreq_MethodAndPath="$OPTARG"
        ;;
    B)
        export B1corr_flag=1 # flag is checked for in create_mask.sh
        export B1_path="$OPTARG"
        ;;
    D)
        export DebugAdditionalInput_flag=1
        export DebugAdditionalInput="$OPTARG"
        ;;
    e)
        export exponential_filter_Hz_flag=1
        export exponential_filter_Hz="$OPTARG"
        ;;
    E)
        export basis_echo_flag=1
        export basis_echo_path="$OPTARG"
        ;;
    f)
        export image_flip_flag=1
        export image_flip_path="$OPTARG"
        ;;
    g)
        export noisedecorrelation_flag=1
        export noisedecorrelation_path="$OPTARG"
        ;;
    G)
        export GradientDelay_flag=1
        export GradientDelay="$OPTARG"
        ;;
    h)
        export hamming_flag=1
        export hamming_factor="$OPTARG"
        ;;
    i)
        export image_normal_flag=1
        export image_normal_path="$OPTARG"
        ;;
    I)
        export InterpolateCSIResolution_flag=1
        export InterpolateCSIResolution="$OPTARG"
        ;;
    j)
        export LCM_ControlPath_flag=1
        export LCM_ControlPath="$OPTARG"
        ;;
    J)
        export control_echo_flag=1
        export control_echo_path="$OPTARG"
        ;;
    k)
        export FirstOrderPhaseModulation_flag=1
        export FID_Truncation_in_ms="$OPTARG"
        ;;
    L)
        export LipidDecon_flag=1
        export LipidDecon_MethodAndNoOfLoops="$OPTARG"
        ;;
    m)
        export mask_flag=1
        export mask_method="$OPTARG"
        ;;
    n)
        export NuisRem_flag=1
        export NuisRem_ControlPath="$OPTARG"
        ;;
    p)
        export FLAIR_flag=1
        export FLAIR_path="$OPTARG"
        ;;
    P)
        export priors_flag=1
        export priors_path="$OPTARG"
        ;;
    r)
        export TwoDCaipParallelImaging_flag=1
        export InPlaneCaipPattern_And_VD_Radius="$OPTARG"
        ;;
    R)
        export SliceParallelImaging_flag=1
        export SliceAliasingPattern="$OPTARG"
        ;;
    s)
        export NonCartTraj_flag=1
        export NonCartTrajFile_path="$OPTARG"
        ;;

    S)
        export julia_reconstruction=1
        export julia_n_threads="$OPTARG"
        julia_mmap=${!OPTIND}
        if [[ -z $julia_mmap ]]; then
            julia_mmap="false"
        fi
        ;;
    t)
        export T1w_flag=1
        export T1w_path="$OPTARG"
        ;;
    T)
        export TimeInterpolation_flag=1
        export TimeInterpolationFactor="$OPTARG"
        ;;
    v)
        export image_VC_flag=1
        export image_VC_path="$OPTARG"
        ;;
    w)
        export WaterReference_flag=1
        export WaterReference_MethodAndFile="$OPTARG"
        ;;
    W)
        export LCM_ControlPath_Water_flag=1
        export LCM_Control_Water_path="$OPTARG"
        ;;
    X)
        export XPACE_motion_correction_flag=1
        export XPACE_motion_correction_path="$OPTARG"
        ;;
    z)
        export ZeroFillMetMaps_flag=1
        export ZeroFillMetMaps="$OPTARG"
        ;;

    # Flags
    F)
        export FirstOrderPhaseCorr_flag=1
        ;;
    K)
        export compiled_matlab_flag=1
        ;;
    l)
        export dont_compute_LCM_flag=1
        ;;
    u)
        export use_phantom_flag=1
        ;;

    ?)
        printf "

Usage: %s

HINT: Things in {} are optional.
mandatory:
-c  [csi file]          Format: DAT, DICOM, or .mat. If a .mat file is passed over, it is expected that everything is already performed like coil combination etc.
                        You can pass over several files of the same type by \'-c \"[csi_path1] [csi_path2] ...\"\'. These files get individually processed and averaged
                        at the end.
-b  [basis files]       Format: .BASIS. Used for LCM fitting (for FID)
-o  [output directory]

optional:
-a  [T1 AntiNoise images]   Format: DICOM. Folder of 3d T1-weighted acquisition containing DICOM files. Used for pre-masking the T1w image to get rid of the noise in air-areas.
-A  [\"Alignment\" Or \"Alignment,Path\" Or \"Overdiscrete,Path\"]  Perform frequency alignment, either based on a given B0-map or based on a dot-product correlation function. The correction can be also done overdiscrete. If a mnc file is given, use this as B0-map, otherwise shift according to water peak of center voxel. For dicom files, provide the folder with the magnitude images, and the phasemap-difference btw the two TEs, e.g. \"Alignment, B0MagPath B0PhaPath\".
-B  [B1 reading]            Path of B1 DICOM data, used for B1 correction.
-D  [DebugAdditionalInput]  A general parameter to provide some additional, not specified input for debug purposes. This should not be used in the stable version of the pipeline, but just if you want to test something quickly.
-e  [LineBroadeningInHz]    Apply an exponential filter to the spectra [Hz].
-E  [basis files]           Format: .BASIS. Used for LCM fitting (spin echo: for fidesi =  fid + echo)
-f  [image FLIP]            Format: DAT or DICOM. Imaging file FLIP (FOV rotated about -180 deg). Used for correcting gradient delays.
-g  [noisedecorr_path]      If this option is used the csi data gets noise decorrelated using noise from passed-over noise file, or if -g \" is given, by noise from the end of the FIDs at the border of the FoV or from the PRESCAN, if available.
-G  [GradDelayPerAngInt_x = [...], GradDelayPerAngInt_y = [...] or GradDelayPerTempInt_x = [...], GradDelayPerTempInt_y = [...]]    Gradient delay in microseconds for CRT sequence. Not used if a mat file is provided with flag \"-s\", and only used for CRT trajectories (so far). The gradient delay can be specified per temporal interleaf (use GradDelayPerTempInt_x/y = ... then), or per angular interleaf (use GradDelayPerAngInt_x/y = ... then). You can specify only one number for x and y, which will then be used for all angular interleaves. If -G is not used no gradient delay is used. If -G \"Default\" or \"\" is used the default values written in InstallProgramPaths.sh are used.
-h  [100]                   Hamming filter CSI data.
-i  [image NORMAL]          Format: DAT or DICOM. The FoV must match that of the CSI file. Used for our coil combination and for creating mask (if no T1 is inputted)
-I  [\"nextpow2\" or \"[x y z]{,kSpace,Ellip}\"]    If nextpow2: Perform zerofilling to the next power of 2 in ROW and COL dimensions (e.g. from 42x42 to 64x64). If vector (e.g. [16 16 1]): Spatially Interpolate to this size. If \",kspace\" is used, perform interpolation in k-space (cut or zerofill in k-space). If additionally \",Ellip\" is used, the k-space after zerofilling/cutting to [x y z] gets elliptically filtered.
-j  [LCM_ControlFile]       ControlFile telling LCModel how to process the data. (for FID) otherwise standard values are assumed. A template file is provided in this package.
-J  [LCM_ControlFile]       ControlFile telling LCModel how to process the data. for ECHO
-L  [LipidRegMethod,RegTerm]    Perform lipid regularization after Bilgic et al. Use either \"L2,[RegTerm]\" or \"L1,Iter\" where RegTerm is a value that penalizes the lipid contamination, and Iter is the number of iterations the L1-regularization should be done. Best method is to try different values, bc unfortunately the data is not normalized, and thus very different values might be needed for different data.
-m  [mask]                  Defines how to create the mask. Options: -m \"bet{,-f +-x.yz -g +-a.bc}\", \"thresh{,lower_threshold=x}\", \"voi\", \"[Path_to_usermade_mask]\". where things in {} are optional, x is a float defining the lower thresold for masking the magnitude. If -m option is not set --> no mask used.
-n  [NuisRemControlFile]    Perform nuisance removal using hsvd according to Chao et al. The control file must specify the number of singular values, the ppm range for water and lipids and the T2's etc. This file must be in MATLAB-format. Please dont write crap in there causing MATLAB to crash or worse...
-p  [FLAIR reading]         path of FLAIR DICOM data.
-P  [prior_knowledge dir]   define directory that contains Extra maps (i.e., 0_pha_map.mnc, 1_pha_map.mnc, shift_map.mnc). THIS WORKS ONLY FOR USING MEGA-OFF PRIOR KNOWLEDGE FOR MEGA-DIFF-FITTING (because 180 deg is added to the 0-order phases)
-r  [InPlaneCaipPattern_And_VD_Radius]  The InPlaneCaipPattern and the VD_Radius as used in ParallelImagingSimReco.m. Example: \"InPlaneCaipPattern = [0 0 0; 0 0 0; 0 0 1]; VD_Radius = 2;\".
-R  [SliceAliasingPattern]
-s  [NonCartTrajFile_path]  Use this option if CSI Data is raw NonCartesian data and pass over trajectory file/path in .m or .mat file (.m file for theoretical \"calculated\" gradients, .mat file for \"measured\" trajectory). For some trajectories (CRT, Antoines rosette/eccentric, egg-trajectory) this is not necessary, as the read-in functions can automatically calculate the trajectory based on the header information. If a measured trajectory is provided with a mat file, it may contain a variable StartingPointAfterLaunchTrack which needs to be a cell with one entry for each angular interleaf, each containing one number saying how many ADC points should be omitted at the beginning in case the measured trajectory was calculated only from a later time point. The file must contain a variable kSpaceTrajectory with subfield .GM (GradientMoment) being a cell with one entry for each angular interleaf, each being a matrix of size [2 ADCPtsPerCircle].
-S  [threads] [mmap]        Use the Julia reconstruction version (less RAM usage, different reconstruction algorithm). [threads=auto] can be auto or a number. [mmap=false] can be \"true\", \"false\" or a path.
-t  [T1 images]             Format: DICOM. Folder of 3d T1-weighted acquisition containing DICOM files. Used for creating mask and for visual purposes. If minc file is given instead of folder, it is treated as the magnitude file.
-T  \"[TruncateFactor ZerofillFactor FillToOrig]\"  Interpolation of FID data in time domain using truncation and zerofilling. TruncateFactor determines how much of the orignal data is left after truncation and must be a value from 0 to 1. ZerofillFactor determines how far the zerofilling happens (relative to the truncated data) and must be larger >1. If FillToOrig is 1, the data is truncated to TruncateFactor and afterwards filled up to the original length (ZerofillFactor is irrelevant in this case). Example: To truncate down to 50 percent and then zerofill to 2x the original size, use [0.5 4 0].
-v  [VC image]              Format: DAT or DICOM. Image of volume or body coil file. Used for sensmap method or for creating mask.
-w  [Water Reference]       Format: DAT or DICOM. LCModel 'Do Water Scaling' or separate water quantification (Water maps are created). The same scan as -c [csi file], but without water suppression.
-X  [XPACE MOTION LOG]      XPACE MOTION LOG

Flags:
-F  If this option is set, the spectra are corrected for the first order phase caused by an acquisition delay of the FID-sequences. You must provide a basis set with an appropriate acquisition delay. DONT USE WITH SPIN ECHO SEQUENCES.
-K	Use compiled MATLAB functions.
        No MATLAB license needed, but the functions must be compiled first (See compile.m)
-l  If this option is set, LCModel is not started, everything else is done normally. Useful for only computing the SNR.
-u  If a phantom was measured. Different settings used for fitting (e.g. some metabolites are omitted)

" $(basename $0) >&2

        exit 2
        ;;
    esac
done

shift $((OPTIND - 1))

###0. GH: measure time elapsed:
START=$(date +%s.%N)

# 1. Create Out Directories
echo -e "\n\n1. Create Directories\n\n"
rm -Rf "$out_path/"*.mat
rm -Rf "$out_path/scalings/"*.mat
mkdir -p "$out_path/maps"
mkdir -p "$out_path/phamaps"
mkdir -p "$out_path/spectra"
mkdir -p "$out_path/scalings"

# Set kSpaceCorrection to Default value written in InstallProgramPaths.sh
if [[ $GradientDelay_flag -eq 1 ]]; then
    if [[ $GradientDelay == "" || $GradientDelay == "Default" || $GradientDelay == "default" ]]; then
        if [ -n "$DefaultGradientDelaysForCRTTrajectory" ]; then
            export GradientDelay_flag=1
            export GradientDelay="$DefaultGradientDelaysForCRTTrajectory"
        fi
    fi
fi

# 2. Write Initial Parameters
echo -e "\n\n2. Write Initial Parameters\n\n"
. ./write_InitialParameters.sh

# read -rp "Stop before Gathering info."
# 3.
############ CREATE TEMPLATES FOR CREATING METABOLIC MAPS LATER ############
echo -e "\n\n3. Gather Information, Create Minc Templates, Prepare Mask Creation.\n\n"
# Run it twice for water referencing (it will automatically process the water ref file the first time, and the other the second time)
if [[ $LCM_ControlPath_Water_flag -eq 1 ]]; then
    run_matlab GetPar_CreateTempl_MaskPart1
fi
run_matlab GetPar_CreateTempl_MaskPart1

# Terminate if there was an error
bash "$tmp_dir/ErrorFile.sh"
if [[ $ErrorInGetPar_CreateTempl -eq 1 ]]; then
    TerminateProgram $DebugFlag
fi

# read -rp "stop before create minc template"
bash "$tmp_dir/CreateMincTemplates.sh"

# read -rp "Stop before creating mask."
## 4.
############# USE MASK OR CREATE MASK OUT OF (IN PRIORITY ORDER): MASK, T1_MAP, IMAGING AC, IMAGING VC, CSI ############
echo -e "\n\n4. CREATE MASK\n\n"
./create_mask.sh

# read -p "Stop before creating B0Map."
## 5.
############# CREATE B0MAP FOR USAGE OF FREQUENCY ALIGNING CSI DATA ############
if [[ $AlignFreq_flag -eq 1 ]] && ! [[ $AlignFreq_MethodAndPath == "" ]]; then
    echo -e "\n\n5. CREATE B0MAP FOR USAGE OF FREQUENCY ALIGNING CSI DATA\n\n"
    ./create_B0Map.sh
fi
if [[ -f "$tmp_dir/mask_brain_hires.mnc" ]]; then
    rm "$tmp_dir/mask_brain_hires.mnc"
fi

# Convert prior knowledge maps to raw if necessary
if [[ $priors_flag -eq 1 ]]; then
    for CurFile in 0_pha_map 1_pha_map shift_map; do
        if [[ ! -f "$priors_path/$CurFile.raw" ]] && [[ -f "$priors_path/$CurFile.mnc" ]]; then
            echo -e "\nconvert $priors_path/$CurFile.mnc"
            minctoraw "$priors_path/$CurFile.mnc" -nonormalize -float >"$priors_path/$CurFile.raw"
        fi
    done
fi

# 6.
###########   Process Data, Prepare LCModel Fitting   ############
echo -e "\n\n\n6. Process Data and Prepare LCModel Processing, first run\n\n"
# Run it twice for water referencing (it will automatically process the water ref file the first time, and the other the second time)
if [[ $WaterReference_flag -eq 1 ]]; then
    echo -e "\nProcess water reference data for water scaling."
    echo -e "\n\nRunning:\n"
    run_mrsi_reconstruction 1
fi
echo -e "\n\n\n6. Process Data and Prepare LCModel Processing, second run\n\n"
date
# If we pass over several IMA or dat files, average them
for ((CurAvg = 1; CurAvg <= NumberOfCSIFiles; CurAvg = CurAvg + 1)); do
    run_mrsi_reconstruction $CurAvg
done

# read -rp "Stop before LCModel fitting"
#7.
########### START LCMODEL PROCESSING OF SINGLE VOXEL DATA ON CPU CORES ############
echo -e "\n\n7. Start LCModel Processing\n\n"
if [[ $dont_compute_LCM_flag -eq 0 ]]; then
    curdir=$(pwd)
    CurrentComputer=$(hostname)

    # Allow either $CurrentComputer or $RunLCModelOn to include a domain name like ".nmr.meduniwien.ac.at"
    Search1=$(echo "$CurrentComputer" | grep -c "${RunLCModelOn}\.")
    Search2=$(echo "$RunLCModelOn" | grep -c "${CurrentComputer}\.")
    if [[ "$RunLCModelOn" == "" ]] || [[ $CurrentComputer == "$RunLCModelOn" ]] || [[ $Search1 -gt 0 ]] || [[ $Search2 -gt 0 ]]; then
        bash "${tmp_dir}/lcm_process_core_parallel.sh"
    else
        rm -fR "$out_path/TempServerDir"
        cp -R "$tmp_dir/" "$out_path/TempServerDir"
        RunFileOnServer=$out_path/TempServerDir/RunLCModel.sh
        cp "$curdir/RunLCModel.sh" "$RunFileOnServer"
        # Run LCModel on different computer,connecting via ssh. You need a key so that you can automatically connect to this computer, without needing to type in the password!
        if [[ "$RunLCModelAs" == "" ]]; then
            ssh -o PasswordAuthentication=no "$RunLCModelOn" "$RunFileOnServer $RunLCModelOn $out_path/TempServerDir"
        else
            ssh -o PasswordAuthentication=no -l "$RunLCModelAs" "$RunLCModelOn" "$RunFileOnServer $RunLCModelOn $out_path/TempServerDir"
        fi
        sleep 40
        rm -fR "$out_path/TempServerDir"
    fi
    sleep 10
fi

# 8.
############ WRITE THE SOURCECODE THAT WAS USED TO OUT-DIR ############
echo -e "\n\n8. Write the used sourcecode to out-dir.\n\n"

# Copy Program where this program lies in
curdir=$(pwd)
ScriptName=$(basename "$curdir")
curlogname=$(basename "$logfile")
BaseNameMatlabFunctions=$(basename "$MatlabFunctionsFolder")

# Archive this script itself
cd ..
tar --exclude='*/tmp*' --transform='s,^,/UsedSourcecode/,' -cf "$out_path/UsedSourcecode_Part1.tar" "$ScriptName"

# Copy the logfile
cp "$logfile" "$out_path/logfile_part1.log"

# Copy MeasurementInfos
cp "$tmp_dir/MeasurementInfos.txt" "$out_path"

# Archive the logfile
cd "$tmp_dir" || exit
tar f "$out_path/UsedSourcecode_Part1.tar" -r "$curlogname"

# Archive the Matlab functions
if ! [[ "$curdir/$BaseNameMatlabFunctions" == "$MatlabFunctionsFolder" ]]; then
    cd "$MatlabFunctionsFolder" || exit
    cd ..
    tar f "$out_path/UsedSourcecode_Part1.tar" -r "$BaseNameMatlabFunctions" --transform='s,^,/UsedSourcecode/,'
fi

# zip everything
gzip "$out_path/UsedSourcecode_Part1.tar"

# Go back to the original folder
cd "$curdir" || exit

# Copy measurement information to out_path (see read_csi_dat_new_v2.m)
cp "$tmp_dir/"*.txt "$out_path"

#9.
############ REMOVE UNECESSARY DATA ############
if [[ $DebugFlag -eq 0 ]]; then
    echo -e "\n\n9. Remove unnecessary data!\n\n"
    find "$out_path" -name '*.RAW' -delete
    find "$out_path" -name '*.control' -delete
fi

#10: GH: finish time measurement
END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
echo -e "\n\n10. Total Processing time of Part 1: " "$DIFF" "\n\n"

TerminateProgram $DebugFlag
