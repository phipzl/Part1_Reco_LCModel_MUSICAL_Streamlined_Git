################################################################################################################
##############################            Readme of Part1_ProcessMRSI           ################################
################################################################################################################


###### Script Summary:
This script prepares CSI data for LCModel fitting. 
It can combine multichannel data with three different methods: The 1stFIDpoint method, the sensmap method and the MUSICAL method.
The CSI data can also be multislice hadamard encoded. However, it can also process single channel or single slice data.
The CSI data can be (but don't have to be) undersampled by 2D-CAIPIRINHA Patterns. 2D-CAIPIRINHA can also be simulated.
The data can be also undersampled with slice-CAIPIRINHA, i.e. several slices can be aliased with a FoV-shift w.r.t. each other. This can also be simulated.
PreProcessing options: 2D-CAIPIRINHA, 1D-(slice)CAIPIRINHA, noise decorrelation of array coil data, coil combination, pre-phasing of data with reference data.
PostProcessing options: Hamming filtering (spatial), zero filling (spatial), Lipid decontamination (according to Bilgic's paper), exponential filter (spectral)



###### How to use this script: 

You need LCModel on a linux computer (preferable Ubuntu 12.04.3 - it's tested there). If you are not in front of this computer,
Connect to it via ssh. In a terminal, go to the the folder containing this script.
Edit the InstallProgramPaths.sh to adapt the script to your paths and needs.
Type "./MRSI_ParallelImaging.sh [options]"
You can save these commands for your datasets in a file like "Run_File_Template.sh" and then run this file via 
"bash Run_Multichannel_Combination.sh". Best practice is to create such a file with the Matlab-script "Part0_GUI/Run_script_creator.m".



The [options] are as follows:

mandatory:
-c	[csi_path]			Format: DAT, IMA or .mat (Siemens raw data or (Siemens) DICOM). Never tested with non-Siemens. The CSI file(s) that should be processed.
						You can pass over several files of the same type by '-c "[csi_path1] [csi_path2] ..."'. These files get individually processed and averaged
						at the end.
						If a .mat file is passed over, it is expected that everything is already performed like coil combination etc.
-b	[basis_paths]		Format: .BASIS. LCModel fits the metabolite spectra included in this file to your measured data. The quality of the basis is MOST crucial
						for the fitting process! See LCModel manual. If several are included, they are used for the different slices. E.g. -b "Basisfile1 Basisfile2"
						will use Basisfile1 for slice 1 of your CSI data and Basisfile2 for slice 2.
-o	[output directory]	To this directory, the output will be written.

optional: 
-i	[image NORMAL]		Format: DAT or IMA (Siemens raw data or (Siemens) DICOM). With this option you can input imaging data of the same slice & FoV as the CSI data. 
						The imaging data is used for
						- MUSICAL coil combination method, if no VC-image (option -v) is input *see Remark
						- the sensmap method, if a VC-image is input *see Remark
						- for phasing the data with the imaging data, if CSI is single-channel.
						- for creating a mask, if no mask (option -m) and T1-image (option -t) is input.  
						- the ACS data for 2D- and slice-CAIPIRINHA (or GRAPPA) reconstruction.
						The FoV must match that of the CSI file, and for MUSICAL coil combination, the echo time must be the same as that of the CSI sequence.
						For spin-echo CSI sequences, the imaging should also be a spin-echo imaging sequence, for FID-CSI-sequences, a GRE sequence.

-f	[image FLIP]		Format: DAT or IMA (Siemens raw data or (Siemens) DICOM). Same as the -i image, but with inversed imaging gradients 
						(can be achieved by rotating the FoV about -180 deg).
						This data is used to correct for gradient delays of GRE sequences. These delays cause linear phases in frequency encoding direction. 

-v	[VC image]	      	Format: DAT or IMA (Siemens raw data or (Siemens) DICOM). Same as -i image, but acquired with volume or body coil. *see Remark.
						This data is used for sensmap coil combination method or for creating a mask (if options -m, -t and -i are not used).

-t	[T1 images] 		Format: DICOM. Folder of 3d T1-weighted acquisition containing DICOM files. Used for creating mask via brain extraction tool BET2.

-m	[mask]				Defines how to create the mask. Options: -m "bet", "thresh", "voi", "[Path_to_usermade_mask]". If not set --> no mask used (i.e. everything is processed).
						bet: Use Brain Extraction Tool of FSL to create a mask out of the magnitude or T1-image.
						thresh: Threshold magnitude or T1. This is useful for phantoms, where BET2 doesnt work. All values above a calculated threshold are set to 1, all below to 0.
						voi: Use VoI to create mask.
						[Path]: Use the minc file given by the path as a mask. 

-h 	[0-100]          	Tells the program to Hamming filter the CSI data and how much in percentage. 100 means "normal" Hamming filter. This is the same as you can set on the 
						Siemens	console.

-r	[MatlabCommand]		Tells the program to simulate 2D-CAIPIRINHA, or if the program cannot automatically find out the 2D-Pattern from the CSI data, how the data was actually 
						undersampled. The string is passed over to Matlab and simply run there. So be sure to input a valid Matlab command.
						Example: -r "InPlaneCaipPattern = [0 0 0; 0 0 0; 0 0 1]; VD_Radius = 2;".

-R	[MatlabCommand]		Same as -r, just for slice-CAIPIRINHA.
						Example: -R "SliceAliasingPattern = [1 4; 2 5;3 6]; FoVShifts_x = [0 0 0 0.5 0.4 0.3]; FoVShifts_y = [0 0 0 0.5 0.4 0.3];"
						This would mean: This would mean: Slice 1 and 4, 2 and 5, and 3 and 6 are aliased,\nwith Slices 1,2,3 not shifted, 
						and Slices 4,5,6 shifted by 0.5 x FoV, 0.4 x FoV and 0.3 x FoV in x- and y-direction.

-g 	[noisedecorr_file]	If this option is used the csi data gets noise decorrelated using noise from passed-over noise file (i.e. a file containing only the noise of all channels, 
						e.g. measured by FID-sequence without gradients and excitation pulses), or if no noisedecorrelation_file is given, by noise from the end of the FIDs at the 
						border of the FoV. In the latter case, just use -g ""

-u	[Nothing]           The phantom_flag, used if a phantom was measured. Different settings used for fitting (e.g. some metabolites are omitted). This option is not really useful,
						we could not really make phantoms work being fitted correctly by LCModel.

-n	[Nothing]           Perform zero-filling to the next power of 2 in ROW and COL dimensions for the CSI data (e.g. zerofill from 42x42 to 64x64)

-l	[Nothing]           If this option is set, LCModel is not started, everything else is done normally. Useful for only computing the SNR.

-j	[LCM_ControlFile]   ControlFile telling LCModel how to process the data. If you want to change the way LCModel processes data, change this file, otherwise standard values are 							assumed. A template file is provided in the package.

-e	[ExpFilterInHz]		Apply an exponential filter to the spectra in Hz. A template file is provided in this package.

-w	[Water Reference]	Format: DAT or DICOM. LCModel 'Do Water Scaling' (W1) or separate water quantification (W2) (Water maps are created). The same scan as -c [csi file], but without water 				suppression.

-W	[Water_Ref Control File] If W2 was used for water quantification a separate water reference control file needs to be input. 	


You have to input the mandatory options.



* Remark:
Which Coil Combination is used for which input?

Your input options determine the coil combination method as follows:
-i	-v	CSI is multichannel	-->	sensmap method
-i		CSI is multichannel	-->	"our" (MUSICAL) method
		CSI is multichannel	-->	1stFIDpoint method
-i		CSI is singlechannel-->	Imaging data is used only for phasing CSI data

The -f option has no influence on this scheme.

###### Water reference processing
Two possible approaches for water reference scan are implemented.

1. Default LCModel version (Do Water Scaling), where a water unsuppressed dataset measured with same parameters as the water suppressed csi is used for referencing. Only one fitting takes place and LCModel handles the water by itself (a blackbox for user, we do not see the fitted water spectra). 
Example of flag:
-w "W1, Path/to/WatRef/file" \

2. A separate water reference file fitting, using a 'special type flag' within LC Model Control file. The water spectra are fitted using ControlWrite.SPTYPE = 'SPTYPE = ''lipid-8'''; which needs to be in the control file. This allows user to see the output of fit, check the water distribution map and apply own corrections. Fitted water data are under /water_spectra folder, metabolite fits are under /spectra as usual. Water distribution maps are stored in /maps together with metabolite maps.
Example of flag:
-w "W1, Path/to/WatRef/file" \
-W Path/to/WaterRef/ControlFile


###### Workflow of Script

1. Create tmp_dir
2. Create logfile and split output to logfile and screen using tee.
3. Find out options.
4. "Install" program, i.e. set paths to LCModel etc.
5. Create output directories.
6. Write initial parameters, i.e. the options passed over to script.
7. Read in the parameters of the CSI data from the headers.
8. Create mask.
9. Start matlab script for reconstruction.
9.1. Read noise, csi and imaging data. Perform Noise decorrelation while reading in data.
9.2. Simulate noise for SNR-computation.
9.3. Perform Parallel Imaging reconstruction.
9.4. Write phase maps of all channels.
9.5. Perform coil combination.
9.6. Perform lipid decontamination according to Bilgic.
9.7. Apply Hamming and Exponential filter.
9.8. Perform same stuff on noise data for SNR computation.
9.9. Write files necessary for LCModel fitting.
9.10. Save processed data in .mat file.
10. Start LCModel processing.
11. Clean up files, copy logfile, copy used sourcecode. 
12. Stop script.



--- Bernhard Strasser, August 30 2012
--- revised: Bernhard Strasser, March 27, 2015
--- revised: Michal Povazan, February 10, 2016
--- revised: Bernhard Strasser, March 4, 2016
  
