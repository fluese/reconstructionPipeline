# MR Reconstruction Pipeline
This is a reconstruction pipeline to reconstruct raw MRI data acquired with a Siemens MR scanner using gradient recalled echo (GRE) sequences. It is built to process raw data which is larger than memory, however, data can be processed entirely in memory and/or parallelized across as many physical CPUs available using MATLAB's paralllel toolbox.

Processing time varies a lot depending on what data you are going to process and what features you enable. It can be roughly a minute for a non-GRAPPA 1mm MPRAGE with SoS to a day for a 250 µm combined with adaptive combine.

## Disclaimer
**Please respect the copywrite of included tools and files! Thanks!**

**External tools used within this reconstruction pipeline**
* ./externalTools/orientatition/ 			-> FreeSurfer (https://surfer.nmr.mgh.harvard.edu/)
* ./externalTools/nifti_tools/				-> MathWorks (https://de.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)
* ./externalTools/mapVBVD/ 					-> Philip Ehses, German Center for Neurodegenerative Diseases, Germany, site Bonn (Personal communication)
* ./externalTools/GRAPPA/ 					-> Miki Lustig (https://people.eecs.berkeley.edu/~mlustig/Software.html)
* ./externalTools/denoise/Coupe/			-> Pierrick Coupé (https://sites.google.com/site/pierrickcoupe/softwares/denoising-for-medical-imaging/mri-denoising)
* ./externalTools/denoise/BM4D/ 			-> Foi & Magioni (http://www.cs.tut.fi/~foi/GCF-BM3D/index.html#ref_software)
* ./externalTools/denoise/awesome/			-> André Pampel, Max Planck Institute for Cognitive and Brain Sciences, Germany, Leipzig (Personal communication)
* ./externalTools/phaseUnwrap/				-> André Pampel, Max Planck Institute for Cognitive and Brain Sciences, Germany, Leipzig (Personal communication)
* ./processingPipeline/adaptiveCombine.m 	-> University of Wuerzburg, Germany


If you are the author of any software used in this pipeline and you are unhappy with either it being included or the way your work is acknowlegded, please get in contact with me using the email below. 

## Instructions
Run 'installReconstructionPipeline.m' first. Choose a folder where data shall be saved to and where SPM12 is located (if installed). Then run 'startReconstruction.m' to reconstruct Siemens raw data.
 
Changes to the default processing can be done either by loading a different preset or by changing the settings manually. Presets are stored in ./presets/ and to manually change settings open 'setOptions.m' and 'setParameters.m', respectively. New presets can be generated using 'createPreset.m' using the currently specified options and parameters in 'setOptions.m' and 'setParameters.m'.

To process a noise dataset (e.g. for decorrelation of channels), use the preset 'noise'. This will change the naming scheme accordingly and run the processing pipeline for as long as needed only.

If you run into errors or anything is unclear, feel free to get in contact with me using the address below. Any kind of feedback is highly appreciated.

## Flow chart
*Reconstruction pipeline*
1. Read in raw data (using mapVBVD)
1. Decorrelate channels (using noise data)
1. 2D GRAPPA reconstruction (using Miki Lustig's code)
1. Denoise data (using BM4D, different NLM, or DnCNN)
1. Gibbs ringing reduction (using 3D Tukey filter)
1. Partial Fourier reconstruction (using zero-padding)
1. Phase unwrapping (using total variation)
1. Combine channels (using root sum-of-squares or adaptive combine)
1. Write NIfTI files (magnitude, phase, or complex)

*Optional additional functions (requires SPM12)*
1. Bias field correction
1. Segmentation
1. Creation of brain mask
1. Distortion correction (not included in public release!)
 
************************************************************************
## Functionalities
* Reconstruction of raw data files that are larger than memory by writing/reading intermediate results to/from disk.
* Reconstruction of data in parallel by splitting data into chunks. Depending on the processing stage data is split across channels or slices. Number of CPUs and chunks to be processed can be user specified or detected automatically.
* Channel decorrelation by estimating noise covariance matrix (requires noise dataset from the same session).
* Channel combination by sum of squares and adaptive combination.
* Channel combination within adaptive combine slice by slice or in 3D.
* Channel compression within adaptive combine.
* Writing of magnitude, phase, and complex valued NIfTI files.
* 2D GRAPPA reconstruction.
* Application of 3D Tukey filter in k space to reduce Gibbs ringing artifact.
* Application of 3D distortion correction (requires gradient coeffient file of Siemens of your specific gradient system).
* Application of 2D phase unwrapping (per slice).
* Application of denoising in complex domain during reconstruction (includes denoising by BM4D, all NLM denoiser's of Coupé et al. as well as using MATLAB's implementation of DnCNNs utilizing the neural network toolbox).
* Application of bias field correction, segmentation, and creation of a brainmask.

## Known bugs
* 3D phase unwrapping and wavelet denoising are not working properly.
* 2D phase unwrapping results in odd results using 250 um data.
* If distortion correction is used with gzipped NIfTIs under windows an error occurs (probably use different tool to load data?).
* installReconstructionPipeline has to be run from root directory of the pipeline.
* Denoising with Coupé denoisers do not work currently.

## Outlook
* Change the way parameters and options are set. To be able to load presets more conveniently. 
* Have a look at 'initializeParallel.m' and clean it up. Get nChunks and stepSize, calculate slice and channel memory even if no autosetup is used.        
* Allow to specify files in case of -nodisplay (currently needs to be written in the according m-files)
* Rework processing status bar in case of -nodisplay. 
* Change maxMemory parameters relative maximum available memory instead of absolute.
* Sanity checks if ~autoSetup for setting up workers, stepSize and nChunks.
* Save output in command window along settings file (using 'diary' ?).
* Implement removal of slice (and phase) oversampling. 
* Display what functions are run in preproc, proc and combine.
* Processing of noise data should be integrated into the processing pipeline without the need of a special preset.
* Write data with correct origin and orientation (currently requires a reference dataset for the correct origin)
* Make use of 'onCleanUp' to close parallel pool if ctrl+c pressed.
* Output data BIDS formated.
* Make use of BART
* Read in ISMRM rawdata format (and H5?).
* Implement NUFFT for radial trajectories.
* Implement SENSE reconstruction.
* Implement reconstruction in units of SNR for quality checks.
* Synthesize missing k space data by POCS or homodyne. Code is available, just needs to be reintegrated into the pipeline.
* Auto-adjust filter strength of Tukey filter depending on resolution.
* Create a simple GUI.
* Implement QA (e.g., mriqc.readthedocs.io/en/latest/iqsm/t1w.html)
* Use niftiread and niftiwrite to read and write NIfTI files (introduced in MATLAB 2017b).
* Move from if then else formalism to switches to make the code cleaner.

## Version
Version 1.0 (05.09.2019)
## Contact
Falk Luesebrink (falk dot luesebrink at med dot ovgu dot de)
