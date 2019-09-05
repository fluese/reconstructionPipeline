-------------------------------------------------------------------

  BM4D software for volumetric data denoising and reconstruction
            Public release ver. 3.2  (30 March 2015)

-------------------------------------------------------------------

Copyright (c) 2010-2015 Tampere University of Technology. 
All rights reserved.
This work should be used for nonprofit purposes only.

Authors:                     Matteo Maggioni
                             Alessandro Foi


BM4D web page:               http://www.cs.tut.fi/~foi/GCF-BM3D


-------------------------------------------------------------------
 Contents
-------------------------------------------------------------------

The package contains these files

*) demo_denoising.m         : denoising demo script
*) demo_reconstruction.m    : reconstruction demo script
*) bm4d.m                   : BM4D volumetric denoising filter [1]
*) helper.m                 : various methods used by the demos
*) ssim_index3d.m           : 3-D SSIM index [4,5]
*) SheppLogan3D.mat         : 3-D Shepp-Logan phantom
*) Transforms.mat           : Default Wavelet transforms [1]
*) t1_icbm_normal_1mm_pn0_rf0.rawb   : BrainWeb T1 phantom [3]

-------------------------------------------------------------------
 Installation & Usage
-------------------------------------------------------------------

Unzip BM4D.zip (contains codes) in a folder that is in the MATLAB 
path. Execute the script "demo_reconstruction.m" to run the
reconstruction demo, or execute the script "demo_denoising.m" to
run a volumetric denoising demo. You can freely modify the 
parameters involved in the filtering at the beginning of each demo.

-------------------------------------------------------------------
 Requirements
-------------------------------------------------------------------

*) MS Windows 64 bit, Linux 64 bit or Mac OS X 64 bit
*) Matlab R2011b or later with installed:
   -- Image Processing Toolbox (only for visualization with "imshow")
   -- Signal Processing Toolbox (only for non-default transforms in BM4D)
   -- Wavelet Toolbox (only for non-default transforms in BM4D)

-------------------------------------------------------------------
 Change log
-------------------------------------------------------------------
v3.2   (30 March 2015)
 ! fixed bug in Rician noise estimation

v3.1.1 (11 December 2014)
 + different transforms can be used within the same cube
 ! fixed bug arising when filtering 2D data
 ! fixed bug in 3rd dimensional inverse transformation

v3.1   (10 November 2014)
 . code optimization to speed up filtering

v3.0   (5 November 2014)
 + introduced low complexity profile in BM4D
 + improved interface of BM4D function
 ! fixed bug in Wiener filtering under Rician noise and unknown sigma
 . removed dependencies with VST package

v2.4   (24 October 2014)
 . faster rician denoising with noise estimation

v2.3   (5 March 2014)
 + handled case of estimated standard deviation sigma=0
 ! minor bug fixes in bm4d function

v2.2.1 (2 March 2014)
 ! introduced error in case of sigma<=0 in demo_denoising

v2.2   (20 September 2013)
 + improved demo_denoising script for Rician spatially varying noise
 . default wavelet transforms do not longer require the wavelet toolbox
 . optimized memory usage

v2.1   (25 July 2013)
 + parametrized thresholding type (hard or soft)
 ! volumetric inputs with depth lower than the depth of the cubes are
   correctly handled, the code scales nicely also for the particular
   case of 2-D inputs

v2.0   (17 April 2012)
 + reconstruction of volumetric phantom data with non-zero phase 
   from noisy and incomplete Fourier-domain (k-space) measurements
 + adaptive denoising for data corrupted by spatially varying noise [2]

v1.0.1 (18 July 2011)
 ! fixed few typos, corrected lambda_thr4D in modified profile

v1.0   (17 July 2011)
 + initial version

-------------------------------------------------------------------
 References
-------------------------------------------------------------------

[1] M. Maggioni, V. Katkovnik, K. Egiazarian, A. Foi, "A Nonlocal 
    Transform-Domain Filter for Volumetric Data Denoising and 
    Reconstruction", IEEE Trans. Image Process., vol. 22, no. 1,
	pp. 119-133, Jan. 2013. doi:10.1109/TIP.2012.2210725

[2] M. Maggioni, A. Foi, "Nonlocal Transform-Domain Denoising of 
    Volumetric Data With Groupwise Adaptive Variance Estimation", 
    Proc. SPIE Electronic Imaging 2012, San Francisco, CA, USA, Jan. 2012

[3] R. Vincent, "Brainweb:  Simulated  brain  database", online at
    http://mouldy.bic.mni.mcgill.ca/brainweb/, 2006.

[4] Z. Wang, A. Bovik, H. Sheikh, E. Simoncelli, "Image quality 
    assessment: from error visibility to structural similarity",
    IEEE Trans. Image Process., vol. 13, no. 4, pp. 600-612, April 2004.

[5] J. V. Manjon, P. Coupe, A. Buades, D. L. Collins, M. Robles, 
    "New methods for MRI denoising based on sparseness and self-similarity", 
    Medical Image Analysis, vol. 16, no. 1, pp. 18-27, January 2012
 
-------------------------------------------------------------------
 Disclaimer
-------------------------------------------------------------------

Any unauthorized use of these routines for industrial or profit-
oriented activities is expressively prohibited. By downloading 
and/or using any of these files, you implicitly agree to all the 
terms of the TUT limited license, as specified in the document
Legal_Notice.txt (included in this package) and online at
http://www.cs.tut.fi/~foi/GCF-BM3D/legal_notice.html

-------------------------------------------------------------------
 Feedback
-------------------------------------------------------------------

If you have any comment, suggestion, or question, please do
contact    Matteo Maggioni   at  matteo.maggioni<at>tut.fi


