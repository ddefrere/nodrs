;+
; NAME: GET_DRS
; 
; PURPOSE:
;   Initializes reduction paramaters (based on output config files)
;
; OUTPUTS:
;   drs, a structure containing the following reduction paramaters:
;     - ADI_MODE   :  Define how the images are processed
;                       - 0: Regular ADI. Default mode.
;                       - 1: PCA.
;     - APER_RAD   :  Radius of the photometric aperture (in pixels)
;     - BCK_IRAD   :  Inner radius of the background ring (in pixels)
;     - BCK_ORAD   :  Outer radius of the background ring
;                  :     - -2: use the whole channel beyond bck_irad from the beam;
;                  :     - -1: force the annulus to have the same number of pixels as the photometric aperture;
;                  :     -  0: use the closest channel edge;
;                  :     -  n: in pixels, the outer radius of the sky annulus (n is a positive integer)
;     - BANDWIDTH  :  [m], central wavelength (superseed filter keywords in the input FITS files)
;     - BCK_CEN    :  Two-element vector with the center position of the background region (relative to XCEN, YCEN)
;     - BCK_LIM    :  Acceptable background range (after nod subtraction)
;     - BCKG_MODE  :  Define background subtraction approach. The background subtraction is done by subtracting to each frame of a given nod a single frame (different each time)
;                     of the corresponding background nod. Use a negative bckg_mode to subtract to each frame the median of the frames of the corresponding background nod.
;                       - 0: No background subtraction
;                       - 1: Background subtraction performed by nod pairs.
;                       - 2: Background subtraction performed using the n closest frames in two adjacent nods, where n is set by N_FRBCK
;                       - 3: Background subtraction performed using the dedicated background frames in the same pointing (i.e., those flagged with an obstype of 3).
;                       - 4: For chopping, not yet implemented
;     - BCKG_SEL   :  Background selection mode when N_FRBCK is not null (0: time, 1: elevation, 2: central value)
;     - BFL_MODE   :  Background floor mode for the background annulus (0: sigma-clipped mean, 1: median). 
;     - BIAS_ESTIM :  Compute the flux in a nearby empty region of the detector
;     - CAL_METHOD :  Interpolation method for the transfer function
;                     - set keyword to 0 for linear interpolation between nearest neighbour (default)
;                     - set keyword to 1 for polynomial interpolation on all calibrator data
;     - CAL_MODE   :  0: calibrate per pointing, 1: calibrate per OB
;     - CHI2_LIM   :  Maximum acceptable chi2 used for null calibration
;     - FLX_MODE   :  Define the method use for flux computation
;                       - 0: aperture photometry (default)
;                       - 1: weighted aperture photometry (only for objects unresolved with one telescope)
;                       - 2: PSF-fitting photometry (not well suited to null data)
;     - FIT_MODE   :  Define the centroid fitting technique (set it negative to fit inverse function, e.g. for AGPM data):
;                       - 0: no stellar centroid fitting
;                       - 1: compute the centroid of a star using a derivative search (using CNTRD function)
;                       - 2: compute the stellar centroid by Gaussian fits to marginal X,Y, sums (using GCNTRD function)
;                       - 3: compute the stellar centroid by Gaussian fit using MPFIT2DPEAK (more time consuming)
;                       - 4: compute the stellar centroid by Lorentzian fit using MPFIT2DPEAK
;                       - 5: compute the stellar centroid by Moffat model fit using MPFIT2DPEAK
;     - IMG_MODE   :  Image combination mode (0: median, 1: mean, 2: resistant mean)
;     - KEEP_RATIO :  Percentage of frames to preserve when computing the weighted average of the null (only used if NULL_MODE = 1) or the visibility (keeping
;                     only frames with the highest corrrelation)
;     - LAMBDA_CEN :  [m], central wavelength (superseed filter keywords in the input FITS files)
;     - LMIRCAM    :  Set if LMIRCam
;     - MIN_FR     :  Minimum number of frames per OB required to run the null computation
;     - MIRAC      :  Set to 1 if MIRAC
;     - N_BTSTRP   :  Number of bootstrap samples to compute the error bar (Bayesian approach used if 0 or 1)
;     - N_BIN      :  Image re-binning size (in pixels)
;     - N_CLIP     :  Size of the output images (in pixels, around beam position)
;     - N_CLIP     :  Size of the output images (in pixels, around beam position)
;     - N_COADD    :  Number of images to coadd before frame derotation
;     - N_FRBCK    :  Maximum Number of frames to preserve in the adjacent nod for background subtraction (see also BCKG_SEL)
;                       - -1: same number of frames as in the current nod
;                       - 0 : all the frames in the adjacent nods
;                       - n : user-defined number of frames in the adjacent nods (must be > 0)
;     - N_FROB     :  Number of frames per OB (superseed OB_MODE)
;     - N_TRANS    :  Number of "bad" transition frames at the beginning of each acquisition (0 by default)
;     - NO_BPM     :  Turn on/off bad pixel correction
;     - NO_DARK    :  Turn on/off dark subraction
;     - NO_FLAT    :  Turn on/off flat fielding
;     - NO_FIND    :  Turn on/off quick beam finding routine (during L0 image calibration)
;     - NOD_COR    :  Set to 1 if nulls appear correlated within a given pointing (but not necessarily from pointing to pointing, e.g. Feb. 2015 beta Leo sequence)
;     - NOMIC      :  If NOMIC
;     - NSC_BFAC   :  Mulitplier on number of histogram bins (nbins_fac*sqrt(Np))
;     - NSC_BINS   :  Bin size for NSC reduction (0: constant, 1: variable)
;     - NSC_CUBE   :  3-element vector containing the number of elements in the chi2 cube along the null, mean phase and rms phase directions respectively (40, 35, 30 if set to 0)
;     - NSC_MODE   :  0: best-fit values, 1: Bayesian results, 2: MEAN of bootstrap samples, 3: MODE of bootstrap samples
;     - NSC_OMIN   :  Minimum number of occurences per bin for the NSC fit
;     - NULL_COR   :  Subtract the null estimated from high-frequency phase noise to each raw null measurements
;     - NULL_LIM   :  Two element vector with the lower and upper limit of the raw null values to keep (e.g., [-0.10,0.10])
;     - NULL_MODE  :  Define the technique used to compute the null
;                       - 0: Default mode. Used the mode of the null distribituion.
;                       - 1: Computed the weighted average of the best NULL_RATIO frames
;                       - 2: Statistical reduction (not yer connected but working separately with the L1 files)
;     - NULL_RANGE ;  Source null range (to scan with NSC) 
;     - NULL_RAD   :  Photometric aperture radius used for null computation (closest to value defined for flux computation, in pixels)
;     - PARA_BIN   ;  Parallactic angle used to bin the images before PCA/ADI processing
;     - PARA_RANGE ;  Range of parallactic angle to preserve before PCA/ADI processing
;     - OB_MODE    :  Define how the frames are grouped for null computation
;                       - 0: Default mode. An OB is defined as all consecutive null files of the same nod position.
;                       - 1: An OB is defined as all consecutive null files of the same pointing (OBSOLETE)
;                       - 2: An OB is defined for each chop/nod sequence (which occurs after 2 nods, not yet implemented)
;                       - 3: An OB is defined as all consecutive null files of the same chop position (not yet implemented)
;     - OFFSET     :  Offset between [XCEN,YCEN] and the center of the image (useful to aligned on a ghost for instance, in pixels)
;     - OVERLAP    :  Force beam overlap (used only with OBSTYPE = 0)
;     - PCA_DIM    :  Dimensions of the ouptut image, on which (s)PCA will be performed [pixels] (0 for full frame)
;     - PCA_KLIP   :  Number of principal components to be used in the (s)PC
;     - PCA_SMART  :  Turn on/off smart PCA
;     - POLYDEG    :  Degree of the polynomial for CAL_METOHD = 1
;     - PRE_CROP   :  Pre-crop the frames at the image reduction stage (4-element vector: x_min, y_min, x_max, y_max)
;     - PRECISION  :  Set to 1 to perform every image operation in double precision (float by default)
;     - R_KEEP     ;  Ratio of frames to preserve (based on correlation coefficient with PSF file)
;     - SIG_FLX    :  Keep only frames where the flux is lower than the mean flux minus + sig_flx*sigma
;     - SIG_POS    :  Keep only frames where the beam position is within sig_pos sigma from the mean beam position
;     - SIG_SLO    :  Keep only frames where the fitted Moffat slope is within sig_slo sigma from the mean beam position
;     - SIG_VIS    :  Keep only frames where the visibility is larger than the mean flux + sig_vis*sigma
;     - SKIP_ADI   :  To turn off ADI processing
;     - SKIP_CAL   :  To turn off image calibration
;     - SKIP_FLX   :  To turn off flux computation
;     - SKIP_NULL  :  To turn off null computation
;     - SKIP_OPEN  ;  Don't read open-loop frames (when appropriate). Not always a good idea (e.g., used for background subtraction)
;     - SKIP_REG   :  Don't register the frames (beam center computed only for tip/tilt information)
;     - SKIP_VIS   :  To turn off visibility computation
;     - SKY_WEIGHT :  Set to 1 in order to weight the number of pixels per column in the background region according to the corresponding number of pixels in the photometric aperture. 
;     - SPLIT_DITH ; Set to 1 to split the OB if different dither patterns are found
;     - SPLIT_NOD  :  Calibrate each nod position separately
;     - SPLIT_TF   :  Set to 1 to split the TF when there is a dead time longer than time_split between two null measurements
;     - SPLIT_HOUR :  UT hour(s) at which the TF will be split
;     - SPLIT_TIME :  Maximum dead time between two null measurements for a single TF (hours)
;     - XCEN       :  Two-element vector with the X coordinate of the star on the detector [pixels] (one element for each side -- superseed the automatic computation)
;     - YCEN       :  Two-element vector with the Y coordinate of the star on the detector [pixels] (one element for each side -- superseed the automatic computation)
;     - X_RANGE    :  Range of X detector positions to preserve during ADI frame selection (in pixels) 
;     - Y_RANGE    :  Range of Y detector positions to preserve during ADI frame selection (in pixels) 
;
; MODIFICATION HISTORY:
;   Version 1.0,  16-FEB-2015, by Denis Defrère, Steward Observatory (ddefrere@email.arizona.edu)
;   Version 1.1,  12-MAR-2015, DD: format output to be compliant with NexSci archiv
;   Version 1.2,  12-DEC-2015, DD: improved comments to reflect better the config files
;   Version 1.3,  28-JUL-2017, DD: added BFL_MODE

PRO GET_DRS, drs, cfg_file

; Read file
drs = READ_CONFIG(cfg_file)

; Append instrument name
struct_add_field, drs, 'instrum', 'tmp'
IF drs.lmircam EQ 1 THEN drs.instrum = 'lmircam'
IF drs.mirac   EQ 1 THEN drs.instrum = 'mirac'
IF drs.nomic   EQ 1 THEN drs.instrum = 'nomic'

; Backward compatibility
IF NOT TAG_EXIST(drs, 'max_time') THEN struct_add_field, drs, 'max_time', 10
IF NOT TAG_EXIST(drs, 'bfl_mode') THEN struct_add_field, drs, 'bfl_mode', 0

END
