;+
; NAME:
;   IMAGE_STDDEV
;
; PURPOSE:
;   This function calculates the local-neighbourhood statistical standard deviation.
;
; CATEGORY:
;   Image Processing
;
; CALLING SEQUENCE:
;    Result = IMAGE_STDDEV(Image, HalfWidth)
;
; DESCRIPTION:
;    Based on the formula for variance:
;      var = (sum of the squares)/n - (square of the sums)/n*n
;   For each array element the variance of the neighbourhood of +- 
;   halfwidth is calculated. The routine avoids any loops and so is fast 
;   and "should" work for any dimension of array.
;
; INPUT VARIABLES:
;   image     : The array of which we calculate the variance. Can be any dimension.
;   halfWidth : The half width of the NEIGHBOURHOOD, indicates we are
;               looking at a neigborhood +/- N from the pixel in each dimension.
;
; KEYWORD PARAMETERS:
;   neighbourhood       : Calculate for the NEIGHBOURHOOD only, not the central pixel.
;   population_estimate : Returns the population estimate of stddev, not the sample stddev
;
; OUTPUT:
;    Returns an array of same dimensions as input array in which each pixel
;    represents the local standard deviation centred at that position.
;
; OPTIONAL OUTPUTS:
;    MEAN_IM: Set to array of local area mean, same dimensionality as input.
;
; RESTRICTIONS:
;    Edges are dealt with by replicating border pixels this is likely to
;    give an underestimate of standard deviation in these regions
;
; EXAMPLE:
;   Example of simple statistical-based filter for removing spike-noise
;     stddev_im =  image_stddev(image,  5, mean=mean_im, /neigh)
;     zim = (image-mean_im)/stddev
;     ids = where(zim gt 3, count)
;     if count gt 0 then image[ids] = mean_im[ids]
;
; MODIFICATION HISTORY:
;   Version 1.0: Olivier Absil (University of Liège), based on the routine IMAGE_VARIANCE by Martin Downing 
;-

FUNCTION IMAGE_STDDEV, image, halfWidth, MEAN_IM=av_im, NEIGHBOURHOOD=NEIGHBOURHOOD, POPULATION_ESTIMATE=POPULATION_ESTIMATE

; full mask size as accepted by SMOOTH()
n = halfWidth*2+1

; sample size
m = n^2
    
; temporary double image copy to prevent overflow
im = DOUBLE(image)
 
; calc average
av_im = SMOOTH(im, n, /EDGE_TRUNCATE)
 
; calc squares image
sq_im = TEMPORARY(im)^2
 
; average squares
asq_im = SMOOTH(sq_im, n, /EDGE_TRUNCATE)

IF KEYWORD_SET(NEIGHBOURHOOD) THEN BEGIN ; remove centre pixel from estimate
  ; calc neighbourhood average (removing centre pixel)
  av_im = (av_im*m - image)/(m-1)
  ; calc neighbourhood average of squares (removing centre pixel)
  asq_im = (asq_im*m - TEMPORARY(sq_im))/(m-1)
  ; adjust sample size
  m = m-1
ENDIF

var_im = TEMPORARY(asq_im) - (av_im^2)
IF KEYWORD_SET(POPULATION_ESTIMATE) THEN var_im = var_im * (DOUBLE(m)/(m-1))

RETURN, SQRT(var_im)

END