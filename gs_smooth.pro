;================================================================
function gs_smooth,n,sigma
;
; prepares a n-element vector of total area=1 and Gaussian shape
; given by sigma, shifted to the origin [bbb.....bbb]
; for Gaussian filter using new=cnv(old,gs_smooth)
;
; when a given FWHM needed, use relation: FWHM = 2.354 sigma
;                                     or: sigma = 0.4248 FWHM
   b = findgen(n)
   c = b
   c = ((b - b(n/2))/sigma)^2/2
; bug fixed Jan.2,91, before not divided by 2
   b(*) = 0.
   d = where(c le 30.) 
   b(d) = exp(-c(d))
   b = shift(b, -n/2)
return, b/total(b)
end
