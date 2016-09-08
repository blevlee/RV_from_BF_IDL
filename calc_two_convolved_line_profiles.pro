function calc_two_convolved_line_profiles, x, p
;1.  Integrity check
;1.1  Check for even velocity spacing in the input x array
x_spacing=x[1:n_elements(x)-1]-x[0:n_elements(x)-2]
x_spacing_median=median(x_spacing)
if x_spacing_median lt 0 then begin
	    print,systime(/utc),'|ERROR|make_bf_gaussfits_temp.two_convolved_line_profiles|The x array should be sorted in increasing order, not decreasing order'
    stop
endif
x_spacing_uniformity_threshold=1.0d-5*x_spacing_median
index_differing_spacings=where(x_spacing lt x_spacing_median-x_spacing_uniformity_threshold or x_spacing gt x_spacing_median+x_spacing_uniformity_threshold,count_differing_spacings)
if count_differing_spacings gt 0 then begin
    print,systime(/utc),'|ERROR|make_bf_gaussfits_temp.two_convolved_line_profiles|The spacing of the x array was not uniform'
    stop
endif
;1.  Set up variables
npoints=n_elements(x)
equatorial_max_rv1=p[0]
limbdarkening_beta1=p[1]
velocity_peak_offset1=p[2]
area_under_profile1=p[3]
equatorial_max_rv2=p[4]
limbdarkening_beta2=p[5]
velocity_peak_offset2=p[6]
area_under_profile2=p[7]*p[3]
instrumental_profile_sigma=p[8]
baseline_offset=p[9]
;vel_halfdomain=(max(x)-min(x))/2.0
;x_midpoint=(max(x)+min(x))/2.0
;;1.2  Set up high-resolution velocity grid for computations.  We will
;;downsample back to the original resolution at the end of this
;;function.
upsample_factor_halfwidth=2L
upsample_factor=upsample_factor_halfwidth*2+1
npoints_hires=upsample_factor*npoints
velocity_hires=dblarr(npoints_hires)
for i=0L,n_elements(velocity_hires)-1 do begin
    lowres_index=i/upsample_factor
    hires_subindex=(i mod upsample_factor)-upsample_factor_halfwidth
    velocity_lowres=x[lowres_index]
    velocity_hires[i]=velocity_lowres+(hires_subindex/double(upsample_factor))*x_spacing_median
    ;print,i,lowres_index,hires_subindex,x[lowres_index],velocity_hires[i]
endfor
vel_halfdomain_hires=(max(velocity_hires)-min(velocity_hires))/2.0
x_midpoint_hires=(max(velocity_hires)+min(velocity_hires))/2.0

;2.  Make star 1's line
;2.1  Convolve delta function with various line broadening effects in
;the star

;;2.1.1  Convolve with rotational broadening
;;!!!I think we need to try keeping the computation at high-resolution
;;until downsampling to the instrumental resolution.  I am unable to
;;fit sharply to the data peak when trying to do all of the steps at
;;final sampling.
;status=calc_rotational_broadening_kernel( $
;                                          n_points_in_kernel=npoints, $
;                                          equatorial_max_rv=equatorial_max_rv1, $
;                                          limbdarkening_beta=limbdarkening_beta1, $
;                                          vel_halfdomain=vel_halfdomain, $ ;km/s
;                                          x_midpoint=x_midpoint, $
;                                          velocity=velocity, $
;                                          normalized_kernel=model_profile $
;                                        )
status=calc_rotational_broadening_kernel( $
                                          n_points_in_kernel=npoints_hires, $
                                          equatorial_max_rv=equatorial_max_rv1, $
                                          limbdarkening_beta=limbdarkening_beta1, $
                                          vel_halfdomain=vel_halfdomain_hires, $ ;km/s
                                          x_midpoint=x_midpoint_hires, $
                                          velocity=velocity_hires_out, $
                                          normalized_kernel=model_profile_hires $
                                        )
;;2.2.  Convolve astrophysical line profile with instrumental broadening
;instrumental_profile=gs_smooth(npoints_hires,instrumental_profile_sigma)
;model_profile=cnv(normalized_kernel,instrumental_profile)
;;2.3.  Apply a velocity shift to the line
shift=velocity_peak_offset1-x_midpoint_hires
model_profile_shifted=interpol(model_profile_hires,velocity_hires_out,velocity_hires_out-shift)
;2.4.  Scale the profile to have the desired area
area=int_tabulated(velocity_hires_out,model_profile_hires,/double)
model_profile_shifted1=model_profile_shifted*area_under_profile1/area

;model_profile_shifted=one_gaussian(x,[0.0,velocity_peak_offset1,equatorial_max_rv1,1.0])
;area=int_tabulated(x,model_profile_shifted,/double)
;model_profile_shifted1=model_profile_shifted*area_under_profile1/area

;oplot,velocity,model_profile_shifted,color=-255

;3.  Make star 2's line
;3.1  Convolve delta function with various line broadening effects in
;the star
;3.1.1  Convolve with rotational broadening
status=calc_rotational_broadening_kernel( $
                                          n_points_in_kernel=npoints_hires, $
                                          equatorial_max_rv=equatorial_max_rv2, $
                                          limbdarkening_beta=limbdarkening_beta2, $
                                          vel_halfdomain=vel_halfdomain_hires, $ ;km/s
                                          x_midpoint=x_midpoint_hires, $
                                          velocity=velocity_hires_out2, $
                                          normalized_kernel=model_profile_hires $
                                        )
;;3.2  Convolve astrophysical line profile with instrumental broadening
;instrumental_profile=gs_smooth(npoints,instrumental_profile_sigma)
;model_profile=cnv(normalized_kernel,instrumental_profile)
;3.3  Apply a velocity shift to the line
shift=velocity_peak_offset2-x_midpoint_hires
model_profile_shifted=interpol(model_profile_hires,velocity_hires_out,velocity_hires_out-shift)
;3.4  Scale the profile to have the desired area
area=int_tabulated(velocity_hires_out,model_profile_hires,/double)
model_profile_shifted2=model_profile_shifted*area_under_profile2/area
;4.  Add star 1 to star 2
model_profile_total_hires=model_profile_shifted1+model_profile_shifted2

;plot,velocity_hires_out,model_profile_shifted1
;oplot,velocity_hires_out,model_profile_shifted1,color=-255
;oplot,velocity_hires_out,model_profile_shifted2,color=-255
;oplot,velocity_hires_out,model_profile_total_hires,color=255*256L

;;5.  Convolve astrophysical line profile with instrumental broadening
instrumental_profile=gs_smooth(npoints_hires,instrumental_profile_sigma*upsample_factor)
model_profile_total_hires=cnv(model_profile_total_hires,instrumental_profile)

;oplot,velocity_hires_out,model_profile_total_hires,color=255
;;;plot,velocity,normalized_kernel
;;;oplot,velocity,model_profile,color=255

;6.  Add the baseline offset
model_profile_total_hires=model_profile_total_hires+baseline_offset

;;7.  Downsample to original resolution
model_profile_total=rebin(model_profile_total_hires,npoints)
;for i=0L,n_elements(velocity_hires)-1 do begin
;    lowres_index=i/upsample_factor
;;    hires_subindex=(i mod upsample_factor)-upsample_factor_halfwidth
;;    velocity_lowres=x[lowres_index]
;;    velocity_hires[i]=velocity_lowres+(hires_subindex/double(upsample_factor))*x_spacing_median
;;    ;print,i,lowres_index,hires_subindex,x[lowres_index],velocity_hires[i]
;endfor
;oplot,x,model_profile_total,color=255+255*256L

;stop
return,model_profile_total
end
