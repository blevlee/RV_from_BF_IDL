function calc_one_convolved_line_profile, x, p
;1.  Integrity check
;1.1  Check for even velocity spacing in the input x array
x_spacing=x[1:n_elements(x)-1]-x[0:n_elements(x)-2]
x_spacing_median=median(x_spacing)
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
instrumental_profile_sigma=p[4]
baseline_offset=p[5]
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
status=calc_rotational_broadening_kernel( $
                                          n_points_in_kernel=npoints_hires, $
                                          equatorial_max_rv=equatorial_max_rv1, $
                                          limbdarkening_beta=limbdarkening_beta1, $
                                          vel_halfdomain=vel_halfdomain_hires, $ ;km/s
                                          x_midpoint=x_midpoint_hires, $
                                          velocity=velocity_hires_out, $
                                          normalized_kernel=model_profile_hires $
                                        )
;;2.2.  Apply a velocity shift to the line
shift=velocity_peak_offset1-x_midpoint_hires
model_profile_shifted=interpol(model_profile_hires,velocity_hires_out,velocity_hires_out-shift)
;;2.3.  Scale the profile to have the desired area
area=int_tabulated(velocity_hires_out,model_profile_hires,/double)
model_profile_total_hires=model_profile_shifted*area_under_profile1/area

;;3.  Convolve astrophysical line profile with instrumental broadening
instrumental_profile=gs_smooth(npoints_hires,instrumental_profile_sigma*upsample_factor)
model_profile_total_hires=cnv(model_profile_total_hires,instrumental_profile)

;;4.  Add the baseline offset
model_profile_total_hires=model_profile_total_hires+baseline_offset

;;5.  Downsample to original resolution
model_profile_total=rebin(model_profile_total_hires,npoints)

return,model_profile_total
end
