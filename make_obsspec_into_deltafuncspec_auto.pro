;;Dependencies:
;    .comp ../read_spectrum.pro
function make_obsspec_into_deltafuncspec_auto, $
                                               lambda=lambda, $
                                               flux=flux, $
                                               linelist_lambda=linelist_lambda, $
                                               linelist_flux=linelist_flux
;;1.  Set up variables
rough_peak_finding_cut=0.9
diff_flagging_threshold=1.5
;;2.  Integrity checks
if ~keyword_set(lambda) then begin
    print,systime(/utc),'|ERROR|make_obsspec_into_deltafuncspec_auto|The lambda array was not defined'
    stop
endif
if ~keyword_set(flux) then begin
    print,systime(/utc),'|ERROR|make_obsspec_into_deltafuncspec_auto|The flux array was not defined'
    stop
endif
n_obs_lambda=n_elements(lambda)
n_obs_flux=n_elements(flux)
if n_obs_lambda ne n_obs_flux then begin
    print,systime(/utc),'|ERROR|make_obsspec_into_deltafuncspec_auto|The flux array was not the same length as the lambda array'
    stop
endif
if n_obs_lambda eq 0 then begin
    print,systime(/utc),'|ERROR|make_obsspec_into_deltafuncspec_auto|The lambda array had no elements'
    stop
endif
if n_obs_flux eq 0 then begin
    print,systime(/utc),'|ERROR|make_obsspec_into_deltafuncspec_auto|The flux array had no elements'
    stop
endif
;;3.  Find points in peaks
index_points_in_peaks=where(flux lt rough_peak_finding_cut)
plot,lambda,flux
oplot,lambda[index_points_in_peaks],flux[index_points_in_peaks],color=255
;;4.  Refine each group into a single line
;;4.1  Detect gaps
;;4.1.1  Determine typical spacing between observed lambdas
lambdadiff=lambda[1:n_obs_lambda-1]-lambda[0:n_obs_lambda-2]
lambdadiff_coeff=robust_linefit(lambda[0:n_obs_lambda-2],lambdadiff,yfit)
print,lambdadiff_coeff
oplot,lambda*lambdadiff_coeff[1]+lambdadiff_coeff[0],color=255
;;4.1.2  Search peaks list for lambda gaps
group_counter=0L
is_first_gap=1
lambda_diffs_within_peakselections=dblarr(n_elements(index_points_in_peaks)-1)
for i=0L,n_elements(index_points_in_peaks)-2 do begin
    lambda_diffs_within_peakselections[i]=lambda[index_points_in_peaks[i+1]]-lambda[index_points_in_peaks[i]]
;;4.1.2.1  Compare current diff with expected diff.  If diff is bigger
;;than threshold*expected diff, flag this one as a break between
;;different peaks.
    expected_diff = lambda[index_points_in_peaks[i]]*lambdadiff_coeff[1]+lambdadiff_coeff[0]

    ;print,'******'
    ;print,i,index_points_in_peaks[i],lambda[index_points_in_peaks[i]]
    ;print,lambda_diffs_within_peakselections[i],diff_flagging_threshold*expected_diff

    if lambda_diffs_within_peakselections[i] gt diff_flagging_threshold*expected_diff then begin
        if is_first_gap eq 1 then begin
            index_gap_start_locations=index_points_in_peaks[i]
            index_peak_start_locations=index_points_in_peaks[0]
            index_peak_end_locations=index_points_in_peaks[i]
            is_first_gap=0
        endif else begin
            index_gap_start_locations=[index_gap_start_locations,index_points_in_peaks[i]]
            index_peak_start_locations=[index_peak_start_locations,index_next_peak_start_location]
            index_peak_end_locations=[index_peak_end_locations,index_points_in_peaks[i]]
        endelse
        index_next_peak_start_location=index_points_in_peaks[i+1]

    ;astrolib
    ;forprint,index_peak_start_locations,index_peak_end_locations

    endif

    ;print,'******'

    ;stop
endfor
plot,lambda[index_points_in_peaks],lambda_diffs_within_peakselections,psym=3
;;4.2  For each peak we found, find the minimum of the line and record it in
;;the linelist.
n_peaks=n_elements(index_peak_start_locations)
linelist_lambda=dblarr(n_peaks)
linelist_flux=dblarr(n_peaks)
for i=0L,n_peaks-1 do begin
    this_peak_starting_lambda=lambda[index_peak_start_locations[i]]
    this_peak_ending_lambda=lambda[index_peak_end_locations[i]]
    ;index_metalevel=where(lambda[index_points_in_peaks] lt this_peak_ending_lambda and lambda[index_points_in_peaks] gt this_peak_starting_lambda)
    ;plot,lambda[index_points_in_peaks[index_metalevel]],
    index_points_in_this_peak=where(lambda ge this_peak_starting_lambda and lambda le this_peak_ending_lambda)
    plot,lambda[index_points_in_this_peak],flux[index_points_in_this_peak],/ynozero,psym=4
    linelist_flux[i]=min(flux[index_points_in_this_peak],subindex_this_peak_lambda)
    linelist_lambda[i]=lambda[index_peak_start_locations[i]+subindex_this_peak_lambda]
endfor
plot,lambda,flux
oplot,linelist_lambda,linelist_flux,psym=4,color=255
;;5.  Scale flux by 10000.0 to match BASS2000 template level
linelist_flux=linelist_flux*10000.0
;;6.  Return success
;forprint,textout='/astro/net/astro-agol/blevlee/DATA/HET/TYC3559_20Nov2012/templates/self/orders/linelist_order'+strtrim(string(order,format='(i3.3)'),2)+'.txt',linelist_lambda,linelist_flux,/nocomment
return,1
end
