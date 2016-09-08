function combine_bf_results, $
                             vel_all_orders=vel_all_orders, $
                             bbb_all_orders=bbb_all_orders, $
                             bb30_all_orders=bb30_all_orders, $
                             sum_bbb=sum_bbb, $
                             sum_bb30=sum_bb30, $
                             sum_vel_axis=sum_vel_axis, $
                             sum_bb_err=sum_bb_err, $
                             window_width_pixels=window_width_pixels, $
                             n_orders=n_orders, $
                             register_peaks=register_peaks, $
                             order_weights=order_weights
;;1.  Set up variables
vel_min=min(vel_all_orders,/nan)
vel_max=max(vel_all_orders,/nan)
sum_vel_axis=(vel_max-vel_min)*(dindgen(window_width_pixels)/double(window_width_pixels)) + vel_min
sum_bbb=dblarr(n_elements(sum_vel_axis))
sum_bb30=dblarr(n_elements(sum_vel_axis))
sum_bb_err=dblarr(n_elements(sum_vel_axis))
sum_bb_err_tempvar=0.0
peak_fitting_domain_halfwidth=3
peak_fitting_domain=peak_fitting_domain_halfwidth+1 ;;number of data points to use for the peak fit in section 2.3
if ~keyword_set(order_weights) then order_weights=replicate(1.0,n_orders)
if n_elements(order_weights) ne n_orders then begin
    print,systime(/utc),'|ERROR|combine_bf_results|The number of elements in order_weights was not equal to n_orders.  Stopping.'
    stop
endif
is_first_order=1L
;;2.  Loop over orders
for i=0L,n_orders-1 do begin
    bbb_tmp=bbb_all_orders[i,*]
    bb30_tmp=bb30_all_orders[i,*]
    vel_tmp=vel_all_orders[i,*]
    bbb_interp=interpol(bbb_tmp,vel_tmp,sum_vel_axis)
    bb30_interp=interpol(bb30_tmp,vel_tmp,sum_vel_axis)
;;2.1 Trim off the ends, where the interpolation was not valid, by setting to NaN
    index_outofrangevels=where(sum_vel_axis gt max(vel_tmp,/nan) or sum_vel_axis lt min(vel_tmp,/nan),count_outofrangevels)
    if count_outofrangevels gt 0 then begin
        bbb_interp[index_outofrangevels]=!values.f_nan
        bb30_interp[index_outofrangevels]=!values.f_nan
    endif
;;2.3 If register_peaks is set, do a peak fit and shift the centre of the peak to the 0th pixel before coadding it:
    if keyword_set(register_peaks) then begin
        max_bf_value=max(bbb_interp,index_max_bf_value,/nan)
;;2.3.1 Check that the peak is far enough away from the edges where the NaN's are
        testarr_edge_proximity=where( $
                                      index_outofrangevels eq index_max_bf_value-peak_fitting_domain or $
                                      index_outofrangevels eq index_max_bf_value+peak_fitting_domain or $
                                      index_max_bf_value lt peak_fitting_domain or $
                                      index_max_bf_value gt n_elements(bbb_interp)-peak_fitting_domain, $
                                      count_edge_proximity $
                                    )
        if count_edge_proximity gt 0 then begin
;;2.3.2 If there was a strange peak maximum detection, do a straight coadd instead of a registered coadd.
            print,systime(/utc),'|WARN|combine_bf_results.pro|The peak of the broadening function on order_loop '+string(i)+' did not satisfy the array-edge avoidance criteria.  Using unregistered BF coadd instead.'
            sum_bbb=sum_bbb+bbb_interp
            sum_bb30=sum_bb30+bb30_interp
        endif else begin
;;2.3.3 Otherwise, proceed to register the peak.
;;2.3.3.1 Extract the data points in the neighbourhood of the peak.
            start_index=index_max_bf_value-peak_fitting_domain_halfwidth
            end_index=index_max_bf_value+peak_fitting_domain_halfwidth
            peakfinding_xaxis=sum_vel_axis[start_index:end_index]
            peakfinding_yaxis=bbb_interp[start_index:end_index]
;;2.3.3.2 Fit a parabola to the neighbourhood
            polydegree=2
            polycoeffs=robust_poly_fit(peakfinding_xaxis,peakfinding_yaxis,polydegree,yfit)
;;2.3.3.3 Compute the location of the peak of the parabolic fit
            directrix=-0.5*polycoeffs[1]/polycoeffs[2]
            ;plot,peakfinding_xaxis,peakfinding_yaxis,psym=4
            ;oplot,peakfinding_xaxis,yfit,color=255
                                ;oplot,[directrix,directrix],[-9e9,9e9],color=-255
;;2.3.3.4 Shift vel_tmp so that velocity=0 is at the directrix
            vel_tmp_shifted=vel_tmp-directrix
;;2.3.3.5 Interpolate onto the common velocity axis
            bbb_interp=interpol(bbb_tmp,vel_tmp_shifted,sum_vel_axis)
            bb30_interp=interpol(bb30_tmp,vel_tmp_shifted,sum_vel_axis)
;;2.3.3.6 Trim off the ends, where the interpolation was not valid, by setting to NaN
            index_outofrangevels=where(sum_vel_axis gt max(vel_tmp_shifted,/nan) or sum_vel_axis lt min(vel_tmp_shifted,/nan),count_outofrangevels)
            if count_outofrangevels gt 0 then begin
                bbb_interp[index_outofrangevels]=!values.f_nan
                bb30_interp[index_outofrangevels]=!values.f_nan
            endif
        endelse
    endif
;;2.4 Normalize the peak height of the profile to 1.0 before coadding
    bbb_interp=bbb_interp/max(bbb_interp,/nan)
    bb30_interp=bb30_interp/max(bb30_interp,/nan)
;;2.5 Coadd the profiles (may or may not have gone through registration)
    sum_bbb=sum_bbb+bbb_interp*order_weights[i]
    sum_bb30=sum_bb30+bb30_interp*order_weights[i]
    sum_bb_err_tempvar=sum_bb_err_tempvar+order_weights[i]
    if is_first_order then begin
        plot,sum_vel_axis,bbb_interp/max(bbb_interp,/nan),xtitle='Peak-registered RV (km/s)',ytitle='Peak-normalized BF',title='BF combination over orders'
        is_first_order=0
    endif else begin
        oplot,sum_vel_axis,bbb_interp/max(bbb_interp,/nan),color=255
    endelse
endfor
oplot,sum_vel_axis,sum_bbb/max(sum_bbb,/nan),color=-255
;;3.  Normalize result by the total of the weights
sum_bbb=sum_bbb/total(order_weights)
sum_bb30=sum_bb30/total(order_weights)
sum_bb_err=sum_bb_err+1.0/sqrt(sum_bb_err_tempvar)
wait,2
;stop
return,1
end
