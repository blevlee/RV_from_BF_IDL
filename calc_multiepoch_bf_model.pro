;+
;; @brief  This is a function to compute a model broadening function
;; for every epoch observed.  The form of the arguments and return
;; value is for use with MPFITFUN.
;; @author Brian L. Lee
;; @date Jan. 24, 2013
;;
;; @dependencies calc_two_convolved_line_profiles
;-
function calc_multiepoch_bf_model, $
                                   x, $
                                   p, $
                                   bf_datavalues=bf_datavalues, $
                                   bf_errs=bf_errs, $
                                   index_epoch_segment_start_array=index_epoch_segment_start_array, $
                                   index_epoch_segment_end_array=index_epoch_segment_end_array, $
                                   n_params_epochspecific=n_params_epochspecific, $
                                   observation_jd=observation_jd, $
                                   phase_at_transit=phase_at_transit, $
                                   plot_to_file=plot_to_file, $
                                   plot_to_screen=plot_to_screen, $
                                   sine_period=sine_period, $
                                   time_of_transit=time_of_transit
;;1.  Integrity checks
if ~keyword_set(index_epoch_segment_start_array) then begin
    print,systime(/utc),'|ERROR|calc_multiepoch_bf_model|The required variable index_epoch_segment_start_array was not successfully received by this function.'
    stop
endif
if ~keyword_set(index_epoch_segment_end_array) then begin
    print,systime(/utc),'|ERROR|calc_multiepoch_bf_model|The required variable index_epoch_segment_end_array was not successfully received by this function.'
    stop
endif
n_seg_starts=n_elements(index_epoch_segment_start_array)
n_seg_ends=n_elements(index_epoch_segment_end_array)
if n_seg_starts ne n_seg_ends then begin
    print,systime(/utc),'|ERROR|calc_multiepoch_bf_model|The number of elements in index_epoch_segment_start_array did not match the number of elements in index_epoch_segment_end_array.'
    stop
endif
n_seg=n_seg_starts
;;2.  Set up variables
;;2.1.  Decode modelling parameters
;;2.1.1.  JD's of the epochs
observation_jd=observation_jd
;;2.1.2.  Parameters already locked by Kepler light curve fitting:
time_of_transit=time_of_transit
phase_at_transit=phase_at_transit
sine_period=sine_period
sine_angular_frequency=2.0*!DPI/sine_period ;;radians/day
;;2.1.4.  Pure epoch-to-epoch parameters
;n_params_epochspecific=4
;instrumental_lsf_sigma=p[0] ;;in index units, not velocity units;
;some epochs may share the same lsf_sigma
;area_under_peak1=p[1] ;;this is the scaling factor for the epoch
;baseline_offset=p[2] ;;expect this to be 0.0
;star1_rv=p[3] ;;normally locked at 0.0 by PARINFO
;;2.1.3.  General parameters to be jointly fit across all epochs
sine_amplitude           = p[n_seg*n_params_epochspecific+0]
sine_offset              = p[n_seg*n_params_epochspecific+1]
star1_limbdarkening      = p[n_seg*n_params_epochspecific+2]
star2_limbdarkening      = p[n_seg*n_params_epochspecific+3]
star1_rotationvelocity   = p[n_seg*n_params_epochspecific+4]
star2_rotationvelocity   = p[n_seg*n_params_epochspecific+5]
fluxratio_star2_to_star1 = p[n_seg*n_params_epochspecific+6]
;;2.1.5.  Epoch-to-epoch parameters derived from the general parameters
star2_rv_allepochs=sine_offset+sine_amplitude*sin(sine_angular_frequency*(observation_jd-time_of_transit)+phase_at_transit) ;;derived from the best-fit orbit parameters and the JD
;;2.2.  
;;3.  Loop over epochs
;;!p.multi=[0,1,3]
;;3.1.  Compute star 2's RV at this phase, given the sine phase, period,
;;semiamplitude, and offset
is_first_seg=1L
for i=0L,n_seg-1 do begin
;;3.1.1.  Extract this segment's data points
    x_seg=x[index_epoch_segment_start_array[i]:index_epoch_segment_end_array[i]]
    y_seg=bf_datavalues[index_epoch_segment_start_array[i]:index_epoch_segment_end_array[i]] ;for plotting and debugging info only-  not used in fit
    y_seg_err=bf_errs[index_epoch_segment_start_array[i]:index_epoch_segment_end_array[i]] ;for plotting and debugging info only-  not used in fit
;;3.1.2.  Extract this segment's parameters
;;p_seg=[p[i*n_params_epochspecific:(i+1)*n_params_epochspecific-1],[star2_rv_allepochs[i],star1_limbdarkening,star2_limbdarkening,star1_rotationvelocity,star2_rotationvelocity,fluxratio_star2_to_star1]]
    p_segspecific=[p[i*n_params_epochspecific:(i+1)*n_params_epochspecific-1]]
    instrumental_lsf_sigma=p_segspecific[0]
    star1_rv=p_segspecific[1]
    star2_rv=star2_rv_allepochs[i]
    area_under_peak1=p_segspecific[2]
    baseline_offset=p_segspecific[3]
;;star1_limbdarkening      = 0.6d
;;star2_limbdarkening      = 0.6d
;;star1_rotationvelocity   = 1.0d
;;star2_rotationvelocity   = 20.0d
;;fluxratio_star2_to_star1 = 0.5d
;;3.1.3.  Sort the segment parameters into the order required by the
;;two-star rotational model function
    p_seg=[ $
            star1_rotationvelocity, $
            star1_limbdarkening, $
            star1_rv, $
            area_under_peak1, $
            star2_rotationvelocity, $
            star2_limbdarkening, $
            star2_rv, $
            fluxratio_star2_to_star1, $
            instrumental_lsf_sigma, $
            baseline_offset $       
          ]
;;3.1.4.  Set up a high-resolution x grid on which to evaluate the
;;model, to be used to make nice smooth plots
    n_points_smoothplot=10000
    x_seg_hires=min(x_seg) + (max(x_seg)-min(x_seg))*(dindgen(n_points_smoothplot)/double(n_points_smoothplot))
;;3.1.4.  Feed x_seg and p_seg to the two-star rotational model
    y_seg_model=calc_two_convolved_line_profiles(x_seg,p_seg)
    ;y_seg_model_hires=calc_two_convolved_line_profiles(x_seg_hires,p_seg)
;;3.1.5.  Also compute the component one-star rotational models for overplotting
    p_seg_star1=[ $
                  p_seg[0], $
                  p_seg[1], $
                  p_seg[2], $
                  p_seg[3], $
                  p_seg[8], $
                  p_seg[9] $
                ]
    p_seg_star2=[ $
                  p_seg[4], $
                  p_seg[5], $
                  p_seg[6], $
                  p_seg[7]*p_seg[3], $
                  p_seg[8], $
                  p_seg[9] $
                ]
    y_seg_model_star1=calc_one_convolved_line_profile(x_seg,p_seg_star1)
    y_seg_model_star2=calc_one_convolved_line_profile(x_seg,p_seg_star2)
    ;y_seg_model_star1_hires=calc_one_convolved_line_profile(x_seg_hires,p_seg_star1)
    ;y_seg_model_star2_hires=calc_one_convolved_line_profile(x_seg_hires,p_seg_star2)
;;3.1.6.  Make plots
    if keyword_set(plot_to_file) then begin
        LoadCT, 39   ; Rainbow+white                                  
        c = { black :  0,  $
              purple :   .12*!d.n_colors, $
              blue :   .25*!d.n_colors, $
              medblue :         .32*!d.n_colors, $
              ltblue :         .40*!d.n_colors, $
              green :  .65*!d.n_colors, $
              yellow :         .75*!d.n_colors, $
              orange :         .80*!d.n_colors, $
              red :    .90*!d.n_colors, $
              white :  !d.n_colors-1   }
        color1=c.purple
        color2=c.red
        color3=c.medblue
        if keyword_set(plot_to_screen) then begin
            print,systime(/utc),'|WARN|calc_multiepoch_bf_model|You asked for plots to screen and to file simultaneously.  This functionality is not implemented.  Plotting to file only.  Setting plot_to_screen=0.'
            plot_to_screen=0
        endif
    endif else begin
        if keyword_set(plot_to_screen) then begin
            color1=255L+255*256L*256L
            color2=255L
            color3=255*256L+255*256L*256L
        endif
    endelse
    if keyword_set(plot_to_screen) or keyword_set(plot_to_file) then begin
        plot, $
          x_seg, $
          y_seg, $
          title='JD='+strtrim(string(observation_jd[i]),2)+', RV_Star_2='+strtrim(string(star2_rv,format='(f10.2)'),2)+' km/s', $
          xtitle='Velocity within line (km/s)', $
          ytitle='Norm. BF', $
          psym=3, $
          symsize=0.5, $
          charsize=1.75, $
          thick=2, $
          xrange=[-40,40], $
          /xs, $
          yrange=[-0.1,1.1], $
          /ys
        oploterror, $
          x_seg, $
          y_seg, $
          y_seg_err, $
          psym=3, $
          symsize=0.5, $
          thick=2, $
          /nohat
        oplot,x_seg,y_seg_model,color=color1
        oplot,x_seg,y_seg_model_star1,color=color2,linestyle=2
        oplot,[star1_rv,star1_rv],[-9e9,9e9],color=color2,linestyle=2
        oplot,x_seg,y_seg_model_star2,color=color3,linestyle=3
        oplot,[star2_rv,star2_rv],[-9e9,9e9],color=color3,linestyle=3
    endif
    if is_first_seg then begin
        y_model=y_seg_model
        is_first_seg=0
    endif else begin
        y_model=[y_model,y_seg_model]
    endelse
;stop
endfor
;stop
return,y_model
end
