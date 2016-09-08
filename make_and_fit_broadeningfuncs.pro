;+
;; @brief make_and_fit_broadeningfuncs is a driver to derive the RV of
;; the second star in TYC 3559 by fitting its small peak in the
;; Rucinski broadening function.
;; @author Brian L. Lee
;; @date Jan. 2013
;;
;; @dependencies BFall_IDL.pro
;-
pro make_and_fit_broadeningfuncs
;;1.  Set up variables
;;In the future, these variables ought to be read in from a parameter
;;file so that they are not hard-coded.
bf_computation_savefilename='all_broadening_functions_auto.sav'
bf_fitting_savefilename='fit_broadening_functions.sav'
if file_test(bf_computation_savefilename) then begin
    restore,bf_computation_savefilename
endif else begin
n_observations=6
window_width_pixels=201
observation_dirnames=[ $
                       '/media/PENDRIVE/astro/net/astro-agol/blevlee/DATA/HET/TYC3559_20Nov2012/for_brian/TYC3559_60k/calspectra/2012-08-13/', $
                       '/media/PENDRIVE/astro/net/astro-agol/blevlee/DATA/HET/TYC3559_20Nov2012/for_brian/TYC3559_60k/calspectra/2012-09-02/', $
                       '/media/PENDRIVE/astro/net/astro-agol/blevlee/DATA/HET/TYC3559_20Nov2012/for_brian/TYC3559_60k/calspectra/lowsnr/2012-08-14/', $
                       '/media/PENDRIVE/astro/net/astro-agol/blevlee/DATA/HET/TYC3559_20Nov2012/for_brian/TYC3559_30k/calspectra/2012-08-09/', $
                       '/media/PENDRIVE/astro/net/astro-agol/blevlee/DATA/HET/TYC3559_20Nov2012/for_brian/TYC3559_30k/calspectra/2012-08-10/', $
                       '/media/PENDRIVE/astro/net/astro-agol/blevlee/DATA/APO3.5/ARCES/tyc3559/textspec/' $
                     ]
observation_filebasenames=[ $
                            'TYC3559_###_b_WLFL_CN_FP.txt', $ ;The ### stands for the number of digits in the order number
                            'TYC3559_###_b_WLFL_CN_FP.txt', $ ;The ### will be replaced in a loop later on
                            'TYC3559_###_b_WLFL_CN_FP.txt', $
                            'TYC3559_###_b_WLFL_CN_FP.txt', $
                            'TYC3559_###_b_WLFL_CN_FP.txt', $
                            'TYC3559-02080-1.0024.ec.###.txt' $
                          ]
observation_jd=[ $
                       2456152.622671d, $
                       2456172.807521d, $
                       2456153.851364d, $
                       2456148.862618d, $
                       2456149.872712d, $
                       2456151.813257d $ ;bad: 2456152.053646d $
                     ]
observation_starting_orders=[ $
                              3, $
                              3, $
                              3, $
                              3, $
                              3, $
                              44 $
                            ]
observation_n_orders=[ $
                       33, $
                       33, $
                       33, $
                       33, $
                       33, $
                       28 $
                     ]
;;2.  Integrity checks
;;2.1.  Check that the requested order files are already present in
;;text format
;;2.1.1.  Create an observation_orders struct to keep track of the
;;info in each order
array_of_recordkeeping_structures=ptrarr(n_observations,/allocate_heap)
for i=0L,n_observations-1 do begin
    struct_for_one_observation={ $
                                 observation_dirnames:observation_dirnames[i], $
                                 observation_filebasenames:observation_filebasenames[i], $
                                 observation_jd:observation_jd[i], $
                                 observation_starting_orders:observation_starting_orders[i], $
                                 observation_n_orders:observation_n_orders[i] $
                               }
    *(array_of_recordkeeping_structures[i])=struct_for_one_observation
endfor
;;2.1.2.  Construct the fully qualified paths to the files and add
;;those filenames to the observation_orders struct
for i=0L,n_observations-1 do begin
;;2.1.2.1.  Unpack the struct
    dirname=(*(array_of_recordkeeping_structures[i])).observation_dirnames ;(0)
    filebasename=(*(array_of_recordkeeping_structures[i])).observation_filebasenames ;(1)
    jd=(*(array_of_recordkeeping_structures[i])).observation_jd
    starting_order=(*(array_of_recordkeeping_structures[i])).observation_starting_orders
    n_orders=(*(array_of_recordkeeping_structures[i])).observation_n_orders
    order_strings=strtrim(string(lindgen(n_orders)+starting_order,format='(i3.3)'),2)
;;2.1.2.2.  Form the fully qualified paths
    order_fullfilenames=replicate(dirname+filebasename,n_orders)
    for j=0L,n_orders-1 do begin
        order_fullfilenames[j]=repstr(order_fullfilenames[j],'###',order_strings[j])
    endfor
;;2.1.2.3.  Put the data back in the struct (now with extras)
    struct_for_one_observation={ $
                                 observation_dirnames:dirname, $
                                 observation_filebasenames:filebasename, $
                                 observation_jd:jd, $
                                 observation_starting_orders:starting_order, $
                                 observation_n_orders:n_orders, $
                                 observation_order_strings:order_strings, $
                                 observation_order_fullfilenames:order_fullfilenames $
                               }
    *(array_of_recordkeeping_structures[i])=struct_for_one_observation
    ;print,(*(array_of_recordkeeping_structures[i])).observation_order_fullfilenames
endfor
;;2.1.3.  Loop through the fully qualified paths to the files
for i=0L,n_observations-1 do begin
    n_orders=(*(array_of_recordkeeping_structures[i])).observation_n_orders
    order_fullfilenames=(*(array_of_recordkeeping_structures[i])).observation_order_fullfilenames
    for j=0L,n_orders-1 do begin
;;2.1.3.1  Check that files are present;  if not, issue an error and stop.
        if file_test(order_fullfilenames[j]) ne 1 then begin
            print,systime(/UTC),'|ERROR|make_and_fit_broadeningfuncs|The data file '+order_fullfilenames[j]+' was not found.  Stopping.'
            stop
        endif else begin
;;2.1.3.2  If files are present, check that they are not corrupted and have the correct format.
        endelse
    endfor
endfor
;;2.1.4.  Garbage collect on any loose blocks of memory after all that
;;pointer manipulation.
heap_gc

;;3.  Loop over number of observations
bf_summed=dblarr(n_observations,window_width_pixels)
err_bf_summed=dblarr(n_observations,window_width_pixels)
velocityaxis_for_bf_summed=dblarr(n_observations,window_width_pixels)
;initindexarray_for_bf_summed=lonarr(n_observations,window_width_pixels)
for i=0L,n_observations-1 do begin
;for i=2,2 do begin
;;3.1.  Unpack the struct
    dirname=(*(array_of_recordkeeping_structures[i])).observation_dirnames
    filebasename=(*(array_of_recordkeeping_structures[i])).observation_filebasenames
    jd=(*(array_of_recordkeeping_structures[i])).observation_jd
    starting_order=(*(array_of_recordkeeping_structures[i])).observation_starting_orders
    n_orders=(*(array_of_recordkeeping_structures[i])).observation_n_orders
    order_strings=(*(array_of_recordkeeping_structures[i])).observation_order_strings
    order_fullfilenames=(*(array_of_recordkeeping_structures[i])).observation_order_fullfilenames
;;3.2.  Make the broadening functions, first pass
;;3.2.1.  Set up arrays to hold the BF's for this epoch
    vel_all_orders=dblarr(n_orders,window_width_pixels)
    bbb_all_orders=dblarr(n_orders,window_width_pixels)
    bb30_all_orders=dblarr(n_orders,window_width_pixels)
    stddev_order_relative=dblarr(n_orders)
    order_linelist_lambda_ptrarr=ptrarr(n_orders,/allocate_heap)
    order_linelist_flux_ptrarr=ptrarr(n_orders,/allocate_heap)
;;3.2.1.  Loop over orders
    for j=0L,n_orders-1 do begin
;;3.2.2.  Load the order file
        print,i,j,order_fullfilenames[j]    
        status=read_spectrum( $
                              spectrum_file=order_fullfilenames[j], $
                              lambda=lambda, $
                              flux=flux $
                            )
;;3.2.3.  If needed, trim junky data off the ends of the arrays and renormalize
        if i eq 5 then begin ;epoch i=5 is the APO spectrum
            start_good_region=200
            end_good_region=n_elements(lambda)-201
            lambda=lambda[start_good_region:end_good_region]
            flux=flux[start_good_region:end_good_region]
            flux=flux/median(flux)
        endif
;;3.2.4.  Make the line list
;;******31Jul2016:  TO DO:  We ought to try inserting the PHOENIX spectrum here to see how that affects the BF.
	lambda_min_data=min(lambda)
	lambda_max_data=max(lambda)
;	stop
;load PHOENIX model spectrum
	status2=read_spectrum(spectrum_file='/media/PENDRIVE/tyc3559/to_brian/lte0550*.txt',lambda=lambda_phoenix,flux=flux_phoenix)
;trim PHOENIX model spectrum since it is too long for doing the BF SVD
	index_lowlambda=where(lambda_phoenix lt lambda_min_data)
	if n_elements(index_lowlambda) ge 2 then begin
		index_lambda_min_data=index_lowlambda[n_elements(index_lowlambda)-1]
	endif else begin
		print,systime(/UTC),'|ERROR|make_and_fit_broadeningfuncs|The PHOENIX model spectrum did not span enough wavelength range (low end).  Stopping.'
	        stop
	endelse
	index_highlambda=where(lambda_phoenix gt lambda_max_data)
	if n_elements(index_highlambda) ge 2 then begin
		index_lambda_max_data=index_highlambda[0]
	endif else begin
		print,systime(/UTC),'|ERROR|make_and_fit_broadeningfuncs|The PHOENIX model spectrum did not span enough wavelength range (high end).  Stopping.'
	        stop
	endelse
	lambda_phoenix=lambda_phoenix[index_lambda_min_data:index_lambda_max_data]
	flux_phoenix=flux_phoenix[index_lambda_min_data:index_lambda_max_data]
;stop
        status=make_obsspec_into_deltafuncspec_auto( $
                                                     lambda=lambda_phoenix, $
                                                     flux=flux_phoenix, $
                                                     linelist_lambda=linelist_lambda, $
                                                     linelist_flux=linelist_flux $
                                                   )
;;3.2.4.1.  Save the line list for use with other spectra off the same
;;instrument which had lower-quality observations.
        *(order_linelist_lambda_ptrarr[j])=linelist_lambda
        *(order_linelist_flux_ptrarr[j])=linelist_flux
        struct_for_one_observation={ $
                                     observation_dirnames:dirname, $
                                     observation_filebasenames:filebasename, $
                                     observation_jd:jd, $
                                     observation_starting_orders:starting_order, $
                                     observation_n_orders:n_orders, $
                                     observation_order_strings:order_strings, $
                                     observation_order_fullfilenames:order_fullfilenames, $
                                     observation_order_linelist_lambda_ptrarr:order_linelist_lambda_ptrarr, $
                                     observation_order_linelist_flux_ptrarr:order_linelist_flux_ptrarr $
                                   }
        *(array_of_recordkeeping_structures[i])=struct_for_one_observation
;;3.2.4.2.  For select epochs i, restore a line list from a previous line list if it would
;;be superior.
;;How to decode the stored linelist pointers:
;;linelist_lambda=*(((*(array_of_recordkeeping_structures[i_master])).observation_order_linelist_lambda_ptrarr)[j])
;;linelist_flux=*(((*(array_of_recordkeeping_structures[i_master])).observation_order_linelist_flux_ptrarr)[j])
        if i ge 1 and i le 4 then begin
            i_master=0
            linelist_lambda=*(((*(array_of_recordkeeping_structures[i_master])).observation_order_linelist_lambda_ptrarr)[j])
            linelist_flux=*(((*(array_of_recordkeeping_structures[i_master])).observation_order_linelist_flux_ptrarr)[j])
        endif
;;3.2.5.  Make the BF for this order and save it in the array of BF's
;;for this epoch
        status=my_bf_script_auto( $
                                  lambda=lambda, $
                                  flux=flux, $
                                  linelist_lambda=linelist_lambda, $
                                  linelist_flux=linelist_flux, $
                                  vel=vel_tmp, $
                                  bbb=bbb_tmp, $
                                  bb30=bb30_tmp, $
                                  ;rv_resolution=0.5, $ ;km/s
                                  window_width_pixels=window_width_pixels, $
                                  wait_time_for_plots=0 $
                                )
        vel_all_orders[j,*]=vel_tmp
        bbb_all_orders[j,*]=bbb_tmp
        bb30_all_orders[j,*]=bb30_tmp
    endfor
;stop
;;3.2.6.  Combine the orderwise BF's:  no weighting
    status=combine_bf_results( $
                               vel_all_orders=vel_all_orders, $
                               bbb_all_orders=bbb_all_orders, $
                               bb30_all_orders=bb30_all_orders, $
                               sum_bbb=sum_bbb, $
                               sum_bb30=sum_bb30, $
                               sum_vel_axis=sum_vel_axis, $
                               window_width_pixels=window_width_pixels, $
                               n_orders=n_orders, $
                               register_peaks=1 $
                             )

;;3.3.  EITHER fit the unweighted summed broadening function OR use the
;;combined BF as a model empirical BF (use the empirical route if the
;;final coadded SNR is good enough)
;;3.3.3.  Create the independent matrix for use in REGRESS:
    x=transpose(sum_bbb)
;;3.3.1.  For empirical model, clean out the NaN's
    index_nan=where(finite(x,/nan),count_nan)
    x[index_nan]=0.0
;;3.3.2.  Create an evenly weighted measure_errors array, but make
;;errors big where there were NaN's
    measure_errors=replicate(1.0,n_elements(x))
    measure_errors[index_nan]=9e9
;;3.3.3.  Loop over orders
    for j=0L,n_orders-1 do begin
;;3.3.3.1.  Evaluate the fractional error in this order relative to the
;;model BF (fit for best amplitude using REGRESS)
;;3.3.3.1.1  Create the dependent matrix for use in REGRESS:
        y=reform(bbb_all_orders[j,*])
        scalingcoeff=regress(x,y,measure_errors=measure_errors,const=const,sigma=sigma)
        scaled_model=sum_bbb*scalingcoeff[0] + const
        plot,sum_vel_axis,y
        oplot,sum_vel_axis,scaled_model,color=255
        stddev_order_relative[j]=stddev(y-scaled_model,/nan)/max(scaled_model,/nan)
    endfor

;;3.4.  Make the broadening functions, second pass
;;3.4.1.  Combine the orderwise BF's, weighting according to the error
;;computed in 3.3.3.1
    status=combine_bf_results( $
                               vel_all_orders=vel_all_orders, $
                               bbb_all_orders=bbb_all_orders, $
                               bb30_all_orders=bb30_all_orders, $
                               sum_bbb=sum_bbb_weighted, $
                               sum_bb30=sum_bb30_weighted, $
                               sum_bb_err=sum_bb_err_weighted, $
                               sum_vel_axis=sum_vel_axis_weighted, $
                               window_width_pixels=window_width_pixels, $
                               n_orders=n_orders, $
                               register_peaks=1, $
                               order_weights=1.0/stddev_order_relative^2 $
                             )
    plot,sum_vel_axis,sum_bbb,psym=-4
    oplot,sum_vel_axis_weighted,sum_bbb_weighted,color=255*256L,psym=-4
    oploterror,sum_vel_axis_weighted,sum_bbb_weighted,sum_bb_err_weighted,color=255
    print,sum_bb_err_weighted
;;3.5  Save the second pass BF for this epoch.
    bf_summed[i,*]=sum_bbb_weighted
    err_bf_summed[i,*]=sum_bb_err_weighted
    velocityaxis_for_bf_summed[i,*]=sum_vel_axis_weighted
;    initindexarray_for_bf_summed[i,*]=lindgen(window_width_pixels)
endfor
save,/all,file=bf_computation_savefilename
endelse

;;4.  Do the joint parameter fit to all epochs.
;;4.1.  Concatenate the BF data in preparation for joint MPFIT
;;4.1.1.  Transpose the data arrays so that they unwrap to 1-D with
;;indexing sequential within each segment
bf_summed=transpose(bf_summed)
err_bf_summed=transpose(err_bf_summed)
velocityaxis_for_bf_summed=transpose(velocityaxis_for_bf_summed)
;;4.1.1.  Mark the indices of the data arrays with a 1-D indexing number before
;;we begin  
initindexarray_for_bf_summed=lindgen(n_observations*window_width_pixels)
initindexarray_for_bf_summed=reform(initindexarray_for_bf_summed,window_width_pixels,n_observations)
;;4.1.1.  Trim out NaN's and excess data far from the line centre
lineprofile_max_velocity=60.0 ;km/s
index_points_to_trim=where(finite(bf_summed,/nan) or abs(velocityaxis_for_bf_summed) gt lineprofile_max_velocity,complement=index_points_to_keep)
bf_summed_trimmed=bf_summed[index_points_to_keep]
err_bf_summed_trimmed=err_bf_summed[index_points_to_keep]
velocityaxis_for_bf_summed_trimmed=velocityaxis_for_bf_summed[index_points_to_keep]
initindexarray_for_bf_summed_trimmed=initindexarray_for_bf_summed[index_points_to_keep]
segment_index_trimmed=initindexarray_for_bf_summed_trimmed/window_width_pixels
;;4.1.2.  Log the starting and ending elements of the epochs
index_epoch_segment_start_array=lindgen(n_observations)
index_epoch_segment_end_array=lindgen(n_observations)
for i=0L,n_observations-1 do begin
    index_trimmed_this_observation=where(segment_index_trimmed eq i,count_trimmed_this_observation)
    if count_trimmed_this_observation lt 5 then begin
        print,systime(/utc),'|ERROR|make_and_fit_broadeningfuncs|Could not find enough valid indices for observation '+strtrim(string(i),2)+'.  Stopping.'
        stop
    endif
    index_epoch_segment_start_array[i]=min(index_trimmed_this_observation)
    index_epoch_segment_end_array[i]=max(index_trimmed_this_observation)
endfor
;;4.2.  Set up starting guesses
;;There's a lot of repetition here, so I should probably make a
;;setting function.
n_params_epochspecific=4
n_params_general=7
p=[dblarr(n_params_epochspecific*n_observations),dblarr(n_params_general)]
;;4.2.0.  Epoch i=0 guesses for the params_epochspecific
i=0
instrumental_lsf_sigma=4.0 ;currently units of velocityaxis tickmarks;  need to change the cnv functionality to make this physical units
star1_rv=0.0
area_under_peak1=10.0
baseline_offset=0.0
p[i*n_params_epochspecific:(i+1)*n_params_epochspecific-1]=[instrumental_lsf_sigma,star1_rv,area_under_peak1,baseline_offset]
;;4.2.1.  Epoch i=1 guesses for the params_epochspecific
i=1
instrumental_lsf_sigma=5.0 ;currently units of velocityaxis tickmarks;  need to change the cnv functionality to make this physical units
star1_rv=0.0
area_under_peak1=9.0
baseline_offset=0.0
p[i*n_params_epochspecific:(i+1)*n_params_epochspecific-1]=[instrumental_lsf_sigma,star1_rv,area_under_peak1,baseline_offset]
;;4.2.2.  Epoch i=2 guesses for the params_epochspecific
i=2
instrumental_lsf_sigma=5.0 ;currently units of velocityaxis tickmarks;  need to change the cnv functionality to make this physical units
star1_rv=0.0
area_under_peak1=9.0
baseline_offset=0.0
p[i*n_params_epochspecific:(i+1)*n_params_epochspecific-1]=[instrumental_lsf_sigma,star1_rv,area_under_peak1,baseline_offset]
;;4.2.3.  Epoch i=3 guesses for the params_epochspecific
i=3
instrumental_lsf_sigma=2.0 ;currently units of velocityaxis tickmarks;  need to change the cnv functionality to make this physical units
star1_rv=0.0
area_under_peak1=14.0
baseline_offset=0.0
p[i*n_params_epochspecific:(i+1)*n_params_epochspecific-1]=[instrumental_lsf_sigma,star1_rv,area_under_peak1,baseline_offset]
;;4.2.4.  Epoch i=4 guesses for the params_epochspecific
i=4
instrumental_lsf_sigma=2.0 ;currently units of velocityaxis tickmarks;  need to change the cnv functionality to make this physical units
star1_rv=0.0
area_under_peak1=13.0
baseline_offset=0.0
p[i*n_params_epochspecific:(i+1)*n_params_epochspecific-1]=[instrumental_lsf_sigma,star1_rv,area_under_peak1,baseline_offset]
;;4.2.5.  Epoch i=5 guesses for the params_epochspecific
i=5
instrumental_lsf_sigma=1.5 ;currently units of velocityaxis tickmarks;  need to change the cnv functionality to make this physical units
star1_rv=0.0
area_under_peak1=15.0
baseline_offset=0.0
p[i*n_params_epochspecific:(i+1)*n_params_epochspecific-1]=[instrumental_lsf_sigma,star1_rv,area_under_peak1,baseline_offset]
;;4.2.6.  General params used over all epochs:
sine_amplitude           = 3.000d ;in the units of the velocityaxis (usually km/s)
sine_offset              = 1.500d ;in the units of the velocityaxis (usually km/s)
star1_limbdarkening      = 0.6d
star2_limbdarkening      = 0.6d
star1_rotationvelocity   = 1.0d ;in the units of the velocityaxis (usually km/s)
star2_rotationvelocity   = 20.0d ;in the units of the velocityaxis (usually km/s)
fluxratio_star2_to_star1 = 0.5d
p[n_observations*n_params_epochspecific:n_observations*n_params_epochspecific+n_params_general-1]=[ $
                                                                                                    sine_amplitude, $
                                                                                                    sine_offset, $
                                                                                                    star1_limbdarkening, $
                                                                                                    star2_limbdarkening, $
                                                                                                    star1_rotationvelocity, $
                                                                                                    star2_rotationvelocity, $
                                                                                                    fluxratio_star2_to_star1 $
                                                                                                  ]
;;4.2.7.  Known parameters from Kepler light curve model:
time_of_transit=2455004.80104d ;;At time of transit, RV_star=0.0 and is about to decrease (planet to head away, star to head towards).
phase_at_transit=!DPI ;;set transit at halfway through the phase in preparation for decreasing RV_star.
sine_period=5.5664943d ;;days
;;4.2.8.  Set up parinfo for MPFIT parameter locking
parinfo=replicate( $
                   { $
                     value:0.D, $
                     limited:[0,0], $
                     limits:[0.D,0], $
                     fixed:0 $
                   }, $
                   n_elements(p) $
                 )
start=p
parinfo[*].value=start
parinfo[0+n_params_epochspecific*lindgen(n_observations)].fixed=0 ;instrument_lsf_sigma
parinfo[1+n_params_epochspecific*lindgen(n_observations)].fixed=1 ;star1_rv
parinfo[2+n_params_epochspecific*lindgen(n_observations)].fixed=0 ;area_under_peak1
parinfo[3+n_params_epochspecific*lindgen(n_observations)].fixed=1 ;baseline_offset
parinfo[n_observations*n_params_epochspecific+2].fixed=1 ;star1_limbdarkening
parinfo[n_observations*n_params_epochspecific+3].fixed=1 ;star2_limbdarkening
;parinfo[n_observations*n_params_epochspecific+4].limited=[1,0] ;star1_rotationvelocity
;parinfo[n_observations*n_params_epochspecific+4].limits=[0.1d,5.0d] ;star1_rotationvelocity
;;4.3.  Run MPFIT
x=velocityaxis_for_bf_summed_trimmed
y=bf_summed_trimmed
yerr=err_bf_summed_trimmed
window,xsize=792,ysize=1024
!p.multi=[0,1,n_observations]
;;4.3.0  DEBUG:  Check behaviour of the model function
y_model_preliminary=calc_multiepoch_bf_model( $
                                  x, $
                                  p, $
                                  bf_datavalues=bf_summed_trimmed, $
                                  bf_errs=err_bf_summed_trimmed, $
                                  index_epoch_segment_start_array=index_epoch_segment_start_array, $
                                  index_epoch_segment_end_array=index_epoch_segment_end_array, $
                                  n_params_epochspecific=n_params_epochspecific, $
                                  observation_jd=observation_jd, $
                                  phase_at_transit=phase_at_transit, $
                                  plot_to_screen=1, $
                                  sine_period=sine_period, $
                                  time_of_transit=time_of_transit $
                                )
;;4.3.1  Run MPFIT.  Use FUNCTARGS to pass extras.
coeffs=mpfitfun( $
                 'calc_multiepoch_bf_model', $
                 x, $
                 y, $
                 yerr, $
                 ;;start, $
                 yfit=yfit, $
                 perror=perror, $
                 dof=dof, $
                 parinfo=parinfo, $
                 functargs={ $
                             bf_datavalues:bf_summed_trimmed, $
                             bf_errs:err_bf_summed_trimmed, $
                             index_epoch_segment_start_array:index_epoch_segment_start_array, $
                             index_epoch_segment_end_array:index_epoch_segment_end_array, $
                             n_params_epochspecific:n_params_epochspecific, $
                             observation_jd:observation_jd, $
                             phase_at_transit:phase_at_transit, $
                             plot_to_screen:1, $
                             sine_period:sine_period, $
                             time_of_transit:time_of_transit $
                           }, $
                 /double $
               )             ;weights=weights,parinfo=parinfo,/double)
;;4.3.2.  Look at the chisq and dof, rescale the error bars to force
;;chisq/dof=1:
chisq_preliminary=total((y-y_model_preliminary)^2/yerr^2)
print,'chisq_prelim_pass1:',chisq_preliminary,' dof:',dof
chisq_final=total((y-yfit)^2/yerr^2)
print,'chisq_final_pass1:',chisq_final,' dof:',dof
error_scaling_factor=sqrt(chisq_final/dof)
yerr=error_scaling_factor*yerr ;err_bf_summed_trimmed
chisq_adj=total((y-yfit)^2/yerr^2)
print,'chisq_adj_pass1:',chisq_adj,' dof:',dof
;;4.3.3.  Refit using the scaled error bars
;;4.3.3.1.  Restore the original starting guesses?
parinfo[*].value=start
;;4.3.3.2.  Run MPFITFUN
coeffs=mpfitfun( $
                 'calc_multiepoch_bf_model', $
                 x, $
                 y, $
                 yerr, $
                 ;;start, $
                 yfit=yfit, $
                 perror=perror, $
                 dof=dof, $
                 parinfo=parinfo, $
                 functargs={ $
                             bf_datavalues:bf_summed_trimmed, $
                             bf_errs:err_bf_summed_trimmed, $
                             index_epoch_segment_start_array:index_epoch_segment_start_array, $
                             index_epoch_segment_end_array:index_epoch_segment_end_array, $
                             n_params_epochspecific:n_params_epochspecific, $
                             observation_jd:observation_jd, $
                             phase_at_transit:phase_at_transit, $
                             plot_to_screen:1, $
                             sine_period:sine_period, $
                             time_of_transit:time_of_transit $
                           }, $
                 /double $
               ) 
;stop
;;4.3.3.3.  Plot the fit from the final coeffs, and overplot the
;;locations of the component peaks that make up the fit.
;;4.3.3.3.1  Set up postscript device if plotting to file
plot_to_file=1
if keyword_set(plot_to_file) then begin
    mydevice=!d.name
    set_plot,'ps'
    device,filename='multiepoch_bf_model.eps',/color,xsize=20.5,ysize=27.5,xoffset=0,yoffset=0;,/landscape
    charsize_legend=0.7
    spacing_legend=0.1
endif
!p.multi=[0,1,6]
y_model_final=calc_multiepoch_bf_model( $
                                  x, $
                                  coeffs, $
                                  bf_datavalues=bf_summed_trimmed, $
                                  bf_errs=yerr, $
                                  index_epoch_segment_start_array=index_epoch_segment_start_array, $
                                  index_epoch_segment_end_array=index_epoch_segment_end_array, $
                                  n_params_epochspecific=n_params_epochspecific, $
                                  observation_jd=observation_jd, $
                                  phase_at_transit=phase_at_transit, $
                                  plot_to_file=plot_to_file, $
                                  plot_to_screen=1, $
                                  sine_period=sine_period, $
                                  time_of_transit=time_of_transit $
                                )
;;4.3.3.4.  Print some results.
chisq_preliminary=total((y-y_model_preliminary)^2/yerr^2)
print,'chisq_prelim_pass2:',chisq_preliminary,' dof:',dof
chisq_final=total((y-yfit)^2/yerr^2)
print,'chisq_final_pass2:',chisq_final,' dof:',dof
error_scaling_factor=sqrt(chisq_final/dof)
yerr=error_scaling_factor*yerr
chisq_adj=total((y-yfit)^2/yerr^2)
print,'chisq_adj_pass2:',chisq_adj,' dof:',dof
for i=0L,n_elements(coeffs)-1 do begin
    print,'Coeff ',i,':',coeffs[i],', err=',perror[i]
endfor
n_coeffs=n_elements(coeffs)
print,'******Key coeffs******'
print,'Star 1 rotational velocity ',coeffs[n_coeffs-3],', err: ',perror[n_coeffs-3]
print,'Star 2 rotational velocity ',coeffs[n_coeffs-2],', err: ',perror[n_coeffs-2]
print,'BF strength ratio star2/star1',coeffs[n_coeffs-1],', err: ',perror[n_coeffs-1]
print,'Star 2 RV semiamplitude ',coeffs[n_coeffs-7],', err: ',perror[n_coeffs-7]
print,'Star 2 RV baseline relative to star 1 ',coeffs[n_coeffs-6],', err: ',perror[n_coeffs-6]
print,'******----------******'
idlastro_legend,[ $
                  'BF strength ratio star2/star1: '+strtrim(string(coeffs[n_coeffs-1],format='(f10.3)'),2)+', err: '+strtrim(string(perror[n_coeffs-1],format='(f10.3)'),2), $
                  'Star 2 RV semiamplitude: '+strtrim(string(coeffs[n_coeffs-7],format='(f10.2)'),2)+', err: '+strtrim(string(perror[n_coeffs-7],format='(f10.2)'),2), $
                  'Star 2 RV baseline relative to star 1: '+strtrim(string(coeffs[n_coeffs-6],format='(f10.2)'),2)+', err: '+strtrim(string(perror[n_coeffs-6],format='(f10.2)'),2), $
                  'Star 1 rotational velocity: '+strtrim(string(coeffs[n_coeffs-3],format='(f10.2)'),2)+', err: '+strtrim(string(perror[n_coeffs-3],format='(f10.2)'),2), $
                  'Star 2 rotational velocity: '+strtrim(string(coeffs[n_coeffs-2],format='(f10.2)'),2)+', err: '+strtrim(string(perror[n_coeffs-2],format='(f10.2)'),2) $
                ], $
                /right, $
                box=0, $
                charsize=charsize_legend, $
                spacing=spacing_legend
if keyword_set(plot_to_file) then begin
    device,/close
    set_plot,mydevice
endif
;;5.  The RV fit was constrained to lie on a sine curve (star2_rv
;;derived from the sine amplitude and offset), so there's
;;nothing further to plot.  If we had fit each peak individually, we
;;would have the additional step here of fitting a sine curve to the
;;individual RV's.

;;6.  Compute msini
;;Assume the host star is 1 solar mass and compute msini of its RV
;;companion:
;;K = (2.0*!DPI*G/P)^(1.0/3.0) * msini/mstar^(2.0/3.0) * 1.0/sqrt(1.0-e^2)
eccentricity=0.0
msun=1.9891e30 ;;kg
mstar=msun
grav_const=6.67e-11 ;;Nm^2/kg^2
period_sec=sine_period*86400.0 ;;s
k=coeffs[n_coeffs-7]*1000.0
msini=k*(2.0*!DPI*grav_const/period_sec)^(-1.0/3.0) * mstar^(2.0/3.0) * sqrt(1.0-eccentricity^2)
mjup=1.89813e27
print,'The msini of the RV companion (assuming its host is 1 solar mass) is: ',strtrim(string(msini/msun),2),' solar masses, or ',strtrim(string(msini/mjup),2),' Jupiter masses.'

stop
end
