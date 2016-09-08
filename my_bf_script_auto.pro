;;Dependencies:
;.comp BFall_IDL.pro
;.comp ../read_spectrum.pro
function my_bf_script_auto, $
                            lambda=lambda, $
                            flux=flux, $
                            linelist_lambda=linelist_lambda, $
                            linelist_flux=linelist_flux, $
                            vel=vel, $
                            bbb=bbb, $
                            bb30=bb30, $
                            window_width_pixels=window_width_pixels, $
                            wait_time_for_plots=wait_time_for_plots
;;1.0  Set up variables
if ~keyword_set(window_width_pixels) then window_width_pixels=201
;;1.1  Read data
wp1=lambda
sp1=flux
ws=linelist_lambda
ss=linelist_flux
;;1.2  Invert spectra to reduce edge effects
ss=1.-ss/10000.0
sp1=1.-sp1
;;1.3  Choose wavelength range to operate on
index_finite_flux=where(finite(sp1),count_finite_flux)
if count_finite_flux lt 10 then print,'There are very few valid flux values in the spectrum.  Stopping.'
if count_finite_flux lt 10 then stop
;user_min=5290;5370.0
;user_max=5390;5430.0
user_min=0.0
user_max=99999.9
lambda_min=max([min(wp1[index_finite_flux]),min(ws),user_min])
lambda_max=min([max(wp1[index_finite_flux]),max(ws),user_max])
;lambda_min=5290.0
;lambda_max=5390.0
print,'***********',lambda_min,lambda_max
;;2.1  Set up parameters for resampling log lambda
;velocity_resolution=0.5 ;km/s
velocity_maximum=2.997924d5*((lambda_max-lambda_min)/lambda_min)
n_resampled_pixels=n_elements(wp1)
if ~keyword_set(velocity_resolution) then velocity_resolution=(velocity_maximum/double(n_resampled_pixels)) ;km/s
r=velocity_resolution/2.997924d5
w1=lambda_min*(1.d0+r)^dindgen(n_resampled_pixels)
;;2.2  Do the resampling
;ssr=interpol(ss,ws,w1)
ssr=calc_nearestneighbourinterpolated_deltafuncspec(linelist_fluxes=ss,linelist_lambda=ws,desired_lambda=w1)
sp1r=interpol(sp1,wp1,w1)
plot,w1,ssr
wait,wait_time_for_plots
;;2.3  Check for bad values at ends
index_finite=where(finite(sp1r),count_finite)
if count_finite eq 0 then begin
    print,systime(/UTC),'|ERROR|my_bf_script|There were no finite values in the program spectrum.'
    stop
endif
;;2.4  Clip bad values
w1=w1[index_finite]
ssr=ssr[index_finite]
sp1r=sp1r[index_finite]
;;2.5  Ensure arrays have an even number of elements
if n_elements(index_finite) mod 2 eq 1 then begin
    w1=w1[0:n_elements(w1)-2]
    ssr=ssr[0:n_elements(ssr)-2]
    sp1r=sp1r[0:n_elements(sp1r)-2]
endif
;;3.1  CCF plot
lag = findgen(window_width_pixels)-floor(window_width_pixels/2)
ccf1=c_correlate(ssr,sp1r,lag)
vel=lag*velocity_resolution
!p.multi=[0,1,3]
plot,vel,ccf1
wait,wait_time_for_plots
;stop
;;4.1
print,systime(/utc),'|my_bf_script.pro|Creating des'
des=map4(ssr,window_width_pixels)
print,systime(/utc),'|my_bf_script.pro|Doing svdc'
svdc,des,ww,u,v,/double
plot_io,ww
;;5.1
sp1rt=sp1r[window_width_pixels/2:n_elements(sp1r)-1-window_width_pixels/2]
broadfunct,sp1rt,ww,u,v,window_width_pixels,bb1
vel=velocity_resolution*(findgen(window_width_pixels)-floor(window_width_pixels/2))
plot,vel,bb1[*,65]
;;6.1
!p.multi=[0,1,1]
sig1=sig(des,bb1,sp1rt)
plot,sig1,yran=[0,1.1*max(sig1)]
;bbb=reform(bb1[*,30])
bbb=reform(bb1[*,window_width_pixels-1])
!p.multi=[0,1,2]
plot,vel,bbb
;;7.1
gs30=gs_smooth(window_width_pixels,3.0)
bb30=cnv(reform(bb1[*,window_width_pixels-1]),gs30)
plot,vel,bb30
;plot,lambda,flux,xrange=[5320,5325]
;plot,lambda,flux,xrange=[5295,5300]
;;plot,vel,bbb,title='TYC 3559, APO3.5+ARCES, BF using BASS2000 solar template',xtitle='RV (km/s) (not pre-corrected for baryvel)',ytitle='Broadening function (bb30)',charsize=1.5,thick=3,charthick=2,xrange=[-150,150],/xs

status=1
return,status
end
