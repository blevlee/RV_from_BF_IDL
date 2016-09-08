function calc_nearestneighbourinterpolated_deltafuncspec, $
                                                          linelist_fluxes=linelist_fluxes, $
                                                          linelist_lambda=linelist_lambda, $
                                                          desired_lambda=desired_lambda
;;1.  Set up variables
;;1.1  Set default flux values to continuum=0 (because of inversion in
;;my_bf_script section 1.2)
desired_flux=dblarr(n_elements(desired_lambda))
;;2.  Loop over the absorption lines in the line list
for i=0L,n_elements(linelist_lambda)-1 do begin
;;2.1  For this line, find the nearest neighbouring wavelength in the
;;desired_lambda
    dummy=min(abs(desired_lambda-linelist_lambda[i]),index_nearest_desired_lambda)
;;2.2  Assign the flux of the line to the nearest neighbouring
;;wavelength in the desired_lambda.
    desired_flux[index_nearest_desired_lambda[0]]=linelist_fluxes[i]
endfor
;;3.  Scale up to a continuum level of 10000.0 counts.
;desired_flux=10000.0*desired_flux
;stop
return,desired_flux
end
