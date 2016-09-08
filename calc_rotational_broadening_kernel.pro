;================================================================

function calc_rotational_broadening_kernel, $
                                            n_points_in_kernel=n_points_in_kernel, $
                                            equatorial_max_rv=equatorial_max_rv, $
                                            limbdarkening_beta=limbdarkening_beta, $
                                            vel_halfdomain=vel_halfdomain, $
                                            x_midpoint=x_midpoint, $
                                            velocity=velocity, $
                                            normalized_kernel=normalized_kernel
;prepares an n-element vector of total area=1 and rotationally
;broadened shape given by equatorial_max_rv and limbdarkening_beta,
;shifted to the origin [bbb.....bbb] for filter using new=cnv(old,calc_rotational_broadening_kernel)
;n_points_in_kernel=n
;equatorial_max_rv=p[0]
;limbdarkening_beta=p[1]
velocity=x_midpoint+2.0*vel_halfdomain*(dindgen(n_points_in_kernel)/double(n_points_in_kernel-1)-0.5)
normalized_coordinate=2.0*vel_halfdomain/equatorial_max_rv*(dindgen(n_points_in_kernel)/double(n_points_in_kernel-1)-0.5)
kernel=dblarr(n_elements(velocity))
kernel_temp=3.0/(3.0+2.0*limbdarkening_beta)*(2.0/!DPI*sqrt(1-normalized_coordinate^2)+limbdarkening_beta/2.0*(1.0-normalized_coordinate^2))
index_valid_region=where(normalized_coordinate le 1.0 and normalized_coordinate ge -1.0)
kernel[index_valid_region]=kernel_temp[index_valid_region]
;dv=replicate(2.0*p[0]*1.0/double(n_points_in_kernel-1),n_elements(kernel))
area=int_tabulated(velocity,kernel,/double) ;total(kernel*dv)
normalized_kernel=kernel/area
;plot,velocity,normalized_kernel,psym=3,/xs
return,1
end
