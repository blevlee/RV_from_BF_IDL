function read_spectrum, $
                        spectrum_file=spectrum_file, $
                        lambda=lambda, $
                        flux=flux, $
                        log_lun=log_lun
;;Read file.
readcol,spectrum_file,lambda,flux,comment='#'
;;Return success.
return,1
end
