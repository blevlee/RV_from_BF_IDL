;================================================================
function cnv,a,b
;               [SMR:  10 Oct 1992]
; general convolution of two vectors or arrays with same dimensions;
; vector or array b must be symmetric at ends:  [bbb....0....bbb]
; and normalized to integral = 1 (for Gaussians, use gs_smooth for 1-D
; and gs2_smooth for 2-D
;
; usage: new_vector=cnv(vector,smoothing_vector)
; 
  return,float(fft(fft(a,-1)*fft(b,-1),+1)*float(n_elements(a)))
end
