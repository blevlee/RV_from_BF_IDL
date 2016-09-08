function sig,des,bb,sp
m = n_elements(bb[0,*])         ; size of vector
sig=fltarr(m)                   ; temporary array
pred = des ## transpose(bb)     ; predicted spectrum
for i=0,m-1 do sig[i]=sqrt(total((pred[i,*]-sp)^2)/m)
return,sig
end
