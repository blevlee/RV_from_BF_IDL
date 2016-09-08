; All BF routines together
; March 2006: Added routines for rotational profiles, Rot*.pro and 
; (hopefully) all needed IDL-astro routines
;
; Order is not exactly alphabetical!
;=============================================================================
pro baryvel, dje, deq, dvelh, dvelb
;+
; NAME:
;       BARYVEL
; PURPOSE:
;       Calculates heliocentric and barycentric velocity components of Earth.
;
; EXPLANATION:
;       BARYVEL takes into account the Earth-Moon motion, and is useful for 
;       radial velocity work to an accuracy of  ~1 m/s.
;
; CALLING SEQUENCE:
;       BARYVEL, dje, deq, dvelh, dvelb
;
; INPUTS:
;       DJE - (scalar) Julian ephemeris date.
;       DEQ - (scalar) epoch of mean equinox of dvelh and dvelb. If deq=0
;               then deq is assumed to be equal to dje.
; OUTPUTS: 
;       DVELH: (vector(3)) heliocentric velocity component. in km/s 
;       DVELB: (vector(3)) barycentric velocity component. in km/s
;
;       The 3-vectors DVELH and DVELB are given in a right-handed coordinate 
;       system with the +X axis toward the Vernal Equinox, and +Z axis 
;       toward the celestial pole.      
;
; PROCEDURE CALLED:
;       Function PREMAT() -- computes precession matrix
;
; NOTES:
;       Algorithm taken from FORTRAN program of Stumpff (1980, A&A Suppl, 41,1)
;       Stumpf claimed an accuracy of 42 cm/s for the velocity.    A 
;       comparison with the JPL FORTRAN planetary ephemeris program PLEPH
;       found agreement to within about 65 cm/s between 1986 and 1994
;
; EXAMPLE:
;       Compute the radial velocity of the Earth toward Altair on 15-Feb-1994
;
;       IDL> jdcnv, 1994, 2, 15, 0, jd          ;==> JD = 2449398.5
;       IDL> baryvel, jd, 2000, vh, vb          
;               ==> vh = [-17.07809, -22.80063, -9.885281]  ;Heliocentric km/s
;               ==> vb = [-17.08083, -22.80471, -9.886582]  ;Barycentric km/s
;
;       IDL> ra = ten(19,50,46.77)*15/!RADEG    ;RA  in radians
;       IDL> dec = ten(08,52,3.5)/!RADEG        ;Dec in radians
;       IDL> v = vb(0)*cos(dec)*cos(ra) + $   ;Project velocity toward star
;               vb(1)*cos(dec)*sin(ra) + vb(2)*sin(dec) 
;
; REVISION HISTORY:
;       Jeff Valenti,  U.C. Berkeley    Translated BARVEL.FOR to IDL.
;       W. Landsman, Cleaned up program sent by Chris McCarthy (SfSU) June 1994
;       Converted to IDL V5.0   W. Landsman   September 1997
;-
 On_Error,2

 if N_params() LT 4 then begin
        print,'Syntax: BARYVEL, dje, deq, dvelh, dvelb'
        print,'    dje - input Julian ephemeris date'
        print,'    deq - input epoch of mean equinox of dvelh and dvelb'
        print,'    dvelh - output vector(3) heliocentric velocity comp in km/s' 
        print,'    dvelb - output vector(3) barycentric velocity comp in km/s'
        return
 endif

;Define constants
  dc2pi = 2*!DPI 
  cc2pi = 2*!PI 
  dc1 = 1.0D0
  dcto = 2415020.0D0
  dcjul = 36525.0D0                     ;days in julian year
  dcbes = 0.313D0
  dctrop = 365.24219572D0               ;days in tropical year (...572 insig)
  dc1900 = 1900.0D0
  AU = 1.4959787D8

;Constants dcfel(i,k) of fast changing elements.
  dcfel = [1.7400353D00, 6.2833195099091D02,  5.2796D-6 $
          ,6.2565836D00, 6.2830194572674D02, -2.6180D-6 $
          ,4.7199666D00, 8.3997091449254D03, -1.9780D-5 $
          ,1.9636505D-1, 8.4334662911720D03, -5.6044D-5 $
          ,4.1547339D00, 5.2993466764997D01,  5.8845D-6 $
          ,4.6524223D00, 2.1354275911213D01,  5.6797D-6 $
          ,4.2620486D00, 7.5025342197656D00,  5.5317D-6 $
          ,1.4740694D00, 3.8377331909193D00,  5.6093D-6 ]
  dcfel = reform(dcfel,3,8)

;constants dceps and ccsel(i,k) of slowly changing elements.
  dceps = [4.093198D-1, -2.271110D-4, -2.860401D-8 ]
  ccsel = [1.675104E-2, -4.179579E-5, -1.260516E-7 $
          ,2.220221E-1,  2.809917E-2,  1.852532E-5 $
          ,1.589963E00,  3.418075E-2,  1.430200E-5 $
          ,2.994089E00,  2.590824E-2,  4.155840E-6 $
          ,8.155457E-1,  2.486352E-2,  6.836840E-6 $
          ,1.735614E00,  1.763719E-2,  6.370440E-6 $
          ,1.968564E00,  1.524020E-2, -2.517152E-6 $
          ,1.282417E00,  8.703393E-3,  2.289292E-5 $
          ,2.280820E00,  1.918010E-2,  4.484520E-6 $
          ,4.833473E-2,  1.641773E-4, -4.654200E-7 $
          ,5.589232E-2, -3.455092E-4, -7.388560E-7 $
          ,4.634443E-2, -2.658234E-5,  7.757000E-8 $
          ,8.997041E-3,  6.329728E-6, -1.939256E-9 $
          ,2.284178E-2, -9.941590E-5,  6.787400E-8 $
          ,4.350267E-2, -6.839749E-5, -2.714956E-7 $
          ,1.348204E-2,  1.091504E-5,  6.903760E-7 $
          ,3.106570E-2, -1.665665E-4, -1.590188E-7 ]
  ccsel = reform(ccsel,3,17)

;Constants of the arguments of the short-period perturbations.
  dcargs = [5.0974222D0, -7.8604195454652D2 $
           ,3.9584962D0, -5.7533848094674D2 $
           ,1.6338070D0, -1.1506769618935D3 $
           ,2.5487111D0, -3.9302097727326D2 $
           ,4.9255514D0, -5.8849265665348D2 $
           ,1.3363463D0, -5.5076098609303D2 $
           ,1.6072053D0, -5.2237501616674D2 $
           ,1.3629480D0, -1.1790629318198D3 $
           ,5.5657014D0, -1.0977134971135D3 $
           ,5.0708205D0, -1.5774000881978D2 $
           ,3.9318944D0,  5.2963464780000D1 $
           ,4.8989497D0,  3.9809289073258D1 $
           ,1.3097446D0,  7.7540959633708D1 $
           ,3.5147141D0,  7.9618578146517D1 $
           ,3.5413158D0, -5.4868336758022D2 ]
  dcargs = reform(dcargs,2,15)

;Amplitudes ccamps(n,k) of the short-period perturbations.
  ccamps = $
    [-2.279594E-5,  1.407414E-5,  8.273188E-6,  1.340565E-5, -2.490817E-7 $
    ,-3.494537E-5,  2.860401E-7,  1.289448E-7,  1.627237E-5, -1.823138E-7 $
    , 6.593466E-7,  1.322572E-5,  9.258695E-6, -4.674248E-7, -3.646275E-7 $
    , 1.140767E-5, -2.049792E-5, -4.747930E-6, -2.638763E-6, -1.245408E-7 $
    , 9.516893E-6, -2.748894E-6, -1.319381E-6, -4.549908E-6, -1.864821E-7 $
    , 7.310990E-6, -1.924710E-6, -8.772849E-7, -3.334143E-6, -1.745256E-7 $
    ,-2.603449E-6,  7.359472E-6,  3.168357E-6,  1.119056E-6, -1.655307E-7 $
    ,-3.228859E-6,  1.308997E-7,  1.013137E-7,  2.403899E-6, -3.736225E-7 $
    , 3.442177E-7,  2.671323E-6,  1.832858E-6, -2.394688E-7, -3.478444E-7 $
    , 8.702406E-6, -8.421214E-6, -1.372341E-6, -1.455234E-6, -4.998479E-8 $
    ,-1.488378E-6, -1.251789E-5,  5.226868E-7, -2.049301E-7,  0.E0 $
    ,-8.043059E-6, -2.991300E-6,  1.473654E-7, -3.154542E-7,  0.E0 $
    , 3.699128E-6, -3.316126E-6,  2.901257E-7,  3.407826E-7,  0.E0 $
    , 2.550120E-6, -1.241123E-6,  9.901116E-8,  2.210482E-7,  0.E0 $
    ,-6.351059E-7,  2.341650E-6,  1.061492E-6,  2.878231E-7,  0.E0 ]
  ccamps = reform(ccamps,5,15)

;Constants csec3 and ccsec(n,k) of the secular perturbations in longitude.
  ccsec3 = -7.757020E-8
  ccsec = [1.289600E-6, 5.550147E-1, 2.076942E00 $
          ,3.102810E-5, 4.035027E00, 3.525565E-1 $
          ,9.124190E-6, 9.990265E-1, 2.622706E00 $
          ,9.793240E-7, 5.508259E00, 1.559103E01 ]
  ccsec = reform(ccsec,3,4)

;Sidereal rates.
  dcsld = 1.990987D-7                   ;sidereal rate in longitude
  ccsgd = 1.990969E-7                   ;sidereal rate in mean anomaly

;Constants used in the calculation of the lunar contribution.
  cckm = 3.122140E-5
  ccmld = 2.661699E-6
  ccfdi = 2.399485E-7

;Constants dcargm(i,k) of the arguments of the perturbations of the motion
; of the moon.
  dcargm = [5.1679830D0,  8.3286911095275D3 $
           ,5.4913150D0, -7.2140632838100D3 $
           ,5.9598530D0,  1.5542754389685D4 ]
  dcargm = reform(dcargm,2,3)

;Amplitudes ccampm(n,k) of the perturbations of the moon.
  ccampm = [ 1.097594E-1, 2.896773E-7, 5.450474E-2,  1.438491E-7 $
           ,-2.223581E-2, 5.083103E-8, 1.002548E-2, -2.291823E-8 $
           , 1.148966E-2, 5.658888E-8, 8.249439E-3,  4.063015E-8 ]
  ccampm = reform(ccampm,4,3)

;ccpamv(k)=a*m*dl,dt (planets), dc1mme=1-mass(earth+moon)
  ccpamv = [8.326827E-11, 1.843484E-11, 1.988712E-12, 1.881276E-12]
  dc1mme = 0.99999696D0

;Time arguments.
  dt = (dje - dcto) / dcjul
  tvec = [1d0, dt, dt*dt]

;Values of all elements for the instant(aneous?) dje.
  temp = (tvec # dcfel) mod dc2pi
  dml = temp[0]
  forbel = temp[1:7]
  g = forbel[0]                         ;old fortran equivalence

  deps = total(tvec*dceps) mod dc2pi
  sorbel = (tvec # ccsel) mod dc2pi
  e = sorbel[0]                         ;old fortran equivalence

;Secular perturbations in longitude.
dummy=cos(2.0)
  sn = sin((tvec[0:1] # ccsec[1:2,*]) mod cc2pi)

;Periodic perturbations of the emb (earth-moon barycenter).
  pertl = total(ccsec[0,*] * sn) + dt*ccsec3*sn[2]
  pertld = 0.0
  pertr = 0.0
  pertrd = 0.0
  for k=0,14 do begin
    a = (dcargs[0,k]+dt*dcargs[1,k]) mod dc2pi
    cosa = cos(a)
    sina = sin(a)
    pertl = pertl + ccamps[0,k]*cosa + ccamps[1,k]*sina
    pertr = pertr + ccamps[2,k]*cosa + ccamps[3,k]*sina
    if k lt 11 then begin
      pertld = pertld + (ccamps[1,k]*cosa-ccamps[0,k]*sina)*ccamps[4,k]
      pertrd = pertrd + (ccamps[3,k]*cosa-ccamps[2,k]*sina)*ccamps[4,k]
    endif
  endfor

;Elliptic part of the motion of the emb.
  phi = (e*e/4d0)*(((8d0/e)-e)*sin(g) +5*sin(2*g) +(13/3d0)*e*sin(3*g))
  f = g + phi
  sinf = sin(f)
  cosf = cos(f)
  dpsi = (dc1 - e*e) / (dc1 + e*cosf)
  phid = 2*e*ccsgd*((1 + 1.5*e*e)*cosf + e*(1.25 - 0.5*sinf*sinf))
  psid = ccsgd*e*sinf / sqrt(dc1 - e*e)

;Perturbed heliocentric motion of the emb.
  d1pdro = dc1+pertr
  drd = d1pdro * (psid + dpsi*pertrd)
  drld = d1pdro*dpsi * (dcsld+phid+pertld)
  dtl = (dml + phi + pertl) mod dc2pi
  dsinls = sin(dtl)
  dcosls = cos(dtl)
  dxhd = drd*dcosls - drld*dsinls
  dyhd = drd*dsinls + drld*dcosls

;Influence of eccentricity, evection and variation on the geocentric
; motion of the moon.
  pertl = 0.0
  pertld = 0.0
  pertp = 0.0
  pertpd = 0.0
  for k = 0,2 do begin
    a = (dcargm[0,k] + dt*dcargm[1,k]) mod dc2pi
    sina = sin(a)
    cosa = cos(a)
    pertl = pertl + ccampm[0,k]*sina
    pertld = pertld + ccampm[1,k]*cosa
    pertp = pertp + ccampm[2,k]*cosa
    pertpd = pertpd - ccampm[3,k]*sina
  endfor

;Heliocentric motion of the earth.
  tl = forbel[1] + pertl
  sinlm = sin(tl)
  coslm = cos(tl)
  sigma = cckm / (1.0 + pertp)
  a = sigma*(ccmld + pertld)
  b = sigma*pertpd
  dxhd = dxhd + a*sinlm + b*coslm
  dyhd = dyhd - a*coslm + b*sinlm
  dzhd= -sigma*ccfdi*cos(forbel[2])

;Barycentric motion of the earth.
  dxbd = dxhd*dc1mme
  dybd = dyhd*dc1mme
  dzbd = dzhd*dc1mme
  for k=0,3 do begin
    plon = forbel[k+3]
    pomg = sorbel[k+1]
    pecc = sorbel[k+9]
    tl = (plon + 2.0*pecc*sin(plon-pomg)) mod cc2pi
    dxbd = dxbd + ccpamv[k]*(sin(tl) + pecc*sin(pomg))
    dybd = dybd - ccpamv[k]*(cos(tl) + pecc*cos(pomg))
    dzbd = dzbd - ccpamv[k]*sorbel[k+13]*cos(plon - sorbel[k+5])

  endfor

;Transition to mean equator of date.
  dcosep = cos(deps)
  dsinep = sin(deps)
  dyahd = dcosep*dyhd - dsinep*dzhd
  dzahd = dsinep*dyhd + dcosep*dzhd
  dyabd = dcosep*dybd - dsinep*dzbd
  dzabd = dsinep*dybd + dcosep*dzbd

;Epoch of mean equinox (deq) of zero implies that we should use
; Julian ephemeris date (dje) as epoch of mean equinox.
  if deq eq 0 then begin
    dvelh = AU * ([dxhd, dyahd, dzahd])
    dvelb = AU * ([dxbd, dyabd, dzabd])
    return
  endif

;General precession from epoch dje to deq.
  deqdat = (dje-dcto-dcbes) / dctrop + dc1900
   prema = premat(deqdat,deq, /FK4)

  dvelh = AU * ( prema # [dxhd, dyahd, dzahd] )
  dvelb = AU * ( prema # [dxbd, dyabd, dzabd] )

  return
  end
;================================================================
pro BFanal3,vel,bf
; simplified analysis of BF's of the contribution of the 3rd
;       component in terms of luminosities
; usage:   BFanal3,i4vel,i415[48,*]

bf=reform(bf)           ; IDL book-keeping
n=n_elements(vel)

set_win
plot,vel,bf

print,'Subtract the baseline:'
print,'   mark two ends of the whole broadened profile'
cursor,x1,y1,/data,/down
oplot,[x1,x1],[y1,y1],psym=6
cursor,x2,y2,/data,/down
oplot,[x2,x2],[y2,y2],psym=6
n1=where(vel ge x1)
n1=n1[0]
n2=where(vel ge x2)
n2=n2[0]
base=interpol([y1,y2],[x1,x2],vel)
oplot,vel,base,line=2
int1=total(bf[n1:n2]-base[n1:n2])

print,'Mark two ends of the 3rd component'
cursor,x3,y3,/data,/down
oplot,[x3,x3],[y3,y3],psym=4
cursor,x4,y4,/data,/down
oplot,[x4,x4],[y4,y4],psym=4
n3=where(vel ge x3)
n3=n3[0]
n4=where(vel ge x4)
n4=n4[0]
base1=interpol([y3,y4],[x3,x4],vel)
oplot,vel[n3:n4],base1[n3:n4],line=1
int2=total(bf[n3:n4]-base1[n3:n4])

print,'Whole profile: ',int1
print,'3rd component: ',int2
print,'Binary:        ',int1-int2
print,'L3/(L1+L2):    ',int2/(int1-int2)

return
end
;================================================================
pro BFanal3x,vel,bf,bf_2
; contribution of the 3rd component in terms of luminosities
;    evaluated from the 3-Gaussian fits
; uses the original BF and the one with the 3rd Gaussian
;    subtracted
; usage:   i=5 & BFanal3x,vel,bf15[i,*],bf15_2[i,*]

bf=reform(bf)           ; IDL book-keeping
n=n_elements(vel)

set_win
plot,vel,bf
oplot,vel,bf_2,line=1

print,'Subtract the baseline:'
print,'   mark two ends of the whole broadened profile'
cursor,x1,y1,/data,/down
oplot,[x1,x1],[y1,y1],psym=6
cursor,x2,y2,/data,/down
oplot,[x2,x2],[y2,y2],psym=6
n1=where(vel ge x1)
n1=n1[0]
n2=where(vel ge x2)
n2=n2[0]
base=interpol([y1,y2],[x1,x2],vel)
oplot,vel,base,line=2
int1=total(bf[n1:n2]-base[n1:n2])

int2=total(bf[n1:n2]-bf_2[n1:n2])  ; simply the difference

print,'Whole profile: ',int1
print,'3rd component: ',int2
print,'Binary:        ',int1-int2
print,'L3/(L1+L2):    ',int2/(int1-int2)

return
end
;================================================================
pro BFpro1,std_fts,w00,n,m,stepV,blank,w1,des,ww,u,v,vel
; processing of the standard (template) spectrum
; ver.March 2005
; input:  std_fts  = name of file with std spec in FITS, with .ext[ension]
;         w00  = starting wave (A) of the log-wave vector, select
;         n    = desired length of the log-wave vector in pix 
;                n must be EVEN, for our spectra typically 1000
;         m    = desired length of the BF, must be ODD number of pixels,
;                typically 111, 121, 131, etc.
;         stepV = step in velocities in the wavelength vector w1
;                for DDO use 11.8 km/s, new 6.5 km/s
;         blank = option to remove a section of the spectrum
;                set [0.,0.] if nothing to remove
; output: w1   = log-wave wavelenghts, for use in BFpro2.pro
;         des  = design array 
;         ww   = singular value vector 
;         u,v  = auxilary arrays
;         vel  = velocity vector, x-axis for broadening functions
;
; usage: 
;
;  old usage without "blank"
;  BFpro1,'rC0017558.fts',5080.0,1016,121,11.8,w1,des,ww,u,v,vel
;               std        w00     n   m stepV w1 des ww u v vel
;
;  new usage
;       March 2004, eg. JY2 centered at 6290A
;  stellar spectrum with the telluric O2 band 6275-6330A removed
;  BFpro1,'K0003000.fits',6165.,1946,131,6.5,[6275.,6330.],w1,des,ww,u,v,vel
;       avoid end problems by using only 6165-6430A
;
;  telluric lines for checking flexure
;  BFpro1,'K0003000.fits',6285.,328,21,6.5,[0.,0.],w1,des,ww,u,v,vel


sts=readfits(std_fts,hs,/silent)     ; reads in FITS 
;parse_head,hs,w0,dw,nn               ; extracts from header
parse_fits_head,hs,nn,w0,dw
                                     ; CRVAL1 => w0
                                     ; CDELT1 => dw
                                     ; NAXIS1 => nn
ans=''
set_win
plot,sts,tit='Whole spectrum plotted. Press any key...'
read,'Press any key...',ans

ss=1.-sts                         ; spectrum inverted, lines are positive
w=w0+dw*findgen(nn)               ; wavelength vector re-created
    
r = stepV/2.997924d5              ; old at 5184A 11.8 km/s/pix, 
                                  ; new at 6290A  6.5 km/s/pix
w1 = w00 * (1.d0+r)^dindgen(n)    ; log-wave vect, equal spacing in vel

; make sure that the rebinned vector is within the original wavelengths
; the program spectra may have their own wavelengths shifted relative 
; to the standard spectrum

print,form='(a,2f10.4)','Original wavelength range: ',w0,w0+dw*(nn-1)
print,form='(a,2f10.4)','Rebinned wavelength range: ',w1[0],w1[n-1]
;read,'Acceptable? (y=any/n)',ans
;if (ans eq 'n') or (ans eq 'N') then goto,E

ssr = interpol(ss,w,w1)              ; spectrum resampled to log-wave
vel=stepV*(findgen(m)-m/2)           ; velocity vector

w=where((w1 lt blank[0]) or (w1 gt blank[1]))  ; telluric band
if w[0] eq -1 then goto,N     ; if nothing to remove
   nw=n_elements(w)
   if (nw mod 2) eq 1 then w=w[0:nw-2] ; checking if size is even 
   
   ssr=ssr[w]                 ; spectrum spliced together
   w1=w1[w]

N:
plot,w1,ssr,tit='Rebinned, positive-spike spectrum'
read,'Acceptable? (y=any/n)',ans

print,'Now be patient, the SVD steps can be slow...'
des = map4(ssr,m)                    ; design array created 
print,'Design array made...'
svdc,des,ww,u,v,/double              ; SVD decomposition

plot_io,ww,tit='!17 Singular values'
print,'BFpro1 finished.'

E:
return
end
;================================================================
pro BFpro2,prg_lst,prg_lst1,w1,ww,u,v,images,spec,bf
; processing of the program spectra, derivation of BF for all spectra
;    as a 2-D array: bf[phase-indx, vel-indx]
; ver March 2005
; this version (May 2002) includes rejection of poor spectra
;    and permits an abort when things go wrong
;
; input:  prg_lst = list of program-star spectrum FITS files
;             (in Win use DOS window: dir *.FTS /B > star.lst)
;         w1   = log-wave vector, calculated with BFpro1.pro
;         ww   = singular values, calculated with BFpro1.pro
;         u,v  = u,v  = auxilary arrays, calculated with BFpro1.pro
; output: 
;         prg_lst1 = new list with some poor spectra rejected
;         images = string array duplicating names, as a check
;         spec = spectra, just in case, not really used
;         bf   = full BF, normally must be smoothed, use BFpro3.pro
;
; usage: BFpro2,'V2357Oph.lst','V2357Oph1.lst',w1,ww,u,v,images,spec,bf
;                  prg list       new list     w1 ww u v images spec bf

ans=''

prgs=strarr(200)                   ; read the list of files
openr,1,prg_lst
i=0
a=''
while(EOF(1) ne 1) do begin
    readf,1,form='(a)',a
    prgs[i]=a
    i = i + 1
endwhile
prgs=prgs[0:i-1]
close,1
nnn=i                                  ; nnn=number of observations

prgs1=prgs                             ; copy of list

m = n_elements(ww)                     ; size of BF vectors
bf=fltarr(nnn,m)                       ; 2-D BF array prepared
n = n_elements(w1)                     ; size of w1 (log-wave)

images=strarr(nnn)                     ; file names for spectra
spec=fltarr(nnn,n)                     ; de-spiked spectra saved
bad=intarr(nnn)                        ; bad[i]=1 => removed spectrum

;-------------------------------------------------------------------
print,'Carriage-return => no change (n) in what follows'
; cycle through all observations including the template at the top

set_win

for i=0,nnn-1 do begin                 ; cycle through spectra

  prg_fts=strtrim(prgs[i],2)
  images[i]=prg_fts                    ; in FITS
  prg=readfits(prg_fts,hp,/silent) 
  print,'Image: ',images[i]
;  parse_head,hp,w0,dw,nn               ; extracts from header
  parse_fits_head,hp,nn,w0,dw          ; CRVAL1 => w0
                                       ; CDELT1 => dw
                                       ; NAXIS1 => nn
  w=w0+dw*findgen(nn)   ; wavelength vector re-created for each
  prg1=prg    ; auxiliary copy

; this step not needed if ends removed by choice of w1 ends
;C2:
;  plot,prg1,tit='!17'+prg_fts 
;  read,'Gross end-trimming? (y/n,<Ret>=n)',ans
;  if (ans eq 'y') or (ans eq 'Y') then begin
;     de_spike2,prg1,prg1,w,nn
;     goto, C2
;  end
;  prg=prg1

C1:
  plot,w,prg1,tit='!17'+prg_fts
  oplot,w1,replicate((min(prg1)+max(prg1))/2,n_elements(w1)),psym=3
  read,'De-spiking? (y/n,<Ret>=n)',ans
  if (ans eq 'y') or (ans eq 'Y') then begin
     de_spike1,prg1,prg1,60 
     goto, C1
  end
  prg=prg1

  read,'Abort sequence(a)! Reject spectrum(y)? (a/y/n,<Ret>=n)',ans
  if (ans eq 'a') or (ans eq 'A') then begin
     print,'ABORTED SEQUENCE'
     goto, E
  end 
  if (ans eq 'y') or (ans eq 'Y') then begin
     print,prg_fts,' rejected!'
     bad[i]=1
  end

  sp = 1. - prg              ; spectrum inverted

  if w1[0] lt w[0] then begin
      print,'Short wavelength end incorrect'
      print,'w1[0]: ',w1[0],'  w[0]: ',w[0],'  w1[0]-w[0]: ',w1[0]-w[0]
;      read,'Accept (y/n)',ans
;      if (ans eq 'n') or (ans eq 'N') then goto,E
  end

  if w1[n-1] gt w[nn-1] then begin
      print,'Long wavelength end incorrect'
      print,'w1[n-1]: ',w1[n-1],'  w[nn-1]: ',w[nn-1], $
                                 '  w1[n-1]-w[nn-1]: ',w1[n-1]-w[nn-1]
;      read,'Accept (y/n)',ans
;      if (ans eq 'n') or (ans eq 'N') then goto,E
  end

  spr = interpol(sp,w,w1)    ; spectrum resampled to log-wave         
  spec[i,*]=1.-spr           ; spectrum saved, just in case
  sprt = trunc(spr,m,n)      ; truncation 
  bf1 = fltarr(m)            ; broadening function declared
  bf1=svsol(u,ww,v,sprt,/double) ; full solution, no singular values removed
  bf1=float(bf1)             ; back to single precision
  bf[i,*]=bf1                ; BF array filled for phase index i

endfor
; end cycle through spectra ----------------------------------------------

E:
w=where(bad eq 0)                 ; good spectra retained
images=images[w]   
spec=spec[w,*]
bf=bf[w,*]
prgs1=prgs[w]

openw,1,prg_lst1                  ; new cleaned list created
printf,1,f='(a)',prgs1
close,1

print,'BFpro2 finished'
return
end
;=======================================================================
pro BFpro2b,std_fts,m,stepV,w1,spec,des,ww,u,v,vel,bf
; derivation of BF for all spectra already extracted and de-spiked
;    with BFpro2.pro; any bad have been already rejected
; use it for cases when a change in the template requires
;    re-determination of BF's or you want to change the size of BF
; the SVD step as in BFpro1.pro is repeated, but for the same
;    wavelength vector w1 as before (with the same rejection of some
;    interval "blank", if applicable)
; 
; Ver March 2006
;
; input:  
;        std_fts = full name of the template spectrum (file with ext)
;        m    = size of the BF in pix
;        stepV = step in km/s
;        w1   = log-wave vector, calculated with BFpro1.pro and kept
;               the same for this processing
;        spec = old spectra processed and stored by BFpro2.pro
; output:
;        des  = design array for the new template 
;        ww   = new singular values, calculated as with BFpro1.pro
;        u,v  = new auxilary arrays, calculated as with BFpro1.pro
;        vel  = new vector of velocities
;        bf   = new full BF, normally must be smoothed, use BFpro3.pro
; Note:  the zero-th element of BF contains the old template against
;        the new template
; 
; usage: BFpro2b,'K0003000.fits',111,6.5,w1,spec,des,ww,u,v,vel,bf_n
;                  std_fits       m  stepV w1 spec des ww u v vel bf 
;        spec -> old spectra, the template spectrum not changed
;        n_BF's -> new BF's
;        des,ww,u,v,vel -> all new; you may want to use different names
;
;
; New template -----------------------------------------------------------

sts=readfits(std_fts,hs,/silent)  ; reads in FITS 
parse_fits_head,hs,nn,w0,dw
                                  ; CRVAL1 => w0
                                  ; CDELT1 => dw
                                  ; NAXIS1 => nn
w=w0+dw*findgen(nn)               ; wavelength vector re-created

ans=''
set_win
plot,w,sts,tit='Whole spectrum plotted. Press any key...'
read,'Acceptable? (y=any/n)',ans
ss=1.-sts                       ; spectrum inverted, lines are positive
ssr = interpol(ss,w,w1)              
;plot,w1,ssr,tit='Rebinned, positive-spike spectrum'

print,'Now be patient, the SVD steps may be slow...'
des = map4(ssr,m)                    ; design array created 
print,'Design array made...'
svdc,des,ww,u,v,/double              ; SVD decomposition

plot_io,ww,tit='!17 Singular values'
vel=stepV*(findgen(m)-m/2)           ; velocity vector

; template done, now cycle through spectra -------------------------------

nnn=n_elements(spec[*,0])    ; number of good spectra 
n=n_elements(spec[0,*])      ; number of pixels in spectra
bf=fltarr(nnn,m)

;spec[0,*]=1.-ssr             ; the template spectrum is NOT
                              ; replaced, all spec kept the same
for i=0,nnn-1 do begin                 
  spr = reform(1.-spec[i,*]) ; positive spikes restored, resampled spec
  sprt = trunc(spr,m,n)      ; truncation 
  bf1 = fltarr(m)            ; broadening function declared
  bf1=svsol(u,ww,v,sprt,/double) ; full solution, no singular values removed
  bf1=float(bf1)             ; back to single precision
  bf[i,*]=bf1                ; BF array filled for phase index i
endfor
; end cycle through spectra ----------------------------------------------

print,'BFpro2b finished'
return
end
;=======================================================================
function BFpro3,bf,sig
; 
; Smoothing of the broadening functions for all phases
;   with arbitrary sigma; adjust sigma to have FWHM corresponding
;   to the projected slit width in pixels (normally 2-3 pix)
;
; input:  bf  = broadening functions calculated with BFpro2.pro 
;         sig = sigma for Gaussian
;                  sigma     FWHM    (both in pixels)
;                   0.5      1.18
;                   0.75     1.68
;                   1.0      2.35
;                   1.5      3.53
;                   2.0      4.7
;                   3.0      7.06
;     generally: FWHM = 2.354 sigma; sigma = 0.4248 FWHM
;
; output = Gaussian smoothed broadening functions with
;              different smoothing measured by dispersion sigma
;         
; usage:  
;          bf30=BFpro3(bf,1.5)

  m = n_elements(bf[*,0])       ; m = number of phases
  n = n_elements(bf[0,*])       ; n = number of pixels in indiv BF

  res=bf

  gs = gs_smooth(n,sig)
  for j=0,m-1 do res[j,*] = cnv(bf[j,*],gs)

return,res
end

;================================================================
pro BFpro4,prg_lst,w00,n,stepV
;
; preparatory processing of all spectra to check wavelength vector sizes;
; run it before any other routine
; ver. March 2005
; input:  prg_fts  = list of files with spectra in FITS format
;         w00  = starting wave (A) of the log-wave vector, select
;         n    = desired length of the log-wave vector  
;                in pixels n must be EVEN, 
;         for old spectra typically 1000, for new typically <2000
;         stepV = step in velocities in the wavelength vector w1
;                for DDO old chip use 11.8 km/s 
;                        new chip     6.5 at 6290A
; output: no output, only screen prinout
;
; the goal is to have only two columns, by adjustment of [w00,n] pair;
; any differences appearing in next columns indicate that the 
;         parameters are incorrect
;
; usage: BFpro4,'V523Cas1.lst',5079.7,1018,11.8
;                for the old chip at 5184A
; new
;        BFpro4,'star.lst',6165.0,1956,6.5
;                new chip at 6190A
;                to avoid shutter problems at ends, use 6165-6430A

r = stepV/2.997924d5                 ; 11.8 km/s/pix
w1 = w00 * (1.d0+r)^dindgen(n)       ; log-wave vect, equal spacing in vel

prgs=strarr(200)                     ; read the list of files
openr,1,prg_lst                      ; for processing, typically 'star.lst'
i=0
a=''
while(EOF(1) ne 1) do begin
    readf,1,form='(a)',a
    prgs[i]=a
    i = i + 1
endwhile
prgs=prgs[0:i-1]
close,1
nnn=i

print,form='(a,3f12.3)','Equal-velocity vector:',w1[0],w1[n-1],w1[1]-w1[0]

for i=0,nnn-1 do begin                  
  prg_fts=strtrim(prgs[i],2)           ; headers of individual spectra
  prg=readfits(prg_fts,hp,/silent)     ; analyzed
  ppp=prg_fts
  parse_head,hp,w0,dw,nn               ; extracts from header
                                       ; CRVAL1 => w0
                                       ; CDELT1 => dw (same as CD1_1)
                                       ; NAXIS1 => nn
  w=w0+dw*findgen(nn)   ; wavelength vector re-created for each spectrum

  ppp=ppp+string(w[0],form='(f12.3)')
  ppp=ppp+string(w[nn-1],form='(f12.3)')
  if w1[0] lt w[0] then $
          ppp=ppp+string(w1[0]-w[0],form='(f12.3)') else $
          ppp=ppp+'         OK '
  if w1[n-1] gt w[nn-1] then $
          ppp=ppp+string(w1[n-1]-w[nn-1],form='(f12.3)') else $
          ppp=ppp+'         OK '
  print,ppp
end

print,'BFpro4 finished.'

return
end
;================================================================
pro broadfunct,sp,w,u,v,m,bb
   bb=fltarr(m,m)
   for i=0,m-1 do begin
      wb=fltarr(m)
      wb[0:i]=w[0:i]
      bb[*,i]=svsol(u,wb,v,sp,/double)
   end
return
end
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
;================================================================
pro crude,bf,vel,velc
; routine for crude measurement of several peaks in BF's in
; search of hidden components, through manual cursor marking
;
; velc will contain geocentric velocities, so convert to helio
; using:      
;     veladd=fltarr(5,n)
;     for i=0,4 do veladd[i,*]=velc[i,*]-hvc[0]+hvc+RVstd
;     where RVstd vecolcity of the template 
;
; Input:
;    bf = array of broadening functions
;    vel = corresponding velocity vector
; Output:
;    velc = updated data in velc array
;
; Use:
;    crude,bf15,vel,velc

  print,'Mark features by using the cursor'
  print,'End gien BF by the cursor left of the window'

  n=n_elements(bf[*,0])  ; how many broadening functions
  velc=replicate(-9999.,5,n)  ; up to 5 components for a given BF
  m=n_elements(vel)

for k=0,n-1 do begin  ; loop through broadening functions
  plot,vel,bf[k,*],tit='!17 BF #'+strtrim(string(k),2)    
                                                       ; plots the BF
  j=0
R:
  cursor,x,y,/down,/data
  if x lt vel[0] then goto,E
  if x gt vel[m-1] then goto,E1
  oplot,[x,x],[-1,2],thick=1,line=1
  velc[j,k]=x
  j=j+1
  goto,R
E:
endfor
E1:
return
end
;================================================================
pro de_spike1,spec1,spec2,wind
; interactive removal of spikes etc problems
; usage : de_spike1,spec1,spec2,100
; input: spec1 - spectrum
;        wind  - size of window around spike
; output: spec2 - de-spiked spectrum

  x=0 & y=0.
  w2=round(wind/2)
  n=n_elements(spec1)
  spec2=spec1

  plot,spec1,tit='Mark approx location of the spike with cursor'
  cursor,x,y,/down
  x = round(x)

  i1 = x-w2>0           ; display window defined
  i2 = x+w2<(n-1)

  plot,spec1,xran=[i1,i2],psym=-1

  print,'Enter cursor 1-st position'
  cursor,x,y,/down,/data
  x1=fix(x)
  oplot,[x1,x1],[spec1(x1),spec1(x1)],psym=2

  print,'Enter cursor 2-nd position'
  cursor,x,y,/down,/data
  x2=fix(x)
  oplot,[x2,x2],[spec1(x2),spec1(x2)],psym=2

  if x2 le x1+1 then begin 
  print,'No change, exiting...' 
  goto,E
  end

  k=(spec1(x2)-spec1(x1))/(x2-x1)
  for j=x1,x2 do spec2(j)=spec2(x1)+k*(j-x1)

  oplot,spec2,line=1

E: return
end
;================================================================
pro de_spike2,spec1,spec2,wave,n
; interactive removal of end problems
; usage : de_spike2,spec1,spec2,w0,nn
; input: spec1 - spectrum
;        wave  - wavelength vector, gets modified
;        n     - its size, gets modified
; output: spec2 - truncated spectrum

   spec2=spec1

   ans=0
   print,'Enter no of pixels to remove, negative for left end'
   read,'Number of pixels?',ans

   if ans eq 0 then goto,E

   if ans lt 0 then begin
      spec2=spec1[abs(ans):*]
      wave=wave[abs(ans):*]
      n=n-abs(ans)
   end

   if ans gt 0 then begin
      spec2=spec1[0:n-1-ans]
      wave=wave[0:n-1-ans]
      n=n-ans
   end 

E:
return
end
;================================================================
pro de_spike3,spec1,spec2,wind
; interactive removal of features from spectra
; usage : de_spike1,spec1,spec2,100
; input: spec1 - spectrum
;        wind  - size of window around the feature
; output: spec2 - cleaned spectrum

  x=0 & y=0.
  w2=round(wind/2)
  n=n_elements(spec1)
  spec2=spec1

  plot,spec1,tit='Mark with cursor approx location of the region to remove' 
  cursor,x,y,/down
  x = round(x)

  i1 = x-w2>0           ; display window defined
  i2 = x+w2<(n-1)

  plot,spec1,xran=[i1,i2],psym=-1

  print,'Enter cursor 1-st position'
  cursor,x,y,/down,/data
  x1=fix(x)
  oplot,[x1,x1],[spec1(x1),spec1(x1)],psym=2

  print,'Enter cursor 2-nd position'
  cursor,x,y,/down,/data
  x2=fix(x)
  oplot,[x2,x2],[spec1(x2),spec1(x2)],psym=2

  if x2 le x1+1 then begin 
  print,'No change, exiting...' 
  goto,E
  end

  k=(spec1(x2)-spec1(x1))/(x2-x1)
  for j=x1,x2 do spec2(j)=spec2(x1)+k*(j-x1)

  oplot,spec2,line=1

E: return
end
;================================================================
function gs_smooth,n,sigma
;
; prepares a n-element vector of total area=1 and Gaussian shape
; given by sigma, shifted to the origin [bbb.....bbb]
; for Gaussian filter using new=cnv(old,gs_smooth)
;
; when a given FWHM needed, use relation: FWHM = 2.354 sigma
;                                     or: sigma = 0.4248 FWHM
   b = findgen(n)
   c = b
   c = ((b - b(n/2))/sigma)^2/2
; bug fixed Jan.2,91, before not divided by 2
   b(*) = 0.
   d = where(c le 30.) 
   b(d) = exp(-c(d))
   b = shift(b, -n/2)
return, b/total(b)
end
;================================================================
function helio_jd,date,ra,dec, B1950 = B1950, TIME_DIFF = time_diff
;+
; NAME:
;      HELIO_JD
; PURPOSE:
;      Convert geocentric (reduced) Julian date to heliocentric Julian date
; EXPLANATION:
;      This procedure corrects for the extra light travel time between the 
;      Earth and the Sun.
;
; CALLING SEQUENCE:
;       jdhelio = HELIO_JD( date, ra, dec, /B1950, /TIME_DIFF)
;
; INPUTS
;       date - reduced Julian date (= JD - 2400000), scalar or vector, MUST
;               be double precision
;       ra,dec - scalars giving right ascension and declination in DEGREES
;               Equinox is J2000 unless the /B1950 keyword is set
;
; OUTPUTS:
;       jdhelio - heliocentric Julian date.   If /TIME_DIFF is set, then
;                 HELIO_JD() instead returns the time difference in seconds
;                 between the geocentric and heliocentric Julian date.
;                 
; OPTIONAL INPUT KEYWORDS 
;       /B1950 - if set, then input coordinates are assumed to be in equinox 
;                B1950 coordinates.
;       /TIME_DIFF - if set, then HELIO_JD() instead the time difference in 
;                seconds between the geocentric and heliocentric Julian date.
           
; EXAMPLE:
;       What is heliocentric julian date of an observation of V402 Cygni
;       (J2000: RA = 20 9 7.8, Dec = 37 09 07) taken June 15, 1973 at 11:40 UT?
;
;       IDL> juldate, [1973,6,15,11,40], jd      ;Get geocentric Julian date
;       IDL> hjd = helio_jd( jd, ten(20,9,7.8)*15., ten(37,9,7) )  
;                                                            
;       ==> hjd = 41848.9881
;
; Wayne Warren (Raytheon ITSS) has compared the results of HELIO_JD with the
; FORTRAN subroutines in the STARLINK SLALIB library (see 
; http://star-www.rl.ac.uk/).    
;                                                  Time Diff (sec)
;      Date               RA(2000)   Dec(2000)  STARLINK      IDL
;
; 1999-10-29T00:00:00.0  21 08 25.  -67 22 00.  -59.0        -59.0
; 1999-10-29T00:00:00.0  02 56 33.4 +00 26 55.  474.1        474.1
; 1940-12-11T06:55:00.0  07 34 41.9 -00 30 42.  366.3        370.2
; 1992-02-29T03:15:56.2  12 56 27.4 +42 10 17.  350.8        350.9
; 2000-03-01T10:26:31.8  14 28 36.7 -20 42 11.  243.7        243.7
; 2100-02-26T09:18:24.2  08 26 51.7 +85 47 28.  104.0        108.8
; PROCEDURES CALLED:
;       bprecess, xyz, zparcheck
;
; REVISION HISTORY:
;       Algorithm from the book Astronomical Photometry by Henden, p. 114
;       Written,   W. Landsman       STX     June, 1989 
;       Make J2000 default equinox, add B1950, /TIME_DIFF keywords, compute
;       variation of the obliquity      W. Landsman   November 1999
;-
 On_error,2
 If N_params() LT 3 then begin
    print,'Syntax -   jdhelio = HELIO_JD( date, ra, dec, /B1950, /TIME_DIFF)'
    print,'      date - reduced Julian date (= JD - 2400000)'
    print,'      Ra and Dec must be in degrees'
 endif

;Because 

 if not keyword_set(B1950) then bprecess,ra,dec,ra1,dec1 else begin
        ra1 = ra
        dec1 = dec
 endelse
 
 radeg = 180.0d/!DPI   
 zparcheck,'HELIO_JD',date,1,[3,4,5],[0,1],'Reduced Julian Date'

 delta_t = (double(date) - 33282.42345905d)/36525.0d
 epsilon_sec = poly( delta_t, [44.836d, -46.8495, -0.00429, 0.00181])
 epsilon = (23.433333d0 + epsilon_sec/3600.0d)/radeg
 ra1 = ra1/radeg
 dec1 = dec1/radeg

 xyz, date, x, y, z

;Find extra distance light must travel in AU, multiply by 1.49598e13 cm/AU,
;and divide by the speed of light, and multiply by 86400 second/year

 time = -499.00522d*( cos(dec1)*cos(ra1)*x + $
                 (tan(epsilon)*sin(dec1) + cos(dec1)*sin(ra1))*y)

 if keyword_set(TIME_DIFF) then return, time else $
      
       return, double(date) + time/86400.0d

 end
;================================================================
pro hjd_ha,prg_lst,lon,lat,jdtime,ha,am
; calculates HA for all observations in FITS pointed by 
;           file prg_lst
; 
; input: prg_lst = text file giving full names of FITS files
;        lon, lat = longitude, latitude of observatory in degr
; output: jdtime = heliocentric JD
;         ha = hour angle velocity correction
;         am = air mass
;
; usage: 
;         hjd_ha,'list.txt',-79.421667d0,43.863333d0,hjd,ha,am
;                           DDO coordinates
;         convention here: long negative = West

prgs=strarr(200)                   ; read the list of files
openr,1,prg_lst
i=0
a=''
while(EOF(1) ne 1) do begin
    readf,1,form='(a)',a
    prgs[i]=a
    i = i + 1
endwhile
prgs=prgs[0:i-1]
close,1
nnn=i

jdt=0.d0

jdtime=dblarr(nnn)     ; JD of mid-obs
ha=dblarr(nnn)         ; HA in degrees                
am=dblarr(nnn)         ; Air Mass

radeg=180.d0/!dpi

; cycle through all observations including the template at the top

for i=0,nnn-1 do begin                  
  prg_fts=strtrim(prgs[i],2)           ; individual spectra
  prg=readfits(prg_fts,hd,/silent)             ;    in FITS 
  parse_head1,hd,jd,t1,t2,epoch,RA,dec ; extracts from header
;  parse_fits_head,hd,nn,w0,dw,jd,t1,t2,ra,dec,equin,obj
;  epoch=equin   ; for compatibility
;  print,'ra (hh), dec (deg): ',ra,dec
  t1=t1/24.                            ; fraction of day
  t2=t2/24.
  if t2 lt t1 then t2=t2+1.
  jdt=jd+(t1+t2)/2                     ; full JD, geoc mid-obs
  jdtime[i]=helio_jd(jdt-2400000.d0,ra*15.,dec) ; (RA,dec) of the epoch
  jdtime[i]=jdtime[i]+2400000.d0                ;    not of time of obs.
; RA in hh, dec in deg

  ct2lst,lst,lon,dummy,jdtime[i]   ; lst=local sidereal time in hours
;  print,'lst (hh): ',lst
  ha[i]=15.*(lst-ra)               ; HA in degrees
  if ha[i] gt 180. then ha[i]=ha[i]-360.
  if ha[i] lt -180. then ha[i]=ha[i]+360.
  cosz = sin(lat/radeg)*sin(dec/radeg)+  $
         cos(lat/radeg)*cos(dec/radeg)*cos(ha[i]/radeg)
  am[i]=1./cosz
endfor

return
end
;=======================================================================
pro hjd_phase,prg_lst,T0,period,jdtime,phase
; calculates heliocentric JD time and phase
;   for all observations in FITS poined by file prg_lst
; input: prg_lst = text file giving full names of FITS files
;        T0      = initial epoch of binary
;        period  = orbital period
; output: jdtime = heliocentric JD
;         phase  = corresponding phase
; usage: hjd_phase,'v401cyg.lst',t0,period,hjd,phase


prgs=strarr(200)                   ; read the list of files
openr,1,prg_lst
i=0
a=''
while(EOF(1) ne 1) do begin
    readf,1,form='(a)',a
    prgs[i]=a
    i = i + 1
endwhile
prgs=prgs[0:i-1]
close,1
nnn=i

jdt=0.d0
jdtime=dblarr(nnn)                     ; JD of mid-obs
phase=dblarr(nnn)

; cycle through all observations including the template at the top

for i=0,nnn-1 do begin                  
  prg_fts=strtrim(prgs[i],2)           ; individual spectra
  prg=readfits(prg_fts,hd)             ;    in FITS 
  parse_head1,hd,jd,t1,t2,epoch,RA,dec ; extracts from header
  t1=t1/24.                            ; fraction of day
  t2=t2/24.
  if t2 lt t1 then t2=t2+1.
  jdt=jd+(t1+t2)/2                     ; full JD, geoc mid-obs
  jdtime[i]=helio_jd(jdt-2400000.d0,ra*15.,dec) ; RA, dec of the epoch
  jdtime[i]=jdtime[i]+2400000.d0                ; not of time of obs.
  phase[i]=(jdtime[i]-t0)/period          
endfor
  phase=phase-floor(phase)

return
end
;================================================================
pro hjd_vel,prg_lst,lon,lat,jdtime,hvc
; calculates heliocentric JD time and velocity correction
;   for all observations in FITS format pointed by file prg_lst
; 
; input: prg_lst = text file giving full names of FITS files
;        lon = longitude of obs, decimal degrees, West -> negative
;        lat = latitude of obs, decimal degrees 
; output: jdtime = heliocentric JD
;         phase  = velocity correction
; usage: 
;       hjd_vel,'v401cyg.lst',-79.421667d0,43.863333d0,hjd,hvc
;                               DDO long    DDO lat

prgs=strarr(200)                   ; read the list of files
openr,1,prg_lst
i=0
a=''
while(EOF(1) ne 1) do begin
    readf,1,form='(a)',a
    prgs[i]=a
    i = i + 1
endwhile
prgs=prgs[0:i-1]
close,1
nnn=i                             ; headers read in

jdt=0.d0
jdtime=dblarr(nnn)                     ; JD of mid-obs
hvc=dblarr(nnn)
vh=dblarr(3)                           ; heliocentric
vb=dblarr(3)                           ; barocentric

radeg=180.d0/!dpi

; cycle through all observations including the template at the top

for i=0,nnn-1 do begin                  
  prg_fts=strtrim(prgs[i],2)           ; individual spectra
  prg=readfits(prg_fts,hd,/silent)             ;    in FITS 
  parse_head1,hd,jd,t1,t2,epoch,RA,dec ; extracts from header
;  parse_fits_head,hd,nn,w0,dw,jd,t1,t2,ra,dec,equin,obj
;  epoch=equin   ; for compatibility
  t1=t1/24.                            ; fraction of day
  t2=t2/24.
  if t2 lt t1 then t2=t2+1.
  jdt=jd+(t1+t2)/2                     ; full JD, geoc time mid-obs
  jdtime[i]=helio_jd(jdt-2400000.d0,ra*15.,dec) ; (RA,dec) of the epoch,
  jdtime[i]=jdtime[i]+2400000.d0                ;    not of time of obs.
; ra in decimal hh, dec in deg

  baryvel,jdtime[i],epoch,vh,vb
  hvc[i] = vb[0]*cos(dec/radeg)*cos(15*ra/radeg) + $  ; vel toward star
       vb[1]*cos(dec/radeg)*sin(15*ra/radeg) + vb[2]*sin(dec/radeg) 
;  hvc[i] = vh[0]*cos(dec/radeg)*cos(15*ra/radeg) + $  ; vel toward star
;       vh[1]*cos(dec/radeg)*sin(15*ra/radeg) + vh[2]*sin(dec/radeg)
;  better use relative to barycentre, vb (not vh relative to Sun) 
  
  ct2lst,lst,lon,dummy,jdtime[i]   ; akward calculation of LST from JDtime
                                   ; LST & RA in decimal hours
  ha=15*(lst-ra)                   ; HA in degrees
  if ha gt  180. then ha=ha-360.
  if ha lt -180. then ha=ha+360.
  rot=-0.465*cos(lat/radeg)*cos(dec/radeg)*sin(ha/radeg)
;print,f='(f20.4,5f12.2)',jdtime[i],15*lst,15*ra,dec,ha,rot
  hvc[i]=hvc[i]+rot                ; Earth rotation correction
endfor

return
end
;
;================================================================
pro ieee_to_host, data, IDLTYPE = idltype
;+
; NAME:
;     IEEE_TO_HOST
; PURPOSE:
;     Translate an IDL variable from IEEE-754 to host representation 
; EXPLANATION:
;     The variable is translated from IEEE-754 (as used, for
;     example, in FITS data ), into the host machine architecture.
;
; CALLING SEQUENCE:
;     IEEE_TO_HOST, data, [ IDLTYPE = , ]
;
; INPUT-OUTPUT PARAMETERS:
;     data - any IDL variable, scalar or vector.   It will be modified by
;             IEEE_TO_HOST to convert from IEEE to host representation.  Byte 
;             and string variables are returned by IEEE_TO_HOST unchanged
;
; OPTIONAL KEYWORD INPUTS:
;     IDLTYPE - scalar integer (1-15) specifying the IDL datatype according
;               to the code given by the SIZE function.     This keyword
;               is usually when DATA is a byte array to be interpreted as
;               another datatype (e.g. FLOAT).
;
; EXAMPLE:
;       A 2880 byte array (named FITARR) from a FITS record is to be 
;       interpreted as floating and converted to the host representaton:
;
;       IDL> IEEE_TO_HOST, fitarr, IDLTYPE = 4     
;
; METHOD:
;       The BYTEORDER procedure is called with the appropriate keyword
;
; MODIFICATION HISTORY:
;      Written, W. Landsman   Hughes/STX   May, 1992
;      Converted to IDL V5.0   W. Landsman   September 1997
;      Under VMS check for IEEE -0.0 values   January 1998
;      VMS now handle -0.0 values under IDL V5.1    July 1998
;      Added new integer datatypes  C. Markwardt/W. Landsman  July 2000
;      Post-V5.1 version, no VMS negative zero check  W. Landsman July 2001
;     
;-
 On_error,2 

 if N_params() EQ 0 then begin
    print,'Syntax - IEEE_TO_HOST, data, [ IDLTYPE = ]'
    return
 endif  

 npts = N_elements( data )
 if npts EQ 0 then $
     message,'ERROR - IDL data variable (first parameter) not defined'

 sz = size(data)
 if not keyword_set( idltype) then idltype = sz[ sz[0]+1]
 
 case idltype of

      1: return                             ;byte

      2: byteorder, data, /NTOHS            ;integer

      3: byteorder, data, /NTOHL            ;long

      4: byteorder, data, /XDRTOF                              ;float

      5: byteorder, data, /XDRTOD                              ;double

      6: byteorder, data, /XDRTOF

      7: return                             ;string

       8: BEGIN				    ;structure

	Ntag = N_tags( data )

	for t=0,Ntag-1 do  begin
          temp = data.(t)
          ieee_to_host, temp
          data.(t) = temp
        endfor 

       END

	9: byteorder, data,/XDRTOD                              ;complex

       12: byteorder, data, /NTOHS

       13: byteorder, data, /NTOHL

       14: if (long(['01'xb,'02'xb,'03'xb,'04'xb],0,1))(0) NE '01020304'x then $
            byteorder, data, /L64swap

       15: if (long(['01'xb,'02'xb,'03'xb,'04'xb],0,1))(0) NE '01020304'x then $  
             byteorder, data, /L64swap

       ELSE: message,'Unrecognized datatype of ' + strtrim(idltype,2)
 ENDCASE

 return
 end 
;================================================================
function is_ieee_big
;+
; NAME:
;	IS_IEEE_BIG
; PURPOSE:
;	Determine if the current machine uses IEEE, big-endian numbers.
; EXPLANATION:
;       (Big endian implies that byteorder XDR conversions are no-ops).
; CALLING SEQUENCE:
;	flag = is_ieee_big()
; INPUT PARAMETERS:
;       None
; RETURNS:
;       1 if the machine appears to be IEEE-compliant, 0 if not.
; COMMON BLOCKS:
;	None.
; SIDE EFFECTS:
;	None
; RESTRICTIONS:
; PROCEDURE:
;       A sample int, long, float and double are converted using
;       byteorder and compared with the original.  If there is no
;       change, the machine is assumed to be IEEE compliant and
;       big-endian.
; MODIFICATION HISTORY:
;       Written 15-April-1996 by T. McGlynn for use in MRDFITS.
;	13-jul-1997	jkf/acc	- added calls to check_math to avoid
;				  underflow messages in V5.0 on Win32 (NT).
;	Converted to IDL V5.0   W. Landsman   September 1997
;-

    itest = 512
    ltest = 102580L
    ftest = 1.23e10
    dtest = 1.23d10
		
		
    it2 = itest
    lt2 = ltest
    ft2 = ftest
    dt2 = dtest
				
    byteorder, it2, /htons
    byteorder, lt2, /htonl
    byteorder, ft2, /ftoxdr
    byteorder, dt2, /dtoxdr

    if itest eq it2  and  ltest eq lt2   and ftest eq ft2  and dtest eq dt2  $
    then begin
    	dum = check_math()
        return, 1
    endif else begin
    	dum = check_math()
        return, 0
    endelse
    end								    
;================================================================
PRO JULDATE, DATE, JD, PROMPT = prompt
;+                                                                  
; NAME:
;     JULDATE
; PURPOSE:                                   
;     Convert from calendar to Reduced Julian Date
;
; EXPLANATION:
;     Julian Day Number is a count of days elapsed since Greenwich mean noon 
;     on 1 January 4713 B.C.  The Julian Date is the Julian day number
;     followed by the fraction of the day elapsed since the preceding noon. 
;
;     This procedure duplicates the functionality of the JULDAY() function in
;     in the standard IDL distribution, but also allows interactive input and
;     gives output as Reduced Julian date (=JD - 2400000.)  
;     (Also note that prior to V5.1 there was a bug in JULDAY() that gave 
;     answers offset by 0.5 days.)
;
; CALLING SEQUENCE:
;     JULDATE, /PROMPT           ;Prompt for calendar Date, print Julian Date
;               or
;     JULDATE, date, jd      
;
; INPUT:
;     DATE -  3 to 6-element vector containing year,month (1-12),day, and 
;              optionally hour, minute, and second all specified as numbers
;              (Universal Time).   Year should be supplied with all digits.
;              Years B.C should be entered as negative numbers (and note that
;              Year 0 did not exist).  If Hour, minute or seconds are not 
;              supplied, they will default to 0. 
;
;  OUTPUT:
;       JD - Reduced Julian date, double precision scalar.  To convert to
;               Julian Date, add 2400000.   JULDATE will print the value of
;               JD at the terminal if less than 2 parameters are supplied, or 
;               if the /PROMPT keyword is set
;      
;  OPTIONAL INPUT KEYWORD:
;       /PROMPT - If this keyword is set and non-zero, then JULDATE will prompt
;               for the calendar date at the terminal.
;
;  RESTRICTIONS:
;       The procedure HELIO_JD can be used after JULDATE, if a heliocentric
;       Julian date is required.
;
;  EXAMPLE:
;       A date of 25-DEC-1981 06:25 UT may be expressed as either
;
;       IDL> juldate, [1981, 12, 25, 6, 25], jd       
;       IDL> juldate, [1981, 12, 25.2673611], jd 
;
;       In either case, one should obtain a Reduced Julian date of 
;       JD = 44963.7673611
;
;  PROCEDURE USED:
;       GETOPT()
;  REVISION HISTORY
;       Adapted from IUE RDAF (S. Parsons)                      8-31-87
;       Algorithm from Sky and Telescope April 1981   
;       Added /PROMPT keyword, W. Landsman    September 1992
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Make negative years correspond to B.C. (no year 0), work for year 1582
;       Disallow 2 digit years.    W. Landsman    March 2000
;-
 On_error,2 

 if ( N_params() EQ 0 ) and (not keyword_set( PROMPT ) ) then begin
     print,'Syntax - JULDATE, date, jd          or JULDATE, /PROMPT'
     print, $
     '  date - 3-6 element vector containing [year,month,day,hour,minute,sec]'
     print,'  jd - output reduced julian date (double precision)'
     return
 endif

 if ( N_elements(date) EQ 0 ) then begin   

    opt = ''                                                          
    rd: read,' Enter Year,Month,Day,Hour, Minute, Seconds (All Numeric): ',opt
    date = getopt( opt, 'F' )

 endif

 case N_elements(date) of      

    6: 
    5: date = [ date, 0.0d]
    4: date = [ date, 0.0d,0.0d]    
    3: date = [ date, 0.0d, 0.0d,0.0d]
    else: message,'Illegal DATE Vector - must have a least 3 elements'

  endcase   

 iy = floor( date[0] ) 
 if iy lt 0 then iy = iy +1  else $
    if iy EQ 0 then message,'ERROR - There is no year 0'                   
 im = fix( date[1] )
 date = double(date)
 day = date[2] + ( date[3] + date[4]/60.0d + date[5]/3600.0d) / 24.0d
;
 if ( im LT 3 ) then begin   ;If month is Jan or Feb, don't include leap day

     iy= iy-1 & im = im+12 

 end

 a = long(iy/100)
 ry = float(iy)

 jd = floor(ry*0.25d) + 365.0d*(ry -1860.d) + fix(30.6001d*(im+1.)) + $
      day  - 105.5d

;Gregorian Calendar starts on Oct. 15, 1582 (= RJD -100830.5)
 if jd GT -100830.5 then jd = jd + 2 - a + floor(a/4)

 if N_params() LT 2 or keyword_set( PROMPT) then begin      
    yr = fix( date[0] )
    print, FORM='(A,I4,A,I3,A,F9.5)',$ 
       ' Year ',yr,'    Month', fix(date[1] ),'    Day', day 
    print, FORM='(A,F15.5)',' Reduced Julian Date:',JD                       
 endif
 
 return                               
 end                                  ; juldate
;================================================================
function map4,x,m
; version for IDL 4.0  (SMR Feb.15'97)
; shifts vectors x vertically within m
; for broadening functions

   t=0
   n=n_elements(x)

   if n gt 32767 then begin
      print,'n - too large' & goto,E
   end
   if (m mod 2) ne 1 then begin
      print,'m - must be odd' & goto,E       
   end
   if (n mod 2) ne 0 then begin
      print,'n must be even' & goto,E
   end

   t = fltarr(m) # fltarr(n-m+1)    ; t(m,n-m)=
                                    ; t(small,large-small) dimensions
   for j=0,m-1 do $
      for i=m/2,n-m/2-1 do t(j,i-m/2)=x(i-j+m/2)
   
E: return,t
end
;================================================================
function one_gs,x,y,a,da
;
; function to fit one Gaussian with zero point (baseline)
; usage: yfit=one_gs(x,y,a,da)
;   input: x,y,a
;   output: yfit,da
; da(4) are diff.corr.'s to parameters a(4) as follows:
;   a(0) - central strength
;   a(1) - central position
;   a(2) - width
;   a(3) - zero point

; In two_gs parameters are:
; a(0) and a(3) - central strengths
; a(1) and a(4) - central positions
; a(2) and a(5) - widths
; a(6) is zero point
;

  n = n_elements(x)
  c = x & e1=x
  c = ((x-a(1))/a(2))^2
  e1(*)=0.
  d = where(c le 30.)
  e1(d) = exp(-c(d))
  res = y - a(0)*e1 - a(3)
  print,'Sum res^2 = ',total(res^2)
;

  t = fltarr(n,4)         ; derivatives
  t(*,0) = e1
  t(*,1) = 2*e1*a(0)*(x-a(1))/a(2)^2
  t(*,2) = t(*,1)*(x-a(1))/a(2)
  t(*,3) = 1.
;
  tt = transpose(t)#t     ; solution
  tr = transpose(t)#res
  svd,tt,w,u,v
  wp = fltarr(4,4)
  mw = max(w)
  for i=0,3 do if w(i) ge mw*1.e-6 then wp(i,i)=1./w(i)
  da = v#wp#(transpose(u)#tr)
  a1 = a + da
;
  c = x & e1=x            ; new Gaussian
  e1(*)=0.
  c = ((x-a1(1))/a1(2))^2
  d = where(c le 30.)
  e1(d) = exp(-c(d))
  res = y - a1(0)*e1 - a1(3)
  print,'Sum new res^2 = ',total(res^2)
return,a1(0)*e1+a1(3)
end
;================================================================
function one_gs1,x,y,a,da
;
; function to fit one Gaussian with zero point (baseline)
; usage: yfit=one_gs1(x,y,a,da)
;   input: x,y,a
;   output: yfit,da
; da(4) are diff.corr.'s to parameters a(4) as follows:
;   a(0) - central strength
;   a(1) - central position
;   a(2) - width, fixed in this version
;   a(3) - zero point
;
; NOTE: this version keeps the width of the Gaussian fixed 

  n = n_elements(x)
  c = x & e1=x
  c = ((x-a(1))/a(2))^2
  e1(*)=0.
  d = where(c le 30.)
  e1(d) = exp(-c(d))
  res = y - a(0)*e1 - a(3)
  print,'Sum res^2 = ',total(res^2)
;

  t = fltarr(n,3)         ; derivatives, width constant
  t(*,0) = e1
  t(*,1) = 2*e1*a(0)*(x-a(1))/a(2)^2
;  t(*,2) = t(*,1)*(x-a(1))/a(2)    ; width not used
;  t(*,3) = 1.
  t(*,2) = 1. 
;
  tt = transpose(t)#t     ; solution
  tr = transpose(t)#res
  svd,tt,w,u,v
;  wp = fltarr(4,4)
  wp = fltarr(3,3)
  mw = max(w)
;  for i=0,3 do if w(i) ge mw*1.e-6 then wp(i,i)=1./w(i)
  for i=0,2 do if w(i) ge mw*1.e-6 then wp(i,i)=1./w(i)
  dda = v#wp#(transpose(u)#tr)
  da = [dda[0:1],0.,dda[2]]
  a1 = a + da
;
  c = x & e1=x            ; new Gaussian
  e1(*)=0.
  c = ((x-a1(1))/a1(2))^2
  d = where(c le 30.)
  e1(d) = exp(-c(d))
  res = y - a1(0)*e1 - a1(3)
  print,'Sum new res^2 = ',total(res^2)
return,a1(0)*e1+a1(3)
end
;================================================================
pro one_star1,star,tp,res
; sine fit to one star
; 
; Input: star - data array [3,*]
;        tp   - 2-el vector with T0,P
; Output: res - 6-el vector:
;               gam, K1 - 2 parameters + dummy 
;               and 2-el vector of errors + dummy
; Usage: one_star1,sv,sv0,sv1
;
  ph = (star[0,*]-tp[0])/tp[1]
  ph = reform(ph)
  ph = ph - floor(ph)
  n = n_elements(ph)

  a = dblarr(2,n)
  a[0,*] = 1.
  a[1,0:n-1]   = -sin(2*!pi*ph)  ; sense as in A-type
;  a[2,n:2*n-1] = +sin(2*!pi*ph) ; no 2nd star
  b = reform(star[1,*])  ; velocities
  w = reform(star[2,*])  ; weights
  v = 0. & cv = 0. ; variance and covariance
  res = solv(a,b,w,v,cv)  ; SVD solution
;  gam = res[0] & K1 = res[1] 
  res=[res,0,sqrt(v),0]  ; zero for place holding of T0

return
end
;================================================================
pro one_star2,star,tp,res,corr
; corrections to the sine fit for one star
; 
; Input: star - data array [5,*]
;        tp   - 2-el vector with T0,P
;        res - 6-el vector of the sine fit
;    consisting of: gam, K1, 0 - 3 parameters 
;                  and 3-el vector of errors
; Output: corr - 3 el vector of diff.corrs 
;          delta of [gam,K1,T0]
; Usage: one_star2,sv,sv0,sv1,sv2
;
  ph = (star[0,*]-tp[0])/tp[1]
  ph = reform(ph)
  ph = ph - floor(ph)
  n = n_elements(ph)

  dif1 = reform(star[1,*])-reform(res[0]-res[1]*sin(2*!pi*ph))
  w = reform(star[2,*]) ; weights
  b = dif1              ; velocity diffs O-C

  a = dblarr(3,n)
  a[0,*] = 1.
  a[1,*]   = -sin(2*!pi*ph)  ; sense as in A-type
  a[2,*]   = +res[1]*2*!pi/tp[1]*cos(2*!pi*ph) ; corr to T0

  v = 0. & cv = 0. ; variance and covariance
  corr = solv(a,b,w,v,cv)  ; SVD solution of corrections
;  gam = corr[0] & K1 = corr[1]
  corr=[corr,sqrt(v)]

  dif1 = reform(star[1,*])-   $
     reform((res[0]+corr[0])-(res[1]+corr[1])*sin(2*!pi*ph))
  print,sqrt(total(dif1^2*w)/total(w))

return
end
;================================================================
pro one_star3,star,tp,res,res1
; bootstrap estimates of errors for one star
; 
; Input: star - data array [5,*]
;        tp   - 2-el vector with T0,P
;        res - 6-el vector of the sine fit
;    consisting of: gam, K1, 0 - 3 parameters 
;                  and 3-el vector of errors
; Output: res1 - 3 +/- 1-sigma ranges and medians of 
;        diff.corr's [gam,K1,T0]
; Usage: one_star3,ah,ah0,ah1,ah3
;
  ph = (star[0,*]-tp[0])/tp[1]
  ph = reform(ph)
  ph = ph - floor(ph)
  n = n_elements(ph)

  dif1 = reform(star[1,*])-reform(res[0]-res[1]*sin(2*!pi*ph))
  w = reform(star[2,*])  ; weights
  b = dif1               ; velocity diffs O-C
  a = dblarr(3,n)
  a[0,*] = 1.
  a[1,*]   = -sin(2*!pi*ph)  ; sense as in A-type
  a[2,*]   = +res[1]*2*!pi/tp[1]*cos(2*!pi*ph) ; corr to T0

  v = 0. & cv = 0. ; variance and covariance

  a1 = a & b1 = b & m = n
  w1 = replicate(1.,m)  ; uniform weights from now on
  for i = 0,m-1 do begin
    a1[*,i] = a[*,i]*w[i]
    b1[i] = b[i]*w[i]
  end

  k = 1000  
  a2 = a1 & b2 = b1
  rang = fltarr(3,k)
  for i = 0,k-1 do begin
      for j = 0,m-1 do begin 
         l = fix(m*randomu(seed))
         a2[*,j] = a1[*,l]
         b2[j] = b1[l]
      endfor     
  rang[*,i] = solv(a2,b2,w1,v,cv)  ; SVD solution of corr's
; order of corrections: gam, K1, T0
  endfor

  for i = 0,2 do rang[i,*]=rang[i,sort(rang[i,*])]
  res1 = dblarr(3,3)
  res1[*,0] = rang[*,fix(0.158*k)]   ; -1 sigma
  res1[*,1] = rang[*,fix(k/2)]       ; median
  res1[*,2] = rang[*,fix(0.841*k)]   ; +1 sigma

return
end
;================================================================
pro parse_fits_head,hd,nn,w0,dw,jd,t1,t2,ra,dec,equin,obj
; 
; parses the FITS header hd to obtain parameters needed for
;    further processing:
;    nn = number of pixels
;    w0 = the starting wavelength,
;    dw = the wavelength increment,
;    jd = JD for the date at midnight
;    t1 = start of observation in decimal UT hours
;    t2 = end of observation in decimal UT hours
;    ra = RA in decimal hours
;    dec = dec in decimal degrees
;    equin = equinox as integer
;    obj = object name
;
; Input: hd
; Output: nn, w0, dw, jd, t1, t2, ra, dec, equin, obj
;
; usage: parse_fits_head,h,nn,w0,dw,jd,t1,t2,ra,dec,equin,obj

; numerical quantities in FITS header ---------------------------------
nn = 0 & w0=0.d0 & dw=0.d0 & equin = 0.

n=n_elements(hd)
;print,strtrim(string(n),2),' lines in the header'

for i=0,n-1 do if (strmid(hd[i],0,6) eq 'NAXIS1') then $
     nn=fix(strmid(hd[i],11,20))
for i=0,n-1 do if (strmid(hd[i],0,6) eq 'CRVAL1') then $
     w0=double(strmid(hd[i],11,30))
for i=0,n-1 do if (strmid(hd[i],0,6) eq 'CDELT1') then $
     dw=double(strmid(hd[i],11,30))
if dw eq 0.0d0 then begin
for i=0,n-1 do if (strmid(hd[i],0,5) eq 'CD1_1') then $
     dw=double(strmid(hd[i],11,30))
end
for i=0,n-1 do if (strmid(hd[i],0,7) eq 'EQUINOX') then $
     equin=float(strmid(hd[i],11,20))

;print,form='(a,i6,3f12.5)','nn, wstart, wend: ', $
;     nn,w0,w0+dw*(nn-1),equin

if nn eq 0 then print,'nn not found'
if w0 eq 0.d0 then print,'w0 not found'
if nn eq 0 then print,'dw not found'

if equin eq 0. then begin 
     for i=0,n-1 do if (strmid(hd[i],0,5) eq 'EPOCH') then $
          equin=float(strmid(hd[i],11,20))
;     print,'EQUINOX not found, EPOCH used instead'
end

; string-based quantities in FITS header -----------------------------

for i=0,n-1 do if (strmid(hd[i],0,8) eq 'DATE-OBS') then $
     st=strmid(hd[i],strpos(hd[i],"'")+1,10)
     if strmid(st,2,1) eq '/' then begin    ; old DDO format
        yr=fix(strmid(st,6,2))+1900
        mo=fix(strmid(st,3,2))
        da=fix(strmid(st,0,2))
     end else begin                     ; new DDO format since Y2K
        yr=fix(strmid(st,0,4))          ; end LCO
        mo=fix(strmid(st,5,2))
        da=fix(strmid(st,8,2))
     end
for i=0,n-1 do if (strmid(hd[i],0,7) eq 'UT-DATE') then $
     st=strmid(hd[i],strpos(hd[i],"'")+1,10)
     if strmid(st,2,1) eq '/' then begin    ; old DDO format
        yr=fix(strmid(st,6,2))+1900
        mo=fix(strmid(st,3,2))
        da=fix(strmid(st,0,2))
     end else begin                     ; new DDO format since Y2K
        yr=fix(strmid(st,0,4))          ; and LCO
        mo=fix(strmid(st,5,2))
        da=fix(strmid(st,8,2))
     end

     ;print,'Date: ',yr,mo,da
     juldate,[yr,mo,da],jd         ; reduced: JD-2400000.0
     jd=2400000.d0+jd              ; coorect for date at midnight

for i=0,n-1 do if (strmid(hd[i],0,8) eq 'TIME-OBS') then $
     st=strmid(hd[i],11,8)                  ; DDO keyword 
     hh=strmid(st,0,2)
     mm=strmid(st,3,2)
     ss=strmid(st,6,2)
     t1=ten(hh,mm,ss)          ; start time in decimal hours

for i=0,n-1 do if (strmid(hd[i],0,8) eq 'UTSTART ') then $
     st=strmid(hd[i],strpos(hd[i],"'")+1,8)   ; LCO keyword
     hh=strmid(st,0,2)
     mm=strmid(st,3,2)
     ss=strmid(st,6,2)
     ;print,'UT start: ',hh,mm,ss
     t1=ten(hh,mm,ss)          ; start time in decimal hours

for i=0,n-1 do if (strmid(hd[i],0,8) eq 'UT-START') then $
     st=strmid(hd[i],strpos(hd[i],"'")+1,8)   ; LCO keyword
     hh=strmid(st,0,2)
     mm=strmid(st,3,2)
     ss=strmid(st,6,2)
     ;print,'UT start: ',hh,mm,ss
     t1=ten(hh,mm,ss)          ; start time in decimal hours

for i=0,n-1 do if (strmid(hd[i],0,8) eq 'ENDTIME ') then $
     st=strmid(hd[i],11,8)                      ; DDO keyword 
     hh=strmid(st,0,2)
     mm=strmid(st,3,2)
     ss=strmid(st,6,2)
     t2=ten(hh,mm,ss)    ; stop time in decimal hours

for i=0,n-1 do if (strmid(hd[i],0,8) eq 'UTEND   ') then $
     st=strmid(hd[i],strpos(hd[i],':')-2,8)     ; LCO keyword
     hh=strmid(st,0,2)
     mm=strmid(st,3,2)
     ss=strmid(st,6,2)
    ;print,'UT end: ',hh,mm,ss
     t2=ten(hh,mm,ss)    ; stop time in decimal hours

for i=0,n-1 do if (strmid(hd[i],0,8) eq 'RA      ') then $
     st=strmid(hd[i],strpos(hd[i],"'")+1,11)   ; general parsing
     hh=strmid(st,0,2)
     mm=strmid(st,3,2)
     ss=strmid(st,6,5)
     ;print,'RA: ',hh,mm,ss
     RA=ten(hh,mm,ss)    ; RA in decimal hours

for i=0,n-1 do if (strmid(hd[i],0,8) eq 'DEC     ') then $
     st=strmid(hd[i],strpos(hd[i],"'")+1,11)   ; general parsing
     stp=strpos(st,":")
     if stp eq 3 then begin 
        hh=strmid(st,0,3)
        mm=strmid(st,4,2)
        ss=strmid(st,7,5)
     end else begin
        hh=strmid(st,0,2)
        mm=strmid(st,3,2)
        ss=strmid(st,6,5)
     end
     ;print,'dec: ',hh,mm,ss
     dec=ten(hh,mm,ss)          ; dec in decimal degrees

for i=0,n-1 do if (strmid(hd[i],0,8) eq 'OBJECT  ') then $
     st=strmid(hd[i],strpos(hd[i],"'")+1,20)   ; general parsing
     stp=strpos(st,"'")
     obj=strmid(st,-stp,stp)
     ;print,'Object: ',obj

return
end
;================================================================
pro parse_head,hd,w0,dw,nn
; 
; parses the FITS header hd to obtain the starting wavelength w0,
;    the wavelength increment dw and the number of pixels nn
; input: hd
; output: w0, dw, nn
;
; usage: parse_head,hs,w0,dw,nn

w0=0.d0
dw=0.d0
nn=0

n=n_elements(hd)
; print,strtrim(string(n),2),' lines in the header'

for i=0,n-1 do if (strmid(hd[i],0,6) eq 'NAXIS1') then $
     nn=fix(strmid(hd[i],11,20))

for i=0,n-1 do if (strmid(hd[i],0,6) eq 'CRVAL1') then $
     w0=double(strmid(hd[i],11,30))

for i=0,n-1 do if (strmid(hd[i],0,6) eq 'CDELT1') then $
     dw=double(strmid(hd[i],11,30))

if dw eq 0.0d0 then begin
for i=0,n-1 do if (strmid(hd[i],0,5) eq 'CD1_1') then $
     dw=double(strmid(hd[i],11,30))
end

;print,form='(a,i6,2f12.5)','nn, wstart, wend: ',nn,w0,w0+dw*(n-1)

if w0 eq 0.d0 then print,'Spectrum parameters not found!'

return
end
;================================================================
pro parse_head1,hd,jd,t1,t2,epoch,RA,dec
; 
; parses the FITS header hd to obtain:
;    Julian Date for the date
;    the starting and ending times t1, t2 in decimal hours
;    the epoch in years
;    RA, dec of the star in decimal hours (RA) and degrees (dec)
; Keywords EPOCH or EQUINOX used the same way for epoch.
;
; input: hd
; output: jd,t1,t2,epoch,RA,dec
;
; usage: parse_head1,hs,jd,t1,t2,ep,RA,dec

jd=0.d0    ; declarations
t1=0.d0
t2=0.d0
epoch=0
RA=0.d0
dec=0.d0

st=''
st1=''

n=n_elements(hd)
; print,strtrim(string(n),2),' lines in the header'

for i=0,n-1 do if (strmid(hd[i],0,8) eq 'DATE-OBS') then $
     st=strmid(hd[i],11,10)
     if strmid(st,2,1) eq '/' then begin    ; old DDO format
       yr=fix(strmid(st,6,2))+1900
       mo=fix(strmid(st,3,2))
       da=fix(strmid(st,0,2))
     end else begin                         ; new DDO format since Y2K
       yr=fix(strmid(st,0,4))
       mo=fix(strmid(st,5,2))
       da=fix(strmid(st,8,2))
     end
     juldate,[yr,mo,da],jd
     jd=2400000.d0+jd

for i=0,n-1 do if (strmid(hd[i],0,8) eq 'TIME-OBS') then $
     st=strmid(hd[i],11,8)
     hh=strmid(st,0,2)
     mm=strmid(st,3,2)
     ss=strmid(st,6,2)
     t1=ten(hh,mm,ss)    ; start time in decimal hours

for i=0,n-1 do if (strmid(hd[i],0,8) eq 'ENDTIME ') then $
     st=strmid(hd[i],11,8)
     hh=strmid(st,0,2)
     mm=strmid(st,3,2)
     ss=strmid(st,6,2)
     t2=ten(hh,mm,ss)    ; stop time in decimal hours

for i=0,n-1 do if (strmid(hd[i],0,8) eq 'RA      ') then $
     st=strmid(hd[i],11,11)
     hh=strmid(st,0,2)
     mm=strmid(st,3,2)
     ss=strmid(st,6,5)
     RA=ten(hh,mm,ss)    ; RA in decimal hours

for i=0,n-1 do if (strmid(hd[i],0,8) eq 'DEC     ') then $
     st=strmid(hd[i],11,12)
     if ((strmid(st,0,1) eq '-') or (strmid(st,0,1) eq '+')) then begin 
        hh=strmid(st,0,3)
        mm=strmid(st,4,2)
        ss=strmid(st,7,5)
     end else begin
        hh=strmid(st,0,2)
        mm=strmid(st,3,2)
        ss=strmid(st,6,5)
     end
; print,hh,' ',mm,' ',ss
     dec=ten(hh,mm,ss)   ; dec in decimal degrees

for i=0,n-1 do if (strmid(hd[i],0,8) eq 'EQUINOX ') then $
     st1=strmid(hd[i],11,20)
     st1=strtrim(st1,2)
     epoch=fix(st1)       ; epoch as integer

if epoch eq 0 then begin 
     for i=0,n-1 do if (strmid(hd[i],0,5) eq 'EPOCH') then $
          equin=float(strmid(hd[i],11,32))
     epoch=fix(equin)
end

return
end
;================================================================
pro parse_head2,hd,jd,t1,t2,epoch,RA,dec
; 
; a modification of parse_head1.pro to read keywords
;    UTSTART, UTEND
;
; parses the FITS header hd to obtain:
;    Julian Date for the date
;    the starting and ending times t1, t2 in decimal hours
;    the epoch in years
;    RA, dec of the star in decimal hours (RA) and degrees (dec)
; input: hd
; output: jd,t1,t2,epoch,RA,dec
;
; usage: parse_head2,hs,jd,t1,t2,ep,RA,dec

jd=0.d0    ; declarations
t1=0.d0
t2=0.d0
epoch=0
RA=0.d0
dec=0.d0

st=''

n=n_elements(hd)
; print,strtrim(string(n),2),' lines in the header'

for i=0,n-1 do if (strmid(hd[i],0,8) eq 'DATE-OBS') then $
     st=strmid(hd[i],11,10)
     if strmid(st,2,1) eq '/' then begin    ; old DDO format
       yr=fix(strmid(st,6,2))+1900
       mo=fix(strmid(st,3,2))
       da=fix(strmid(st,0,2))
     end else begin                         ; new DDO format since Y2K
       yr=fix(strmid(st,0,4))
       mo=fix(strmid(st,5,2))
       da=fix(strmid(st,8,2))
     end
     juldate,[yr,mo,da],jd
     jd=2400000.d0+jd

for i=0,n-1 do if (strmid(hd[i],0,8) eq 'TIME-OBS') then $
     st=strmid(hd[i],11,8)
     hh=strmid(st,0,2)
     mm=strmid(st,3,2)
     ss=strmid(st,6,2)
     t1=ten(hh,mm,ss)    ; start time in decimal hours

for i=0,n-1 do if (strmid(hd[i],0,8) eq 'UTSTART ') then $
     st=strmid(hd[i],strpos(hd[i],"'")+1,8)   ; LCO keyword
     hh=strmid(st,0,2)
     mm=strmid(st,3,2)
     ss=strmid(st,6,2)
     t1=ten(hh,mm,ss)          ; start time in decimal hours

for i=0,n-1 do if (strmid(hd[i],0,8) eq 'ENDTIME ') then $
     st=strmid(hd[i],11,8)
     hh=strmid(st,0,2)
     mm=strmid(st,3,2)
     ss=strmid(st,6,2)
     t2=ten(hh,mm,ss)    ; stop time in decimal hours

for i=0,n-1 do if (strmid(hd[i],0,8) eq 'UTEND   ') then $
     st=strmid(hd[i],strpos(hd[i],':')-2,8)     ; LCO keyword
     hh=strmid(st,0,2)
     mm=strmid(st,3,2)
     ss=strmid(st,6,2)
     t2=ten(hh,mm,ss)    ; stop time in decimal hours

for i=0,n-1 do if (strmid(hd[i],0,8) eq 'RA      ') then $
     st=strmid(hd[i],11,11)
     hh=strmid(st,0,2)
     mm=strmid(st,3,2)
     ss=strmid(st,6,5)
     RA=ten(hh,mm,ss)    ; RA in decimal hours

for i=0,n-1 do if (strmid(hd[i],0,8) eq 'DEC     ') then $
     st=strmid(hd[i],11,12)
     if strmid(st,0,1) eq '-' then begin 
        hh=strmid(st,0,3)
        mm=strmid(st,4,2)
        ss=strmid(st,7,5)
     end else begin
        hh=strmid(st,0,2)
        mm=strmid(st,3,2)
        ss=strmid(st,6,5)
     end
; print,hh,' ',mm,' ',ss
     dec=ten(hh,mm,ss)   ; dec in decimal degrees

for i=0,n-1 do if (strmid(hd[i],0,8) eq 'EQUINOX ') then $
     st=strmid(hd[i],11,20)
     st=strtrim(st,2)
     epoch=fix(st)       ; epoch as integer

return
end
;================================================================
pro plot_1star,star,tp,res,title
; plot of the sine fit for one star
; 
; Input: star - data array [3,*]
;        tp   - 2-el vector with T0,P
;        res - 6-el vector:
;               gam, K1, 0 - parameters
;                 or with added corrections      
;               (3-el vector of errors disregarded)
; Usage: plot_1star,sv,sv0,sv1,'SV Equ'   ; no corrections
;        plot_1star,sv,sv0,sv1+sv2,'SV Equ' ; with corr's
;
  
  ph = (star[0,*]-tp[0]-res[2])/tp[1]  
;    res[2] is the correction to T0
  ph = reform(ph)
  ph = ph - floor(ph)
  n = n_elements(ph)
  amp = abs(res[0])+abs(res[1])         ; rough range in y-coord
;  amp=350.

  good1=where(star[2,*] gt 0.)   ; weight > 0.
;  good2=where(star[4,*] gt 0.)
  full1=where(star[2,*] gt 0.5)  ; full weight
;  full2=where(star[4,*] gt 0.5)
  bad1=where(star[2,*] eq 0.)
;  bad2=where(star[4,*] eq 0.)

  set_ps
  device,/port
  plotsym,0,1.25             ; 1st star - circles
  plot,ph[good1],star[1,good1],psym=8,    $
      xran=[-0.1,+1.1],        $
      yran=[-amp,amp]*2,     $
      tit='!17 '+title
  oplot,ph[good1]-1,star[1,good1],psym=8
  oplot,ph[good1]+1,star[1,good1],psym=8
  plotsym,0,1.25,/fill
  oplot,ph[full1],star[1,full1],psym=8
  oplot,ph[full1]-1,star[1,full1],psym=8
  oplot,ph[full1]+1,star[1,full1],psym=8

;  plotsym,4,1.25              ; 2nd star - triangles
;  oplot,ph[good2],star[3,good2],psym=8
;  oplot,ph[good2]-1,star[3,good2],psym=8
;  oplot,ph[good2]+1,star[3,good2],psym=8
;  oplot,[-0.5,+1.5],[res[0],res[0]]
;  plotsym,4,1.25,/fill
;  oplot,ph[full2],star[3,full2],psym=8
;  oplot,ph[full2]-1,star[3,full2],psym=8
;  oplot,ph[full2]+1,star[3,full2],psym=8

  if (bad1[0] ne -1) then begin   ; marks of bad obs
    k = n_elements(bad1)
    for i=0,k-1 do oplot,[ph[bad1[i]],ph[bad1[i]]],[-1.,-0.9]*amp
  end
;  if (bad2[0] ne -1) then begin
;    k = n_elements(bad2)
;    for i=0,k-1 do oplot,[ph[bad2[i]],ph[bad2[i]]],[-1.,-0.9]*amp
;  end

  oplot,[-0.5,+1.5],[res[0],res[0]]

  x = findgen(200)/200
  oplot,x,res[0]-sin(2*!pi*x)*res[1]  ; sense as in A-type
  oplot,x-1,res[0]-sin(2*!pi*x)*res[1]
  oplot,x+1,res[0]-sin(2*!pi*x)*res[1]
;  oplot,x,res[0]+sin(2*!pi*x)*res[2]
;  oplot,x-1,res[0]+sin(2*!pi*x)*res[2]
;  oplot,x+1,res[0]+sin(2*!pi*x)*res[2]
  device,/close

return
end
;================================================================
pro plot_2star,star,tp,res,title
; plot of the sine fit to two stars
; 
; Input: star - data array [5,*]
;        tp   - 2-el vector with T0,P
;        res - 8-el vector:
;               gam, K1, K2, 0 - parameters
;                 or with added corrections      
;               (4-el vector of errors disregarded)
; Usage: plot_2star,ah,ah0,ah1,'AH Aur'   ; no corrections
;        plot_2star,ah,ah0,ah1+ah2,'AH Aur' ; with corr's
;
  
  star=double(star)
  tp=double(tp)

  ph = (star[0,*]-tp[0]-res[3])/tp[1]  
;    res[3] is the correction to T0
  ph = reform(ph)
  ph = ph - floor(ph)
  n = n_elements(ph)
  amp = res[0]+res[1]         ; rough range in y-coord
  amp = amp > res[0]+res[2]

  good1=where(star[2,*] gt 0.)   ; weight > 0.
  good2=where(star[4,*] gt 0.)
  full1=where(star[2,*] gt 0.5)  ; full weight
  full2=where(star[4,*] gt 0.5)
  bad1=where(star[2,*] eq 0.)
  bad2=where(star[4,*] eq 0.)

  set_ps
  device,/port
  plotsym,0,1.25             ; 1st star - circles
  plot,ph[good1],star[1,good1],psym=8,    $
      xran=[-0.1,+1.1],        $
;      yran=[-amp,amp]*1.1,     $
      yran=[-350,350],         $
      xmin=-1,ymin=-1,         $
      tit='!17 '+title
  oplot,ph[good1]-1,star[1,good1],psym=8
  oplot,ph[good1]+1,star[1,good1],psym=8
  plotsym,0,1.25,/fill
  oplot,ph[full1],star[1,full1],psym=8
  oplot,ph[full1]-1,star[1,full1],psym=8
  oplot,ph[full1]+1,star[1,full1],psym=8

  plotsym,4,1.25              ; 2nd star - triangles
  oplot,ph[good2],star[3,good2],psym=8
  oplot,ph[good2]-1,star[3,good2],psym=8
  oplot,ph[good2]+1,star[3,good2],psym=8
  oplot,[-0.5,+1.5],[res[0],res[0]]
  plotsym,4,1.25,/fill
  oplot,ph[full2],star[3,full2],psym=8
  oplot,ph[full2]-1,star[3,full2],psym=8
  oplot,ph[full2]+1,star[3,full2],psym=8

  if (bad1[0] ne -1) then begin   ; marks of bad obs
    k = n_elements(bad1)
    for i=0,k-1 do oplot,[ph[bad1[i]],ph[bad1[i]]],[-1.,-0.9]*amp
  end
  if (bad2[0] ne -1) then begin
    k = n_elements(bad2)
    for i=0,k-1 do oplot,[ph[bad2[i]],ph[bad2[i]]],[-1.,-0.9]*amp
  end

  if (bad1[0] ne -1) then begin   ; marks of bad obs
    k = n_elements(bad1)
    for i=0,k-1 do oplot,[ph[bad1[i]],ph[bad1[i]]]-1,[-1.,-0.9]*amp
  end
  if (bad2[0] ne -1) then begin
    k = n_elements(bad2)
    for i=0,k-1 do oplot,[ph[bad2[i]],ph[bad2[i]]]-1,[-1.,-0.9]*amp
  end

  if (bad1[0] ne -1) then begin   ; marks of bad obs
    k = n_elements(bad1)
    for i=0,k-1 do oplot,[ph[bad1[i]],ph[bad1[i]]]+1,[-1.,-0.9]*amp
  end
  if (bad2[0] ne -1) then begin
    k = n_elements(bad2)
    for i=0,k-1 do oplot,[ph[bad2[i]],ph[bad2[i]]]+1,[-1.,-0.9]*amp
  end


  oplot,[-0.5,+1.5],[res[0],res[0]]

  x = findgen(200)/200
  oplot,x,res[0]-sin(2*!pi*x)*res[1]  ; sense as in A-type
  oplot,x-1,res[0]-sin(2*!pi*x)*res[1]
  oplot,x+1,res[0]-sin(2*!pi*x)*res[1]
  oplot,x,res[0]+sin(2*!pi*x)*res[2]
  oplot,x-1,res[0]+sin(2*!pi*x)*res[2]
  oplot,x+1,res[0]+sin(2*!pi*x)*res[2]
  device,/close

return
end
;================================================================
pro plot_vel,phase,hvc,rvpm,rvsm,wp,ws,star
; 
; usage: plot_vel,phase,hvc,rvpm1,rvsm1,weightp,weights,'V401 Cyg'
;
; plots measured velocities corrected to the template


set_ps

n=n_elements(phase)

tit='!17'+star

plotsym,0,1.0,/fill
w=where(wp eq 1.)
plot,phase[w],rvpm[w]+hvc[w]-hvc[0],   $
     tit=tit,xran=[0,1],yran=[-500,+500],psym=8
plotsym,0,1.0
w=where(wp eq 0.5)
if w[0] ne -1 then oplot,phase[w],rvpm[w]+hvc[w]-hvc[0],psym=8

plotsym,4,1.0,/fill
w=where(ws eq 1.)
if w[0] ne -1 then oplot,phase[w],rvsm[w]+hvc[w]-hvc[0],psym=8
plotsym,4,1.0
w=where(ws eq 0.5)
if w[0] ne -1 then oplot,phase[w],rvsm[w]+hvc[w]-hvc[0],psym=8

device,/close

return
end
;================================================================
pro plotBFphase,bf,phase,vel,hvc,nr1,nr2,star
; 
; usage: plotBFphase,bf20,phase,vel,hvc,2,10,'V401 Cyg'
;  
; plots in relative heliocentric velocities the selected BF's
;    for comparison; select by the index #: nr1, nr2

set_ps

n=n_elements(phase)
if nr1 lt 2 then begin print,'nr1 must be >= 2' & goto,E & end
if nr2 gt n-1 then begin print,'nr2 must be <=',n-1 & goto,E & end

tit='!17'+star+': BF for '+strtrim(string(nr1),2)+' - '+strtrim(string(nr2),2)

plot,vel+hvc[1]-hvc[0],bf[1,*],tit=tit,xran=[-700,700],yran=[-0.2,1.2]

for i=nr1,nr2 do oplot,vel+hvc[i]-hvc[0],bf[i,*]

oplot,vel,bf[0,*],line=2

device,/close

E: 
return
end
;================================================================
function premat, equinox1, equinox2, FK4 = FK4
;+
; NAME:
;       PREMAT
; PURPOSE:
;       Return the precession matrix needed to go from EQUINOX1 to EQUINOX2.  
; EXPLANTION:
;       This matrix is used by the procedures PRECESS and BARYVEL to precess 
;       astronomical coordinates
;
; CALLING SEQUENCE:
;       matrix = PREMAT( equinox1, equinox2, [ /FK4 ] )
;
; INPUTS:
;       EQUINOX1 - Original equinox of coordinates, numeric scalar.  
;       EQUINOX2 - Equinox of precessed coordinates.
;
; OUTPUT:
;      matrix - double precision 3 x 3 precession matrix, used to precess
;               equatorial rectangular coordinates
;
; OPTIONAL INPUT KEYWORDS:
;       /FK4   - If this keyword is set, the FK4 (B1950.0) system precession
;               angles are used to compute the precession matrix.   The 
;               default is to use FK5 (J2000.0) precession angles
;
; EXAMPLES:
;       Return the precession matrix from 1950.0 to 1975.0 in the FK4 system
;
;       IDL> matrix = PREMAT( 1950.0, 1975.0, /FK4)
;
; PROCEDURE:
;       FK4 constants from "Computational Spherical Astronomy" by Taff (1983), 
;       p. 24. (FK4). FK5 constants from "Astronomical Almanac Explanatory
;       Supplement 1992, page 104 Table 3.211.1.
;
; REVISION HISTORY
;       Written, Wayne Landsman, HSTX Corporation, June 1994
;       Converted to IDL V5.0   W. Landsman   September 1997
;-    
  On_error,2                                           ;Return to caller

  npar = N_params()

   if ( npar LT 2 ) then begin 

     print,'Syntax - PREMAT, equinox1, equinox2, /FK4]'
     return,-1 

  endif 

  deg_to_rad = !DPI/180.0d
  sec_to_rad = deg_to_rad/3600.d0

   t = 0.001d0*( equinox2 - equinox1)

 if not keyword_set( FK4 )  then begin  
           st = 0.001d0*( equinox1 - 2000.d0)
;  Compute 3 rotation angles
           A = sec_to_rad * T * (23062.181D0 + ST*(139.656D0 +0.0139D0*ST) $
            + T*(30.188D0 - 0.344D0*ST+17.998D0*T))

           B = sec_to_rad * T * T * (79.280D0 + 0.410D0*ST + 0.205D0*T) + A

        C = sec_to_rad * T * (20043.109D0 - ST*(85.33D0 + 0.217D0*ST) $
              + T*(-42.665D0 - 0.217D0*ST -41.833D0*T))

 endif else begin  

           st = 0.001d0*( equinox1 - 1900.d0)
;  Compute 3 rotation angles

           A = sec_to_rad * T * (23042.53D0 + ST*(139.75D0 +0.06D0*ST) $
            + T*(30.23D0 - 0.27D0*ST+18.0D0*T))

           B = sec_to_rad * T * T * (79.27D0 + 0.66D0*ST + 0.32D0*T) + A

           C = sec_to_rad * T * (20046.85D0 - ST*(85.33D0 + 0.37D0*ST) $
              + T*(-42.67D0 - 0.37D0*ST -41.8D0*T))

 endelse  

  sina = sin(a) &  sinb = sin(b)  & sinc = sin(c)
  cosa = cos(a) &  cosb = cos(b)  & cosc = cos(c)

  r = dblarr(3,3)
  r[0,0] = [ cosa*cosb*cosc-sina*sinb, sina*cosb+cosa*sinb*cosc,  cosa*sinc]
  r[0,1] = [-cosa*sinb-sina*cosb*cosc, cosa*cosb-sina*sinb*cosc, -sina*sinc]
  r[0,2] = [-cosb*sinc, -sinb*sinc, cosc]

  return,r
  end
;================================================================
pro read_1star,file,star
; reading in of obs file: one star only visible
; use: read_1star,'SV_Equ.prn',SV
;
; input file must have format: time,RV,weight


openr,1,file
star=dblarr(5,200)

d0=1.d0 & d1=d0 & d2=d0
i=0
while not eof(1) do begin
   readf,1,form='(f12.0,2f7.0)',d0,d1,d2
   star(0,i)=d0 & star(1,i)=d1 & star(2,i)=d2
   i = i+1
   end
close,1
star=star(*,0:i-1)

return
end
;================================================================
pro read_2star,file,star
; reads in of obs files for two stars
; use: read_2star,'AH_Aur.prn',AH_Aur

openr,1,file
star=dblarr(5,200)

d0=1.d0 & d1=d0 & d2=d0 & d3=d0 & d4=d0
i=0
while not eof(1) do begin
   readf,1,form='(f12.0,4f7.0)',d0,d1,d2,d3,d4
   star(0,i)=d0 & star(1,i)=d1 & star(2,i)=d2
   star(3,i)=d3 & star(4,i)=d4
   i = i+1
   end
close,1
star=star(*,0:i-1)

return
end
;================================================================
pro read_3star,file,star
; reads in of obs files for 3rd star from Morbey's style solutions
;      
; use: read_3star,'SWLyn3.prn',SW3

; format: JD, Vobs, Vpred, phase, weight
;       f12.0, 4f7.0

openr,1,file
star=dblarr(5,200)

d0=1.d0 & d1=d0 & d2=d0 & d3=d0 & d4=d0
i=0
while not eof(1) do begin
   readf,1,form='(f12.0,4f7.0)',d0,d1,d2,d3,d4
   star(0,i)=d0 & star(1,i)=d1 & star(2,i)=d2
   star(3,i)=d3 & star(4,i)=d4
   i = i+1
   end
close,1
star=star(*,0:i-1)

return
end
;================================================================
; MOTE: Readfits.pro uses many auxiliary routines
;
;+
; NAME:
;       READFITS
; PURPOSE:
;       Read a FITS file into IDL data and header variables. 
; EXPLANATION:
;       Under Unix, READFITS() can also read gzip or Unix compressed FITS files.
;       See http://idlastro.gsfc.nasa.gov/fitsio.html for other ways of
;       reading FITS files with IDL.
;
; CALLING SEQUENCE:
;       Result = READFITS( Filename,[ Header, heap, /NOSCALE, EXTEN_NO=,
;                       NSLICE=, /SILENT , NaNVALUE =, STARTROW =, NUMROW = ,
;                       /No_Unsigned ] )
;
; INPUTS:
;       FILENAME = Scalar string containing the name of the FITS file  
;                 (including extension) to be read.   If the filename has
;                  a *.gz extension, it will be treated as a gzip compressed
;                  file.   If it has a .Z extension, it will be treated as a
;                  Unix compressed file.
;
; OUTPUTS:
;       Result = FITS data array constructed from designated record.
;                If the specified file was not found, then Result = -1
;
; OPTIONAL OUTPUT:
;       Header = String array containing the header from the FITS file.
;       heap = For extensions, the optional heap area following the main
;              data array (e.g. for variable length binary extensions).
;
; OPTIONAL INPUT KEYWORDS:
;
;       EXTEN_NO - scalar integer specify the FITS extension to read.  For
;               example, specify EXTEN = 1 or /EXTEN to read the first 
;               FITS extension.    Extensions are read using recursive
;               calls to READFITS.
;
;       NaNVALUE - This scalar is only needed on architectures (such as VMS
;               prior to IDL V5.1) that do not recognize the IEEE "not a number"
;               (NaN) convention.   It specifies the value to translate any IEEE
;               NAN values in the FITS data array.  
;   
;       /NOSCALE - If present and non-zero, then the ouput data will not be
;                scaled using the optional BSCALE and BZERO keywords in the 
;                FITS header.   Default is to scale.
;
;       /NO_UNSIGNED - By default, if theIDL Version is 5.2 or greater, and the
;               header indicates an unsigned integer (BITPIX = 16, BZERO=2^15,
;               BSCALE=1) then FITS_READ will output an IDL unsigned integer 
;               data type (UINT).   But if /NO_UNSIGNED is set, or the IDL 
;               version is before 5.2, then the data is converted to type LONG.  
;
;       NSLICE - An integer scalar specifying which N-1 dimensional slice of a 
;                N-dimensional array to read.   For example, if the primary 
;                image of a file 'wfpc.fits' contains a 800 x 800 x 4 array, 
;                then 
;
;                 IDL> im = readfits('wfpc.fits',h, nslice=2)
;                           is equivalent to 
;                 IDL> im = readfits('wfpc.fits',h)
;                 IDL> im = im[*,*,2]
;                 but the use of the NSLICE keyword is much more efficient.
;
;       NUMROW -  Scalar non-negative integer specifying the number of rows 
;                 of the image or table to read.   Useful when one does not 
;                 want to read the entire image or table.
;
;       POINT_LUN  -  Position (in bytes) in the FITS file at which to start
;                 reading.   Useful if READFITS is called by another procedure
;                 which needs to directly read a FITS extension.    Should 
;                 always be a multiple of 2880, and not be used with EXTEN_NO
;                 keyword.
;
;       /SILENT - Normally, READFITS will display the size the array at the
;                 terminal.  The SILENT keyword will suppress this
;
;        STARTROW - Non-negative integer scalar specifying the row
;               of the image or extension table at which to begin reading. 
;               Useful when one does not want to read the entire table.
;
; EXAMPLE:
;       Read a FITS file test.fits into an IDL image array, IM and FITS 
;       header array, H.   Do not scale the data with BSCALE and BZERO.
;
;              IDL> im = READFITS( 'test.fits', h, /NOSCALE)
;
;       If the file contain a FITS extension, it could be read with
;
;              IDL> tab = READFITS( 'test.fits', htab, /EXTEN )
;
;       The function TBGET() can be used for further processing of a binary 
;       table, and FTGET() for an ASCII table.
;       To read only rows 100-149 of the FITS extension,
;
;              IDL> tab = READFITS( 'test.fits', htab, /EXTEN, 
;                                   STARTR=100, NUMR = 50 )
;
;       To read in a file that has been compressed:
;
;              IDL> tab = READFITS('test.fits.gz',h)
;
; ERROR HANDLING:
;       If an error is encountered reading the FITS file, then 
;               (1) the system variable !ERROR is set (via the MESSAGE facility)
;               (2) the error message is displayed (unless /SILENT is set),
;                   and the message is also stored in !ERR_STRING
;               (3) READFITS returns with a value of -1
; RESTRICTIONS:
;       (1) Cannot handle random group FITS
;
; NOTES:
;       (1) If data is stored as integer (BITPIX = 16 or 32), and BSCALE
;       and/or BZERO keywords are present, then the output array is scaled to 
;       floating point (unless /NOSCALE is present) using the values of BSCALE
;       and BZERO.   In the header, the values of BSCALE and BZERO are then 
;       reset to 1. and 0., while the original values are written into the 
;       new keywords O_BSCALE and O_BZERO.     If the BLANK keyword was
;       present, then any input integer values equal to BLANK in the input
;       integer image are unchanged by BSCALE or BZERO
;       
;       (2) The use of the NSLICE keyword is incompatible with the NUMROW
;       or STARTROW keywords.
;
;       (3) READFITS() underwent a substantial rewrite in October 1998 to 
;       eliminate recursive calls, and improve efficiency when reading
;       extensions.
;            1. The NUMROW and STARTROW keywords can now be used when reading
;              a primary image (extension = 0).
;            2. There is no error check for moving past the end of file when
;               reading the data array.
;
;       (4) On some Unix shells, one may get a "Broken pipe" message if reading
;        a compressed file, and not reading to the end of the file (i.e. the
;        decompression has not gone to completion).     This is an informative
;        message only, and should not affect the output of READFITS.   
; PROCEDURES USED:
;       Functions:   IS_IEEE_BIG(), SXPAR(), WHERENAN()
;       Procedures:  IEEE_TO_HOST, SXADDPAR, SXDELPAR
;
; MODIFICATION HISTORY:
;       Original Version written in 1988, W.B. Landsman   Raytheon STX
;       Revision History prior to June 1997 removed          
;       Recognize BSCALE, BZERO in IMAGE extensions             WBL Jun-97
;       Added NSLICE keyword                                    WBL Jul-97
;       Added ability to read heap area after extensions        WBL Aug-97      
;       Suppress *all* nonfatal messages with /SILENT           WBL Dec-97
;       Converted to IDL V5.0                                   WBL Dec-1997
;       Fix NaN assignment for int data       C. Gehman/JPL     Mar-98
;       Fix bug with NaNvalue = 0.0           C. Gehman/JPL     Mar-98
;       Major rewrite to eliminate recursive calls when reading extensions
;                  W.B. Landsman  Raytheon STX                    October 1998
;       Add /binary modifier needed for Windows  W. Landsman    April 1999
;       Read unsigned datatypes, added /no_unsigned   W. Landsman December 1999
;       Output BZERO = 0 for unsigned data types   W. Landsman   January 2000
;       Open with /swap_if_little_endian if since V5.1 W. Landsman February 2000
;       Fixed logic error when using NSLICE keyword W. Landsman March 2000
;       Fixed byte swapping problem for compressed files on little endian 
;             machines introduced in Feb 2000     W. Landsman       April 2000
;-
;
function READFITS, filename, header, heap, NOSCALE = noscale, $
                   SILENT = silent, EXTEN_NO = exten_no, NUMROW = numrow, $
                   POINTLUN = pointlun, STARTROW = startrow, $
                   NaNvalue = NaNvalue, NSLICE = nslice, $
                   NO_UNSIGNED = no_unsigned


  On_error,2                    ;Return to user

; Check for filename input

   if N_params() LT 1 then begin                
      print,'Syntax - im = READFITS( filename, [ h, heap, /NOSCALE, /SILENT,'
      print,'                 NaNValue = ,EXTEN_NO =, STARTROW = , NUMROW='
      print,'                 NSLICE = , /No_UnSigned]'
      return, -1
   endif

; Set default keyword values

   silent = keyword_set( SILENT )
   if N_elements(exten_no) EQ 0 then exten_no = 0

;  Check if this is a compressed file.
                
    len = strlen(filename)
    if len gt 3 then tail = strmid(filename, len-3, 3) else tail = ' '
    ucmprs = ''
    if strlowcase(strmid(tail,1,2)) eq '.z' then  $
                  ucmprs = 'uncompress -c '
    if strlowcase(tail) eq '.gz'  then ucmprs = 'gzip -cd '
    gzip = ucmprs NE ''

;  Go to the start of the file.

   if !VERSION.RELEASE GE '5.1' then $
   openr, unit, filename, ERROR=error,/get_lun,/BLOCK,/binary, $
                          /swap_if_little_endian else $
   openr, unit, filename, ERROR=error,/get_lun,/BLOCK,/binary
   if error NE 0 then begin
        message,/con,' ERROR - Unable to locate file ' + filename
        return, -1
   endif

;  Handle compressed files.   On some Unix machines, users might wish to force
;  use of /bin/sh in the line spawn, ucmprs+filename, unit=unit,/sh

        if gzip then begin
                free_lun, unit
                if (!version.os_family EQ 'unix') then begin
                        spawn, ucmprs+filename, unit=unit
                endif else begin
                        print, 'READFITS: Only Unix IDL supports piped spawns'
                        print, '         File must be uncompressed manually'
                        return, -1                      
                endelse
                
        endif else begin
   endelse
         
  if keyword_set(POINTLUN) then begin
       if gzip then  readu,unit,bytarr(pointlun,/nozero) $
               else point_lun, unit, pointlun
  endif

  for ext = 0, exten_no do begin
               
;  Read the next header, and get the number of bytes taken up by the data.

       block = string(replicate(32b,80,36))
       w = [-1]
       header = ' '
       
       while w[0] EQ -1 do begin
          
       if EOF(unit) then begin 
            message,/ CON, $
               'EOF encountered attempting to read extension ' + strtrim(ext,2)
            free_lun,unit
            return,-1
       endif

      readu, unit, block
      w = where(strlen(block) NE 80, Nbad)
      if (Nbad GT 0) then begin
           message,'Warning-Invalid characters in header',/INF,NoPrint=Silent
           block[w] = string(replicate(32b, 80))
      endif
      w = where(strmid(block, 0, 8) eq 'END     ', Nend)
      if Nend EQ 0 then header = [header, block] $
                   else header = [header, block[0:w[0]]]
      endwhile

      header = header[1:*]
      if (ext EQ 0 ) and (keyword_set(pointlun) EQ 0) then $
             if strmid( header[0], 0, 8)  NE 'SIMPLE  ' then begin
              message,/CON, $
                 'ERROR - Header does not contain required SIMPLE keyword'
                free_lun, unit
                return, -1
      endif

                
; Get parameters that determine size of data region.
                
       bitpix =  sxpar(header,'BITPIX')
       naxis  = sxpar(header,'NAXIS')
       gcount = sxpar(header,'GCOUNT') > 1
       pcount = sxpar(header,'PCOUNT')
                
       if naxis GT 0 then begin 
            dims = sxpar( header,'NAXIS*')           ;Read dimensions
            ndata = dims[0]
            if naxis GT 1 then for i = 2, naxis do ndata = ndata*dims[i-1]
                        
                endif else ndata = 0
                
                nbytes = (abs(bitpix) / 8) * gcount * (pcount + ndata)

;  Move to the next extension header in the file.

      if ext LT exten_no then begin
                nrec = long((nbytes + 2879) / 2880)
                if nrec GT 0 then begin     
                if gzip then begin 
                        buf = bytarr(nrec*2880L,/nozero)
                        readu,unit,buf 
                        endif else  begin 
                        point_lun, -unit,curr_pos
                        point_lun, unit,curr_pos + nrec*2880L
                endelse
                endif
       endif
       endfor

 case BITPIX of 
           8:   IDL_type = 1          ; Byte
          16:   IDL_type = 2          ; Integer*2
          32:   IDL_type = 3          ; Integer*4
         -32:   IDL_type = 4          ; Real*4
         -64:   IDL_type = 5          ; Real*8
        else:   begin
                message,/CON, 'ERROR - Illegal value of BITPIX (= ' +  $
                strtrim(bitpix,2) + ') in FITS header'
                free_lun,unit
                return, -1
                end
  endcase     

; Check for dummy extension header

 if Naxis GT 0 then begin 
        Nax = sxpar( header, 'NAXIS*' )   ;Read NAXES
        ndata = nax[0]
        if naxis GT 1 then for i = 2, naxis do ndata = ndata*nax[i-1]

  endif else ndata = 0

  nbytes = (abs(bitpix)/8) * gcount * (pcount + ndata)
 
  if nbytes EQ 0 then begin
        if not SILENT then message, $
                "FITS header has NAXIS or NAXISi = 0,  no data array read",/CON
        free_lun, unit
        return,-1
 endif

; Check for FITS extensions, GROUPS

 groups = sxpar( header, 'GROUPS' ) 
 if groups then message,NoPrint=Silent, $
           'WARNING - FITS file contains random GROUPS', /INF

; If an extension, did user specify row to start reading, or number of rows
; to read?

   if not keyword_set(STARTROW) then startrow = 0
   if naxis GE 2 then nrow = nax[1] else nrow = ndata
   if not keyword_set(NUMROW) then numrow = nrow
   if exten_no GT 0 then begin
        xtension = strtrim( sxpar( header, 'XTENSION' , Count = N_ext),2)
        if N_ext EQ 0 then message, /INF, NoPRINT = Silent, $
                'WARNING - Header missing XTENSION keyword'
   endif 
      
   if (exten_no GT 0) and ((startrow NE 0) or (numrow NE nrow)) then begin
        if startrow GE nax[1] then begin
           message,'ERROR - Specified starting row ' + strtrim(startrow,2) + $
          ' but only ' + strtrim(nax[1],2) + ' rows in extension',/CON
           free_lun,unit
           return,-1
        endif 
        nax[1] = nax[1] - startrow    
        nax[1] = nax[1] < numrow
        sxaddpar, header, 'NAXIS2', nax[1]
        if gzip then begin
                if startrow GT 0 then begin
                        tmp=bytarr(startrow*nax[0])
                        readu,unit,tmp
                endif 
        endif else begin 
              point_lun, -unit, pointlun          ;Current position
              point_lun, unit, pointlun + startrow*nax[0]
    endelse
    endif else if (N_elements(NSLICE) EQ 1) then begin
        lastdim = nax[naxis-1]
        if nslice GE lastdim then message,/CON, $
        'ERROR - Value of NSLICE must be less than ' + strtrim(lastdim,2)
        nax = nax[0:naxis-2]
        sxdelpar,header,'NAXIS' + strtrim(naxis,2)
        naxis = naxis-1
        sxaddpar,header,'NAXIS',naxis
        ndata = ndata/lastdim
        nskip = nslice*ndata*abs(bitpix/8) 
        if gzip then begin 
            if (Ndata GT 0) then  begin 
                  buf = bytarr(nskip,/nozero)
                  readu,unit,buf   
             endif
        endif else begin 
                   point_lun, -unit, currpoint          ;Current position
                   point_lun, unit, currpoint + nskip
        endelse
  endif


  if not (SILENT) then begin   ;Print size of array being read

         if exten_no GT 0 then message, $
                     'Reading FITS extension of type ' + xtension, /INF
         snax = strtrim(NAX,2)
         st = snax[0]
         if Naxis GT 1 then for I=1,NAXIS-1 do st = st + ' by '+SNAX[I] $
                            else st = st + ' element'
         st = 'Now reading ' + st + ' array'
         if (exten_no GT 0) and (pcount GT 0) then st = st + ' + heap area'
         message,/INF,st   
   endif

; Read Data in a single I/O call

    data = make_array( DIM = nax, TYPE = IDL_type, /NOZERO)

     readu, unit, data
    if (exten_no GT 0) and (pcount GT 0) then begin
        theap = sxpar(header,'THEAP')
        skip = theap - N_elements(data)
        if skip GT 0 then begin 
                temp = bytarr(skip,/nozero)
                readu, unit, skip
        endif
        heap = bytarr(pcount*gcount*abs(bitpix)/8)
        readu, unit, heap
    endif

    free_lun, unit

; If necessary, replace NaN values, and convert to host byte ordering
; Byte ordering is needed if 
        
   check_NaN = (bitpix LT 0) and ( N_elements(NaNvalue) GT 0 )
   if check_NaN then NaNpts = whereNaN( data, Count)
   if !VERSION.RELEASE LT '5.1' or gzip then $ 
          if not is_ieee_big() then ieee_to_host, data
   if check_NaN then $
        if ( Count GT 0 ) then data[ NaNpts] = NaNvalue

; Scale data unless it is an extension, or /NOSCALE is set
; Use "TEMPORARY" function to speed processing.  

   do_scale = not keyword_set( NOSCALE )
   if (do_scale and (exten_no GT 0)) then do_scale = xtension EQ 'IMAGE' 
   if do_scale then begin

          Nblank = 0
          if bitpix GT 0 then begin
                blank = sxpar( header, 'BLANK', Count = N_blank) 
                if N_blank GT 0 then $ 
                        blankval = where( data EQ blank, Nblank)
          endif

          Bscale = float( sxpar( header, 'BSCALE' , Count = N_bscale))
          Bzero = float( sxpar(header, 'BZERO', Count = N_Bzero ))
 
; Check for unsigned integer (BZERO = 2^15) or unsigned long (BZERO = 2^31)
; Write uint(32768) rather than 32768U to allow compilation prior to V5.2

          if not keyword_set(No_Unsigned) and $
                 !VERSION.RELEASE GE '5.2' then begin
            no_bscale = (Bscale EQ 1) or (N_bscale EQ 0)
            unsgn_int = (bitpix EQ 16) and (Bzero EQ 32768) and no_bscale
            unsgn_lng = (bitpix EQ 32) and (Bzero EQ 2147483648) and no_bscale
            unsgn = unsgn_int or unsgn_lng
           endif else unsgn = 0

          if unsgn then begin
                 sxaddpar, header, 'BZERO', 0
                 sxaddpar, header, 'O_BZERO', bzero, $
                          'Original Data is unsigned Integer'
                   if unsgn_int then $ 
                        data =  uint(data) - uint(32768) else $
                   if unsgn_lng then  data = ulong(data) - ulong(2147483648) 
                
          endif else begin
 
          if N_Bscale GT 0  then $ 
               if ( Bscale NE 1. ) then begin
                   data = temporary(data) * Bscale 
                   sxaddpar, header, 'BSCALE', 1.
                   sxaddpar, header, 'O_BSCALE', Bscale,' Original BSCALE Value'
               endif

         if N_Bzero GT 0  then $
               if (Bzero NE 0) then begin
                     data = temporary( data ) + Bzero
                     sxaddpar, header, 'BZERO', 0.
                     sxaddpar, header, 'O_BZERO', Bzero,' Original BZERO Value'
               endif
        
        endelse

        if (Nblank GT 0) and ((N_bscale GT 0) or (N_Bzero GT 0)) then $
                data[blankval] = blank

        endif

; Return array

        return, data    
 end 
;=============================================================================
pro rectif_fits,file
; cursor interactive rectification of FITS files
; utilizes rectif_simple.pro
; assumes that the extension is "fits"
;
; input:  file  = FITS file
; output: new file written to the same directory, added "r" to the name
;
; useage:
;    rectif_fits,'C0040424s.fits'
; this produces file: C0040424sr.fits
;
; - use: set_win 
; - enlarge the screen window for better definition
; - work with the cursor from the left-most edge of the spectrum, 
;       and end with setting the cursor to the right of the right frame

     im=readfits(file,h,/silent)
     rectif_simple,im,im1

     file1=strcompress(file)     ; name of new file
     k=strpos(file1,'.fits')     ; if extension .fts, then edit these 2 lines
     file1=strmid(file1,0,k)+'r.fits'

     writefits, file1,im1,h
     print,'New file written:    ',file1
return
end
;=============================================================================
pro rectif_simple,s,s1
; rectification of spectra
; continuum at the cursor points
; 
;    s  = input spectrum
;    s1 = output spectrum, rectified
; use:
;    rectif_simple,s,s1
;
; for FITS rectification (input & output in FITS) 
;    use rectif_fits.pro

x=0. & y=0.
ans=''
xx=fltarr(100)
yy=fltarr(100)

n=n_elements(s)
x1=findgen(n)
y1=fltarr(n)

plot,s,tit='Mark continuum left->right, exit beyond right edge'

k=0
R:
   ; measuring loop
cursor,x,y,/down,/data
  x = round(x)
  if x ge n then goto,E  
           ; get out of the loop out of right frame
  xx[k]=x
  yy[k]=y
;  print,k,xx[k],yy[k]
  oplot,xx[0:k],yy[0:k]
  k=k+1
if k ge 99 then goto,E
goto,R    

E:
xx=xx[0:k-1]
yy=yy[0:k-1]

xx=[0,xx]       ; ends added
yy=[yy[0],yy]
xx=[xx,n-1]
yy=[yy,yy[k]]

;print,xx,yy
y1=interpol(yy,xx,x1)

read,'Result (press any key): ',ans
s1=s/y1
plot,s1

return
end
;=======================================================================
function Rot_one,x,bf,a,da
; Fits a rotational profile to a given BF for one star; gives 
; corrections da to the parameters a. Must be iterated: a=a+da.
; Input:
;    x = velocity (X- axis)
;    bf = broadening function (Y-axis)
;    a = four parameters a[0:3]: at input approx values
;        a[0] = strength, 
;        a[1] = position, 
;        a[2] = width, 
;        a[3] = baseline level.
; Output:
;    da = corrections to a;
;    y = fit, very aproximate, 
;  
; Usage:
;    y=Rot_one(x,bf,a,da)
;        
; assumed limb darkening u=0.75 (not crucial, can be changed)
;

u=0.75

m=n_elements(x)

c1=1-u                ; auxiliary quantities
c2=!pi*u/4

h1=fltarr(m)
w=where(abs(x-a[1]) le a[2])
h1[w]=1.-((x[w]-a[1])/a[2])^2 
nw=n_elements(w)      ; make sure h1 is not zero at ends
if(h1[w[nw-1]] eq 0.) then w=w[0:nw-2]
if(h1[w[0]] eq 0.) then w=w[1:*]

p0=fltarr(m) & p1=p0 & p2=p0  ; declare derivatives and set to zero
p0=(c1*sqrt(h1)+c2*h1)/(c1+c2)   ; profile shape
p1[w]=2.*a[0]/a[2]^2*(x[w]-a[1])*(c1/2/sqrt(h1[w])+c2)/(c1+c2) 
p2[w]=p1[w]*(x[w]-a[1])/a[2] 
dev=bf-(a[0]*p0+a[3])            ; deviations

t=fltarr(4,m)   ; design matrix for LSQ solution
t[0,*]=p0
t[1,*]=p1
t[2,*]=p2
t[3,*]=1.

svdc,t,w,u,v,/double    ; LSQ solution
da=svsol(u,w,v,dev,/double)

a1=a+da 
return,a1[0]*p0+a1[3]
end
;=======================================================================
function Rot_two,x,bf,a,da
; Fits a rotational profile to a given BF for two stars; gives 
; corrections da to the parameters a. Must be iterated: a=a+da.
; Input:
;    x = velocity (X- axis)
;    bf = broadening function (Y-axis)
;    a = seven parameters a[0:6]: at input approx values
;        a[0] = strength 1,
;        a[1] = position 1, 
;        a[2] = width 1, 
;        a[3] = strength 2, 
;        a[4] = position 2, 
;        a[5] = width 2, 
;        a[6] = baseline level.
; Output:
;    da = corrections to a;
;    y = fit, very aproximate, 
; 
; Usage:
;    y=Rot_Two(x,bf,a,da)
;        
; assumed limb darkening u=0.75 (not crucial, can be changed or changed
;    into a parameter)
;

u=0.75

m=n_elements(x)

c1=1-u                ; auxiliary quantities
c2=!pi*u/4

h1=fltarr(m)
w1=where(abs(x-a[1]) le a[2])
h1[w1]=1.-((x[w1]-a[1])/a[2])^2 
nw=n_elements(w1)      ; make sure h1 is not zero at ends
if(h1[w1[nw-1]] eq 0.) then w1=w1[0:nw-2]
if(h1[w1[0]] eq 0.) then w1=w1[1:*]

h2=fltarr(m)
w2=where(abs(x-a[4]) le a[5])
h2[w2]=1.-((x[w2]-a[4])/a[5])^2 
nw=n_elements(w2)      ; make sure h2 is not zero at ends
if(h2[w2[nw-1]] eq 0.) then w2=w2[0:nw-2]
if(h2[w2[0]] eq 0.) then w2=w2[1:*]

p0=fltarr(m) & p1=p0 & p2=p0  ; declare derivatives comp.#1
p0=(c1*sqrt(h1)+c2*h1)/(c1+c2)   ; shape & strength #1
p1[w1]=2.*a[0]/a[2]^2*(x[w1]-a[1])*(c1/2/sqrt(h1[w1])+c2)/(c1+c2) ; pos. #1  
p2[w1]=p1[w1]*(x[w1]-a[1])/a[2]  ; width #1

p3=fltarr(m) & p4=p3 & p5=p3  ; declare derivatives comp.#2
p3=(c1*sqrt(h2)+c2*h2)/(c1+c2)   ; shape & strength #2
p4[w2]=2.*a[3]/a[5]^2*(x[w2]-a[4])*(c1/2/sqrt(h2[w2])+c2)/(c1+c2)  ; pos. #2 
p5[w2]=p4[w2]*(x[w2]-a[4])/a[5]  ; width #2

p=(a[0]*p0)>(a[3]*p3)+a[6] ; where overlap, profile is upper envelope 
dev=bf-p                   ; deviations

w3=where((p2 gt 0.) and (p5 gt 0.))  ; overlap part
if w3[0] ne -1 then begin 
     p0[w3]=0. &  p3[w3]=0.
     p1[w3]=0. &  p4[w3]=0. 
     p2[w3]=0. &  p5[w3]=0. 
end

t=fltarr(7,m)   ; design matrix for LSQ solution
t[0,*]=p0
t[1,*]=p1
t[2,*]=p2
t[3,*]=p3
t[4,*]=p4
t[5,*]=p5
t[6,*]=1.

svdc,t,w,u,v,/double    ; LSQ solution
da=svsol(u,w,v,dev,/double)

a1=a+da 
return,(a1[0]*p0>a1[3]*p3)+a1[6]
end
;========================================================================
function Rot_two1,x,bf,a,da
; Fits a rotational profile to a given BF for two stars; gives 
; corrections da to the parameters a. Must be iterated: a=a+da.
; 
; Version with fixed widths for more stable solutions
;
; Input:
;    x = velocity (X- axis)
;    bf = broadening function (Y-axis)
;    a = seven parameters a[0:6]: at input approx values
;        a[0] = strength 1,
;        a[1] = position 1, 
;        a[2] = width 1,     FIXED
;        a[3] = strength 2, 
;        a[4] = position 2, 
;        a[5] = width 2,     FIXED
;        a[6] = baseline level.
; Output:
;    da = corrections to a;
;    y = fit, very aproximate, 
; 
; Usage:
;    y=Rot_Two(x,bf,a,da)
;        
; assumed limb darkening u=0.75 (not crucial, can be changed)
;

u=0.75

m=n_elements(x)

c1=1-u                ; auxiliary quantities
c2=!pi*u/4

h1=fltarr(m)
w1=where(abs(x-a[1]) le a[2])
h1[w1]=1.-((x[w1]-a[1])/a[2])^2 
nw=n_elements(w1)      ; make sure h1 is not zero at ends
if(h1[w1[nw-1]] eq 0.) then w1=w1[0:nw-2]
if(h1[w1[0]] eq 0.) then w1=w1[1:*]

h2=fltarr(m)
w2=where(abs(x-a[4]) le a[5])
h2[w2]=1.-((x[w2]-a[4])/a[5])^2 
nw=n_elements(w2)      ; make sure h2 is not zero at ends
if(h2[w2[nw-1]] eq 0.) then w2=w2[0:nw-2]
if(h2[w2[0]] eq 0.) then w2=w2[1:*]

p0=fltarr(m) & p1=p0 & p2=p0  ; declare derivatives comp.#1
p0=(c1*sqrt(h1)+c2*h1)/(c1+c2)   ; profile shape
p1[w1]=2.*a[0]/a[2]^2*(x[w1]-a[1])*(c1/2/sqrt(h1[w1])+c2)/(c1+c2) 
p2[w1]=p1[w1]*(x[w1]-a[1])/a[2]
 
p3=fltarr(m) & p4=p3 & p5=p3  ; declare derivatives comp.#2
p3=(c1*sqrt(h2)+c2*h2)/(c1+c2)   ; profile shape
p4[w2]=2.*a[3]/a[5]^2*(x[w2]-a[4])*(c1/2/sqrt(h2[w2])+c2)/(c1+c2) 
p5[w2]=p4[w2]*(x[w2]-a[4])/a[5]
 
p=(a[0]*p0)>(a[3]*p3)+a[6]   ; if overlap, the upper envelope 
dev=bf-p                     ; deviations

w3=where((p2 gt 0.) and (p5 gt 0.))  ; overlap part
if w3[0] ne -1 then begin 
     p0[w3]=0. &  p3[w3]=0.
     p1[w3]=0. &  p4[w3]=0. 
     p2[w3]=0. &  p5[w3]=0. 
end

t=fltarr(5,m)   ; reduced size with fixed a[2] and a[5]; p2 & p5 not used
t[0,*]=p0
t[1,*]=p1
t[2,*]=p3
t[3,*]=p4
t[4,*]=1.

svdc,t,w,u,v,/double    ; LSQ solution
da1=svsol(u,w,v,dev,/double)            ; 5 unknowns
da=[da1[0:1],0.,da1[2:3],0.,da1[4]]     ; 7 corrections

a1=a+da
return,(a1[0]*p0>a1[3]*p3)+a1[6]
end
;========================================================================
function solv,design,right,weight,var,covar
;                               !! variance and covariance
; SVD least-squares solution 
; result given as value of the function
; design=fltarr(small,large) right=fltarr(large) (=dimensions)
; weight multiplies each row of the design array
; 
  right=reform(right)
  n=n_elements(design(*,0))
  m=n_elements(design(0,*))
  design1=design & right1=right
  for i=0,m-1 do begin 
    design1(*,i)=design(*,i)*weight(i)
    right1(i)=right(i)*weight(i)
  endfor
;
    svd,transpose(design1),w,u,v
    wp=fltarr(n,n) & for i=0,n-1 do wp(i,i)=1./w(i)
    x=v#wp#(transpose(u)#right1)
;
  var=fltarr(n)
  for j=0,n-1 do for i=0,n-1 do $
     var(j)=var(j)+v(j,i)^2*wp(i,i)
;
  covar = fltarr(n,n)
  for j=0,n-1 do for l=0,j do begin
      s = 0.
      for k=0,n-1 do s = s + wp(k,k) * v(j,k) * v(k,l)
      covar(j,l) = s
      covar(l,j) = s
  endfor
;
return,x
end
;================================================================
function SXPAR, hdr, name, abort, COUNT=matches, COMMENT = comments, $
                                  NoContinue = NoContinue, SILENT = silent
;+
; NAME:
;      SXPAR
; PURPOSE:
;      Obtain the value of a parameter in a FITS header
;
; CALLING SEQUENCE:
;      result = SXPAR( Hdr, Name, [ Abort, COUNT=, COMMENT =, /NoCONTINUE  ])   
;
; INPUTS:
;      Hdr =  FITS header array, (e.g. as returned by READFITS) 
;             string array, each element should have a length of 80 characters      
;
;      Name = String name of the parameter to return.   If Name is of the
;             form 'keyword*' then an array is returned containing values of
;             keywordN where N is an integer.  The value of keywordN will be
;             placed in RESULT(N-1).  The data type of RESULT will be the
;             type of the first valid match of keywordN found.
;
; OPTIONAL INPUTS:
;       ABORT - string specifying that SXPAR should do a RETALL
;               if a parameter is not found.  ABORT should contain
;               a string to be printed if the keyword parameter is not found.
;               If not supplied, SXPAR will return quietly with COUNT = 0
;               (and !ERR = -1) if a keyword is not found.
;
; OPTIONAL INPUT KEYWORDS: 
;       /NOCONTINUE = If set, then continuation lines will not be read, even
;                 if present in the header
;       /SILENT - Set this keyword to suppress warning messages about duplicate
;                 keywords in the FITS header.
;
; OPTIONAL OUTPUT KEYWORDS:
;       COUNT - Optional keyword to return a value equal to the number of 
;               parameters found by SXPAR, integer scalar
;
;       COMMENT - Array of comments associated with the returned values
;
; OUTPUTS:
;       Function value = value of parameter in header.
;               If parameter is double precision, floating, long or string,
;               the result is of that type.  Apostrophes are stripped
;               from strings.  If the parameter is logical, 1b is
;               returned for T, and 0b is returned for F.
;               If Name was of form 'keyword*' then a vector of values
;               are returned.
;
; SIDE EFFECTS:
;       !ERR is set to -1 if parameter not found, 0 for a scalar
;       value returned.  If a vector is returned it is set to the
;       number of keyword matches found.    The use of !ERR is deprecated, and
;       instead the COUNT keyword is preferred
;
;       If a keyword (except HISTORY or COMMENT) occurs more than once in a 
;       header, a warning is given, and the *last* occurence is used.
;
; EXAMPLES:
;       Given a FITS header, h, return the values of all the NAXISi values
;       into a vector.    Then place the history records into a string vector.
;
;       IDL> naxisi = sxpar( h ,'NAXIS*')         ; Extract NAXISi value
;       IDL> history = sxpar( h, 'HISTORY' )      ; Extract HISTORY records
;
; PROCEDURE:
;       The first 8 chacters of each element of Hdr are searched for a 
;       match to Name.  The value from the last 20 characters is returned.  
;       An error occurs if there is no parameter with the given name.
;
;       If a numeric value has no decimal point it is returned as type
;       LONG.   If it contains more than 8 numerals, or contains the 
;       characters 'D' or 'E', then it is returned as type DOUBLE.  Otherwise
;       it is returned as type FLOAT.    Very large integer values, outside
;       the range of valid LONG, are returned as DOUBLE.
;
;       If the value is too long for one line, it may be continued on to the
;       the next input card, using the OGIP CONTINUE convention.  For more info,
;       http://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/ofwg_recomm/r13.html
;
;       Complex numbers are recognized as two numbers separated by one or more
;       space characters.
;
;       If a numeric value has no decimal point (or E or D) it is returned as
;       type LONG.  If it contains more than 8 numerals, or contains the
;       character 'D', then it is returned as type DOUBLE.  Otherwise it is
;       returned as type FLOAT.    If an integer is too large to be stored as
;       type LONG, then it is returned as DOUBLE.
;
; NOTES:
;       The functions SXPAR() and FXPAR() are nearly identical, although
;       FXPAR() has slightly more sophisticated parsing.   There is no
;       particular reason for having two nearly identical procedures, but
;       both are too widely used to drop either one.
;
; PROCEDURES CALLED:
;       GETTOK(), VALID_NUM()
; MODIFICATION HISTORY:
;       DMS, May, 1983, STPAR Written.
;       D. Lindler Jan 90 added ABORT input parameter
;       J. Isensee Jul,90 added COUNT keyword
;       W. Thompson, Feb. 1992, added support for FITS complex values.
;       W. Thompson, May 1992, corrected problem with HISTORY/COMMENT/blank
;               keywords, and complex value error correction.
;       W. Landsman, November 1994, fix case where NAME is an empty string 
;       W. Landsman, March 1995,  Added COMMENT keyword, ability to read
;               values longer than 20 character
;       W. Landsman, July 1995, Removed /NOZERO from MAKE_ARRAY call
;       T. Beck May 1998, Return logical as type BYTE
;       W. Landsman May 1998, Make sure integer values are within range of LONG
;       Converted to IDL V5.0, May 1998
;       W. Landsman Feb 1998, Recognize CONTINUE convention 
;       W. Landsman Oct 1999, Recognize numbers such as 1E-10 as floating point
;       W. Landsman Jan 2000, Only accept integer N values when name = keywordN
;       W. Landsman Dec 2001, Optional /SILENT keyword to suppress warnings
;       W. Landsman/D. Finkbeiner  Mar 2002  Make sure extracted vectors 
;             of mixed data type are returned with the highest type.
;-
;----------------------------------------------------------------------
 if N_params() LT 2 then begin
     print,'Syntax -     result =  sxpar( hdr, name, [abort])'
     print,'   Input Keywords:    /NOCONTINUE, /SILENT'
     print,'   Output Keywords:   COUNT=,  COMMENT= '
     return, -1
 endif 

 VALUE = 0
 if N_params() LE 2 then begin
      abort_return = 0
      abort = 'FITS Header'
 end else abort_return = 1
 if abort_return then On_error,1 else On_error,2

;       Check for valid header

  s = size(hdr)         ;Check header for proper attributes.
  if ( s[0] NE 1 ) or ( s[2] NE 7 ) then $
           message,'FITS Header (first parameter) must be a string array'

  nam = strtrim( strupcase(name) )      ;Copy name, make upper case     


;  Determine if NAME is of form 'keyword*'.  If so, then strip off the '*', and
;  set the VECTOR flag.  One must consider the possibility that NAM is an empty
;  string.

   namelength1 = (strlen(nam) - 1 ) > 1         
   if strpos( nam, '*' ) EQ namelength1 then begin    
            nam = strmid( nam, 0, namelength1)  
            vector = 1                  ;Flag for vector output  
            name_length = strlen(nam)   ;Length of name 
            num_length = 8 - name_length        ;Max length of number portion  
            if num_length LE 0 then  $ 
                  message, 'Keyword length must be 8 characters or less'

;  Otherwise, extend NAME with blanks to eight characters.

    endif else begin  
                while strlen(nam) LT 8 do nam = nam + ' ' ;Make 8 chars long
                vector = 0      
    endelse


;  If of the form 'keyword*', then find all instances of 'keyword' followed by
;  a number.  Store the positions of the located keywords in NFOUND, and the
;  value of the number field in NUMBER.

        histnam = (nam eq 'HISTORY ') or (nam eq 'COMMENT ') or (nam eq '') 
        if N_elements(start) EQ 0 then start = -1l
        start = long(start[0])
        if (not vector) and (start GE 0) then begin
            if N_elements(precheck)  EQ 0 then precheck = 5
            if N_elements(postcheck) EQ 0 then postcheck = 20
            nheader = N_elements(hdr)
            mn = (start - precheck)  > 0
            mx = (start + postcheck) < nheader-1
            keyword = strmid(hdr[mn:mx], 0, 8)
        endif else begin
            restart:
            start   = -1l
            keyword = strmid( hdr, 0, 8)
        endelse

        if vector then begin
            nfound = where(strpos(keyword,nam) GE 0, matches)
            if ( matches gt 0 ) then begin
                numst= strmid( hdr[nfound], name_length, num_length)
                number = replicate(-1, matches)
                for i = 0, matches-1 do         $
                    if VALID_NUM( numst[i], num,/INTEGER) then number[i] = num
                igood = where(number GE 0, matches)
                if matches GT 0 then begin
                    nfound = nfound[igood]
                    number = number[igood]
                endif
            endif

;  Otherwise, find all the instances of the requested keyword.  If more than
;  one is found, and NAME is not one of the special cases, then print an error
;  message.

        endif else begin
            nfound = where(keyword EQ nam, matches)
            if (matches EQ 0) and (start GE 0) then goto, RESTART
            if (start GE 0) then nfound = nfound + mn
            if (matches GT 1) and (not histnam) then        $
                if not keyword_set(silent) then $
                message,/informational, 'Warning - keyword ' +   $
                nam + ' located more than once in ' + abort
            if (matches GT 0) then start = nfound[matches-1]
        endelse


; Process string parameter 

 if matches GT 0 then begin
  line = hdr[nfound]
  svalue = strtrim( strmid(line,9,71),2)
  if histnam then $
        value = strtrim(strmid(line,8,71),2) else for i = 0,matches-1 do begin
      if ( strmid(svalue[i],0,1) EQ "'" ) then begin   ;Is it a string?
                  test = strmid( svalue[i],1,strlen( svalue[i] )-1)
                  next_char = 0
                  off = 0
                  value = '' 
          NEXT_APOST:
                  endap = strpos(test, "'", next_char)      ;Ending apostrophe  
                  if endap LT 0 then $ 
                            MESSAGE,'Value of '+name+' invalid in '+abort
                  value = value + strmid( test, next_char, endap-next_char )  

;  Test to see if the next character is also an apostrophe.  If so, then the
;  string isn't completed yet.  Apostrophes in the text string are signalled as
;  two apostrophes in a row.

                 if strmid( test, endap+1, 1) EQ "'" then begin    
                    value = value + "'"
                    next_char = endap+2         
                    goto, NEXT_APOST
                 endif      

; Extract the comment, if any
                
                slash = strpos( test, "/", endap )
                if slash LT 0 then comment = '' else    $
                        comment = strmid( test, slash+1, strlen(test)-slash-1 )

; This is a string that could be continued on the next line.  Check this
; possibility with the following four criteria: *1) Ends with '&'
; (2) Next line is CONTINUE  (3) LONGSTRN keyword is present (recursive call to
; SXPAR) 4. /NOCONTINE is not set

    if not keyword_set(nocontinue) then begin
                off = off + 1
                val = strtrim(value,2)

                if (strlen(val) gt 0) and $
                  (strmid(val, strlen(val)-1, 1) EQ '&') and $
                  (strmid(hdr[nfound[i]+off],0,8) EQ 'CONTINUE') then begin
                   if (size(sxpar(hdr, 'LONGSTRN',/NoCONTINUE)))[1] EQ 7 then begin                    
                  value = strmid(val, 0, strlen(val)-1)
                  test = hdr[nfound[i]+off]
                  test = strmid(test, 8, strlen(test)-8)
                  test = strtrim(test, 2)
                  if strmid(test, 0, 1) NE "'" then message, $
                    'ERROR: Invalidly CONTINUEd string in '+ abort
                  next_char = 1
                  GOTO, NEXT_APOST
                ENDIF
               ENDIF
    ENDIF


; Process non-string value  

          endif else begin

                test = svalue[i]
                slash = strpos( test, "/" )
                if slash GT 0 then begin
                        comment = strmid( test, slash+1, strlen(test)-slash-1 )
                        test = strmid( test, 0, slash )
                end else comment = ''

; Find the first word in TEST.  Is it a logical value ('T' or 'F')

                test2 = test
                value = gettok(test2,' ')
               if ( value EQ 'T' ) then value = 1b else $
               if ( value EQ 'F' ) then value = 0b else begin

;  Test to see if a complex number.  It's  a complex number if the value and
;  the next word, if any, are both valid values.

                if strlen(test2) EQ 0 then goto, NOT_COMPLEX
                value2 = gettok( test2, ' ') 
                if value2 EQ '' then goto, NOT_COMPLEX
                On_ioerror, NOT_COMPLEX
                value2 = float(value2)
                value = complex(value,value2)
                goto, GOT_VALUE

;  Not a complex number.  Decide if it is a floating point, double precision,
;  or integer number.

NOT_COMPLEX:
                On_IOerror, GOT_VALUE
                  if (strpos(value,'.') GE 0) or (strpos(value,'E') GT 0) $
                  or (strpos(value,'D') GE 0) then begin  ;Floating or double?
                      if ( strpos(value,'D') GT 0 ) or $  ;Double?
                         ( strlen(value) GE 8 ) then value = double(value) $
                                                else value = float(value)
                       endif else begin                   ;Long integer
                            lmax = 2.0d^31 - 1.0d
                            lmin = -2.0d31
                            value = double(value)
                            if (value GE lmin) and (value LE lmax) then $
                                value = long(value)
                       endelse

GOT_VALUE:
                On_IOerror, NULL
                endelse
             endelse; if c eq apost

;  Add to vector if required

         if vector then begin
               if ( i EQ 0 ) then begin
                     maxnum = max(number)
                     dtype = size(value,/type)
                     result = make_array( maxnum, TYPE = dtype )
                     comments = strarr( maxnum )
               endif 
               if size(value,/type) GT dtype then begin   ;Do we need to recast?
                    result = result + 0*value
                    dtype = size(value,/type)
               endif
               result[ number[i]-1 ] =  value
               comments[ number[i]-1 ] = comment
          endif else $
                comments = comment
  endfor

  if vector then begin
         !ERR = matches     
         return, result
  endif else !ERR = 0

endif  else  begin    
     if abort_return then message,'Keyword '+nam+' not found in '+abort
     !ERR = -1
endelse     

return, value       

END                 
;================================================================
function three_gs2,x,y,a,da
;
; function to fit three Gaussians
; usage: yfit=three_gs2(x,y,a,da)
; input: x,y,a
; output: yfit,da
; da are diff.corr.'s to parameters a
; a(0), a(3) and a(6) - central strengths
; a(1), a(4) and a(7) - central positions
; a(2), a(5) and a(8) - widths
; a(9) is zero point
;
; widths a[2], a[5] are kept constant to determine
;    the set for the 3rd Gaussian best: a[6,7,8]

  n = n_elements(x)
  c = x & e1=x & e2=x & e3=x
  
  c = ((x-a(1))/a(2))^2
  e1[*]=0.
  d = where(c le 30.)
  if (d[0] ne -1) then e1(d) = exp(-c(d))

  c = ((x-a(4))/a(5))^2
  e2[*]=0.
  d = where(c le 30.)
  if (d[0] ne -1) then e2(d) = exp(-c(d))

  c = ((x-a(7))/a(8))^2
  e3[*]=0.
  d = where(c le 30.)
  if (d[0] ne -1) then e3(d) = exp(-c(d))

  res = y - a(0)*e1 - a(3)*e2 - a(6)*e3-a(9)
  print,'Sum res^2 = ',total(res^2)
;

;  t = fltarr(n,10)  
  t = fltarr(n,8) 
        ; 0=str1, 1=pos1, 2=str2, 3=pos2, 4=str3, 5=pos3, 6=wid3, 7=off
  t(*,0) = e1
  t(*,2) = e2
  t(*,4) = e3
  t(*,1) = 2*e1*a(0)*(x-a(1))/a(2)^2
  t(*,3) = 2*e2*a(3)*(x-a(4))/a(5)^2
  t(*,5) = 2*e3*a(6)*(x-a(7))/a(8)^2
;  t(*,2) = t(*,1)*(x-a(1))/a(2)
;  t(*,5) = t(*,4)*(x-a(4))/a(5)
  t(*,6) = t(*,5)*(x-a(7))/a(8)
  t(*,7) = 1.
;
  tt = transpose(t)#t
  tr = transpose(t)#res
  svd,tt,w,u,v
  wp = fltarr(8,8)
  mw = max(w)
  for i=0,7 do if w(i) ge mw*1.e-10 then wp(i,i)=1./w(i)
  dda = v#wp#(transpose(u)#tr)
  da = [dda(0:1),0.,dda(2:3),0.,dda(4:6),dda(7)] ; no changes to widths
  a1 = a + da                                    ; of comp. 1 & 2
;
  c = ((x-a(1))/a(2))^2
  e1[*]=0.
  d = where(c le 30.)
  if (d[0] ne -1) then e1(d) = exp(-c(d))

  c = ((x-a(4))/a(5))^2
  e2[*]=0.
  d = where(c le 30.)
  if (d[0] ne -1) then e2(d) = exp(-c(d))

  c = ((x-a(7))/a(8))^2
  e3[*]=0.
  d = where(c le 30.)
  if (d[0] ne -1) then e3(d) = exp(-c(d))

  res = y - a1(0)*e1 - a1(3)*e2 - a1(6)*e3-a1(9)
  print,'Sum new res^2 = ',total(res^2)
return,a1(0)*e1+a1(3)*e2+a1(6)*e3+a1(9)
end
;================================================================
function trunc,spec,m,n
; must be: m-odd, n-even
   return,spec[m/2:n-m/2-1]
end
;================================================================
function two_gs,x,y,a,da
;
; function to fit two Gaussians
; usage: yfit=two_gs(x,y,a,da)
; input: x,y,a
; output: yfit,da
; da are diff.corr.'s to parameters a
; a(0) and a(3) - central strengths
; a(1) and a(4) - central positions
; a(2) and a(5) - widths
; a(6) is zero point
;
  n = n_elements(x)
  c = x & e1=x & e2=x
  c = ((x-a(1))/a(2))^2
  e1(*)=0. & e2(*)=0.
  d = where(c le 30.)
  e1(d) = exp(-c(d))
  c = ((x-a(4))/a(5))^2
  d = where(c le 30.)
  e2(d) = exp(-c(d))
  res = y - a(0)*e1 - a(3)*e2-a(6)
  print,'Sum res^2 = ',total(res^2)
;
  t = fltarr(n,7)
  t(*,0) = e1                             ; strenght
  t(*,3) = e2
  t(*,1) = 2*e1*a(0)*(x-a(1))/a(2)^2      ; position
  t(*,4) = 2*e2*a(3)*(x-a(4))/a(5)^2
  t(*,2) = t(*,1)*(x-a(1))/a(2)           ; width
  t(*,5) = t(*,4)*(x-a(4))/a(5)
  t(*,6) = 1.
;
  tt = transpose(t)#t
  tr = transpose(t)#res
  svd,tt,w,u,v
  wp = fltarr(7,7)
  mw = max(w)
  for i=0,6 do if w(i) ge mw*1.e-6 then wp(i,i)=1./w(i)
  da = v#wp#(transpose(u)#tr)
  a1 = a + da
;
  c = x & e1=x & e2=x
  e1(*)=0. & e2(*)=0.
  c = ((x-a1(1))/a1(2))^2
  d = where(c le 30.)
  e1(d) = exp(-c(d))
  c = ((x-a1(4))/a1(5))^2
  d = where(c le 30.)
  e2(d) = exp(-c(d))
  res = y - a1(0)*e1 - a1(3)*e2-a1(6)
  print,'Sum new res^2 = ',total(res^2)
return,a1(0)*e1+a1(3)*e2+a1(6)
end
;================================================================
function two_gs1,x,y,a,da
;
; function to fit two Gaussians
; usage: yfit=two_gs1(x,y,a,da)
; input: x,y,a
; output: yfit,da
; da are diff.corr.'s to parameters a
; a(0) and a(3) - central strengths
; a(1) and a(4) - central positions
; a(2) and a(5) - widths                 [a(5) kept constant]
; a(6) is zero point
;
; NOTE: this version keeps the width of the 2nd Gaussian fixed 
; but inputs and outputs are the same as for "two_gs"
;
  n = n_elements(x)
  c = x & e1=x & e2=x
  c = ((x-a(1))/a(2))^2
  e1(*)=0. & e2(*)=0.
  d = where(c le 30.)
  e1(d) = exp(-c(d))
  c = ((x-a(4))/a(5))^2
  d = where(c le 30.)
  e2(d) = exp(-c(d))
  res = y - a(0)*e1 - a(3)*e2-a(6)
  print,'Sum res^2 = ',total(res^2)
;
  t = fltarr(n,6)
  t(*,0) = e1
  t(*,3) = e2
  t(*,1) = 2*e1*a(0)*(x-a(1))/a(2)^2
  t(*,4) = 2*e2*a(3)*(x-a(4))/a(5)^2
  t(*,2) = t(*,1)*(x-a(1))/a(2)   ; to solve for width of 1st Gaussian
  t(*,5) = 1.
;
  tt = transpose(t)#t
  tr = transpose(t)#res
  svd,tt,w,u,v
  wp = fltarr(6,6)
  mw = max(w)
  for i=0,5 do if w(i) ge mw*1.e-6 then wp(i,i)=1./w(i)
  dda = v#wp#(transpose(u)#tr)
  da = [dda(0:4),0.,dda(5)]   ; width of 2nd Gaussian unchanged
  a1 = a + da
;
  c = x & e1=x & e2=x
  e1(*)=0. & e2(*)=0.
  c = ((x-a1(1))/a1(2))^2
  d = where(c le 30.)
  e1(d) = exp(-c(d))
  c = ((x-a1(4))/a1(5))^2
  d = where(c le 30.)
  e2(d) = exp(-c(d))
  res = y - a1(0)*e1 - a1(3)*e2-a1(6)
  print,'Sum new res^2 = ',total(res^2)
return,a1(0)*e1+a1(3)*e2+a1(6)
end
;================================================================
function two_gs2,x,y,a,da
;
; function to fit two Gaussians
; usage: yfit=two_gs2(x,y,a,da)
; input: x,y,a (values of x,y for function and 7 parameters a)
; output: yfit,da (fit and corrections to a)
; a(0) and a(3) - central strengths
; a(1) and a(4) - central positions
; a(2) and a(5) - widths       [both kept constant]
; a(6) is zero point
;
; NOTE: this version keeps widths of both Gaussians fixed 
;    as full solutions normally unstable

  n = n_elements(x)
  c = x & e1=x & e2=x
  e1(*)=0. & e2(*)=0.
  c = ((x-a(1))/a(2))^2
  d = where(c le 30.)
  e1(d) = exp(-c(d))
  c = ((x-a(4))/a(5))^2
  d = where(c le 30.)
  e2(d) = exp(-c(d))
  res = y - a(0)*e1 - a(3)*e2-a(6)
  print,'Sum res^2 = ',total(res^2)
;
  t = fltarr(n,5)  ; widths not determined
  t(*,0) = e1
  t(*,2) = e2
  t(*,1) = 2*e1*a(0)*(x-a(1))/a(2)^2
  t(*,3) = 2*e2*a(3)*(x-a(4))/a(5)^2
  t(*,4) = 1.
;
  tt = transpose(t)#t   ; solution
  tr = transpose(t)#res
  svd,tt,w,u,v
  wp = fltarr(5,5) 
  mw = max(w)
  for i=0,4 do if w(i) ge mw*1.e-6 then wp(i,i)=1./w(i)
  dda = v#wp#(transpose(u)#tr)          ; corrections computed
  da = [dda(0:1),0.,dda(2:3),0.,dda(4)] ; no changes to widths
  a1 = a + da
;
  c = x & e1=x & e2=x  ; initialize
  e1(*)=0. & e2(*)=0.
  c = ((x-a1(1))/a1(2))^2
  d = where(c le 30.)
  e1(d) = exp(-c(d))   ; first Gaussian
  c = ((x-a1(4))/a1(5))^2
  d = where(c le 30.)
  e2(d) = exp(-c(d))   ; second Gaussian
  res = y - a1(0)*e1 - a1(3)*e2-a1(6)  ; residuals
  print,'Sum new res^2 = ',total(res^2)
return,a1(0)*e1+a1(3)*e2+a1(6)
end
;================================================================
pro two_star1,star,tp,res
; sine fit to two stars
; 
; Input: star - data array [5,*]: JD, vel1, weigh1, vel2, weigh2 (in cols)
;        tp   - 2-el vector with T0,P
; Output: res - 8-el vector:
;            gam, K1, K2 - 3 parameters + dummy 
;            and 3-el vector of errors + dummy
; Usage: two_star1,ah,ah0,ah1
;

  star = double(star)
  tp=double(tp)

  ph = (star[0,*]-tp[0])/tp[1]
  ph = reform(ph)
  ph = ph - floor(ph)
  n = n_elements(ph)

  a = dblarr(3,2*n)
  a[0,*] = 1.
  a[1,0:n-1]   = -sin(2*!pi*ph)  ; sense as in A-type
  a[2,n:2*n-1] = +sin(2*!pi*ph)
  b = [reform(star[1,*]),reform(star[3,*])]  ; velocities
  w = [reform(star[2,*]),reform(star[4,*])]  ; weights
  v = 0. & cv = 0. ; variance and covariance
  res = solv(a,b,w,v,cv)  ; SVD solution
;  gam = res[0] & K1 = res[1] & K2 = res[2]
  res=[res,0,sqrt(v),0]  ; zero for place holding

return
end
;================================================================
pro two_star2,star,tp,res,corr
; corrections to the sine fit for two stars
; 
; Input: star - data array [5,*]
;        tp   - 2-el vector with T0,P
;        res - 8-el vector of the sine fit
;    consisting of: gam, K1, K2, 0 - 4 parameters 
;                  and 4-el vector of errors
; Output: corr - 4 el vector of diff.corrs 
;          delta of [gam,K1,K2,T0]
; Usage: two_star2,ah,ah0,ah1,ah2
;

  star=double(star)
  tp=double(tp)

  ph = (star[0,*]-tp[0])/tp[1]
  ph = reform(ph)
  ph = ph - floor(ph)
  n = n_elements(ph)

  dif1 = reform(star[1,*])-reform(res[0]-res[1]*sin(2*!pi*ph))
  dif2 = reform(star[3,*])-reform(res[0]+res[2]*sin(2*!pi*ph))
  w = [reform(star[2,*]),reform(star[4,*])]  ; weights
  b = [dif1,dif2]       ; velocity diffs O-C

  a = dblarr(4,2*n)
  a[0,*] = 1.
  a[1,0:n-1]   = -sin(2*!pi*ph)  ; sense as in A-type
  a[2,n:2*n-1] = +sin(2*!pi*ph)
  a[3,0:n-1]   = +res[1]*2*!pi/tp[1]*cos(2*!pi*ph) ; corr to T0
  a[3,n:2*n-1] = -res[2]*2*!pi/tp[1]*cos(2*!pi*ph)

  v = 0. & cv = 0. ; variance and covariance
  corr = solv(a,b,w,v,cv)  ; SVD solution of corrections
;  gam = corr[0] & K1 = corr[1] & K2 = corr[2]
  corr=[corr,sqrt(v)]

  dif1 = reform(star[1,*])-  $
      reform((res[0]+corr[0])-(res[1]+corr[1])*sin(2*!pi*ph))
  dif2 = reform(star[3,*])-  $
      reform((res[0]+corr[0])+(res[2]+corr[2])*sin(2*!pi*ph))
  w1 = reform(star[2,*])
  w2 = reform(star[4,*])
  print,'std dev one obs'
  print,sqrt(total(dif1^2*w1)/total(w1))
  print,sqrt(total(dif2^2*w2)/total(w2))

return
end
;================================================================
pro two_star3,star,tp,res,res1
; error ranges for the sine fit for two stars by bootstrap
; 
; Input: star - data array [5,*]
;        tp   - 2-el vector with T0,P
;        res - 8-el vector of the sine fit
;    consisting of: gam, K1, K2, 0 - 4 parameters 
;                  and 4-el vector of errors
; Output: res1 - 4 +/- 1-sigma ranges and medians of 
;        diff.corr's [gam,K1,K2,T0]
; Usage: two_star3,ah,ah0,ah1,ah3
;

  star=double(star)
  tp=double(tp)

  ph = (star[0,*]-tp[0])/tp[1]
  ph = reform(ph)
  ph = ph - floor(ph)
  n = n_elements(ph)

  dif1 = reform(star[1,*])-reform(res[0]-res[1]*sin(2*!pi*ph))
  dif2 = reform(star[3,*])-reform(res[0]+res[2]*sin(2*!pi*ph))
  w = [reform(star[2,*]),reform(star[4,*])]  ; weights
  b = [dif1,dif2]       ; velocity diffs O-C
  a = dblarr(4,2*n)
  a[0,*] = 1.
  a[1,0:n-1]   = -sin(2*!pi*ph)  ; sense as in A-type
  a[2,n:2*n-1] = +sin(2*!pi*ph)
  a[3,0:n-1]   = +res[1]*2*!pi/tp[1]*cos(2*!pi*ph) ; corr to T0
  a[3,n:2*n-1] = -res[2]*2*!pi/tp[1]*cos(2*!pi*ph)

  v = 0. & cv = 0. ; variance and covariance

  a1 = a & b1 = b & m = 2*n
  w1 = replicate(1.,m)  ; uniform weights from now on
  for i = 0,m-1 do begin
    a1[*,i] = a[*,i]*w[i]
    b1[i] = b[i]*w[i]
  end

  k = 1000                 ; number of re-samplings, can be changed
  a2 = a1 & b2 = b1
  rang = dblarr(4,k)
  for i = 0,k-1 do begin
      for j = 0,m-1 do begin 
         l = fix(m*randomu(seed))
         a2[*,j] = a1[*,l]
         b2[j] = b1[l]
      endfor     
  rang[*,i] = solv(a2,b2,w1,v,cv)  ; SVD solution of corr's
; order of corrections: gam, K1, K2, T0
  endfor

  for i = 0,3 do rang[i,*]=rang[i,sort(rang[i,*])]
  res1 = dblarr(4,3)
  res1[*,0] = rang[*,fix(0.158*k)]   ; -1 sigma
  res1[*,1] = rang[*,fix(k/2)]       ; median
  res1[*,2] = rang[*,fix(0.841*k)]   ; +1 sigma

return
end
;=================================================================
; below: two routines for display, not really related to the BF's
;=================================================================
pro set_win,dummy
;
; routine to set plot system variables for Windows
;
  !p.font=-1
  !p.charthick=1.
  !p.charsize=1.
  !x.charsize=1.
  !y.charsize=1.
  !z.charsize=1.
  !p.thick=1.
  !x.thick=1.
  !y.thick=1.
  !z.thick=1.
; 
; just resets, good for any plot configuration
;
  !p.multi=[0,0,0]
  !p.position=[0.15,0.10,0.9,0.92]
  !x.ticks=0
  !y.ticks=0
  !z.ticks=0
  !x.tickname='' ; vectors of strings for ticks, must be (!x.ticks+1)
  !y.tickname='' ;                                   of them
  !z.tickname=''
  !x.tickv=[0,0] ; same but as values, not strings
  !y.tickv=[0,0]
  !z.tickv=[0,0]
  !x.range=[0,0]
  !y.range=[0,0]
  !z.range=[0,0]
  !x.style=1
  !y.style=1
  !z.style=1

;  !x.margin=[10,3]
;  !y.margin=[4,2]

;  !p.title='!17'
;  !x.title='!17'
;  !y.title='!17'
;  !z.title='!17'
;
  set_plot,'win'
return
end
;========================================================================
pro set_ps,dummy
;
; routine to set plot system variables for PostScript printer
;
  !p.font=-1
  !p.charthick=1.2
  !p.charsize=1.2
  !x.charsize=1.2
  !y.charsize=1.2
  !z.charsize=1.2
  !p.thick=5.
  !x.thick=5.
  !y.thick=5.
  !z.thick=5.
; 
; just resets, good for any plot configuration
;
  !p.multi=[0,0,0]
  !p.position=[0.15,0.15,0.85,0.95]
  !x.ticks=0
  !y.ticks=0
  !z.ticks=0
  !x.tickname='' ; vectors of strings for ticks, must be (!x.ticks+1)
  !y.tickname='' ;                                   of them
  !z.tickname=''
  !x.tickv=[0,0] ; same but as values, not strings
  !y.tickv=[0,0]
  !z.tickv=[0,0]
  !x.range=[0,0]
  !y.range=[0,0]
  !z.range=[0,0]
  !x.style=1
  !y.style=1
  !z.style=1
;  !p.title='!17'
;  !x.title='!17'
;  !y.title='!17'
;  !z.title='!17'
;
;
print,'After plotting enter: '
print,'  device,/close '
print,'  (unless closed by calling routine)'
print,'For landscape use: device,/land '
;print,'   useful: !p.position=[0.1,0.05,0.9,0.95] '
  set_plot,'ps'
return
end



;=====================================================================
; Below routines to show BFs as a 2-D image
;=====================================================================
pro BFimage1,bf,hvc,vel,RVstd,phase,star,bf_2d,ph_2d
; 2-D BF's in hel system
; produces a 2-D image of sorted BF's in phases, as they are,
; in perhaps an unequal distribution in phase
; 
; input: bf,hvc,vel,RVstd,phase,star
;        bf - broadening function, say bf15
;        hvc - heliocentric velocity correction
;        vel - velocity vector
;        RVstd - radial vel of the standard used for BF
;        phase - phase vector
;        name of the star for idl.ps image
;
; output: bf_2d - 2-dimensional image with sorted phases
;        ph_2d - corresponding phases
; 
; also   idl.ps image. NOT IN EQUAL SPACING OF PHASE.
;
; use:
;      BFimage1,bf15,hvc,vel,+38.97,phase,'V395 And',bf_2d,ph_2d
;

  k1=n_elements(bf[*,0])  ; number of phases
  k2=n_elements(bf[0,*])  ; number of velocity bins  

  bf_2d=fltarr(k1,k2)

  for i=0,k1-1 do begin
;    junk=vel-rv3m[i]    ; shifted by measured vels of 3rd
     junk=vel+hvc[0]-hvc[i]-RVstd  ; shifted to hel Vstd
     bf_2d[i,*]=interpol(bf[i,*],vel,junk)
  end
  bf_2d=bf_2d[1:*,*]      ; remove std

  kphase=sort(phase[1:*])
  sphase=phase[1:*]
  sphase=sphase[kphase]

  bf_2d=bf_2d[kphase,*]
  ph_2d=sphase           ; sorted phase
;---------------------------------------------------------
set_ps

device,xsize=6.0,ysize=4.0,/inch,bits=8,/port,xoffset=1.25

pos=[0.17,0.15,0.98,0.93]

xsize=(pos[2]-pos[0])*!d.x_Vsize
ysize=(pos[3]-pos[1])*!d.y_Vsize
xstart=pos[0]*!d.x_Vsize
ystart=pos[1]*!d.y_Vsize

tv,bytscl(bf_2d),xstart,ystart,XSize=xsize,Ysize=ysize
plot,findgen(k1),findgen(k2),/nodata,/noerase,pos=pos, $
   tit='!17 '+star,xtit='!17 Bin number in phase', $
   ytit='!17 Bin number in velocity'

device,/close

return
end
;=================================================================
pro BFimage2,bf_2d,vel,phase,del_ph,star,bf_2x,key
; 2-D BF's in hel system in equal intervals of phase
;   if necessary, filling missing observations (key=0)
;
;   run BFimage1.pro first to create bf_2d
; 
; input: bf_2d,hvc,vel,phase,del_ph,star,key,bf_2x
;        bf_2d - broadening function array created using
;                                      BFimage1.pro
;        vel - velocity vector
;        phase - phase vector
;        del_ph - increment (bin size) in phase
;        name of the star for idl.ps image
;        if key=1, then skip filling the missing bins 
;                  currently key=0 does not work
;
; output: bf_2x - 2-dim image with sorted equi-distant phases
;         idl.ps image

; use:
;   BFimage2,bf_2d,vel,phase,0.02,'V395 And',bf_2x,1

  k1=n_elements(bf_2d[*,0])  ; number of phases
  k2=n_elements(bf_2d[0,*])  ; number of velocity bins  
  k3=fix((1.0+1.e-5)/del_ph) ; number of phase intervals
  if abs(k3*del_ph - 1.0) ge 1.e-5 then begin
      print,'Aborting: Select del_ph for integer division of 1.0'
      goto,EE
  end

  kphase=sort(phase[1:*])  ; template not counted
  sphase=phase[1:*]
  sphase=sphase[kphase]

  bf_2x=fltarr(k3,k2)

  ph_bin=del_ph/2.+del_ph*findgen(k3)
  ph_flag=intarr(k3)

for i=0,k3-1 do begin
   k=where((sphase gt ph_bin[i]-del_ph/2.) and $
                    (sphase le ph_bin[i]+del_ph/2))
   if k[0] eq -1 then begin      ; no data i the phase interval
;       bf_2x[i,*]=1. 
       bf_2x[i,*]=max(bf_2d)/2.  ; reasonable grey shade
       print,i,ph_bin[i]
       ph_flag[i]=1 
   end else begin
       n = n_elements(k)
       print,i,ph_bin[i],k
       for j=0,n-1 do bf_2x[i,*]=bf_2x[i,*]+bf_2d[k[j],*]
       bf_2x[i,*]=bf_2x[i,*]/n   ; simple averaging in the bin
   end
end

; ---------------------------------------------------------------

print,'Total bins: ',k3
k=where(ph_flag eq 1)
if k[0] eq -1 then goto,E

print,'Missing: ',k
n = n_elements(k)

; skip interpolation 
;                     - for some reason interpolation does not work
if key eq 1 then goto,E

i1=0 & i2=n
if k[0] eq 0 then i1=1
if k[n-1] eq k3-1 then i2=k3-1

;for i=0,n-1 do bf_2x[k[i],*]=(bf_2x[k[i]-1,*]+bf_2x[k[i]+1,*])/2.
for i=i1,i2-1 do bf_2x[k[i],*]=(bf_2x[k[i]-1,*]+bf_2x[k[i]+1,*])/2.

E:
; ---------------------------------------------------------------

  k1=n_elements(bf_2x[*,0])  ; number of phases
  k2=n_elements(bf_2x[0,*])  ; number of velocity bins  

set_ps

device,xsize=6.0,ysize=4.0,/inch,bits=8,/port,xoffset=1.25

pos=[0.17,0.15,0.98,0.93]

xsize=(pos[2]-pos[0])*!d.x_Vsize
ysize=(pos[3]-pos[1])*!d.y_Vsize
xstart=pos[0]*!d.x_Vsize
ystart=pos[1]*!d.y_Vsize

tv,bytscl(bf_2x),xstart,ystart,XSize=xsize,Ysize=ysize
;plot,findgen(k1)/k1,findgen(k2),/nodata,/noerase,pos=pos, $
plot,findgen(k1)/k1,vel,/nodata,/noerase,pos=pos, $
   tit='!17'+star+' (bin ='+string(del_ph,f='(f6.3)')+')', $
   xstyle=1,ystyle=5, $
   xtit='!17 Phase'

axis,yaxis=0,ystyle=1,yrange=[vel[0],vel[k2-1]], $
   ytit='!17 Velocity (km/s)'
axis,yaxis=1,ystyle=1,yrange=[vel[0],vel[k2-1]], $
   ychars=0.001

device,/close

EE:
return
end
;=================================================================
pro BFimage3,bf_2d,vel,phase,fwhm,ph_y,star,bf_2y,no_y
;
; 2-D BF's in hel system in equal intervals of phase
; smoothing along phase direction and new BFs given for
; phases set in vector ph_y
;
; input: bf_2d,vel,phase,fwhm,ph_y,star
;
;        bf_2d - broadening function array (helio) created using
;                   BFimage1.pro, not in equal spacing in phase
;        vel - velocity vector
;        phase - phase vector
;        fwhm - width of the smoothing Gaussian
;        ph_y - vector of phases to interpolate into
;            eg. 0.01 + 0.02*findgen(50)   [0.01 (0.02) 0.99]
;                0.025*findgen(40)         [0.00 (0.025) 0.975]
;                0.05*findgen(20)          [0.00 (0.05) 0.95]
;            or can be just a vector of required phases
;        star - name of the star for idl.ps image
; output:
;        bf_2y - smooth BFs for phases in ph_y
;        no_y - number of original BF data per FWHM
;
; use:
; BFimage3,bf_2d,vel,phase,0.025,0.025*findgen(40),'V395 And',bf_2y,no_y
; or
; BFimage3,bf_2d,vel,phase,0.025,ph_y,'V395 And',bf_2y,no_y
;
; Routine modified by Theo Pribulla, February 12, 2009

    k1=n_elements(bf_2d[*,0])  ; number of phases, one less than phases
    k2=n_elements(bf_2d[0,*])  ; number of velocity bins
    k3=n_elements(ph_y)        ; number of new phases
    no_y=fltarr(k3)            ; number of obs per FWHM element
    sphase=fltarr(3*k1)
    ans=''

    kphase=sort(phase[1:*]) & ssphase=phase[1:*] &  ssphase=ssphase[kphase]
    sphase[0:k1-1]=ssphase-1.0
    sphase[k1:2*k1-1]=ssphase
    sphase[2*k1:3*k1-1]=ssphase+1.0

    ; print,sphase

    bf_2y=fltarr(k3,k2)               ; new smoothed image
    cut=fltarr(3*k1)
    x=(findgen(3001)-1000)/1000.      ; for smoothing: [-1,2]

    sig=fwhm/2.354e-3                 ; sigma in pixels for 0.001 step in x
    gs=gs_smooth(3001,sig)            ; Gaussian kernel

for i=0,k2-1 do begin
   cut[0:k1-1]=reform(bf_2d[*,i])     ; horizontal cut for a given velocity
   cut[k1:2*k1-1]=reform(bf_2d[*,i])
   cut[2*k1:3*k1-1]=reform(bf_2d[*,i])
   cut1=interpol(cut,sphase,x)        ; rebin to 3001 phases & smoothing
   cut2=cnv(cut1,gs)
   bf_2y[*,i]=interpol(cut2,x,ph_y)
endfor

for i=0,k3-1 do begin
   w=where((sphase gt ph_y[i]-fwhm/2.) and (sphase le ph_y[i]+fwhm/2.))
   if w[0] ne -1 then no_y[i]=n_elements(w)
endfor

for i=0,k3-1 do print,f='(f7.3,f7.1)',ph_y[i],no_y[i]
; ---------------------------------------------------------------
set_ps

device,xsize=6.0,ysize=4.0,/inch,bits=8,/port,xoffset=1.25
pos=[0.17,0.15,0.98,0.93]

xsize=(pos[2]-pos[0])*!d.x_Vsize
ysize=(pos[3]-pos[1])*!d.y_Vsize
xstart=pos[0]*!d.x_Vsize
ystart=pos[1]*!d.y_Vsize

tv,bytscl(bf_2y),xstart,ystart,XSize=xsize,Ysize=ysize
plot,ph_y,vel,/nodata,/noerase,pos=pos, $
   tit='!17'+star+' (FWHM ='+string(fwhm,f='(f6.3)')+')', $
   xstyle=1,ystyle=5, $
   xtit='!17 Phase'

axis,yaxis=0,ystyle=1,yrange=[vel[0],vel[k2-1]], $
   ytit='!17 Velocity (km/s)'
axis,yaxis=1,ystyle=1,yrange=[vel[0],vel[k2-1]], $
   ychars=0.001

device,/close

return
end
;===================================================================
pro BFimage4,bf_2x,vel,star
; contour plot of the final BF image
; uses the equally-rebinned image
;
; usage:
;       BFimage4,bf_2x,vel,'V395 And'

k1=n_elements(bf_2x[*,0])

set_ps
device,bits=8,color=0,/port

contour,1.-bf_2x,findgen(k1)/k1,vel,        $
    xran=[0,1],                    $
    yran=[-450,+450],              $
    xstyl=1,ystyl=1,/cell_fill,    $
;    lev=[0,0.05,0.5,0.8, 0.94,0.96,0.98],  $
    lev=[0.01,0.05,0.1,0.4,0.6,0.7,0.75,0.78,0.8,0.82,0.85,0.87,0.9],  $
    tit='!17 '+star,  $
    xtit='!17 Phase',ytit='!17 Velocity (km/s)'

;x=findgen(100)/100.
;oplot,x,-160.*sin(2*!pi*x)+30,line=2,color=255
;oplot,x,+130.*sin(2*!pi*x)+30,line=2,color=255

device,/close

return
end
;====================================================================
pro BFimage5,bf_2y,vel,fwhm,ph_y,star,bf_2z
;
; utilizes existing BFimage bf_2y rebinned to equal intervals in
; phase in helio system; a version of BFimage3.pro
; 
; input: bf_2y,vel,fwhm,ph_y,star
; 
;        bf_2y - broadening function array (helio) created using
;                                      BFimage2.pro
;        vel - velocity vector
;        fwhm - width of the smoothing Gaussian
;        ph_y - vector of phases
;            eg. 0.01 + 0.02*findgen(50)   [0.01 (0.02) 0.99]
;                0.025*findgen(40)         [0.00 (0.025) 0.975]
;                0.05*findgen(20)          [0.00 (0.05) 0.95]
;            or can be just a vector of required phases
;        star - name of the star for idl.ps image
; output:
;        bf_2z - smooth BFs for phases in ph_y
;
; use:
; BFimage5,bf_2y,vel,0.025,0.025*findgen(40),'V395 And',bf_2z

    k1=n_elements(bf_2y[*,0])  ; number of phases, one less than phases
    k2=n_elements(bf_2y[0,*])  ; number of velocity bins  

;    kphase=sort(phase[1:*]) & sphase=phase[1:*] &  sphase=sphase[kphase]

    bf_2z=bf_2y                ; new smoothed image
    x=findgen(1001)/1000.      ; for smoothing: [0,1]

    sig=fwhm/2.354e-3          ; sigma in pixels for 0.001 step in x
    gs=gs_smooth(1001,sig)     ; Gaussian kernel

for i=0,k2-1 do begin
   cut=reform(bf_2y[*,i])     ; horizontal cut for a given velocity
   help,cut,ph_y
   cut1=interpol(cut,ph_y,x)   ; rebin to 1001 & smoothing
   cut2=cnv(cut1,gs)
   bf_2z[*,i]=interpol(cut2,x,ph_y)
endfor

; ---------------------------------------------------------------
set_ps

device,xsize=6.0,ysize=4.0,/inch,bits=8,/port,xoffset=1.25
pos=[0.17,0.15,0.98,0.93]

xsize=(pos[2]-pos[0])*!d.x_Vsize
ysize=(pos[3]-pos[1])*!d.y_Vsize
xstart=pos[0]*!d.x_Vsize
ystart=pos[1]*!d.y_Vsize

tv,bytscl(bf_2z),xstart,ystart,XSize=xsize,Ysize=ysize
plot,ph_y,vel,/nodata,/noerase,pos=pos, $
   tit='!17'+star+' (FWHM ='+string(fwhm,f='(f6.3)')+')', $
   xstyle=1,ystyle=5, $
   xtit='!17 Phase'

axis,yaxis=0,ystyle=1,yrange=[vel[0],vel[k2-1]], $
   ytit='!17 Velocity (km/s)'
axis,yaxis=1,ystyle=1,yrange=[vel[0],vel[k2-1]], $
   ychars=0.001

device,/close

return
end
;================================================================
;+
; NAME: 
;     VALID_NUM
; PURPOSE:               
;     Check if a string is a valid number representation.
; EXPLANATION:              
;     The input string is parsed for characters that may possibly
;     form a valid number.  It is more robust than simply checking
;     for an IDL conversion error because that allows strings such
;     as '22.3qwert' to be returned as the valid number 22.3
;     See also the original NUM_CHK which returns the status in 
;     the opposite sense.
;
; CALLING SEQUENCE: 
;     IDL> status = valid_num(string  [,value]  [,/integer])
;    
; Inputs      : string  -  the string to be tested
;               
; Opt. Inputs : None
;               
; Outputs     : The function returns 1 for valid, 0 for invalid number
;               
; Opt. Outputs: value	- The value the string decodes to.  This will be
;			  returned as a double precision number unless /INTEGER
;			  is present, in which case a long integer is returned.
;               
; Keywords    : Integer   -  if present code checks specifically for an integer.
;
; Calls       : None
;               
; Restrictions: None
;               
; Category    : Utilities, Numerical
;               
; Prev. Hist. : Small changes from NUM_CHK by Andrew Bowen, 
;                                             Tessella Support Services, 8/3/93
;
; Written     : CDS version by C D Pike, RAL, 24-May-93
;               
; Modified    : Version 1, C D Pike, RAL, 24-May-93
;		Version 2, William Thompson, GSFC, 14 October 1994
;			Added optional output parameter VALUE to allow
;			VALID_NUM to replace STRNUMBER in FITS routines.
;
; Version     : Version 1  24-May-93
;	Converted to IDL V5.0   W. Landsman   September 1997
;-            

FUNCTION valid_num, string, value, INTEGER=integer

		;**** Set defaults for keyword ****
  IF NOT (KEYWORD_SET(integer)) THEN integer=0

		;**** arrays of legal characters ****
  numbers 	= '0123456789'
  signs 	= '+-'
  decimal 	= '.'
  exponents 	= 'ED'

		;**** trim leading and trailing blanks/compress white ****
		;**** space and convert any exponents to uppercase.   ****
  numstr = strupcase(strtrim(strcompress(string),2))

		;**** length of input string ****
  len = strlen(numstr)

  ok = 1

  if integer eq 0 then stage = 1 else stage = 6

  for i = 0, len-1 do begin

    char = strmid(numstr,i,1)

		;**** the parsing steps 1 to 8 are for floating   ****
		;**** point, steps 6 to 8, which test for a legal ****
		;**** exponent, can be used to check for integers ****

;**** The parsing structure is as follows.  Each character in the ****
;**** string is checked against the valid list at the current     ****
;**** stage.  If no match is found an error is reported.  When a  ****
;**** match is found the stage number is updated as indicated     ****
;**** ready for the next character.  The valid end points are     ****
;**** indicated in the diagram.					  ****
;
;Stage	1		2		3		4
;
;Valid	sign	--> 2	dec-pt	--> 3	digit	--> 5	dec-pt	--> 5
;  "	dec-pt	--> 3	digit	--> 4			digit	--> 4
;  "	digit	--> 4					exp't	--> 6
;  "							END
;
;Stage	5		6		7		8
;
;Valid	digit	--> 5	sign	--> 7	digit	--> 8	digit	-->8
;  "	exp't	--> 6	digit	--> 8			END
;  "	END
;

    CASE stage OF

      1 : begin
        if 		strpos(signs,char) ge 0 	then stage = 2 $
	else if 	decimal eq char 		then stage = 3 $
	else if 	strpos(numbers,char) ge 0 	then stage = 4 $
	else 		ok = 0
      end

      2 : begin
	if	 	decimal eq char 		then stage = 3 $
	else if 	strpos(numbers,char) ge 0 	then stage = 4 $
	else 		ok = 0
      end

      3 : begin
	if	 	strpos(numbers,char) ge 0 	then stage = 5 $
	else 		ok = 0
      end

      4 : begin
	if	 	decimal eq char 		then stage = 5 $
	else if 	strpos(numbers,char) ge 0 	then stage = 4 $
	else if		strpos(exponents,char) ge 0	then stage = 6 $
	else 		ok = 0
      end

      5 : begin
	if	 	strpos(numbers,char) ge 0 	then stage = 5 $
	else if		strpos(exponents,char) ge 0	then stage = 6 $
	else 		ok = 0
      end

      6 : begin
        if 		strpos(signs,char) ge 0 	then stage = 7 $
	else if 	strpos(numbers,char) ge 0 	then stage = 8 $
	else 		ok = 0
      end

      7 : begin
	if	 	strpos(numbers,char) ge 0 	then stage = 8 $
	else 		ok = 0
      end

      8 : begin
	if	 	strpos(numbers,char) ge 0 	then stage = 8 $
	else 		ok = 0
      end

    ENDCASE

  end

		;**** check that the string terminated legally ****
		;**** i.e in stages 4, 5 or 8                  ****
  if (stage ne 4) and (stage ne 5) and (stage ne 8) then ok = 0

		;**** If requested, then form the value. ****

  if (n_params() eq 2) and ok then begin
	if keyword_set(integer) then value = long(string) else	$
		value = double(string)
  endif

		;**** return error status to the caller ****
  RETURN, ok


END
;=====================================================================
function gettok,st,char
;+
; NAME:
;	GETTOK                                    
; PURPOSE:
;	Retrieve the first part of the string up to a specified character
; EXPLANATION:
;	GET TOKen - Retrieve first part of string until the character char 
;	is encountered.
;
; CALLING SEQUENCE:
;	token = gettok( st, char )
;
; INPUT:
;	char - character separating tokens, scalar string
;
; INPUT-OUTPUT:
;	st - (scalar) string to get token from (on output token is removed)
;
; OUTPUT:
;	token - scalar string value is returned 
;
; EXAMPLE:
;	If ST is 'abc=999' then gettok(ST,'=') would return
;	'abc' and ST would be left as '999' 
;
; NOTES:
;       A version of GETTOK that accepts vector strings is available for users 
;       of IDL V5.3 or later from  http://idlastro.gsfc.nasa.gov/ftp/v53/
; HISTORY
;	version 1  by D. Lindler APR,86
;	Remove leading blanks    W. Landsman (from JKF)    Aug. 1991
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;----------------------------------------------------------------------
  On_error,2                           ;Return to caller

; if char is a blank treat tabs as blanks

  tab = string(9b)
  while strpos(st,tab) GE 0 do begin    ;Search for tabs
	pos = strpos(st,tab)
	strput,st,' ',pos
  endwhile

  st = strtrim(st,1)              ;Remove leading blanks

; find character in string

  pos = strpos(st,char)
  if pos EQ -1 then begin         ;char not found?
	token = st
 	st = ''
 	return, token
  endif

; extract token

 token = strmid(st,0,pos)
 len = strlen(st)
 if pos EQ (len-1) then st = '' else st = strmid(st,pos+1,len-pos-1)

;  Return the result.

 return,token
 end
;=====================================================================

; updated October 15, 2005
; End of routines
