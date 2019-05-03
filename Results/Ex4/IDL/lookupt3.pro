function lookupT3, x

; generate look-up-table of modelled matrix T3(beta1, aoi,epsilon) with 
; aoi   = angle of incidence, 
; beta1 = angle in a plane perpendicular to scattering plane.
;         by the Unisotropy the scattering signal get's averaged within -beta1 and beta1. -> rotation of T3 matrix and average. 
; the lookuptable shows beta1(H,alpha). H,alpha calculated from T3(beta1,..) (modelled)

; Calculation of T3(beta1,..)
; calculate fresnel coefficients. epsilon = epsilon'-i*epsilon'' (complex dielectric constant)
; assume epsilon'' is zero. (correct at least for ice.. epsilon'' ~ 1e-1..1e-3)

;Neps    = 50
;epsilon = (findgen(Neps)/Neps)^2*(40-2)+2
;beta1   = (findgen(Neps)+1)/Neps*(!pi/2)


epsilon = transpose([x(0,*)])
beta1   = transpose([x(1,*)])
Neps    = n_elements(epsilon)
; print, epsilon
; print, beta1
H0      = scope_varfetch('H0_fetch',Level=-1)
alpha0  = scope_varfetch('alpha0_fetch',Level=-1)
theta_aoi = scope_varfetch('theta_aoi',Level=-1)

Rs  = ( cos(theta_aoi) - sqrt(epsilon - sin(theta_aoi)^2) ) / ( cos(theta_aoi) + sqrt(epsilon - sin(theta_aoi)^2) )  
Rp  = ((epsilon - 1)*(sin(theta_aoi)^2 - epsilon*(1+sin(theta_aoi)^2))) / (epsilon*cos(theta_aoi) + sqrt(epsilon - sin(theta_aoi)^2) )^2
Cf1 = abs(Rs+Rp)^2  
Cf2 = (Rs+Rp)*conj(Rs-Rp)
Cf3 = 0.5*abs(Rs - Rp)^2

; calculate non-zero matrix elements of modelled T3 matrix.
T3m = fltarr(3,3,Neps,Neps)
T3m[0,0,*,*] = transpose(Cf1) ## (fltarr(Neps)+1)
T3m[1,0,*,*] = transpose(Cf2) ## (sin(2*beta1)/(2*beta1))
T3m[0,1,*,*] = T3m[1,0,*,*]
T3m[1,1,*,*] = transpose(Cf3) ## (1+sin(4*beta1)/(4*beta1))
T3m[2,2,*,*] = transpose(Cf3) ## (1-sin(4*beta1)/(4*beta1))

; calculate eigenvalues and vectors of T3m
T3m_eigv  = complexarr([Neps,Neps,3])
T3m_eigvc = complexarr([Neps,Neps,3,3])
for i=0,Neps-1 do begin
    for j=0,Neps-1 do begin
       ; note: the eigenvectors are returned as row-vectors!!
       T3m_eigv[i,j,*]=LA_EIGENQL(reform(T3m[*,*,i,j]),eigenvectors=eigvec)
       T3m_eigvc[i,j,*,*] = eigvec
    endfor
endfor
; calculate Entropy, Anisotropy and alpha for T3m
pT3m = fltarr([Neps,Neps,3])
pT3m[*,*,0] = T3m_eigv[*,*,0]/total(T3m_eigv,3)
pT3m[*,*,1] = T3m_eigv[*,*,1]/total(T3m_eigv,3)
pT3m[*,*,2] = T3m_eigv[*,*,2]/total(T3m_eigv,3) 
HT3m = -total(pT3m*alog(pT3m),3)/alog(3)
;AT3m = real((T3m_eigv[*,*,1]-T3m_eigv[*,*,0])/(T3m_eigv[*,*,1]+T3m_eigv[*,*,0]))
alphaT3m = fltarr([Neps,Neps,3])
alphaT3m[*,*,0] = acos(abs(T3m_eigvc[*,*,0,0]))
alphaT3m[*,*,1] = acos(abs(T3m_eigvc[*,*,0,1]))
alphaT3m[*,*,2] = acos(abs(T3m_eigvc[*,*,0,2]))
alphaT3m_mean = total(alphaT3m*pT3m,3)
; print, 'HT3m = ',  HT3m
; print, 'alpha = ', alphaT3m_mean
return, [HT3m-H0, alphaT3m_mean-alpha0]
end