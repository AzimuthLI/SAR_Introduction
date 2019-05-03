function Xbragg, epsilon, beta1, theta_aoi,anisotropy=anisotropy
; input: vector [epsilon, beta1] of length Neps, Nbeta
; output: 2Neps x Nbeta array:[H(epsilon,beta1,theta), alpha(epsilon,beta1,theta)]
;        or Neps x Nbeta array:A(epsilon,beta1,theta)
 
; generate a model of the matrix T3(beta1, aoi,epsilon) with the parameters 
; aoi   = angle of incidence, 
; beta1 = angle in a plane perpendicular to scattering plane.
;         by the Unisotropy the scattering signal get's averaged within -beta1 and beta1. -> rotation of T3 matrix and average. 
;
; and calculates H and alpha or the Anisotropy calculated from T3(beta1,..) (modelled)

; Calculation of T3(beta1,..)
; calculate fresnel coefficients. epsilon = epsilon'-i*epsilon'' (complex dielectric constant)
; assume epsilon'' is zero. (correct at least for ice.. epsilon'' ~ 1e-1..1e-3)

Neps    = n_elements(epsilon)
Nbeta   = n_elements(beta1)

; generate parameters for matrix elements
Rs  = complexarr(Neps)
Rp  = complexarr(Neps)
Rs  = ( cos(theta_aoi) - sqrt(complex(epsilon - sin(theta_aoi)^2,0)) ) / ( cos(theta_aoi) + sqrt(complex(epsilon - sin(theta_aoi)^2,0)) )
Rp  = ((epsilon - 1)*(sin(theta_aoi)^2 - epsilon*(1+sin(theta_aoi)^2))) / (epsilon*cos(theta_aoi) + sqrt(complex(epsilon - sin(theta_aoi)^2,0)) )^2
Cf1 = abs(Rs+Rp)^2  
Cf2 = (Rs+Rp)*conj(Rs-Rp)
Cf3 = 0.5*abs(Rs - Rp)^2

; calculate non-zero matrix elements of modelled T3 matrix.
T3m = complexarr(3,3,Neps,Nbeta)
T3m[0,0,*,*] = transpose(Cf1) ## (fltarr(Nbeta)+1)
T3m[1,0,*,*] = transpose(Cf2) ## (sinc(2*beta1))
T3m[0,1,*,*] = T3m[1,0,*,*]
T3m[1,1,*,*] = transpose(Cf3) ## (1+sinc(4*beta1))
T3m[2,2,*,*] = transpose(Cf3) ## (1-sinc(4*beta1))

; calculate eigenvalues and vectors of T3m
T3m_eigv  = complexarr([Neps,Nbeta,3])
T3m_eigvc = complexarr([Neps,Nbeta,3,3])
for i=0,Neps-1 do begin
    for j=0,Nbeta-1 do begin
       ; note: the eigenvectors are returned as row-vectors!!
       T3m_eigv[i,j,*]=LA_EIGENQL(reform(T3m[*,*,i,j]),eigenvectors=eigvec)
       T3m_eigvc[i,j,*,*] = eigvec
    endfor
endfor
;stop

; calculate pseudo-probabilities p_i 
pT3m = fltarr([Neps,Nbeta,3])
pT3m[*,*,0] = abs(T3m_eigv[*,*,0])/total(T3m_eigv,3)
pT3m[*,*,1] = abs(T3m_eigv[*,*,1])/total(T3m_eigv,3)
pT3m[*,*,2] = abs(T3m_eigv[*,*,2])/total(T3m_eigv,3)

; calculate Anisotropy 
if keyword_set(anisotropy) then begin
  A = fltarr(Neps,Nbeta)+1
  p2p3 = reform(pT3m[*,*,1]+pT3m[*,*,0])
  p2p3_1 = reform(pT3m[*,*,1])
  p2p3_0 = reform(pT3m[*,*,0])
  ;stop
  AI = where(p2p3 ge 1e-6)
  A[AI] = (p2p3_1[AI]-p2p3_0[AI])/p2p3[AI]
  res = A 
  ;stop;
; calculate Entropy and alpha 
endif else begin
  HT3m = fltarr(Neps,Nbeta)
  for j=0,2 do begin
    pI  = where(pT3m[*,*,j] gt 0,/NULL)
    if n_elements(pI) gt 0 then begin
      pT3ms = reform(pT3m[*,*,j]) 
      HT3m[pI] += -pT3ms[pI]*alog(pT3ms[pI])/alog(3)
    end  
  end  
    
  alphaT3m = fltarr([Neps,Nbeta,3])
;  if total(where(abs(T3m_eigvc) gt 1)) gt 0 then print, 'eig_vec gt 1!!!'
  
  alphaT3m[*,*,0] = acos(abs(T3m_eigvc[*,*,0,0]))
  alphaT3m[*,*,1] = acos(abs(T3m_eigvc[*,*,0,1]))
  alphaT3m[*,*,2] = acos(abs(T3m_eigvc[*,*,0,2]))
  alphaT3m_mean = total(alphaT3m*pT3m,3)
  
  res = [HT3m, alphaT3m_mean]

endelse

return, res 
end