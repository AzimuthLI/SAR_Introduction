  pro polarimetry

;dataset = 'i_al_af_1206'
dataset = 'i_1206'
suffix = '' ; '_corr'

crop         = 0
fullpol      = 0
calc_sigma   = 0 
calc_coh     = 0
calc_ev_ana  = 0
calc_Xbragg  = 0
load_Xbragg  = 0
calc_Freeman = 1


; load data for different polarizations (HH, VV, HV,VH)
; this are elements of the matrix S = [[HH, HV], [VH, VV]]
print, 'load files...'
t0 = systime(1)
hh   = loadxdr(dataset + '_hh'+suffix+'.dat',/flipud)
vv   = loadxdr(dataset + '_vv'+suffix+'.dat',/flipud)
rla  = loadxdr(dataset + '_rla_xdr_float.dat',/float,/flipud)
if file_test(dataset + '_hv'+suffix+'.dat') and file_test(dataset + '_hv'+suffix+'.dat') then begin
  hv   = loadxdr(dataset + '_hv'+suffix+'.dat',/flipud)
  vh   = loadxdr(dataset + '_vh'+suffix+'.dat',/flipud)
end else xx   = loadxdr(dataset + '_xx'+suffix+'.dat',/flipud) 

; reduce hv +vh to xx if not fully polarimetric calculations is wanted.
if keyword_set(hv) and keyword_set(hv) then xx = (hv+vh)/2 

  
if crop then begin
  ddim_new = [400,800]
  hh = hh[0:ddim_new[0],0:ddim_new[1]]
  vv = vv[0:ddim_new[0],0:ddim_new[1]]
  if fullpol then begin
  hv = hv[0:ddim_new[0],0:ddim_new[1]]
  vh = vh[0:ddim_new[0],0:ddim_new[1]]
  end
  xx = xx[0:ddim_new[0],0:ddim_new[1]]
  rla = rla[0:ddim_new[0],0:ddim_new[1]]
end
ddim = size(hh,/dimensions)
fscale = sin(rla)*1e-6
fscl = sqrt(sin(rla)*1e-6)

if calc_sigma then begin
  ; calculate absolute powers
  print, '['+str(systime(1)-t0,'%5.1f s')+'] calculate sigma...' 
  ; use float() to convert complex -> float. abs(x)^2 = x*conj(x) 
  sigma0_hh = nlscale(float(hh*conj(hh))*fscale,'cliph_W->dB',range=[-40,5]) 
  sigma0_vv = nlscale(float(vv*conj(vv))*fscale,'cliph_W->dB',range=[-40,5])
  sigma0_xx = nlscale(float(xx*conj(xx))*fscale,'cliph_W->dB',range=[-40,5])
end

if calc_coh then begin
  print, '['+str(systime(1)-t0,'%5.1f s')+'] calculate coherences...'
  ; calculate polarimetric coherences
  c_hhvv = coherence(hh,vv,5)
  c_hhxx = coherence(hh,xx,5)
  c_vvxx = coherence(vv,xx,5)
  c_llrr = coherence((hh+complex(0,1)*vv)/sqrt(2),(hh-complex(0,1)*vv)/sqrt(2),5)
end

  print, '['+str(systime(1)-t0,'%5.1f s')+'] calculate C3...'
  ; calculate (monostatic) covariance matrix C_3 by Omega-taget vector:
  Omega = complexarr([ddim,3])
  Omega[*,*,0] = hh * fscl
  Omega[*,*,1] = sqrt(2)*xx * fscl
  Omega[*,*,2] = vv * fscl
  
  C3 = mprod(Omega,Omega,Navg=5)

print, '['+str(systime(1)-t0,'%5.1f s')+'] calculate T...'
; calculate pauli scattering vector k and lexicographic covariance matrix
pauli_k = complexarr([ddim,3+fullpol])
pauli_k[*,*,0] = (hh + vv)/sqrt(2) * fscl
pauli_k[*,*,1] = (hh - vv)/sqrt(2) * fscl
if fullpol then begin
  pauli_k[*,*,2] = (hv+vh)/sqrt(2) * fscl
  pauli_k[*,*,3] = (hv-vh)/sqrt(2) * fscl
end else pauli_k[*,*,2] = 2/sqrt(2) * xx * fscl
;stop
T = mprod(pauli_k,pauli_k,Navg=5)


print, '['+str(systime(1)-t0,'%5.1f s')+'] calculate total powers...'
; calculate the total power of the matrix S, C and T
if fullpol then S_total = fscale*float(hh*conj(hh) + vv*conj(vv) + hv*conj(hv)+vh*conj(vh)) $
else            S_total = fscale*float(hh*conj(hh) + vv*conj(vv) + 2*xx*conj(xx))

C_total = fltarr(ddim)
C3dim = (size(C3, /dim))[2]   
for i=0,C3dim-1 do C_total +=float(C3[*,*,i,i])

T_total = fltarr(ddim)
Tdim = (size(T, /dim))[2]
for i=0,Tdim-1 do T_total +=float(T[*,*,i,i])

; calculate Eigenvectors of T.
; unfortunately there is no way to do it simultaniously, 
; so each pixel has to be calculated separatly... :(
print, '['+str(systime(1)-t0,'%5.1f s')+'] calculating eigenvectors numerically...'
eig_val = complexarr([ddim,Tdim])
eig_vec = complexarr([ddim,Tdim,Tdim])

t_ev = systime(1)
for i=0,ddim[0]-1 do begin
  if (t_ev+10 lt systime(1)) then begin 
    t_ev = systime(1)
    print, '['+str(systime(1)-t0,'%5.1f s')+'] ' + str(float(i)/ddim[0]*100) + '%'
  end
  
  for j=0,ddim[1]-1 do begin
     ; note: the eigenvectors are returned as row-vectors!!
     ; note: eigenvalues are returned in ascending order: use reverse to get an descending order     
     eig_val[i,j,*]=reverse(LA_EIGENQL(reform(T[i,j,*,*]),eigenvectors=eigvec))
     eig_vec[i,j,*,*] = reverse(eigvec,2)
  endfor
endfor

if calc_ev_ana then begin 
  print, '['+str(systime(1)-t0,'%5.1f s')+'] calculate eigenvalues analytically...'
  ; eigenvalues analytical; the analytical eigenvectors don't work somehow... 
  eig_val_analyt = eigenvaluesT3(T)
end


; add noise filtering using 4th eigenvalues here!
ev_num_uc = eig_val[*,*,0:2] ; save uncorrected eigenvalues for plotting.
if fullpol then for j=0,2 do eig_val[*,*,j] = eig_val[*,*,j]-eig_val[*,*,3] 

print, '['+str(systime(1)-t0,'%5.1f s')+'] calculate pseudoprobabilities...'
p = fltarr([ddim,3])
p[*,*,0] = eig_val[*,*,0]/total(eig_val,3)
p[*,*,1] = eig_val[*,*,1]/total(eig_val,3)
p[*,*,2] = eig_val[*,*,2]/total(eig_val,3) 

print, '['+str(systime(1)-t0,'%5.1f s')+'] calculate Entropy H...'
H = -total(p*alog(p),3)/alog(3)

print, '['+str(systime(1)-t0,'%5.1f s')+'] calculate Anisotropy A...'
A = real((eig_val[*,*,1]-eig_val[*,*,2])/(eig_val[*,*,1]+eig_val[*,*,2]))

print, '['+str(systime(1)-t0,'%5.1f s')+'] calculate angles alpha1, alpha2, alpha3...'
; the first component of the eigenvectors is can be parametrisized as cos(alpha)*exp(i*phi).
; This can be used to retrive alpha.
alpha = fltarr([ddim,3])

alpha[*,*,0] = acos(abs(eig_vec[*,*,0,0]))
alpha[*,*,1] = acos(abs(eig_vec[*,*,0,1]))
alpha[*,*,2] = acos(abs(eig_vec[*,*,0,2])) 
; below is wrong. Ioles code is also wrong.
;alpha[*,*,0] = acos(abs(eig_vec[*,*,0,0]))
;alpha[*,*,1] = acos(abs(eig_vec[*,*,1,0]))
;alpha[*,*,2] = acos(abs(eig_vec[*,*,2,0])) 


; calculate averaged alpha
alpha_mean = total(alpha*p,3)

; calculate 2D Histogramms
print, '['+str(systime(1)-t0,'%5.1f s')+'] calculate 2D-Histograms...'  
Nbins = 100
h2d_Halpha = HIST_2D(H,alpha_mean,bin2=!pi/2/Nbins,bin1=1.0/Nbins,min1=0,min2=0,max2=!pi/2, max1=1)
h2d_HA     = HIST_2D(H,A,bin1=1.0/Nbins,bin2=1.0/Nbins,min1=0,min2=0,max1=1, max2=1)
h2d_ax= findgen(Nbins+1)/(Nbins+1)

; calculate X-Bragg model
; preselection of datapoints which allow inversion

; calculates surface roughness ks = 1-A (according to paper...)

; generate look-up-table of modelled matrix T(beta1, aoi,epsilon) with 
; aoi   = angle of incidence, 
; beta1 = angle in a plane perpendicular to scattering plane.
;         by the Unisotropy the scattering signal get's averaged within -beta1 and beta1. -> rotation of T matrix and average. 
; the lookuptable shows beta1(H,alpha). H,alpha calculated from T(beta1,..) (modelled)

; Calculation of T(beta1,..)
; calculate fresnel coefficients. epsilon = epsilon'-i*epsilon'' (complex dielectric constant)
; assume epsilon'' is zero. (correct at least for ice.. epsilon'' ~ 1e-1..1e-3)
 
print, '['+str(systime(1)-t0,'%5.1f s')+'] Invert X-Bragg model: epsilon, beta1 ...'

xbragg_bd = xbragg_boundary(rla)
  
  ; input: x 2xN array [epsilon, beta1] of column height N
; output: 2NxN array:[H(epsilon,beta1), alpha(epsilon,beta1)]
if calc_Xbragg then begin   
  ;stop
  epsilon = fltarr(ddim)
  beta1   = fltarr(ddim)
  inv = 0
  t_xbragg = systime(1)
  rla_min = min(rla)
  rla_max = max(rla)
  for i=0,ddim[0]-1 do begin
    for j=0,ddim[1]-1 do begin
      
      ; every 10 seconds print status and save epsilon and beta
      if (t_xbragg+10 lt systime(1)) then begin
        t_xbragg = systime(1)
        print, '['+str(systime(1)-t0,'%7.1f s')+'] ' + str(float(i*ddim[1]+j)/(ddim[0]*ddim[1])*100) + '%'
        save, epsilon, filename=dataset+'_epsilon.tmp'
        save, beta1, filename=dataset+'_beta1.tmp'  
      end

      ; define parameters of function lookupT which are fetched 
      ; by the scope_varfetch function inside lookupT3
      H0_fetch        = H[i,j]  
      alpha0_fetch    = alpha_mean[i,j]
      theta_aoi_fetch = rla[i,j]
      boundary_angleIndex = ceil((theta_aoi_fetch-rla_min)/(rla_max-rla_min)*99)
      
      ; preselect H for inversion
      Hlow = xbragg_bd.H_low[boundary_angleIndex,*]
      Hlow = Hlow(where(finite(Hlow)))
      HlowI = findfirst(Hlow ge H0_fetch)
  
      Hhigh = xbragg_bd.H_high[boundary_angleIndex,*]
      Hhigh = Hhigh(where(finite(Hhigh)))
      HhighI = findfirst(Hhigh ge H0_fetch)
      
      if (HlowI ne !NULL) and (HhighI ne !NULL) then begin
        ; preselect alpha for inversion             
        if (alpha0_fetch lt xbragg_bd.alpha_high[boundary_angleIndex,HhighI]) and (alpha0_fetch gt xbragg_bd.alpha_low[boundary_angleIndex,HlowI]) then begin
          inv = inv + 1
          skip_newton = 0
          CATCH, Error_status
          if Error_status ne 0 then begin
            PRINT, 'Error caught from newton: Errorstatus: ', Error_status + string(13b) + 'Error message: ', !ERROR_STATE.MSG
            CATCH, /CANCEL
            skip_newton = 1
          endif
          if (skip_newton eq 0) then res = NEWTON([10,!pi/4], 'Xbraggwrapper', CHECK=converged, ITMAX=100) else res = [0,0] ;[!VALUES.F_NAN,!VALUES.F_NAN]
        end else res = [0,0]
      end else  res = [0,0] 
      epsilon[i,j] = res[0]
      beta1[i,j] = res[1]
    endfor
  endfor
  CATCH, /CANCEL
  print, inv, ' pixel inverted'
end else if load_Xbragg then begin
  print, '['+str(systime(1)-t0,'%7.1f s')+'] Invert X-Bragg model: epsilon, beta1 loaded from harddisk, not calculated! ...'
  if file_test(dataset+'_eps.dat',/READ,/REGULAR) then restore, dataset+'_eps.dat', description=epsilon else print, 'File '''+dataset+'_epsilon.dat'' not found.'
  if file_test(dataset+'_beta.dat',/READ,/REGULAR) then restore, dataset+'_beta.dat', description=beta1 else print, 'File '''+dataset+'_beta1.dat'' not found.'         
  epsilon=eps
  delvar, eps
end
if (calc_Xbragg or load_Xbragg) then begin
  ; set NaN-results of Newton to zero.
  epsilon[where(~finite(epsilon))] = 0
  beta1[where(~finite(beta1))]     = 0
  
  print, '['+str(systime(1)-t0,'%7.1f s')+'] Invert X-Bragg model: Anisotropy ...'
  
  ; calculate modelled anisotropy
  A_Xbragg    = fltarr(ddim)
  epsNotZeroI = where(epsilon ne 0)
  
  for j=0,n_elements(epsNotZeroI)-1 do A_Xbragg[epsNotZeroI[j]] = Xbragg([epsilon[epsNotZeroI[j]]],beta1[epsNotZeroI[j]],rla[epsNotZeroI[j]],/anisotropy)
  print, '['+str(systime(1)-t0,'%7.1f s')+'] Invert X-Bragg model: Water content ...'
  
  ; calculate water content with a linear estimation of the paper "Electromagnetic determination of soil water content: Measurements in coaxial transmission lines" by Topp, Davis, Annan, WATER RESOURCES RESEARCH, VOL. 16, NO. 3, PP. 574-582, 1980
  mv = fltarr(ddim)
  epsgt2 = where(epsilon gt 2.5)
  ; empirical function for inversion
  mv[epsgt2] = 0.1*sqrt(epsilon[epsgt2]-2.5)-0.07   
  ; correct/clip water contend function
  mv[where(mv gt 1)] = 1
  mv[where(mv lt 0)] = 0
end   

if calc_Freeman then begin
  print, '['+str(systime(1)-t0,'%7.1f s')+'] Freeman-Durden decomposition...'
  print, 'use the four-component Yamaguchi decomposition based on Coherency Matrix T...'
  ; using section IV of the paper A Four Component Decomposition of POLSAR Images based on the Coherency Matrix T, 
  ; published by Yamaguchi, Yajima and Yamada: IEEE Geosc. and Remote Sensing Lett., vol.3,no.3, 2006
   
  ; helical scattering coefficient
  ; if fc > 0: right-circular
  ; if fc < 0:  left-circular
  fc = 2*imag(T[*,*,2,1])
  
  ; volume scattering coefficient
  fv = (4*abs(T[*,*,2,2]) - 2*abs(fc))
  
  fs = fltarr(ddim)
  fd = fltarr(ddim)
  alpha_square_fd = fltarr(ddim)
  beta_square_fd  = fltarr(ddim)
  
  ; calculate coefficients A, B, C
  A_fd = float(T[*,*,1,1] - T[*,*,2,2])
  B_fd = float(T[*,*,0,0] - 2*T[*,*,2,2]) + abs(fc)
  C_fd = abs(T[*,*,1,0])
  
  ; select dominant surface scattering
  surfaceIpos = where(real(smooth(hh*conj(vv),2,/edge_truncate)) gt 0,/NULL)
  ; calculate coefficients for dominant surface scatterint  
       fs[surfaceIpos] = B_fd[surfaceIpos]
  beta_square_fd[surfaceIpos] = abs(C_fd[surfaceIpos]/B_fd[surfaceIpos])^2
       fd[surfaceIpos] = A_fd[surfaceIpos] - abs(C_fd[surfaceIpos])^2/B_fd[surfaceIpos]
  
  ; indicator for double bounce indicates no surface scattering 
  surfaceIneg = where(real(smooth(hh*conj(vv),5,/edge_truncate)) lt 0,/NULL)
  ; calculate coefficients for dominant double-bounce scattering
        fs[surfaceIneg] = B_fd[surfaceIneg] - abs(C_fd[surfaceIneg])^2/A_fd[surfaceIneg]
  alpha_square_fd[surfaceIneg] = abs(C_fd[surfaceIneg]/A_fd[surfaceIneg])^2
        fd[surfaceIneg] = A_fd[surfaceIneg]
        
  ; finally, calculate scattering powers:
  P_fd_surf  = abs(fs)*(1+ beta_square_fd)
  P_fd_dihed = fd*(1+alpha_square_fd)
  P_fd_vol   = fv
  P_fd_helix = fc
  ;undefine, fc, fv, fs, fd, alpha_fd, beta_fd   


  print, '['+str(systime(1)-t0,'%7.1f s')+'] Freeman-Durden decomposition...'
  print, 'use the four-component Yamaguchi decomposition based on Covariance Matrix C...'
  ; using section IV of the paper A Four Component Decomposition of POLSAR Images based on the Coherency Matrix T, 
  ; published by Yamaguchi, Yajima and Yamada: IEEE Geosc. and Remote Sensing Lett., vol.3,no.3, 2006
  ; here, also the T-matrix is involved, as it serves with several quantities which are not required to calculate separately.
  Pc_fd = 2*abs(imag(T[*,*,2,1]))
  Pc_sign = (imag(T[*,*,2,1]) gt 0)*2-1
  c_square = reform(0.5*C3[*,*,1,1])
  P_volSelector = 10/alog(10)*alog(C3[*,*,2,2]/C3[*,*,0,0]); C22 = <|b|^2>, C00 = <|a|^2>
  lt_2db = where(P_volSelector lt -2,/NULL)
  in_2db = where((P_volSelector le 2) and (P_volSelector ge -2),/NULL)
  gt_2db = where(P_volSelector gt 2,/NULL)  
  
  S_fd  = fltarr(ddim)
  D_fd  = fltarr(ddim)
  C_fd  = fltarr(ddim)
  Pv_fd = fltarr(ddim)
  Pd_fd = fltarr(ddim)
  Ps_fd = fltarr(ddim)
  Pv_fd[lt_2db] = 15.0/2*c_square[lt_2db]-15.0/8*Pc_fd[lt_2db]
  Pv_fd[in_2db] =      8*c_square[in_2db]-2*Pc_fd[in_2db]
  Pv_fd[gt_2db] = 15.0/2*c_square[gt_2db]-15.0/8*Pc_fd[gt_2db]
  
  ; remove helix scattering
  Pv_fd[where(Pv_fd lt 0,/NULL)] = 0
    
  S_fd[lt_2db] = (T[*,*,0,0])[lt_2db]-0.5*Pv_fd[lt_2db]
  S_fd[in_2db] = (T[*,*,0,0])[in_2db]-4*c_square[in_2db]+Pc_fd[in_2db] 
  S_fd[gt_2db] = (T[*,*,0,0])[gt_2db]-0.5*Pv_fd[gt_2db]
  
  D_fd[lt_2db] = (T[*,*,1,1])[lt_2db]-7.0/4*c_square[lt_2db]-1.0/16*Pc_fd[lt_2db]
  D_fd[in_2db] = (T[*,*,1,1])[in_2db]-2*c_square[in_2db] 
  D_fd[gt_2db] = (T[*,*,1,1])[gt_2db]-7.0/4*c_square[gt_2db]-1.0/16*Pc_fd[gt_2db]
  
  C_fd[lt_2db] = (T[*,*,1,0])[lt_2db]-1.0/6*Pv_fd[lt_2db]
  C_fd[in_2db] = (T[*,*,1,0])[in_2db] 
  C_fd[gt_2db] = (T[*,*,1,0])[gt_2db]+1.0/6*Pv_fd[gt_2db]

  ; select high surface and dihedral scatterers 
  PsdI = where(Pv_fd + Pc_fd lt C_total,/NULL)
  PvcI = where(Pv_fd + Pc_fd ge C_total,/NULL)
  
  ; select double bounce and surface scatterers
  C0 = (C3[*,*,2,0])[PsdI] - c_square[PsdI] + 0.5*Pc_fd[PsdI]
  PdbI = PsdI[where(real(C0) lt 0,/NULL)]
  PsfI = PsdI[where(real(C0) ge 0,/NULL)]
  
  Ps_fd[PdbI] = S_fd[PdbI] - abs(C_fd[PdbI])^2/D_fd[PdbI]
  Pd_fd[PdbI] = D_fd[PdbI] + abs(C_fd[PdbI])^2/D_fd[PdbI]
  
  Ps_fd[PsfI] = S_fd[PsfI] + abs(C_fd[PsfI])^2/D_fd[PsfI]
  Pd_fd[PsfI] = D_fd[PsfI] - abs(C_fd[PsfI])^2/D_fd[PsfI]    
  
  ; select indices with four component of scattering types
  fd_4C  = PsdI[where((Ps_fd[PsdI] gt 0) and (Pd_fd[PsdI] gt 0),/NULL)]  
  ; three components, with no dihedral scatterers, but surface scatterers
  fd_3Cs = PsdI[where((Ps_fd[PsdI] gt 0) and (Pd_fd[PsdI] lt 0),/NULL)]
  Pd_fd[fd_3Cs] = 0
  Ps_fd[fd_3Cs] = C_total[fd_3Cs] - Pv_fd[fd_3Cs] - Pc_fd[fd_3Cs]
  
  ; three components, with no surface scatterers, but dihedral scatterers
  fd_3Cd = PsdI[where((Ps_fd[PsdI] lt 0) and (Pd_fd[PsdI] gt 0),/NULL)]
  Ps_fd[fd_3Cd] = 0
  Pd_fd[fd_3Cd] = C_total[fd_3Cd] - Pv_fd[fd_3Cd] - Pc_fd[fd_3Cd]
  
  ; two components: No surface or dihedral scatterers: Ps = Pd = 0
  Pv_fd[PvcI] = C_total[PvcI] - Pc_fd[PvcI]   
end  
  

print, '['+str(systime(1)-t0,'%7.1f s')+'] plotting results...'
@pol_plotdata


stop
end