pro polinsar
@loadparameters_and_data
coh_win =5; coherency window. 1..10 

; calculate coherence between picture pairs
c_hh12 = coherence(hh1,hh2,coh_win)
c_hh13 = coherence(hh1,hh3,coh_win)
c_vv12 = coherence(vv1,vv2,coh_win)
c_vv13 = coherence(vv1,vv3,coh_win)
c_xx12 = coherence(xx1,xx2,coh_win)
c_xx13 = coherence(xx1,xx3,coh_win)

; define pauli vector
k1 = sqrt(0.5)*[ [[hh1+vv1]], [[hh1-vv1]], [[2*xx1]] ]
k2 = sqrt(0.5)*[ [[hh2+vv2]], [[hh2-vv2]], [[2*xx2]] ]
k3 = sqrt(0.5)*[ [[hh3+vv3]], [[hh3-vv3]], [[2*xx3]] ]

; calculate pauli coherences
c_k112 = coherence(k1[*,*,0],k2[*,*,0],coh_win)
c_k113 = coherence(k1[*,*,0],k3[*,*,0],coh_win)
c_k212 = coherence(k1[*,*,1],k2[*,*,1],coh_win)
c_k213 = coherence(k1[*,*,1],k3[*,*,1],coh_win)

; flat earth correction, analytically
; an additional phase is added to each dataset so that there are 
; no big phase differents in the interferograms any more.  
hh2f = rm_flatearth(hh2, k, H, B_12, theta, grng_res,0.228*lambda)
vv2f = rm_flatearth(vv2, k, H, B_12, theta, grng_res,0.228*lambda)
xx2f = rm_flatearth(xx2, k, H, B_12, theta, grng_res,0.228*lambda)

hh3f = rm_flatearth(hh3, k, H, B_13, theta, grng_res,0.006*lambda) 
vv3f = rm_flatearth(vv3, k, H, B_13, theta, grng_res,0.006*lambda) 
xx3f = rm_flatearth(xx3, k, H, B_13, theta, grng_res,0.006*lambda)

hh3f2 = rm_flatearth(hh3, k, H, B_13-B_12, theta, grng_res,0.277*lambda)
vv3f2 = rm_flatearth(vv3, k, H, B_13-B_12, theta, grng_res,0.277*lambda)
xx3f2 = rm_flatearth(xx3, k, H, B_13-B_12, theta, grng_res,0.277*lambda)


; calculate coherence between picture pairs
c_hh12f = coherence(hh1,hh2f, coh_win)
c_hh13f = coherence(hh1,hh3f, coh_win)
c_vv12f = coherence(vv1,vv2f, coh_win)

c_vv13f = coherence(vv1,vv3f, coh_win)
c_xx12f = coherence(xx1,xx2f, coh_win)
c_xx13f = coherence(xx1,xx3f, coh_win)

c_hh23f = coherence(hh2,hh3f2,coh_win)
c_vv23f = coherence(vv2,vv3f2,coh_win)
c_xx23f = coherence(xx2,xx3f2,coh_win)


; define pauli vector
k2f = sqrt(0.5)*[ [[hh2f+vv2f]], [[hh2f-vv2f]], [[2*xx2f]] ]
k3f = sqrt(0.5)*[ [[hh3f+vv3f]], [[hh3f-vv3f]], [[2*xx3f]] ]
k20 = sqrt(0.5)*[ [[hh2+vv2]], [[hh2-vv2]], [[2*xx2]] ]
k23f= sqrt(0.5)*[ [[hh3f2+vv3f2]], [[hh3f2-vv3f2]], [[2*xx3f2]] ]

; calculate pauli coherences
c_k112f = coherence( k1[*,*,0], k2f[*,*,0],coh_win)
c_k113f = coherence( k1[*,*,0], k3f[*,*,0],coh_win)
c_k212f = coherence( k1[*,*,1], k2f[*,*,1],coh_win)
c_k213f = coherence( k1[*,*,1], k3f[*,*,1],coh_win)
c_k123f = coherence(k20[*,*,1],k23f[*,*,1],coh_win)
c_k223f = coherence(k20[*,*,1],k23f[*,*,1],coh_win)


; calculate vertical wave-vector, which is the interference-wave-vector
; of the two Beams from each side of the baseline 
kz12     = k/sin(theta) * (B_12)*cos(theta)^2/H
kz13     = k/sin(theta) * (B_13)*cos(theta)^2/H

; set parameters for lookup-table   
hv_min    =  0   ; (m) mimimal height of trees
hv_max12    = !pi/abs(kz12); 30   ; (m) maximal height of trees: set to height sensitivity range of interferometry.
hv_max13    = 1.6*!pi/abs(kz13); 30   ; (m) maximal height of trees
sigma_min = 0.00 ; (dB/m) mimimal extinction from trees
sigma_max = 2.00 ; (dB/m) maximal extinction from trees
hv_N      = 50   ; number of cells of lookup table
sigma_N   = 50   ; number of cells of lookup table
; calculate height and extinction axis for look-up table
hv12     = findgen(hv_N)/(hv_N-1)*(hv_max12-hv_min) + hv_min
hv13     = findgen(hv_N)/(hv_N-1)*(hv_max13-hv_min) + hv_min
sigma  = (findgen(sigma_N)/(sigma_N-1)*(sigma_max-sigma_min) + sigma_min)/8.686

; calculate look-up-tables
c_lut12  = lookup_tab(hv12, sigma, kz12, theta)
c_lut13  = lookup_tab(hv13, sigma, kz13, theta)
c_lut23  = c_lut12


for l = 0,2 do begin
  ; generate five-element coherency vector
  case l of
    0: begin 
         c5 = [[[c_hh12f]],[[c_vv12f]],[[c_xx12f]],[[c_k112f]],[[c_k212f]]]
         clut = c_lut12
         hv = hv12
       end
    1: begin
         c5 = [[[c_hh13f]],[[c_vv13f]],[[c_xx13f]],[[c_k113f]],[[c_k213f]]]
         clut = c_lut13
         hv = hv13
       end
    2: begin
         c5 = [[[c_hh23f]],[[c_vv23f]],[[c_xx23f]],[[c_k123f]],[[c_k223f]]]
         clut = c_lut23
         hv = hv12 ; the same, as hv23 would be.
       end
  end
  
  ground_phase = fltarr(ddim)     
  gamma_vol    = complexarr(ddim)
  h_veg        = fltarr(ddim)
  extin        = fltarr(ddim)
  
  for i=0, ddim[0]-1 do for j=0,ddim[1]-1 do begin
    ; calculate orientation of 5 (or more) points in complex unit circle 
    phi = cloudorientation(reform(c5[i,j,*]),significance=significance)
    if significance gt 0.2 then begin ; significance varies between 0 (circular distribution) and 1 (straight line)
      c5r = reform(c5[i,j,*])*exp(-complex(0,1)*phi); rotate by negative phi, because if phi < 0 the ellipse is rotated clockwise
      ; as the coordinatesystem has been rotated about phi the distances to x=0 on the straight line are just the real part
      mu = real(c5r)
      ; the section on the y-axis is the mean() of the imaginary component of c5r  
      ; find the center of xx-components vs. the non-xx components to get a direction where the ground phase is.
      ; the ground phase in the not rotated coordinate system is the angle of intersection of the horizontal line with the circle + phi.
      if (mean(mu[2]) gt mean(mu[[0,1]])) $
        then ground_phase[i,j] = acos(mean(imag(c5r))) + phi +!pi/2 $
        else ground_phase[i,j] = asin(mean(imag(c5r))) + phi
      ground_phase[i,j] = ((ground_phase[i,j] + !pi) mod (2*!pi)) - !pi ; restrict ground-phase to -pi..pi
    end else begin
      ; if the ground-phase could not be determined properly, assume a ground phase of zero.
      ; this is true, while it is sure, that the ground-phase is zero: This is not the case for any arbitrary topography.
      ground_phase[i,j] = 0 ;!VALUES.F_NAN
    end    

    ; find left-most and right-most point on the line
    mumax = max(mu,muImax,SUBSCRIPT_MIN=muImin)
    ; depending, where the groundpoint is, choose the most distant coherence, which is normally the XX-phase
    ; note, that the cohence is selected in the not-rotated coordinate system.
    ; As the lookup-table is calculated for zero ground-phase, 
    ; we have to rotate the estimated coherence c5[..] by the groundphase in clockwise direction, 
    ; towards the ground-phase point.
    ; ????? Minus/or Plus in exp? 
    if (mean(mu[2]) gt mean(mu[[0,1]])) $
      then gamma_vol[i,j] = c5[i,j,muImax]*exp( complex(0,1)*ground_phase[i,j]) $
      else gamma_vol[i,j] = c5[i,j,muImin]*exp(-complex(0,1)*ground_phase[i,j])
          
  
    ; inversion
    if finite(ground_phase[i,j]) then begin
      error      = min(abs(gamma_vol[i,j]-clut),minI) ; proximum of lookup-table with measured data
      min_pos    = array_indices(clut,minI)
      h_veg[i,j] = hv[min_pos[0]]
      extin[i,j] = sigma[min_pos[1]]*8.686; Neper to dB/m           
    end
  
  end   
  
  case l of
    0: begin 
         gp12        = ground_phase
         gamma_vol12 = gamma_vol
         h_veg12     = h_veg
         extin12     = extin                   
       end
    1: begin 
         gp13        = ground_phase 
         gamma_vol13 = gamma_vol
         h_veg13     = h_veg
         extin13     = extin                   
       end
    2: begin 
         gp23        = ground_phase 
         gamma_vol23 = gamma_vol
         h_veg23     = h_veg
         extin23     = extin                            
       end
  end
end

@plotdata
stop
end