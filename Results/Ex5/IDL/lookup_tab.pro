function lookup_tab , hv, sigmaa, kz, theta
; generation of look-up table of coherence
; 
; hv     : height
; sigmaa : extinction coefficent
; kz     : vertical wave number
; theta  : incidence angle

 
j = complex(0,1)
LupT = complexarr(N_elements(hv),N_elements(sigmaa))

for i=0, n_elements(hv)-1 do for l=0, n_elements(sigmaa)-1 do begin
  ex_eff = 2*sigmaa[l]/cos(theta)
  case 1 of
  (hv[i] ne 0) and (sigmaa[l] ne 0): LupT[i,l] = ex_eff/(exp(hv[i]*ex_eff)-1) * (exp(hv[i]*(j*kz+ex_eff))-1)/(j*kz+ex_eff)
  (hv[i] ne 0) and (sigmaa[l] eq 0):  LupT[i,l] = (exp(j*hv[i]*kz)-1)/(j*hv[i]*kz)
  (hv[i] eq 0): LupT[i,l] = 1
  end  
;  LupT[i,l] = ((hv[i] ne 0) and (sigmaa[l] ne 0) ? ex_eff/(exp(hv[i]*ex_eff)-1) * (exp(hv[i]*(j*kz+ex_eff))-1)/(j*kz+ex_eff) : 1)
  ;LupT[i,l] = ((hv[i] ne 0) and (sigmaa[l] eq 0) ? ex_eff/(exp(hv[i]*ex_eff)-1) * (exp(hv[i]*(j*kz+ex_eff))-1)/(j*kz+ex_eff) : 1) 
endfor

return, LupT

end 