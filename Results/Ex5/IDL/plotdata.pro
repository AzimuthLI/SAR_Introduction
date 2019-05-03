
hue=hue_table(0) ; create circular color table
wid = objarr(10)  & wi=0
wti = objarr(10)  & wt=0
img = objarr(100) & im=0 
cba = objarr(100) & cb=0
leg = objarr(10)  & lg=0
plt = objarr(100) & pl=0
 
win = [0,1,1,1,1]
histos = 1

if win[0] then begin
  wid[wi++] = window(window_title = sprintf('PolInSAR: Source data: %i x %i px',ddim),dim=wdim, loc=wloc) 
  wti[wt++] = text(0.25,0.96,'PolInSAR: Source data',alignment=0.5,font_size=12)
  
  img[im++] = image(clip(db(abs(hh1*pscl)^2),rg=[-40,0]),pos=imgpos(661,scale=[1,1]),title='$HH_1$',/current)
  img[im++] = image(clip(db(abs(vv1*pscl)^2),rg=[-40,0]),pos=imgpos(662,scale=[1,1]),title='$VV_1$',/current)
  img[im++] = image(clip(db(abs(xx1*pscl)^2),rg=[-40,0]),pos=imgpos(663,scale=[1,1]),title='$XX_1$',/current)
  
  img[im++]   = image(clip(db(abs(hh2*pscl)^2),rg=[-40,0]),pos=imgpos(667),title='$HH_2$',/current)
  img[im++]   = image(clip(db(abs(vv2*pscl)^2),rg=[-40,0]),pos=imgpos(668),title='$VV_2$',/current)
  img[im++]   = image(clip(db(abs(xx2*pscl)^2),rg=[-40,0]),pos=imgpos(669),title='$XX_2$',/current)
  
  img[im++]   = image(clip(db(abs(hh3*pscl)^2),rg=[-40,0]),pos=imgpos(66,13),title='$HH_3$',/current)
  img[im++]   = image(clip(db(abs(vv3*pscl)^2),rg=[-40,0]),pos=imgpos(66,14),title='$VV_3$',/current)
  img[im++]   = image(clip(db(abs(xx3*pscl)^2),rg=[-40,0]),pos=imgpos(66,15),title='$XX_3$',/current)
  for j=1,9 do cba[cb++] = colbar(img[im-j],'dB',rg=[-40,0],scl=0.7)
  
  
  w_title = text(0.75,0.96,'Coherences',alignment=0.5,font_size=12)
  img[im++]  = image(clip(abs(c_hh12),rg=[0,1],/ext),pos=imgpos(664),title='$\gamma_{(1,2), HH}$',/current)
  img[im++]  = image(clip(abs(c_vv12),rg=[0,1],/ext),pos=imgpos(665),title='$\gamma_{(1,2), VV}$',/current)
  img[im++]  = image(clip(abs(c_xx12),rg=[0,1],/ext),pos=imgpos(666),title='$\gamma_{(1,2), XX}$',/current)
  
  img[im++]  = image(clip(abs(c_hh13),rg=[0,1],/ext),pos=imgpos(66,16),title='$\gamma_{(1,3), HH}$',/current)
  img[im++]  = image(clip(abs(c_vv13),rg=[0,1],/ext),pos=imgpos(66,17),title='$\gamma_{(1,3), VV}$',/current)
  img[im++]  = image(clip(abs(c_xx13),rg=[0,1],/ext),pos=imgpos(66,18),title='$\gamma_{(1,3), XX}$',/current)
  for j=1,6 do cba[cb++] = colbar(img[im-j],scl=0.7)
  
  img[im++]  = image(atan(c_hh12,/p),pos=imgpos(66,10),title='$\gamma_{(1,2), HH}$',/current,rgb_table=hue)
  img[im++]  = image(atan(c_vv12,/p),pos=imgpos(66,11),title='$\gamma_{(1,2), VV}$',/current,rgb_table=hue)
  img[im++]  = image(atan(c_xx12,/p),pos=imgpos(66,12),title='$\gamma_{(1,2), XX}$',/current,rgb_table=hue)
  
  img[im++]  = image(atan(c_hh13,/p),pos=imgpos(66,22),title='$\gamma_{(1,3), HH}$',/current,rgb_table=hue)
  img[im++]  = image(atan(c_vv13,/p),pos=imgpos(66,23),title='$\gamma_{(1,3), VV}$',/current,rgb_table=hue)
  img[im++]  = image(atan(c_xx13,/p),pos=imgpos(66,24),title='$\gamma_{(1,3), XX}$',/current,rgb_table=hue)
  for j=1,6 do cba[cb++] = colbar(img[im-j],rg=[-!pi,!pi],scl=0.7,tickstr=['$-\pi$','$\pi$'],minor=1)
  
  ; plot pauli vector and lexicographic vector 
  img_k1 = fltarr([ddim,3])
  img_k2 = fltarr([ddim,3])
  img_k3 = fltarr([ddim,3])
  img_l1 = fltarr([ddim,3])
  img_l2 = fltarr([ddim,3])
  img_l3 = fltarr([ddim,3])
  pI = [1,2,0]
  for j=0,2 do begin
    img_k1[*,*,j] = clip(db(abs(k1[*,*,pI[j]]*pscl)^2),rg=[-40,5])
    img_k2[*,*,j] = clip(db(abs(k2[*,*,pI[j]]*pscl)^2),rg=[-40,5])
    img_k3[*,*,j] = clip(db(abs(k3[*,*,pI[j]]*pscl)^2),rg=[-40,5])
  end
    img_l1[*,*,1] = clip(db(abs(xx1*pscl)^2),rg=[-40,5])
    img_l2[*,*,1] = clip(db(abs(xx2*pscl)^2),rg=[-40,5])
    img_l3[*,*,1] = clip(db(abs(xx3*pscl)^2),rg=[-40,5])
    img_l1[*,*,0] = clip(db(abs(hh1*pscl)^2),rg=[-40,5])
    img_l2[*,*,0] = clip(db(abs(hh2*pscl)^2),rg=[-40,5])
    img_l3[*,*,0] = clip(db(abs(hh3*pscl)^2),rg=[-40,5])
    img_l1[*,*,2] = clip(db(abs(vv1*pscl)^2),rg=[-40,5])
    img_l2[*,*,2] = clip(db(abs(vv2*pscl)^2),rg=[-40,5])
    img_l3[*,*,2] = clip(db(abs(vv3*pscl)^2),rg=[-40,5])  
  
  ; pauli base
  img[im++]  = image(img_k1,pos=imgpos(66,19),title='$Pauli k_1$',/current)
  img[im++]  = image(img_k2,pos=imgpos(66,25),title='$Pauli k_2$',/current)
  img[im++]  = image(img_k3,pos=imgpos(66,31),title='$Pauli k_3$',/current)
  plt[pl]    = plothisto(img_k1[*,*,0],'-r',pos=imgpos(66,21,scale=[1,0.7]),hmax=hmax1,/clipedges)
  plt[pl]    = plothisto(img_k1[*,*,1],'-g',overplot=plt[pl],hmax=hmax2,/clipedges)
  plt[pl]    = plothisto(img_k1[*,*,2],'-b',overplot=plt[pl],hmax=hmax3,/clipedges,title='$histo k_1$')
  plt[pl].yrange = [0,max([hmax1,hmax2,hmax3])] & plt[pl].xrange = [-40,5] & pl++
  
  ; lexicographic base
  img[im++]  = image(img_l1,pos=imgpos(66,20),title='$Lexicogr. img1$',/current)
  img[im++]  = image(img_l2,pos=imgpos(66,26),title='$Lexicogr. img2$',/current)
  img[im++]  = image(img_l3,pos=imgpos(66,32),title='$Lexicogr. img3$',/current)
  
  ; plot coherences of pauli vectors
  
  img[im++]  = image(clip(abs(c_k112),rg=[0,1],/ext),pos=imgpos(66,27),title='$\gamma_{(1,2), HH+VV}$',/current)
  img[im++]  = image(clip(abs(c_k113),rg=[0,1],/ext),pos=imgpos(66,28),title='$\gamma_{(1,3), HH+VV}$',/current)
  img[im++]  = image(clip(abs(c_k212),rg=[0,1],/ext),pos=imgpos(66,29),title='$\gamma_{(1,2), HH-VV}$',/current)
  img[im++]  = image(clip(abs(c_k213),rg=[0,1],/ext),pos=imgpos(66,30),title='$\gamma_{(1,3), HH-VV}$',/current)
  for j=1,4 do cba[cb++] = colbar(img[im-j],scl=0.7)
  
  img[im++]  = image(atan(c_k112,/p),pos=imgpos(66,33),title='$\gamma_{(1,2), HH+VV}$',/current,rgb_table=hue)
  img[im++]  = image(atan(c_k113,/p),pos=imgpos(66,34),title='$\gamma_{(1,3), HH+VV}$',/current,rgb_table=hue)
  img[im++]  = image(atan(c_k212,/p),pos=imgpos(66,35),title='$\gamma_{(1,2), HH-VV}$',/current,rgb_table=hue)
  img[im++]  = image(atan(c_k213,/p),pos=imgpos(66,36),title='$\gamma_{(1,3), HH-VV}$',/current,rgb_table=hue)
  for j=1,4 do cba[cb++] = colbar(img[im-j],rg=[-!pi,!pi],scl=0.7,tickstr=['$-\pi$','$\pi$'],minor=1)
end

; Plot flat earth / baseline corrected data
if win[1] then begin
  ; plot flat_earth corrected coherences
  wid[wi++] = window(window_title = sprintf('PolInSAR: Source data: %i x %i px',ddim),dim=wdim, loc=wloc) 
  wti[wt++] = text(0.33,0.96,'Coherences (flat earth corrected)',alignment=0.5,font_size=12)
  
  ; amplitudes of coherences
  img[im++]  = image(clip(abs(c_hh12f),rg=[0,1],/ext),pos=imgpos(661),title='$|\gamma_{(1,2), HH}|$',/current)
  img[im++]  = image(clip(abs(c_vv12f),rg=[0,1],/ext),pos=imgpos(662),title='$|\gamma_{(1,2), VV}|$',/current)
  img[im++]  = image(clip(abs(c_xx12f),rg=[0,1],/ext),pos=imgpos(663),title='$|\gamma_{(1,2), XX}|$',/current)
  
  img[im++]  = image(clip(abs(c_hh13f),rg=[0,1],/ext),pos=imgpos(66,13),title='$|\gamma_{(1,3), HH}|$',/current)
  img[im++]  = image(clip(abs(c_vv13f),rg=[0,1],/ext),pos=imgpos(66,14),title='$|\gamma_{(1,3), VV}|$',/current)
  img[im++]  = image(clip(abs(c_xx13f),rg=[0,1],/ext),pos=imgpos(66,15),title='$|\gamma_{(1,3), XX}|$',/current)
  
  img[im++]  = image(clip(abs(c_k112),rg=[0,1],/ext),pos=imgpos(66,25),title='$|\gamma_{(1,2), HH+VV}|$',/current)
  img[im++]  = image(clip(abs(c_k113),rg=[0,1],/ext),pos=imgpos(66,26),title='$|\gamma_{(1,3), HH+VV}|$',/current)
  img[im++]  = image(clip(abs(c_k212),rg=[0,1],/ext),pos=imgpos(66,27),title='$|\gamma_{(1,2), HH-VV}|$',/current)
  img[im++]  = image(clip(abs(c_k213),rg=[0,1],/ext),pos=imgpos(66,28),title='$|\gamma_{(1,3), HH-VV}|$',/current)
  for j=1,10 do cba[cb++] = colbar(img[im-j],scl=0.7)
  
  
  img[im++]  = image(clip(atan(c_hh12f,/p),rg=[-!pi,!pi],/ext),pos=imgpos(66,7),title='$\phi(\gamma_{(1,2), HH)$',/current,rgb_table=hue)
  img[im++]  = image(clip(atan(c_vv12f,/p),rg=[-!pi,!pi],/ext),pos=imgpos(66,8),title='$\phi(\gamma_{(1,2), VV})$',/current,rgb_table=hue)
  img[im++]  = image(clip(atan(c_xx12f,/p),rg=[-!pi,!pi],/ext),pos=imgpos(66,9),title='$\phi(\gamma_{(1,2), XX})$',/current,rgb_table=hue)
  
  img[im++]  = image(clip(atan(c_hh13f,/p),rg=[-!pi,!pi],/ext),pos=imgpos(66,19),title='$\phi(\gamma_{(1,3), HH})$',/current,rgb_table=hue)
  img[im++]  = image(clip(atan(c_vv13f,/p),rg=[-!pi,!pi],/ext),pos=imgpos(66,20),title='$\phi(\gamma_{(1,3), VV})$',/current,rgb_table=hue)
  img[im++]  = image(clip(atan(c_xx13f,/p),rg=[-!pi,!pi],/ext),pos=imgpos(66,21),title='$\phi(\gamma_{(1,3), XX})$',/current,rgb_table=hue)
  
  img[im++]  = image(clip(atan(c_k112f,/p),rg=[-!pi,!pi],/ext),pos=imgpos(66,31),title='$\phi(\gamma_{(1,2), HH+VV})$',/current,rgb_table=hue)
  img[im++]  = image(clip(atan(c_k113f,/p),rg=[-!pi,!pi],/ext),pos=imgpos(66,32),title='$\phi(\gamma_{(1,3), HH+VV})$',/current,rgb_table=hue)
  img[im++]  = image(clip(atan(c_k212f,/p),rg=[-!pi,!pi],/ext),pos=imgpos(66,33),title='$\phi(\gamma_{(1,2), HH-VV})$',/current,rgb_table=hue)
  img[im++]  = image(clip(atan(c_k213f,/p),rg=[-!pi,!pi],/ext),pos=imgpos(66,34),title='$\phi(\gamma_{(1,3), HH-VV})$',/current,rgb_table=hue)
  for j=1,10 do cba[cb++] = colbar(img[im-j],rg=[-!pi,!pi],scl=0.7,tickstr=['$-\pi$','$\pi$'],minor=1)
  ;stop
  if histos eq 1 then begin
  plt[pl]    = plothisto(abs(c_hh12),':g',range=[0,1],/clipedges,pos=imgpos(66,4,scale=[1,0.7]),hmax=hmax0,title='$histo: |\gamma_{(1,2)}|$')
  plt[pl]    = plothisto(abs(c_vv12),':b',range=[0,1],/clipedges,overplot=plt[pl],hmax=hmax1)
  plt[pl]    = plothisto(abs(c_xx12),':r',range=[0,1],/clipedges,overplot=plt[pl],hmax=hmax2)
  plt[pl]    = plothisto(abs(c_hh12f),'g',range=[0,1],/clipedges,overplot=plt[pl],hmax=hmax3)
  plt[pl]    = plothisto(abs(c_vv12f),'b',range=[0,1],/clipedges,overplot=plt[pl],hmax=hmax4)
  plt[pl]    = plothisto(abs(c_xx12f),'r',range=[0,1],/clipedges,overplot=plt[pl],hmax=hmax5)
  plt[pl].xrange = [0,1] & plt[pl].yrange = [0, 0.5*max([hmax1, hmax2, hmax3, hmax4, hmax5])] & plt[pl].xtickvalues=[0,1] &  pl++
  
  plt[pl]    = plothisto(atan(c_hh12,/p),':g',pos=imgpos(66,10,scale=[1,0.7]),hmax=hmax0,title='$histo: \phi(\gamma_{(1,2)})$')
  plt[pl]    = plothisto(atan(c_vv12,/p),':b',overplot=plt[pl],hmax=hmax1)
  plt[pl]    = plothisto(atan(c_xx12,/p),':r',overplot=plt[pl],hmax=hmax2)
  plt[pl]    = plothisto(atan(c_hh12f,/p),'g',overplot=plt[pl],hmax=hmax3)
  plt[pl]    = plothisto(atan(c_vv12f,/p),'b',overplot=plt[pl],hmax=hmax4)
  plt[pl]    = plothisto(atan(c_xx12f,/p),'r',overplot=plt[pl],hmax=hmax5)
  plt[pl].yrange = [0,0.2*max([hmax1, hmax2, hmax3, hmax4, hmax5])] & plt[pl].xrange=[-!pi,!pi]
  plt[pl].xtickvalues=[-!pi,!pi] & plt[pl].xtickname=['$-\pi$','$\pi$']& pl++
  
  plt[pl]    = plothisto(abs(c_hh13),':g',range=[0,1],/clipedges,pos=imgpos(66,16,scale=[1,0.7]),hmax=hmax0,title='$histo: |\gamma_{(1,3)}|$')
  plt[pl]    = plothisto(abs(c_vv13),':b',range=[0,1],/clipedges,overplot=plt[pl],hmax=hmax1)
  plt[pl]    = plothisto(abs(c_xx13),':r',range=[0,1],/clipedges,overplot=plt[pl],hmax=hmax2)
  plt[pl]    = plothisto(abs(c_hh13f),'g',range=[0,1],/clipedges,overplot=plt[pl],hmax=hmax3)
  plt[pl]    = plothisto(abs(c_vv13f),'b',range=[0,1],/clipedges,overplot=plt[pl],hmax=hmax4)
  plt[pl]    = plothisto(abs(c_xx13f),'r',range=[0,1],/clipedges,overplot=plt[pl],hmax=hmax5)
  plt[pl].xrange = [0,1] & plt[pl].yrange = [0,0.1*max([hmax1, hmax2, hmax3, hmax4, hmax5])] & plt[pl].xtickvalues=[0,1] & pl++
  
  plt[pl]    = plothisto(atan(c_hh13,/p),':g',pos=imgpos(66,22,scale=[1,0.7]),hmax=hmax0,title='$histo: \phi(\gamma_{(1,3)})$')
  plt[pl]    = plothisto(atan(c_vv13,/p),':b',overplot=plt[pl],hmax=hmax1)
  plt[pl]    = plothisto(atan(c_xx13,/p),':r',overplot=plt[pl],hmax=hmax2)
  plt[pl]    = plothisto(atan(c_hh13f,/p),'g',overplot=plt[pl],hmax=hmax3)
  plt[pl]    = plothisto(atan(c_vv13f,/p),'b',overplot=plt[pl],hmax=hmax4)
  plt[pl]    = plothisto(atan(c_xx13f,/p),'r',overplot=plt[pl],hmax=hmax5)
  plt[pl].yrange = [0,0.33*max([hmax1, hmax2, hmax3, hmax4, hmax5])] & plt[pl].xrange=[-!pi,!pi]
  plt[pl].xtickvalues=[-!pi,!pi] & plt[pl].xtickname=['$-\pi$','$\pi$']& pl++
  
  plt[pl]    = plothisto(abs(c_k112),':b',pos=imgpos(66,23,scale=[1,0.7]),hmax=hmax0,xtickvalues=[0,1],title='$histo: |\gamma_{k(1,2)}|$')
  plt[pl]    = plothisto(abs(c_k212),':r',overplot=plt[pl],hmax=hmax1)
  plt[pl]    = plothisto(abs(c_k113),':b',overplot=plt[pl],hmax=hmax2)
  plt[pl]    = plothisto(abs(c_k213),':r',overplot=plt[pl],hmax=hmax3)
  plt[pl]    = plothisto(abs(c_k112f),'b',overplot=plt[pl],hmax=hmax4)
  plt[pl]    = plothisto(abs(c_k212f),'r',overplot=plt[pl],hmax=hmax5)
  plt[pl]    = plothisto(abs(c_k113f),'b',overplot=plt[pl],hmax=hmax6)
  plt[pl]    = plothisto(abs(c_k213f),'r',overplot=plt[pl],hmax=hmax7)
  plt[pl].xrange = [0,1] & plt[pl].yrange = [0,0.07*max([hmax1, hmax2, hmax3, hmax4, hmax5,hmax6,hmax7])] & plt[pl].xtickvalues=[0,1] & pl++
  
  plt[pl]    = plothisto(atan(c_k112,/p),':b',pos=imgpos(66,24,scale=[1,0.7]),hmax=hmax0,title='$histo: \phi(\gamma_{k(1,2)})$')
  plt[pl]    = plothisto(atan(c_k212,/p),':r',overplot=plt[pl],hmax=hmax1)
  plt[pl]    = plothisto(atan(c_k113,/p),':b',overplot=plt[pl],hmax=hmax2)
  plt[pl]    = plothisto(atan(c_k213,/p),':r',overplot=plt[pl],hmax=hmax3)
  plt[pl]    = plothisto(atan(c_k112f,/p),'b',overplot=plt[pl],hmax=hmax4)
  plt[pl]    = plothisto(atan(c_k212f,/p),'r',overplot=plt[pl],hmax=hmax5)
  plt[pl]    = plothisto(atan(c_k113f,/p),'b',overplot=plt[pl],hmax=hmax6)
  plt[pl]    = plothisto(atan(c_k213f,/p),'r',overplot=plt[pl],hmax=hmax7)
  plt[pl].yrange = [0,0.6*max([hmax1, hmax2, hmax3, hmax4, hmax5,hmax6,hmax7])] & plt[pl].xrange=[-!pi,!pi]
  plt[pl].xtickvalues=[-!pi,!pi] & plt[pl].xtickname=['$-\pi$','$\pi$']& pl++
  
  ;plt[pl]    = plotucirc(c_hh12, pos=imgpos(66, 5,scale=[1,0.7]),/current,title='$\gamma_{(1,2), HH}$','')
  plt[pl++]  = plotucirc(c_hh12f,'g',pos=imgpos(66, 5,scale=[1,0.7]),/current,title='$\gamma_{(1,2), HH}$')
  plt[pl++]  = plotucirc(c_vv12f,'b',pos=imgpos(66, 6,scale=[1,0.7]),/current,title='$\gamma_{(1,2), VV}$')
  plt[pl++]  = plotucirc(c_xx12f,'r',pos=imgpos(66,11,scale=[1,0.7]),/current,title='$\gamma_{(1,2), XX}$')
  
  plt[pl++]  = plotucirc(c_hh13f,'g',pos=imgpos(66,17,scale=[1,0.7]),/current,title='$\gamma_{(1,3), HH}$')
  plt[pl++]  = plotucirc(c_vv13f,'b',pos=imgpos(66,18,scale=[1,0.7]),/current,title='$\gamma_{(1,3), VV}$')
  plt[pl++]  = plotucirc(c_xx13f,'r',pos=imgpos(66,12,scale=[1,0.7]),/current,title='$\gamma_{(1,3), XX}$')
  
  plt[pl++]  = plotucirc(c_k112f,'b',pos=imgpos(66,29,scale=[1,0.7]),/current,title='$\gamma_{(1,2), HH+VV}$')
  plt[pl++]  = plotucirc(c_k113f,'r',pos=imgpos(66,30,scale=[1,0.7]),/current,title='$\gamma_{(1,3), HH+VV}$')
  plt[pl++]  = plotucirc(c_k212f,'b',pos=imgpos(66,35,scale=[1,0.7]),/current,title='$\gamma_{(1,2), HH-VV}$')
  plt[pl++]  = plotucirc(c_k213f,'r',pos=imgpos(66,36,scale=[1,0.7]),/current,title='$\gamma_{(1,3), HH-VV}$')
  end
end

if (win[2] eq 1) then begin
  wid[wi++] = window(window_title = sprintf('PolInSAR: Source data: %i x %i px, H=%f',ddim,H),dim=wdim, loc=wloc)
  img[im++]  = image(clip(atan(c_hh12f,/p),rg=[-!pi,!pi],/ext),pos=imgpos(44,1),title='$\gamma_{(1,2), HH}$',/current,rgb_table=hue)
  img[im++]  = image(clip(atan(c_vv12f,/p),rg=[-!pi,!pi],/ext),pos=imgpos(44,2),title='$\gamma_{(1,2), VV}$',/current,rgb_table=hue)
  img[im++]  = image(clip(atan(c_xx12f,/p),rg=[-!pi,!pi],/ext),pos=imgpos(44,3),title='$\gamma_{(1,2), XX}$',/current,rgb_table=hue)
  
  img[im++]  = image(clip(atan(c_hh13f,/p),rg=[-!pi,!pi],/ext),pos=imgpos(44,5),title='$\gamma_{(1,3), HH}$',/current,rgb_table=hue)
  img[im++]  = image(clip(atan(c_vv13f,/p),rg=[-!pi,!pi],/ext),pos=imgpos(44,6),title='$\gamma_{(1,3), VV}$',/current,rgb_table=hue)
  img[im++]  = image(clip(atan(c_xx13f,/p),rg=[-!pi,!pi],/ext),pos=imgpos(44,7),title='$\gamma_{(1,3), XX}$',/current,rgb_table=hue)

  img[im++]  = image(clip(atan(c_hh23f,/p),rg=[-!pi,!pi],/ext),pos=imgpos(44,9),title='$\gamma_{(2,3), HH}$',/current,rgb_table=hue)
  img[im++]  = image(clip(atan(c_vv23f,/p),rg=[-!pi,!pi],/ext),pos=imgpos(44,10),title='$\gamma_{(2,3), VV}$',/current,rgb_table=hue)
  img[im++]  = image(clip(atan(c_xx23f,/p),rg=[-!pi,!pi],/ext),pos=imgpos(44,11),title='$\gamma_{(2,3), XX}$',/current,rgb_table=hue)
  ;for j=1,9 do cba[cb++] = colbar(img[im-j],scl=1)
  
  plt[pl]    = plothisto(atan(c_hh12f,/p),'r',range=[-!pi,!pi],pos=imgpos(14,4,scale=[1,0.7]),hmax=hmax0,title='$histo: \phi(\gamma_{k(i,j)})$')
  plt[pl]    = plothisto(atan(c_vv12f,/p),'g',range=[-!pi,!pi],overplot=plt[pl],hmax=hmax1)
  plt[pl]    = plothisto(atan(c_xx12f,/p),'b',overplot=plt[pl],hmax=hmax2)
  plt[pl]    = plothisto(atan(c_hh13f,/p),':r',overplot=plt[pl],hmax=hmax3)
  plt[pl]    = plothisto(atan(c_vv13f,/p),':g',overplot=plt[pl],hmax=hmax4)
  plt[pl]    = plothisto(atan(c_xx13f,/p),':b',overplot=plt[pl],hmax=hmax5)
  plt[pl]    = plothisto(atan(c_hh23f,/p),'--r',overplot=plt[pl],hmax=hmax6)
  plt[pl]    = plothisto(atan(c_vv23f,/p),'--g',overplot=plt[pl],hmax=hmax7)
  plt[pl]    = plothisto(atan(c_xx23f,/p),'--b',overplot=plt[pl],hmax=hmax7)
  plt[pl].yrange = [0,7000]
  plt[pl].xrange = [-0.15,0.15]
  pl++
end

if (win[3] eq 1) then begin
  wid[wi++] = window(window_title ='Ground phases',dim=wdim, loc=wloc) 
  wti[wt++] = text(0.5,0.96,'PolInSAR: Ground phases',alignment=0.5,font_size=12)
  gp12[where(~finite(gp12))] = 0;!pi
  gp13[where(~finite(gp13))] = 0;!pi
  gp23[where(~finite(gp23))] = 0;!pi
  img[im++]  = image(clip(gp12,rg=[-!pi,!pi],/ext), pos=imgpos(55,1),title='$Grnd phase 1-2$',/current,rgb_table=hue)
  img[im++]  = image(clip(gp13,rg=[-!pi,!pi],/ext), pos=imgpos(55,2),title='$Grnd phase 1-3$',/current,rgb_table=hue)
  img[im++]  = image(clip(gp23,rg=[-!pi,!pi],/ext),pos=imgpos(55,3),title='$Grnd phase 2-3$',/current,rgb_table=hue)
  for j=1,3 do cba[cb++] = colbar(img[im-j],rg=[-!pi,!pi],scl=0.7,tickstr=['$-\pi$','$\pi$'],minor=1)
  plt[pl]    = plothisto(gp12,'r',pos=imgpos(55,4,scale=[1,0.63]),hmax=hmax1,title='histogram')
  plt[pl]    = plothisto(gp13,'b',o=plt[pl],hmax=hmax2)
  plt[pl]    = plothisto(gp23,'k',o=plt[pl],hmax=hmax3)
  plt[pl].xrange = [-!pi,!pi]
  plt[pl++].yrange = [0,max([hmax1, hmax2, hmax3])]
  
  img[im++]  = image(abs(gamma_vol12), pos=imgpos(55,6),title='$|\gamma_{vol,1-2}|$',/current,rgb_table=0)
  img[im++]  = image(abs(gamma_vol13), pos=imgpos(55,7),title='$|\gamma_{vol,1-3}|$',/current,rgb_table=0)
  img[im++]  = image(abs(gamma_vol23), pos=imgpos(55,8),title='$|\gamma_{vol,2-3}|$',/current,rgb_table=0)
  for j=1,3 do cba[cb++] = colbar(img[im-j],rg=[0,1],scl=0.7,minor=0)
  
  img[im++]  = image(atan(gamma_vol12,/p), pos=imgpos(55,11),title='$\phi(\gamma_{vol,1-2})$',/current,rgb_table=hue)
  img[im++]  = image(atan(gamma_vol13,/p), pos=imgpos(55,12),title='$\phi(\gamma_{vol,1-3})$',/current,rgb_table=hue)
  img[im++]  = image(atan(gamma_vol23,/p), pos=imgpos(55,13),title='$\phi(\gamma_{vol,2-3})$',/current,rgb_table=hue)
  for j=1,3 do cba[cb++] = colbar(img[im-j],rg=[-!pi,!pi],scl=0.7,tickstr=['$-\pi$','$\pi$'],minor=1)
  
  sel = [20,125,25,140]
  hveg12c = h_veg12[sel[0]:sel[1],sel[2]:sel[3]]
  hveg13c = h_veg13[sel[0]:sel[1],sel[2]:sel[3]]
  hveg23c = h_veg23[sel[0]:sel[1],sel[2]:sel[3]]
  extin12c = extin12[sel[0]:sel[1],sel[2]:sel[3]]
  extin13c = extin13[sel[0]:sel[1],sel[2]:sel[3]]
  extin23c = extin23[sel[0]:sel[1],sel[2]:sel[3]]
  gamma_vol12c = gamma_vol12[sel[0]:sel[1],sel[2]:sel[3]]
  gamma_vol13c = gamma_vol13[sel[0]:sel[1],sel[2]:sel[3]]
  gamma_vol23c = gamma_vol23[sel[0]:sel[1],sel[2]:sel[3]]
  gammaN = n_elements(gamma_vol12c)  
  
  plt[pl++]    = plothisto(atan(gamma_vol12c,/p),'r',pos=imgpos(55,14,scale=[1,0.63]),hmax=hmax1,title='histogram',name='$B_{12}$')
  plt[pl++]    = plothisto(atan(gamma_vol13c,/p),'b',o=plt[pl-2],hmax=hmax2,name='$B_{13}$')
  plt[pl]      = plothisto(atan(gamma_vol23c,/p),'k',o=plt[pl-2],hmax=hmax3,name='$B_{23}$')
  plt[pl].xrange   = [-!pi,!pi]
  plt[pl++].yrange = [0,max([hmax1, hmax2, hmax3])]
  leg[lg++]  = legend(target=[plt[pl+[-1,-2,-3]]],position=[0.8,0.5],font_size=10,sample_width=0.04) 
  

  img[im++]  = image(clip(h_veg12,rg=[0,25],/ext), pos=imgpos(55,16),title='$h_{veg,1-2}$',/current,rgb_table=13)
  img[im++]  = image(clip(h_veg13,rg=[0,25],/ext), pos=imgpos(55,17),title='$h_{veg,1-3}$',/current,rgb_table=13)
  img[im++]  = image(clip(h_veg23,rg=[0,25],/ext), pos=imgpos(55,18),title='$h_{veg,2-3}$',/current,rgb_table=13)
  for j=1,3 do cba[cb++] = colbar(img[im-j],'m',rg=[0,25],scl=0.7,minor=1)
  
 
  
  plt[pl]    = plothisto(hveg12c,'r',pos=imgpos(55,19,scale=[1,0.63]),hmax=hmax1,nbins=(hv_N+1),xtitle='height h (m)',title='histogram')
  plt[pl]    = plothisto(hveg13c,'b',o=plt[pl],hmax=hmax2,nbins=(hv_N+1))
  plt[pl]    = plothisto(hveg23c,'k',o=plt[pl],hmax=hmax3,nbins=(hv_N+1))
  plt[pl].xrange   = [0,25]
  plt[pl++].yrange = [0,max([hmax1, hmax2, hmax3])]
  t = text(0.81,0.34,'Average height',font_size=10)
  t = text(0.81,0.315,'$B_{12}:$  '+str(mean(hveg12c),'%4.1f')+' m',font_size=10)
  t = text(0.81,0.30, '$B_{13}:$  '+str(mean(hveg13c(where(hveg13c le 18))),'%4.1f')+' m',font_size=10)
  t = text(0.81,0.285,'$B_{23}:$  '+str(mean(hveg23c),'%4.1f')+' m',font_size=10)
    
  img[im++]  = image(extin12, pos=imgpos(55,21),title='$\sigma_{veg,1-2}$',/current,rgb_table=13)
  img[im++]  = image(extin13, pos=imgpos(55,22),title='$\sigma_{veg,1-3}$',/current,rgb_table=13)
  img[im++]  = image(extin23, pos=imgpos(55,23),title='$\sigma_{veg,2-3}$',/current,rgb_table=13)
  for j=1,3 do cba[cb++] = colbar(img[im-j],'dB',rg=[0,2],scl=0.7,minor=1)
  
  for j=1,6 do plt[pl++] = plot([sel[0],sel[1],sel[1],sel[0],sel[0]],[sel[2],sel[2],sel[3],sel[3],sel[2]],'w',axis_style=0,o=img[im-j])
  for j=7,9 do plt[pl++] = plot([sel[0],sel[1],sel[1],sel[0],sel[0]],[sel[2],sel[2],sel[3],sel[3],sel[2]],'k',axis_style=0,o=img[im-j])
  
  plt[pl]    = plothisto(extin12c,'r',pos=imgpos(55,24,scale=[1,0.63]),hmax=hmax1,nbins=(sigma_N+1),xtitle='extinction \sigma (dB)',title='histogram')
  plt[pl]    = plothisto(extin13c,'b',o=plt[pl],hmax=hmax2,nbins=(sigma_N+1))
  plt[pl]    = plothisto(extin23c,'k',o=plt[pl],hmax=hmax3,nbins=(sigma_N+1))
  plt[pl].xrange   = [0,2]
  plt[pl++].yrange = [0,0.07*max([hmax1, hmax2, hmax3])]
  t = text(0.81,0.15,'Average extinction',font_size=10)
  t = text(0.81,0.125,'$B_{12}:$  '+str(mean(extin12c),'%4.2f')+' dB',font_size=10)
  t = text(0.81,0.11, '$B_{13}:$  '+str(mean(extin13c),'%4.2f')+' dB',font_size=10)
  t = text(0.81,0.095,'$B_{23}:$  '+str(mean(extin23c),'%4.2f')+' dB',font_size=10)

end

if win[4] eq 1 then begin
  wid[wi++] = window(window_title ='Look-up table',dim=wdim, loc=wloc) 
  wti[wt++] = text(0.5,0.96,'PolInSAR: Look-up table',alignment=0.5,font_size=12)
  
  img[im++]  = image(clip(abs(c_lut12),rg=[0,1],/ext),             hv12,sigma*8.686,pos=imgpos(34,1,scale=0.7), title='$|\gamma_{\nu,1-2}(h_v,\sigma)|$',     xtitle='$height h_{\nu} (m)$',ytitle='$extinction \sigma  (dB)$',/current, rgb_table=13, aspect_ratio=hv_max12/sigma_max, axis_style=1, xrange=[hv_min,hv_max12], yrange=[sigma_min,sigma_max], xmajor=3,ymajor=3,xminor=0,yminor=0)
  img[im++]  = image(clip(atan(c_lut12,/phase),rg=[-!pi,!pi],/ext),hv12,sigma*8.686,pos=imgpos(34,2,scale=0.7), title='$\phi(\gamma_{\nu,1-2}(h_v,\sigma))$', xtitle='$height h_{\nu} (m)$',ytitle='$extinction \sigma  (dB)$',/current, rgb_table=hue, aspect_ratio=hv_max12/sigma_max, axis_style=1, xrange=[hv_min,hv_max12], yrange=[sigma_min,sigma_max], xmajor=3,ymajor=3,xminor=0,yminor=0)
  img[im++]  = image(clip(abs(c_lut13),rg=[0,1],/ext),             hv13,sigma*8.686,pos=imgpos(34,4,scale=0.7), title='$|\gamma_{\nu,1-3}(h_v,\sigma)|$',     xtitle='$height h_{\nu} (m)$',ytitle='$extinction \sigma  (dB)$',/current, rgb_table=13, aspect_ratio=hv_max13/sigma_max, axis_style=1, xrange=[hv_min,hv_max13], yrange=[sigma_min,sigma_max], xmajor=3,ymajor=3,xminor=0,yminor=0)
  img[im++]  = image(clip(atan(c_lut13,/phase),rg=[-!pi,!pi],/ext),hv13,sigma*8.686,pos=imgpos(34,5,scale=0.7), title='$\phi(\gamma_{\nu,1-3}(h_v,\sigma))$', xtitle='$height h_{\nu} (m)$',ytitle='$extinction \sigma  (dB)$',/current, rgb_table=hue, aspect_ratio=hv_max13/sigma_max, axis_style=1, xrange=[hv_min,hv_max13], yrange=[sigma_min,sigma_max], xmajor=3,ymajor=3,xminor=0,yminor=0)
  for j=2,4,2 do c = colbar(img[im-j],orient=1,fmt='%3.1f',minor=1)
  for j=1,3,2 do c = colbar(img[im-j],orient=1,fmt='%4.1f',rg=[-!pi,!pi],tickstr=['$-\pi$','$\pi$'],minor=1)
  
  phi = findgen(90)/89*2*!pi
  
  plt[pl]   = plot(reform(real(gamma_vol12c),gammaN),reform(imag(gamma_vol12c),gammaN),'.',color='#A0dfff',pos=imgpos(347,scale=0.8),aspect_ratio=1, title='$lut_{12}$',/c)
  plt[pl]   = plot(reform(real(c_lut12),hv_N*sigma_N),reform(imag(c_lut12),hv_N*sigma_N),'.k',/o)
  plt[pl++] = plot(cos(phi),sin(phi),/o,axis_style=0)
  
  plt[pl]   = plot(reform(real(gamma_vol13c),gammaN),reform(imag(gamma_vol13c),gammaN),'.',color='#A0dfff',pos=imgpos(348,scale=0.8),aspect_ratio=1, title='$lut_{13}$',/c) 
  plt[pl]   = plot(reform(real(c_lut13),hv_N*sigma_N),reform(imag(c_lut13),hv_N*sigma_N),'.k',/o)  
  plt[pl++] = plot(cos(phi),sin(phi),/o,axis_style=0)

  plt[pl]   = plot(reform(real(gamma_vol23c),gammaN),reform(imag(gamma_vol23c),gammaN),'.',color='#A0dfff',pos=imgpos(349,scale=0.8),aspect_ratio=1, title='$lut_{13}$',/c) 
  plt[pl]   = plot(reform(real(c_lut23),hv_N*sigma_N),reform(imag(c_lut23),hv_N*sigma_N),'.k',/o)  
  plt[pl++] = plot(cos(phi),sin(phi),/o,axis_style=0)

  cnt1 = contour(smooth(h_veg12,10,/edge_wrap),planar=0,rgb_table=8,/fill,c_value=findgen(50)/50*15-1,rgb_indices=[findgen(50)/50*255],axis_style=0,aspect_z=2,pos=imgpos(34,10,move=[-0.7,0.2],scale=2.1),/c)
  cnt1 = contour(smooth(h_veg13,10,/edge_wrap),planar=0,rgb_table=8,/fill,c_value=findgen(50)/50*15-1,rgb_indices=[findgen(50)/50*255],axis_style=0,aspect_z=2,pos=imgpos(34,11,move=[-0.8,-0.1],scale=2.1),/c)
  cnt1 = contour(smooth(h_veg23,10,/edge_wrap),planar=0,rgb_table=8,/fill,c_value=findgen(50)/50*15-1,rgb_indices=[findgen(50)/50*255],axis_style=0,aspect_z=2,pos=imgpos(34,12,move=[-0.9,-0.5],scale=2.1),/c)
end