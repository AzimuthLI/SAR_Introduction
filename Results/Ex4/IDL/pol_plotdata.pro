scrsize = get_screen_size()
wscr = scrsize[0]
hscr = scrsize[1]
showwin = [0,0,0, 0,0,0, 0,0,0, 1] 

if showwin[0] then begin
; part 1: plot absolute power, colorbars and histogramm.
w_data = window(window_title = 'data, 5px smoothed',dimensions=[wscr*0.36,hscr*0.8], location=[wscr*2/3,0])
tx = text(0.5,0.97,'data and pauli vector, 5px smoothed',alignment=0.5,font_size=12)

if calc_sigma then begin
  ; plot image HH 
  img1 = image(smooth(sigma0_hh,5,/edge_truncate),position=imgpos(3,2,1,move=[0,-0.1]),/current, Aspect_ratio=0.5,title='$\sigma_{0,HH}$')
  chh = colbar(img1, 'dB')
  h_hh = plothisto(sigma0_hh,range=[-40,5],pos=imgpos(3,2,1,scale=[1,0.15],move=[0,0.42]),/clipedge)
  
  ; plot image VV  
  img2 = image(smooth(sigma0_vv,5,/edge_truncate),position=imgpos(3,2,2,move=[0,-0.1]),/current, Aspect_ratio=0.5,title='$\sigma_{0,VV}$')
  cvv  = colbar(img2, 'dB')
  h_vv = plothisto(sigma0_vv,range=[-40,5],pos=imgpos(3,2,2,scale=[1,0.15],move=[0,0.42]),/clipedge)
  
  ; plot image XX
  img3 = image(smooth(sigma0_xx,5,/edge_truncate),position=imgpos(3,2,3,move=[0,-0.1]),/current, Aspect_ratio=0.5,title='$\sigma_{0,XX}$')
  cxx  = colbar(img3, 'dB')
  h_xx = plothisto(sigma0_xx,range=[-40,5],pos=imgpos(3,2,3,scale=[1,0.15],move=[0,0.42]),/clipedge)
end

; plot Pauli vector
imgdata4 = fltarr(ddim[0],ddim[1],3)
imgdata4[*,*,0] = smooth(nlscale(abs(pauli_k[*,*,1]),'cliph_W->dB',range=[-20,5]),5,/edge_truncate)
imgdata4[*,*,1] = smooth(nlscale(abs(pauli_k[*,*,2]),'cliph_W->dB',range=[-20,5]),5,/edge_truncate)
imgdata4[*,*,2] = smooth(nlscale(abs(pauli_k[*,*,0]),'cliph_W->dB',range=[-20,5]),5,/edge_truncate) 
imgdata4a = fltarr(ddim[0],ddim[1],3)
imgdata4a[*,*,0] = smooth(nlscale(abs(hh)*fscl,'cliph_W->dB',range=[-20,5]),5,/edge_truncate)
imgdata4a[*,*,1] = smooth(nlscale(abs(xx)*fscl,'cliph_W->dB',range=[-20,5]),5,/edge_truncate)
imgdata4a[*,*,2] = smooth(nlscale(abs(vv)*fscl,'cliph_W->dB',range=[-20,5]),5,/edge_truncate) 

img4 = image(imgdata4,position=imgpos(3,2,4),/current, Aspect_ratio=0.5,title='Pauli scattering vector')
h_R = plothisto(nlscale(abs(pauli_k[*,*,1]),'cliph_W->dB',range=[-20,5]),'-r2',pos=imgpos(3,2,5,scale=[1,0.65],move=[0,0.0]),/clipedge,name='HH - VV  (dihedral)')
h_G = plothisto(nlscale(abs(pauli_k[*,*,2]),'cliph_W->dB',range=[-20,5]),'-g2',/clipedge,overplot=h_R,name='2 HV        (volume)')
h_B = plothisto(nlscale(abs(pauli_k[*,*,0]),'cliph_W->dB',range=[-20,5]),'-b2',/clipedge,overplot=h_R,name='HH + VV (surface)')

img4a = image(imgdata4a,position=imgpos(3,2,6),/current, Aspect_ratio=0.5,title='Lexicographic vector')
h_Rl = plothisto(nlscale(abs(hh)*fscl,'cliph_W->dB',range=[-20,5]),'--r1',/clipedge,overplot=h_R,name='HH')
h_Gl = plothisto(nlscale(abs(xx)*fscl,'cliph_W->dB',range=[-20,5]),'--g1',/clipedge,overplot=h_R,name='XX')
h_Bl = plothisto(nlscale(abs(vv)*fscl,'cliph_W->dB',range=[-20,5]),'--b1',/clipedge,overplot=h_R,name='VV')
h_R.title = 'histogram of Pauli vector (-) !Clexicographic vector (--)'
l = legend(target=[h_R,h_G,h_B],position=[0.35,0.075],font_size=9,sample_width=0.04)
l = legend(target=[h_Rl,h_Gl,h_Bl],position=[0.58,0.075],font_size=9,sample_width=0.04)
undefine, imgdata4
end

if showwin[1] and calc_coh then begin
; part 2: Coherences
w_data = window(window_title = 'coherences',dimensions=[wscr*.36,hscr*0.8], location=[wscr*0.666,0])
tx = text(0.5,0.97,'5 px Coherences',alignment=0.5,font_size=12)

img6data = smooth(abs(c_hhvv),5,/edge_truncate)
img7data = smooth(abs(c_hhxx),5,/edge_truncate)
img8data = smooth(abs(c_vvxx),5,/edge_truncate)
img8adata = smooth(abs(c_llrr),5,/edge_truncate)
img6   = image(img6data,pos=imgpos(4,2,1,move=[0,-0.05]),/current, Aspect_ratio=0.5,title='$|\gamma_{HH,VV}|$')
chhvva = colbar(img6, '')
hc_hvp = plothisto(abs(c_hhvv),pos=imgpos(4,2,1,scale=[1,0.15],move=[0,0.4]),xtickvalues=[0,0.5,1],xtickname=['0','0.5','1'])

img7   = image(img7data,pos=imgpos(4,2,2,move=[0,-0.05]),/current, Aspect_ratio=0.5,title='$|\gamma_{HH,XX}|$')
chhxxa = colbar(img7, '')
hc_hvp = plothisto(abs(c_hhxx),pos=imgpos(4,2,2,scale=[1,0.15],move=[0,0.4]),xtickvalues=[0,0.5,1],xtickname=['0','0.5','1'])

img8   = image(img8data,pos=imgpos(4,2,3,move=[0,-0.05]),/current, Aspect_ratio=0.5,title='$|\gamma_{VV,XX}|$')
chhxxa = colbar(img8, '')
hc_hvp = plothisto(abs(c_vvxx),pos=imgpos(4,2,3,scale=[1,0.15],move=[0,0.4]),xtickvalues=[0,0.5,1],xtickname=['0','0.5','1'])

img8   = image(img8adata,pos=imgpos(4,2,4,move=[0,-0.05]),/current, Aspect_ratio=0.5,title='$|\gamma_{LL,RR}|$')
chhxxa = colbar(img8, '')
hc_hvp = plothisto(abs(c_llrr),pos=imgpos(4,2,4,scale=[1,0.15],move=[0,0.4]),xtickvalues=[0,0.5,1],xtickname=['0','0.5','1'])

undefine, img6data, img7data, img8data

img10  = image(smooth(atan(c_hhvv,/phase),5,/edge_truncate),position=imgpos(4,2,5,move=[0,0]),/current, Aspect_ratio=0.5,title='$\phi(\gamma_{HH,VV})$')
chhvvp = colbar(img10, '',tickstr = ['-!9p!x','!9p!x'])
hc_hvp = plothisto(atan(c_hhvv,/phase),range=[-!pi,!pi],pos=imgpos(4,2,5,scale=[1,0.15],move=[0,0.45]),xtickvalues=[-!pi,0.5,!pi],xtickname=['-!9p','0','!9p'])

img11  = image(smooth(atan(c_hhxx,/phase),5,/edge_truncate),position=imgpos(4,2,6,move=[0,0]),/current, Aspect_ratio=0.5,title='$\phi(\gamma_{HH,XX})$')
chhxxp = colbar(img11, '',tickstr = ['-!9p!x','!9p!x'])
hc_hxp = plothisto(atan(c_hhxx,/phase),range=[-!pi,!pi],pos=imgpos(4,2,6,scale=[1,0.15],move=[0,0.45]),xtickvalues=[-!pi,0.5,!pi],xtickname=['-!9p','0','!9p'])

img12  = image(smooth(atan(c_vvxx,/phase),5,/edge_truncate),position=imgpos(4,2,7,move=[0,0]),/current, Aspect_ratio=0.5,title='$\phi(\gamma_{VV,XX})$')
cvvxxp = colbar(img12, '',tickstr = ['-!9p!x','!9p!x'])
hc_vhp = plothisto(atan(c_vvxx,/phase),range=[-!pi,!pi],pos=imgpos(4,2,7,scale=[1,0.15],move=[0,0.45]),xtickvalues=[-!pi,0.5,!pi],xtickname=['-!9p','0','!9p'])

img12a  = image(smooth(atan(c_llrr,/phase),5,/edge_truncate),position=imgpos(4,2,8,move=[0,0]),/current, Aspect_ratio=0.5,title='$\phi(\gamma_{LL,RR})$')
cllrrp = colbar(img12a, '',tickstr = ['-!9p!x','!9p!x'])
hc_lrp = plothisto(atan(c_llrr,/phase),range=[-!pi,!pi],pos=imgpos(4,2,8,scale=[1,0.15],move=[0,0.45]),xtickvalues=[-!pi,0.5,!pi],xtickname=['-!9p','0','!9p'])

end

if showwin[2] then begin
; part 3: Covariance Matrix C3

w_data = window(window_title = 'C - matrix',dimensions=[wscr*.36,hscr*0.8], location=[wscr*0.666,0])
tx = text(0.5,0.97,'Covariance matrix C',alignment=0.5,font_size=12)

img13data = nlscale(abs(C3[*,*,0,0]),'cliph_W->dB',range=[-30,0])
img14data = nlscale(abs(C3[*,*,1,1]),'cliph_W->dB',range=[-30,0])
img15data = nlscale(abs(C3[*,*,2,2]),'cliph_W->dB',range=[-30,0])
img16data = nlscale(abs(C3[*,*,0,2]),'cliph_W->dB',range=[-40,0])
img17data = nlscale(abs(C3[*,*,1,2]),'cliph_W->dB',range=[-40,0])
img18data = nlscale(abs(C3[*,*,0,1]),'cliph_W->dB',range=[-40,0])
img19data = atan(C3[*,*,0,2],/phase)
img20data = atan(C3[*,*,1,2],/phase)
img21data = atan(C3[*,*,0,1],/phase)

img13  = image(img13data,position=imgpos(4,3,1),/current, Aspect_ratio=0.5,title='$C_{11}$')
cimg13 = colbar(img13, 'dB')
img14  = image(img14data,position=imgpos(4,3,2),/current, Aspect_ratio=0.5,title='$C_{22}$')
cimg14 = colbar(img14, 'dB')
img15  = image(img15data,position=imgpos(4,3,3),/current, Aspect_ratio=0.5,title='$C_{33}$')
cimg15 = colbar(img15, 'dB')

himg13 = plothisto(img13data,pos=imgpos(4,3,4,scale=[1,0.5],move=[0,0.125]),/clipedges, name='$|C_{11}|$',title='histogram')
himg14 = plothisto(img14data,/clipedges,overplot=himg13,color='red', name='$|C_{22}|$')
himg15 = plothisto(img15data,/clipedges,overplot=himg13,color='blue', name='$|C_{33}|$',xtickvalues=[-30,-20, -10, 0])
l = legend(target=[himg13,himg14,himg15],position=[0.825,0.73],font_size=9)


img16  = image(img16data,position=imgpos(4,3,5),/current, Aspect_ratio=0.5,title='$|C_{13}|$')
cimg16 = colbar(img16, 'dB')
img17  = image(img17data,position=imgpos(4,3,6),/current, Aspect_ratio=0.5,title='$|C_{23}|$')
cimg17 = colbar(img17, 'dB')
img18  = image(img18data,position=imgpos(4,3,7),/current, Aspect_ratio=0.5,title='$|C_{12}|$')
cimg18 = colbar(img18, 'dB')

himg16 = plothisto(img16data,pos=imgpos(4,3,8,scale=[1,0.5],move=[0,0.125]),/clipedges, name='$|C_{13}|$',title='histogram')
himg17 = plothisto(img17data,/clipedges,overplot=himg16,color='red', name='$|C_{23}|$')
himg18 = plothisto(img18data,/clipedges,overplot=himg16,color='blue', name='$|C_{12}|$',yrange=[0,4e4])
l = legend(target=[himg16,himg17,himg18],pos=[0.825,0.41],font_size=9)


img19  = image(img19data,position=imgpos(4,3,9),/current, Aspect_ratio=0.5,title='$\phi(C_{13})$')
cimg19 = colbar(img19, '', range=[-!pi,!pi],tickstr = ['-!9p','!9p'])
img20  = image(img20data,position=imgpos(4,3,10),/current, Aspect_ratio=0.5,title='$\phi(C_{23})$')
cimg20 = colbar(img20, '', range=[-!pi,!pi],tickstr = ['-!9p','!9p'])
img21  = image(img21data,position=imgpos(4,3,11),/current, Aspect_ratio=0.5,title='$\phi(C_{12})$')
cimg21 = colbar(img21, '', range=[-!pi,!pi],tickstr = ['-!9p','!9p'])

himg19 = plothisto(img19data,pos=imgpos(4,3,12,scale=[1,0.5],move=[0,0.125]), name='$|C_{13}|$',title='histogram',/norm)
himg20 = plothisto(img20data,overplot=himg19,color='red', name='$|C_{23}|$',/norm)
himg21 = plothisto(img21data,overplot=himg19,color='blue', name='$|C_{12}|$',xtickvalues=[-!pi,0,!pi],xtickname = ['-!9p','0','!9p'],/norm)
l = legend(target=[himg19,himg20,himg21],position=[0.825,0.1],font_size=9)

undefine, img13data, img14data, img15data, img16data, img17data, img18data, img19data, img20data, img21data
end

if showwin[3] then begin
; part 4: Coherency Matrix T

w_data = window(window_title = 'T - matrix',dimensions=[wscr*.36,hscr*0.8], location=[wscr*0.666,0])
tx = text(0.5,0.97,'Coherency Matrix T',alignment=0.5,font_size=12)

img22data = nlscale(abs(T[*,*,0,0]),'cliph_W->dB',range=[-33,5])
img23data = nlscale(abs(T[*,*,1,1]),'cliph_W->dB',range=[-33,5])
img24data = nlscale(abs(T[*,*,2,2]),'cliph_W->dB',range=[-33,5])
if fullpol then img24adata = nlscale(abs(T[*,*,3,3]),'cliph_W->dB',range=[-45,5])
img25data = nlscale(abs(T[*,*,0,2]),'cliph_W->dB',range=[-40,0])
img26data = nlscale(abs(T[*,*,1,2]),'cliph_W->dB',range=[-40,0])
img27data = nlscale(abs(T[*,*,0,1]),'cliph_W->dB',range=[-40,0])
img28data = atan(T[*,*,0,2],/phase)
img29data = atan(T[*,*,1,2],/phase)
img30data = atan(T[*,*,0,1],/phase)


img22  = image(img22data,position=imgpos(4,3,1),/current, Aspect_ratio=0.5,title='$T_{11}$')
cimg22 = colbar(img22, 'dB')
img23  = image(img23data,position=imgpos(4,3,2),/current, Aspect_ratio=0.5,title='$T_{22}$')
cimg23 = colbar(img23, 'dB')
img24  = image(img24data,position=imgpos(4,3,3),/current, Aspect_ratio=0.5,title='$T_{33}$')
cimg24 = colbar(img24, 'dB')
;stop
himg22 = plothisto(img22data,pos=imgpos(4,3,4,scale=[1,0.5],move=[0,0.125]),/clipedges, name='$|T_{11}|$',title='histogram',hmax=h22max)
himg23 = plothisto(img23data,/clipedges,overplot=himg22,color='red', name='$|T_{22}|$',hmax=h23max)
himg24 = plothisto(img24data,/clipedges,overplot=himg22,color='blue', name='$|T_{33}|$',xtickvalues=[-30,-20, -10, 0],hmax=h24max)
if fullpol then himg24a = plothisto(img24adata,/clipedges,overplot=himg22,color='green', name='$|T_{44}|$',xtickvalues=[-40,-20, 0],hmax=h24amax)
himg22.yrange = [0,max([h22max,h23max,h24max,(fullpol ? h24amax : [])])]
l = legend(target=[himg22,himg23,himg24,(fullpol ? himg24a : [])],position=[0.825,0.73],font_size=9)

img25  = image(img25data,position=imgpos(4,3,5),/current, Aspect_ratio=0.5,title='$|T_{13}|$')
cimg25 = colbar(img25, 'dB')
img26  = image(img26data,position=imgpos(4,3,6),/current, Aspect_ratio=0.5,title='$|T_{23}|$')
cimg26 = colbar(img26, 'dB')
img27  = image(img27data,position=imgpos(4,3,7),/current, Aspect_ratio=0.5,title='$|T_{12}|$')
cimg27 = colbar(img27, 'dB')

himg25 = plothisto(img25data,pos=imgpos(4,3,8,scale=[1,0.5],move=[0,0.125]),/clipedges, name='$|T_{13}|$',title='histogram',hmax=h25max)
himg26 = plothisto(img26data,/clipedges,overplot=himg25,color='red', name='$|T_{23}|$',hmax=h26max)
himg27 = plothisto(img27data,/clipedges,overplot=himg25,color='blue', name='$|T_{12}|$',yrange=[0,3.5e4],hmax=h27max)
himg25.yrange = [0,1.05*max([h25max,h26max,h27max])]
l = legend(target=[himg25,himg26,himg27],position=[0.825,0.41],font_size=9)


img28  = image(img28data,position=imgpos(4,3,9),/current, Aspect_ratio=0.5,title='$\phi(T_{13})$')
cimg28 = colbar(img28, '', range=[-!pi,!pi],tickstr = ['-!9p!x','!9p!x'])
img29  = image(img29data,position=imgpos(4,3,10),/current, Aspect_ratio=0.5,title='$\phi(T_{23})$')
cimg29 = colbar(img29, '', range=[-!pi,!pi],tickstr = ['-!9p!x','!9p!x'])
img30  = image(img30data,position=imgpos(4,3,11),/current, Aspect_ratio=0.5,title='$\phi(T_{12})$')
cimg30 = colbar(img30, '', range=[-!pi,!pi],tickstr = ['-!9p!x','!9p!x'])

himg28 = plothisto(img28data,pos=imgpos(4,3,12,scale=[1,0.5],move=[0,0.125]), name='$|T_{13}|$',title='histogram')
himg29 = plothisto(img29data,overplot=himg28,color='red', name='$|T_{23}|$')
himg30 = plothisto(img30data,overplot=himg28,color='blue', name='$|T_{12}|$',xtickvalues=[-!pi,0,!pi],xtickname = ['-!9p','0','!9p'])
l = legend(target=[himg28,himg29,himg30],position=[0.825,0.1],font_size=9)

undefine, imgdata22, imgdata23, imgdata24, imgdata25, imgdata26, imgdata27, imgdata28, imgdata29, imgdata30

end

if showwin[4] then begin
; part 5: Total Power
w_data = window(window_title = 'Total Power',dimensions=[wscr*.36,hscr*0.8], location=[wscr*0.666,0])
tx = text(0.5,0.97,'Total power (5 px smoothed)',alignment=0.5,font_size=12)

img31data = nlscale(smooth(S_total,5,/edge_truncate),'cliph_W->dB',range= [-20,5])
if calc_C3 then img32data = nlscale(C_total,'cliph_W->dB',range= [-20,5])
img33data = nlscale(T_total,'cliph_W->dB',range= [-20,5])

img31  = image(img31data,position=imgpos(3,2,1),/current, Aspect_ratio=0.5,title='$total power ||S||_2^2$')
cimg31 = colbar(img31, 'dB')
if calc_C3 then img32  = image(img32data,position=imgpos(3,2,2),/current, Aspect_ratio=0.5,title='$total power ||C_3||_2^2$')
if calc_C3 then cimg32 = colbar(img32, 'dB')
img33  = image(img33data,position=imgpos(3,2,3),/current, Aspect_ratio=0.5,title='$total power ||T_3||_2^2$')
cimg33 = colbar(img33, 'dB')

himg31 = plothisto(img31data,pos=imgpos(1,2,2),/clipedges, name='$||S||_2^2$',title='histogram',thick=3)
if calc_C3 then himg32 = plothisto(img32data,/clipedges, name='$||C_3||_2^2$',overplot=himg31,color='red',thick=2)
himg33 = plothisto(img33data,/clipedges, name='$||T_3||_2^2$',overplot=himg31,color='#00ff00',linestyle='--')
l = legend(target=[himg31, def(himg32), himg33],location=[0.45,-0.15])

tl    = strarr(7)
tl[0] = 'Note:!C!C'
tl[1] = '$as Span(S) = Tr(T_3) = Tr(C_3) !C!C'
tl[2] = '                     = |S_{HH}|^2 + |S_{VV}|^2 + 2|S_{XX}|^2$!C!C'
tl[3] = 'the total power in each basis must be the same.!C!C'
tl[4] = 'But: the histograms depend strongly if you average!C'
tl[5] = 'the data before or after the Watt to dB conversion!C'
tl[6] = 'as this is a nonlinear conversion!'
tx = text(0.23,0.03, tl[0]+tl[1]+tl[2]+tl[3]+tl[4]+tl[5]+tl[6],font_size=10)
undefine, img22data, img23data, img24data, img25data, img26data, img27data, img28data, img29data, img30data
end

if showwin[5] then begin
; part 6: Eigenvalues
w_data = window(window_title = 'Eigenvectors and values',dimensions=[wscr*.36,hscr*0.8], location=[wscr*0.666,0])
tx = text(0.5,0.97,'Eigenvectors and values',alignment=0.5,font_size=12)

img34data = nlscale(abs(ev_num_uc),'cliph_W->dB',range=[-45,5])
if calc_ev_ana then img35data = nlscale(abs(eig_val_analyt),'cliph_W->dB',range=[-45,5])
if calc_ev_ana then img37data = nlscale(abs(ev_num_uc-eig_val_analyt)/(0.5*abs(ev_num_uc+eig_val_analyt)),'cliph_W->dB',range=[-76,0])
img34  = image(img34data,position=imgpos(3,2,1),/current, Aspect_ratio=0.5,title='$Eigenvalues of T_3 !C (numeric)$')
if calc_ev_ana then img35  = image(img35data,position=imgpos(3,2,2),/current, Aspect_ratio=0.5,title='$Eigenvalues of T_3 !C (analytical)$')
if calc_ev_ana then img36  = image(img37data,position=imgpos(3,2,4),/current, Aspect_ratio=0.5,title='$rel. diff. analytic vs. numeric !C\Sigma_i(|\lambda_{ana,i} - \lambda_{num,i}|/0.5*|\lambda_{ana,i}+\lambda_{num,i}|)$',rgb_table=13)
;cimg36 = colbar(img36, 'dB', range=[-70,-35])
if fullpol then img36a  = image(nlscale(abs(eig_val[*,*,3]),'cliph_W->dB',range=[-45,-20]),position=imgpos(3,2,6),/current, Aspect_ratio=0.5,title='$4th eigenvalue \lambda_4$',rgb_table=1)
;cimg36a = colbar(img36a, 'dB', range=[-45,-30])

himg34 = plothisto(img34data[*,*,0],pos=imgpos(3,2,3,scale=[1,0.35],move=[0,0.15]),/norm, /clipedges, name='$\lambda_{num,1}$',color='red',thick=2,xtitle='dB',title='Histogram')
himg35 = plothisto(img34data[*,*,1],/clipedges,/norm,  name='$\lambda_{num,2}$',overplot=himg34,color='green',thick=2)
himg36 = plothisto(img34data[*,*,2],/clipedges,/norm,  name='$\lambda_{num,3}$',overplot=himg34,color='blue',thick=2)
if fullpol then himg36a= plothisto(nlscale(abs(eig_val[*,*,3]),'cliph_W->dB',range=[-45,-20]),/clipedges,/norm,  name='$\lambda_{num,4}$',overplot=himg34,color='black',thick=2)
if calc_ev_ana then begin
    himg37 = plothisto(img35data[*,*,0],/clipedges,/norm,  name='$\lambda_{ana,1}$',overplot=himg34,color='#ff8800',line_style=':',thick=0.5)
    himg38 = plothisto(img35data[*,*,1],/clipedges,/norm,  name='$\lambda_{ana,2}$',overplot=himg34,color='#00ff00',line_style=':',thick=0.5)
    himg39 = plothisto(img35data[*,*,2],/clipedges,/norm,  name='$\lambda_{ana,3}$',overplot=himg34,color='#66bbff',line_style=':',thick=0.5)
    
    himg40 = plothisto(img37data[*,*,0],pos=imgpos(3,2,5,scale=[1,0.65]),/norm,/clipedges, name='$(\lambda_{num,1} - \lambda_{ana,1}):\Sigma\lambda_{j,1}$',color='red',thick=2,xtitle='dB',title='Histogram !C relative error: analytic vs. numeric')
    himg41 = plothisto(img37data[*,*,1],/clipedges,/norm, name='$(\lambda_{num,2} - \lambda_{ana,2}):\Sigma\lambda_{j,2}$',overplot=himg40,color='green',thick=2)
    himg42 = plothisto(img37data[*,*,2],/clipedges,/norm, name='$(\lambda_{num,3} - \lambda_{ana,3}):\Sigma\lambda_{j,3}$',overplot=himg40,color='blue',thick=2)
end  
if ~fullpol then himg42.xrange = [-75,-35]
l_tg = [himg34, himg35, himg36]
if fullpol and calc_ev_ana then l_tg = [l_tg, himg36a]
l = legend(target=l_tg,position=[0.70,0.67],sample_width=0.04)
if calc_ev_ana then l = legend(target=[himg37, himg38, himg39],position=[0.83,0.67],sample_width=0.04)
if calc_ev_ana then l = legend(target=[himg40, himg41, himg42],position=[0.4,0.53],sample_width=0.04)
end

if showwin[6] then begin
  ; part 7: Entropy H, Anisotropy A, scattering parameter alpha, and alpha_mean
  w_data = window(window_title = 'H, A, alpha, <alpha>',dimensions=[wscr*.36,hscr*0.8], location=[wscr*0.666,0])
  tx = text(0.5,0.97,'Entropy, Anisotropy and scattering parameter <$\alpha$>',alignment=0.5,font_size=12)

  img40data = fltarr([ddim,3])
  img40data[*,*,0] = smooth(alpha[*,*,0],5,/edge_truncate)
  img40data[*,*,1] = smooth(alpha[*,*,1],5,/edge_truncate)
  img40data[*,*,2] = smooth(alpha[*,*,2],5,/edge_truncate)
  
  img38  = image(smooth(H,5,/edge_truncate),position=imgpos(3,2,1),/current, Aspect_ratio=0.5,title='$Entropy H$')
  img39  = image(smooth(A,5,/edge_truncate),position=imgpos(3,2,2),/current, Aspect_ratio=0.5,title='$Anisotropy A$')
  img40  = image(img40data,position=imgpos(3,2,4),/current, Aspect_ratio=0.5,title='$angle [\alpha_1,\alpha_2,\alpha_3]$')
  img41  = image(smooth(alpha_mean,5,/edge_truncate),position=imgpos(3,2,5),/current, Aspect_ratio=0.5,title='$average: <\alpha>$')
  cimg41 = colbar(img41, '', range=[0,!pi/2],tickstr = ['0','!9p!x/2'])
  cimg38 = colbar(img38, '', range=[0,1],tickstr = ['0','1'])
  cimg39 = colbar(img39, '', range=[0,1],tickstr = ['0','1'])
  
  himg43 = plothisto(H,pos=imgpos(3,2,3,scale=[1,0.65]),/clipedges, name='$H$',color='black',thick=2,title='Histogram')
  himg44 = plothisto(A,/clipedges, name='$A$',color='green',thick=2,overplot=himg43)
  himg45 = plothisto(alpha[*,*,0]*180/!pi,range=[0,90],pos=imgpos(3,2,6,scale=[1,0.65]),/clipedges, name='$\alpha_1$',color='red',thick=2,title='scattering parameter',xtickvalues=[0, 20, 40,50,70,90],hmax=h45max)  
  himg46 = plothisto(alpha[*,*,1]*180/!pi,range=[0,90],/clipedges, name='$\alpha_2$',color='green',thick=2,overplot=himg45,hmax=h46max)
  himg47 = plothisto(alpha[*,*,2]*180/!pi,range=[0,90],/clipedges, name='$\alpha_3$',color='blue',thick=2,overplot=himg45,hmax=h47max)
  himg48 = plothisto(alpha_mean*180/!pi,range=[0,90],/clipedges, name='$<\alpha>$',color='#aaaaaa',thick=2,overplot=himg45,hmax=h48max)
  himg45.yrange = [0,max([h45max,h46max,h47max,h48max])]
  tx = text(0.72,0.05,'surface',orientation=-45)
  tx = text(0.8,0.05,'dipol',orientation=-45)
  tx = text(0.88,0.05,'dihedral',orientation=-45)
  l = legend(target=[himg43, himg44],position=[0.88,0.85],sample_width=0.04)
  l = legend(target=[himg45, himg46, himg47,himg48],position=[0.685,0.035],sample_width=0.04,orientation=1,font_size=8)
  
end

if showwin[7] then begin
  ; part 8: Entropy H, Anisotropy A, scattering parameter alpha, and alpha_mean
  w_data = window(window_title = '2D Histograms',dimensions=[wscr*.36,hscr*0.8], location=[wscr*0.666,0])
  tx = text(0.5,0.97,'2D-Histograms',alignment=0.5,font_size=12)
  
  img42  = image(max(h2d_Halpha) - h2d_Halpha,h2d_ax,h2d_ax*90,xmajor=1,position=imgpos(1,2,1,scale=[0.85,0.85]),/current, aspect_ratio=0.7/90,title='$entropy-alpha, H-\alpha$',axis_style=2,xtickvalues=[0,0.5,1],ytickvalues=[0,30,60,90],yminor=3,xminor=4,yticklen=0.03,xticklen=0.03,xtitle='Entropy H',ytitle='<alpha> $\alpha$')
  img43  = image(max(h2d_HA) - h2d_HA,h2d_ax,h2d_ax,    position=imgpos(1,2,2,scale=[0.87,0.85]),/current,  aspect_ratio=0.7, title='$entropy-anisotropy, H-A$',axis_style=2,xtickvalues=[0,0.5,1],ytickvalues=[0,0.5,1],yminor=4,xminor=4,yticklen=0.03,xticklen=0.03,xtitle='Entropy H',ytitle = 'Anisotropy A')
  img42l = plot([0.5,0.5],[0,90],overplot=img42,color='#888888')
  img42l = plot([0,0.5],[42.5,42.5],overplot=img42,color='#888888')
  img42l = plot([0,0.5],[47.5,47.5],overplot=img42,color='#888888')
  
  img42l = plot([0.5,1],[40,40],overplot=img42,color='#888888')
  img42l = plot([0.5,0.9],[50,50],overplot=img42,color='#888888')
  img42l = plot([0.9,1],[55,55],overplot=img42,color='#888888')
  img42l = plot([0.9,0.9],[0,90],overplot=img42,color='#888888')
  
  img422 = plot([0.5,0.5],[0,1],overplot=img43,color='#888888')
  img422 = plot([0,1],[0.5,0.5],overplot=img43,color='#888888')
  img422 = plot([0.9,0.9],[0,1],overplot=img43,color='#888888')
  tx = text(0.28,0.69,'Bragg')
  tx = text(0.28,0.715,'dipole')
  tx = text(0.28,0.8,'dihedral')
  tx = text(0.55,0.83,'double reflection')
  tx = text(0.51,0.715,'anisotropic')
  tx = text(0.6,0.647,'random')
  tx = text(0.82,0.72,'random-!Canisotropic')
  tx = text(0.811,0.81,'complex')
  Halphacurve = halphaline()
  img42l = plot(Halphacurve(*,0),Halphacurve[*,1],overplot=img42,color='#888888')
  HAcurve = haline()
  img422 = plot(HAcurve(*,0),HAcurve[*,1],overplot=img43,color='#888888')
end

if showwin[8] and (calc_Xbragg or load_Xbragg) then begin
  ; part 9: plot absolute power, colorbars and histogramm.
  w_data = window(window_title = 'extended Bragg Model',dimensions=[wscr*0.36,hscr*0.8], location=[wscr*2/3,0])
  tx = text(0.5,0.97,string(format='(%"X-BRAGG: Extended Bragg Model: Invertion possible for %4.1f%% of total pixels.")',100.0*total(epsilon gt 0)/n_elements(epsilon)),alignment=0.5,font_size=12)
  
  ; plot lookup-table for inversion
  Nlookup = 20
  epsilon_bd = (findgen(Nlookup)/(Nlookup-1))^2*18+2
  beta1_bd   = findgen(Nlookup)/(Nlookup-1)*!pi/2
  res0 = Xbragg(epsilon_bd,beta1_bd,min(rla))
  res1 = Xbragg(epsilon_bd,beta1_bd,max(rla))
  A0   = Xbragg(epsilon_bd,beta1_bd,min(rla),/anisotropy)
  A1   = Xbragg(epsilon_bd,beta1_bd,max(rla),/anisotropy)
  Hlookup0     = res0[0:Nlookup-1,*]
  alphalookup0 = res0[Nlookup:Nlookup*2-1,*]
  Hlookup1     = res1[0:Nlookup-1,*]
  alphalookup1 = res1[Nlookup:Nlookup*2-1,*]
  ; plot data
  p0 = plot([0,0.7],[0,33],/nodata,xtitle='Entropy H',ytitle='$scattering angle \alpha$',title = '$H(\epsilon'',\beta_1) vs. \alpha(\epsilon'',\beta_1,\theta)$',$
       position=imgpos(2,2,1,scale=[0.75,0.45],move=[0.1,0.2]),/current,xminor=0,yminor=0,xticklen=0.02,yticklen=0.02)
  for j=0,Nlookup-1 do p0 = plot(Hlookup0[*,j],180*alphalookup0[*,j]/!pi,'r1+-',color=255*[float(j)/Nlookup, 0, 1-float(j)/Nlookup],sym_size=0.5,vert_color=indgen(Nlookup),rgb_table=13,overplot=p0)
  for j=0,Nlookup-1 do p0 = plot(Hlookup1[*,j],180*alphalookup1[*,j]/!pi,'r1+-',color=255*[float(j)/Nlookup, 0, 1-float(j)/Nlookup],sym_size=0.5,vert_color=indgen(Nlookup),rgb_table=13,overplot=p0)
  tx = text(0.183,3,'$\uparrow\epsilon'' = 2$',/DATA)
  tx = text(0.646,20.5,'$\epsilon'' = 100\rightarrow$',/DATA,baseline=[0,1],updir=[-1,0],alignment=1,vertical_alignment=0.45)
  tx = text(0,30,'$\downarrow\beta_1 = 0 $',/DATA,color=[120,20,80])
  tx = text(0.46,30,'$\beta_1 = \pi/2\downarrow$',/DATA,color='r')
  tx = text(0.3,20,'$\theta = 57°$',/DATA)
  tx = text(0.06,2,'$\theta = 24°$',/DATA,baseline=[0,1],updir=[-1,0])
  p0.yrange = [0,33] 
  p0.xrange = [0,0.7]
  
  ; create color table and save it in array ctab
  LOADCT, 13, NCOLORS=100, RGB_TABLE=ctab
  p1 = plot([0,0.7],[0,33],/nodata,xtitle='Entropy H',ytitle='$scattering angle \alpha$',title='Boundarys of Inversion',$
       position=imgpos(2,2,2,scale=[0.75,0.45],move=[0.1,0.2]),/current,xminor=0,yminor=0,xticklen=0.02,yticklen=0.02)
  j = 0.0
  while j le 99 do begin
    p1a = plot([transpose(xbragg_bd.H_low[floor(j),*]),reverse(transpose(xbragg_bd.H_high[floor(j),*]))],[transpose(xbragg_bd.alpha_low[floor(j),*]),reverse(transpose(xbragg_bd.alpha_high[floor(j),*]))]*180/!pi,color=transpose(ctab[floor(j),*]),overplot=p1)
    j += 100.0/12    
  end
  p1.xrange = [0,0.8]
  p1.yrange = [0,33]
  
  tx = text(0.64,0.72,'$\theta = 24°$')
  tx = text(0.79,0.76,'$\theta = 41°$',color='#00eeee')
  tx = text(0.88,0.82,'$\theta = 24°$',color='#ff6600')
  
  p2 = plot(beta1_bd*180/!pi,A0[*,0],'b2',xtitle='Width of $\beta_1$ distribution (deg)',ytitle='Anisotropy A',title = '$A(\epsilon'',\beta_1) vs. \alpha(\epsilon'',\beta_1,\theta)$',$
       position=imgpos(2,2,1,scale=[0.75,0.45],move=[0.1,-0.45]),/current,xminor=0,yminor=0,xticklen=0.02,yticklen=0.02)
  p2 = plot(beta1_bd*180/!pi,A0[*,Nlookup-1],overplot=p2,'b')  
  p2 = plot(beta1_bd*180/!pi,A1[*,0],overplot=p2,'r')
  p2 = plot(beta1_bd*180/!pi,A1[*,Nlookup-1],overplot=p2,'r1')
  
  p3 = plot(epsilon_bd, 0.1*sqrt(epsilon_bd-2.5)-0.07,'b2',xtitle='dielectric constant $\epsilon_r$',ytitle='Soil moisture content $m_v$',title = '$m_v(\epsilon'')$',$
       position=imgpos(2,2,2,scale=[0.75,0.45],move=[0.1,-0.45]),/current,xminor=0,yminor=0,xticklen=0.02,yticklen=0.02)
  ks = 1-A_xbragg
  ks[where(ks eq 1)] = 0
  img44 = image(enhancepx(enhancepx(ks)),position=imgpos(3,3,7,scale=[0.9,0.9],move=[0,0.05]),/current, aspect_ratio=0.5,title='$Surface roughness ks = 1 - A$',rgb_table=13)  
  img45 = image(enhancepx(enhancepx(mv)),position=imgpos(3,3,8,scale=[0.9,0.9],move=[0,0.05]),/current, aspect_ratio=0.5,title='$Surface moisture contend m_v$',rgb_table=13)
  cimg44 = colbar(img44, '', range=[0,1])
  cimg45 = colbar(img45, '', range=[0,1])
  himg48 = plothisto(ks,/clipedges, name='$ks$',color='#00ff00',thick=1,position=imgpos(3,3,9,scale=[0.9,0.9],move=[0,0.05]),/current,title='histogram ks, $m_v$',hmax=h48max)  
  himg49 = plothisto(mv,/clipedges, name='$m_v$',color='#0000ff',thick=1,overplot=himg48,hmax=h49max)
  himg48.yrange = [0, max(h48max,h49max)]
  l = legend(target=[himg48, himg49],position=[0.87,0.29],sample_width=0.04)  
end

if showwin[9] and calc_Freeman then begin
  ; part 9: plot absolute power, colorbars and histogramm.
  w_data = window(window_title = '4-component Freeman-Durden decomposition: T-Matrix based by Yamaguchi',dimensions=[wscr*0.36,hscr*0.8], location=[wscr*2/3,0])
  tx = text(0.5,0.94,'4-component Freeman-Durden decomposition: !C T-Matrix vs. C -Matrix based algorithm,  by Yamaguchi et al.',alignment=0.5,font_size=12)

  img50data = fltarr([ddim,3])
  img50data[*,*,0] = nlscale(abs(P_fd_dihed),'cliph_W->dB',range=[-40,10])
  img50data[*,*,1] = nlscale(abs(P_fd_vol),'cliph_W->dB',range=[-40,10])
  img50data[*,*,2] = nlscale(abs(P_fd_surf),'cliph_W->dB',range=[-40,10])
  img51data = nlscale(abs(P_fd_helix),'cliph_W->dB',range=[-45,0])
  img50   = image(img50data,pos=imgpos(3,2,1,scale=[1,0.9]),/current, aspect_ratio=0.5,title='$T-based: P_d, P_v, P_s$') 
  img51   = image(img51data,pos=imgpos(3,2,2,scale=[1,0.9]),/current, aspect_ratio=0.5,title='T-based: Helicity',rgb_table=1)
  c = colbar(img51,'dB',range=[-45,0])
;  himg50a = plothisto(img50data[*,*,0],'r2',/clipedges, name='$P_{dihedral}$',position=imgpos(3,2,3,scale=[1,0.665],move=[0,-0.005]),/current,title='histogram P_d, P_v, P_s',hmax=h50amax)  
  himg50a1 = plothisto(nlscale(P_fd_dihed[where(P_fd_dihed gt 0)],'cliph_W->dB'),'r2',/clipedges, name='$P_{dihedral}$',position=imgpos(3,2,3,scale=[1,0.665],move=[0,-0.005]),/current,title='T-based: histogram P_d, P_v, P_s',hmax=h50a1max)
    
  ;himg50b = plothisto(img50data[*,*,1],'g2',/clipedges, name='$P_{volume}$',overplot=himg50a,hmax=h50bmax)
  ;himg50c = plothisto(img50data[*,*,2],'b2',/clipedges, name='$P_{surface}$',overplot=himg50a,hmax=h50cmax)
  ;himg50a = plothisto(nlscale(P_fd_surf[where(P_fd_dihed gt 0)],'cliph_W->dB',range=[-40,0]),'g2',/clipedges, name='$P_{helix(+)}$',overplot=himg50a,hmax=h50a1max)
  himg50a2 = plothisto(nlscale(-P_fd_dihed[where(P_fd_dihed lt 0)],'cliph_W->dB',range=[-40,10]),'--r2',/clipedges, name='$P_{dihed(-)}$',overplot=himg50a1,hmax=h50a2max)  
  himg50b1 = plothisto(nlscale(P_fd_vol[where(P_fd_vol gt 0)],'cliph_W->dB',range=[-40,10]),'g2',/clipedges, name='$P_{vol(+)}$',overplot=himg50a1,hmax=h50b1max)
  himg50b2 = plothisto(nlscale(-P_fd_vol[where(P_fd_vol lt 0)],'cliph_W->dB',range=[-40,10]),'--g2',/clipedges, name='$P_{vol(-)}$',overplot=himg50a1,hmax=h50b2max)
  himg50c1 = plothisto(nlscale(P_fd_surf[where(P_fd_surf gt 0)],'cliph_W->dB',range=[-40,10]),'b2',/clipedges, name='$P_{surf(+)}$',overplot=himg50a1,hmax=h50c1max)
  ;himg50c = plothisto(nlscale(-P_fd_surf[where(P_fd_surf lt 0)],'cliph_W->dB',range=[-40,0]),'--b2',/clipedges, name='$P_{helix(+)}$',overplot=himg50a,hmax=h50c2max)
  himg50d1 = plothisto(nlscale(P_fd_helix[where(P_fd_helix gt 0)],'cliph_W->dB',range=[-40,10]),'k1',/clipedges, name='$P_{helix(+)}$',overplot=himg50a1,hmax=h50d1max)
  himg50d2 = plothisto(nlscale(-P_fd_helix[where(P_fd_helix lt 0)],'cliph_W->dB',range=[-40,10]),'--k1',/clipedges, name='$P_{helix(-)}$',overplot=himg50a1,hmax=h50d2max)
  himg50a1.yrange = [0,max([h50a1max, h50a2max, h50b1max, h50b2max, h50c1max,h50d1max,h50d2max])]
  l = legend(target=[himg50a1, himg50a2, himg50b1, himg50b2, himg50c1],position=[0.68,0.865],sample_width=0.04)
  l = legend(target=[himg50d1, himg50d2],position=[0.87,0.865],sample_width=0.04)

  img52data = fltarr([ddim,3])
  img52data[*,*,0] = nlscale(Pd_fd,'cliph_W->dB',range=[-40,10])
  img52data[*,*,1] = nlscale(Pv_fd,'cliph_W->dB',range=[-40,10])
  img52data[*,*,2] = nlscale(Ps_fd,'cliph_W->dB',range=[-40,10])
  img53data = nlscale(Pc_fd,'cliph_W->dB',range=[-45,0])
  img52   = image(img52data,pos=imgpos(3,2,4,scale=[1,0.9],move=[0,0.15]),/current, aspect_ratio=0.5,title='$C-based: P_d, P_v, P_s$') 
  img53   = image(img53data,pos=imgpos(3,2,5,scale=[1,0.9],move=[0,0.15]),/current, aspect_ratio=0.5,title='C-based: Helicity',rgb_table=1)
  c = colbar(img53,'dB',range=[-45,0])
  himg52a = plothisto(img52data[*,*,0],'r2',/clipedges, name='$P_{dihedral}$',position=imgpos(3,2,6,scale=[1,0.665],move=[0,0.145]),/current,title='C-based: histogram P_d, P_v, P_s',hmax=h52amax)
  himg52b = plothisto(img52data[*,*,1],'g2',/clipedges, name='$P_{volume}$',overplot=himg52a,hmax=h52bmax)
  himg52c = plothisto(img52data[*,*,2],'b2',/clipedges, name='$P_{surface}$',overplot=himg52a,hmax=h52cmax)
  himg52d1 = plothisto(img53data[where(Pc_sign ge 0)],'k1',/clipedges, name='$P_{helix}(+)$',overplot=himg52a,hmax=h52d1max)
  himg52d2 = plothisto(img53data[where(Pc_sign lt 0)],'--k1',/clipedges, name='$P_{helix}(-)$',overplot=himg52a,hmax=h52d2max)
  himg52a.yrange = [0,max([h52amax, h52bmax,h52cmax,h52d1max,h52d2max])]
  l = legend(target=[himg52a, himg52b, himg52c],position=[0.675,0.127],sample_width=0.04)
  l = legend(target=[himg52d1, himg52d2],position=[0.846,0.127],sample_width=0.04)

    
end

