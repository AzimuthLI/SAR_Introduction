; get screen dimensions
scrsize = get_screen_size()
wscr = scrsize[0]
hscr = scrsize[1]

; window for raw data
w_data = window(window_title = 'Raw data (5px smoothed)',dimensions=[wscr*1/3,hscr*0.8], location=[wscr*2/3,0])
t = text(0.5,0.97,'Raw data (5px smoothed)',alignment=0.5,font_size=12)
w_dataID1 = igetcurrent()

imdata1 = smooth(nlscale(abs(data1),'clipihisto'),5,/edge_truncate)
imdata2 = smooth(nlscale(abs(data2),'clipihisto'),5,/edge_truncate)
imdata3 = fltarr(3,ddim[0],ddim[1])
imdata3[0,*,*] = (imdata1-min(imdata1))/(max(imdata1)-min(imdata1))
imdata3[1,*,*] = (imdata2-min(imdata2))/(max(imdata2)-min(imdata2))
imdata3[2,*,*] = (imdata3[0,*,*]+imdata3[1,*,*])/2

; plot images and superposition
img1 = image(imdata1,position=imgpos(3,2,1),/current, Aspect_ratio=0.5,title='Image 1')
img2 = image(imdata1,position=imgpos(3,2,2),/current, Aspect_ratio=0.5,title='Image 2')
img3 = image(imdata3,position=imgpos(3,2,3),/current, Aspect_ratio=0.5,title='Superposition 1+2')

; plot correlation function
corr_shift = nlscale(shift(abs(corr),ddim[0]/2,ddim[1]/2),'minmax')
xpx = indgen(ddim(0))-ddim(0)/2
ypx = indgen(ddim(1))-ddim(1)/2
img4 = image(corr_shift,xpx, ypx,position=imgpos(3,2,4),/current, Aspect_ratio=0.5,rgb_table=5,title = '2D correlation function',axis_style=2)
p4a = polyline([0,0],[-50,+50],color='white',thick=1,target=img4,/DATA)
p4b = polyline([0-20,+20],[0,0],color='white',thick=1,target=img4,/DATA)
p4c = plot([arrayindices[0]],[arrayindices[1]],symbol='o',color='white',/CURRENT,/OVERPLOT,sym_size=0.5)

; plot superposition of matched images
im1 = smooth(nlscale(abs(image1),'clipihisto'),5,/edge_truncate)
im2 = smooth(nlscale(abs(image2),'clipihisto'),5,/edge_truncate)
im3 = fltarr(3,ddim(0)-abs(rg_shift),ddim(1)-abs(az_shift))
im3[0,*,*] = (im1-min(im1))/(max(im1)-min(im1))
im3[1,*,*] = (im2-min(im2))/(max(im2)-min(im2))
im3[2,*,*] = (im3[0,*,*]+im3[1,*,*])/2
img6 = image(im3,position=imgpos(3,2,6),/current, Aspect_ratio=0.5,title='Superposition 1+2, shifted')

; Open window for interferograms
w_data = window(window_title = 'Interferogram',dimensions=[wscr*1/3,hscr*0.8], location=[wscr*2/3,0])
t = text(0.5,0.98,'Interferogram (5px smoothed)',alignment=0.5,font_size=10)
w_dataID2 = igetcurrent()

; plot interferogram without flat earth removal
imdata7 = smooth(nlscale(abs(Interf),'log'),5,/edge_truncate)
imdata8 = smooth(nlscale(atan(Interf,/phase),'clipihisto'),5,/edge_truncate)
img7 = image(imdata7,position=imgpos(3,2,1),/current, Aspect_ratio=0.5,title='Intensity of interferogram')
img8 = image(imdata8,position=imgpos(3,2,4),/current, Aspect_ratio=0.5,title='Phase of interferogram')

; plot interferogram with flat earth removed
imdata9  = smooth(nlscale(abs(Interf_flat1),'log'),5,/edge_truncate)
imdata10 = smooth(nlscale(atan(Interf_flat1,/phase),'clipihisto'),5,/edge_truncate)
img9     = image(imdata9,position=imgpos(3,2,2),/current, Aspect_ratio=0.5,title='Intensity of !C flat interferogram')
img10    = image(imdata10,position=imgpos(3,2,5),/current, Aspect_ratio=0.5,title='Phase of !C flat interferogram')


; Open window for spectra, etc.
w_data = window(window_title = 'Interferogram spectra',dimensions=[wscr*0.36,hscr*0.8], location=[wscr*2/3,0])
w_dataID3 = igetcurrent()

; plot spectra of shifted images (from flat earth removal)
p11 = plot(xpixels,rg_spec1,position=imgpos(3,4,3,scale=[0.8,0.75],move=[0.1,0.1]),/current, title='Spectra of shifted !C images',name = 'image 1')
p12 = plot(xpixels,rg_spec2, /overplot, color='red', name = 'image 2')
l11 = legend(target=[p11, p12])

; plot spectra of range filtered images
p11a = plot(xpixels,shift(total(abs(fft(img1filtered )),2), idim[0]/2),position=imgpos(3,4,6,scale=[0.8,0.75],move=[0.1,0]),/current, title='Spectra of shifted !C images',name = 'image 1')
p12a = plot(xpixels,shift(total(abs(fft(img2_flat_fil)),2), idim[0]/2), /overplot, color='red', name = 'image 2')

; plot spectrum of interferogram
imdata11 = smooth(nlscale(Interf_fft,'log'),5,/edge_truncate)
img11     = image(imdata11,position=imgpos(3,2,1),/current, Aspect_ratio=0.5,title='Spectrum of !C interferogram',rgb_table=coltab())

; plot azimut integrated spectrum of interferogram
imdata12 = total(imdata11,2)
p13      = plot(xpixels,imdata12,position=imgpos(3,4,2,scale=[0.75,0.8],move=[0.02,0.1]),/current, title='Spectrum of interferogram',name = 'R(k)')

; plot (linear) phase of spectrum of interferogram and show derrivation of linearity
imdata13 = fltarr(idim)
imdata14 = fltarr(idim)
for j=0,idim[1]-1 do imdata13[*,j] = phunwrap(atan(Interf_fft[*,j],/phase))
mean_phase = total(imdata13,2)/idim[1]
for j=0,idim[1]-1 do imdata14[*,j] = imdata13[*,j] - mean_phase
img13    = image(imdata13,position=imgpos(3,2,4,scale=[1,0.9],move=[0,-0.05]),/current, Aspect_ratio=0.5,title='Phase of !C interferogram',rgb_table=coltab())
p14      = plot(xpixels,mean_phase,position=imgpos(3,4,8,scale=[0.7,0.7],move=[0.1,0]),/current, title='Avg rg_phase of Interferogram',name = 'R(k)')
p15      = plot(xpixels,total(imdata14,2)/idim[1],position=imgpos(3,4,11,scale=[0.7,0.7],move=[0.1,0]),/current, title='Derrivation from !C mean rg_phase',name = 'R(k)')

; plot spectra of raw images and show filters for range filtering
p16 = plot(xpixels,rg_spec1,position=imgpos(3,4,5,scale=[0.9,0.75],move=[0.05,0.0]),/current, title='Frequency Spectrum $I_1, I_2$ !C and range filters',name = 'image 1')
p17 = plot(xpixels,rg_spec2a, /overplot, color='red', name = 'image 2')

signalbw1 = fwhm(xpixels,rg_spec1,POS = fwhm1x, YVAL = fwhm1y)
signalbw2 = fwhm(xpixels,rg_spec1,POS = fwhm2x, YVAL = fwhm2y)
p4a = polyline([fwhm1x[0],fwhm1x[0]],[0,fwhm1y[0]],color='blue',thick=1,target=p16,/DATA)
p4a = polyline([fwhm1x[1],fwhm1x[1]],[0,fwhm1y[1]],color='blue',thick=1,target=p16,/DATA)
p4a = polyline([fwhm2x[0],fwhm2x[0]],[0,fwhm2y[0]],color='blue',thick=1,target=p16,/DATA)
p4a = polyline([fwhm2x[1],fwhm2x[1]],[0,fwhm2y[1]],color='blue',thick=1,target=p16,/DATA)
print, 'signalbw1 = ' + str(signalbw1)
print, 'signalbw2 = ' + str(signalbw2)

p18 = plot(xpixels,filter1*fwhm1y[0]*0.5,/overplot,color='black')
p19 = plot(xpixels,filter2*fwhm1y[0]*0.6,/overplot,color='red')

; goback to window 2 and add range filtered interferogram
isetcurrent, w_dataID2
img15data = smooth(nlscale(atan(Interf_flat_fil,/phase),'clipihisto'),5,/edge_truncate)
img15    = image(img15data,position=imgpos(3,2,6),/current, Aspect_ratio=0.5,title='Phase of flat & range-filtered !C interferogram')
img16data = smooth(nlscale(Interf_flat_fil,'log'),5,/edge_truncate)
img16    = image(img16data,position=imgpos(3,2,3),/current, Aspect_ratio=0.5,title='Phase of flat & range-filtered !C interferogram')

; Open window for coherence plots
w_data = window(window_title = 'Coherence',dimensions=[wscr*0.36,hscr*0.8], location=[wscr*2/3,0])
w_dataID4 = igetcurrent()

;plot coherences of shifted and range filtered images
img17 = image(smooth(coherence,5,/edge_truncate),rgb_table=1,position=imgpos(2,1,1,scale=[1,0.57],move=[0,0.2]),/current, Aspect_ratio=0.5,title='Coherence shifted images') 
img18 = image(smooth(coherencef,5,/edge_truncate),rgb_table=1,position=imgpos(2,1,2,scale=[1,0.57],move=[0,0.2]),/current, Aspect_ratio=0.5,title='Coherence filtered images')
p17 = plot(coh_histo,position=imgpos(1,3,3,scale=[0.65,1],move=[-0.14,0.1]),/current, title='Coherence histogram',name = 'shifted images')
p18 = plot(cohf_histo, name = 'shifted & filtered images',/overplot, color='red')

; plot coherence function and normalized difference of coherences of the two images 
img19data = (smooth(coherencef,10,/edge_truncate)-smooth(coherence,10,/edge_truncate))/smooth(coherencef,10,/edge_truncate)
img = image(img19data,position=imgpos(1,3,3,scale=[1,1],move=[0.35,0.1]),rgb_table=1,aspect_ratio=0.5,/current,title='Normalized difference !C of coherences')