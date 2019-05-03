p4a = plot(f/1e6, abs(shift(s_fft,round(N/2))), xrange=[0, 40], color='#ff0000', /CURRENT, position=[x0g[0], y0g[1],x1g[0],y1g[1]])
p4a.xtitle='Frequency !9n !3 (MHz)'
p4a.ytitle='Amplitude Spectrum !C of Signal (a.u.)'
p4a.XTICKVALUES = [0, 20, 40]
p4a.YTICKVALUES = [0, 500, 1000] 
p4a.xminor = 1
p4a.yminor = 0
p4a.xticklen=0.025
p4a.yticklen=0.025

fw  = fwhm(f/1e6,abs(shift(s_fft,round(N/2))),Pos=l4x,yval=l4y)
t4 = text(mean(l4x),l4y[0]*1.05,'$bw = $'+string(format='(%"%4.1f")',fw)+'$ MHz$',/DATA,alignment=0.5, font_size=9,target=p4a)
l4 = marrow(l4x,l4y,width=1,height=40,/DATA,target=p4a)
p4a.xrange = [0, 40]
