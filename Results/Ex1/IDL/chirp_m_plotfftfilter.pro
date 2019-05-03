; plot frequency spectrum of signal 
p5a = plot(f/1e6, abs(h1_fft), color='#00ff00', Linestyle=0,/CURRENT, position=[x0g[1], y0g[1],x1g[1],y1g[1]], name='filter $\ith\rm(\nu)$')
p5b = plot(f/1e6, abs(h2_fft), color='#0000ff', name='filter $\ith_1\rm(\nu)$',/OVERPLOT)
p5c = plot(f/1e6, abs(shift(h_fft,N2)),  color='#ff0000', name='filter $\ith_1\rm(\nu)$',/OVERPLOT)

; add arrow and text
fw = fwhm(f/1e6,abs(h1_fft),Pos=l5x,yval=l5y)
t5 = text(mean(l5x),l5y[1],'$bw=$'+string(format='(%"%4.1f")',fw)+'$ MHz$',/DATA,alignment=0.5, font_size=10,target=p5a)
l5 = marrow(l5x,[l5y[1],l5y[1]],width=1,height=40,/DATA,target=p5a)

p5a.xtitle='Frequency !9n !3 (MHz)'
p5a.ytitle='Amplitude Spectrum !C of Filter Function (a.u.)'
p5a.XTICKVALUES = [0, 20, 40]
p5a.YTICKVALUES = [0, 500, 1000] 
p5a.xminor = 1
p5a.yminor = 0
p5a.xticklen=0.025
p5a.yticklen=0.025
legend5 = legend(target=[p5a, p5b, p5c],position=[lxp[1],lyp[1]-0.12], font_size=8, sample_width=0.03, HORIZONTAL_ALIGNMENT='right', vertical_spacing=0.005, shadow=0)
p5a.xrange = [-20, 40]