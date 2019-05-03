p5a = plot(f/1e6, abs(shift(h_fft,round(N/2))), color='#0000ff', Linestyle=0, xrange=[0, 40],/CURRENT, position=[x0g[1], y0g[1],x1g[1],y1g[1]])
p5a.xtitle='Frequency !9n !3 (MHz)'
p5a.ytitle='Amplitude Spectrum !C of Filter Function (a.u.)'
p5a.XTICKVALUES = [0, 20, 40]
p5a.YTICKVALUES = [0, 500, 1000] 
p5a.xminor = 1
p5a.yminor = 0
p5a.xticklen=0.025
p5a.yticklen=0.025
