p6a = plot(f/1e6, abs(shift(cs_fft,round(N/2))), color='#bb00bb', Linestyle=0, xrange=[0, 40],/CURRENT, position=[x0g[2], y0g[1],x1g[2],y1g[1]], name = '$\itcs\rm(\it\nu\rm)$')
p6a.xtitle='Frequency !9n !3 (MHz)'
p6a.ytitle='Amplitude Spectrum !C of Compressed Signal (a.u.)'
p6a.XTICKVALUES = [0, 20, 40]
p6a.YTICKVALUES = [0, 5e5, 1e6] 
p6a.xminor = 1
p6a.yminor = 0
p6a.xticklen=0.025
p6a.yticklen=0.025

