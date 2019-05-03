p6a = plot(f/1e6, abs(shift(cs_fft,round(N/2))), color='#bb00bb', Linestyle=0,/CURRENT, position=[x0g[2], y0g[1],x1g[2],y1g[1]], name = '$\itcs\rm(\it\nu\rm)$')
p6b = plot(f/1e6, abs(cs1_fft), color='#00ff00', name='$cs_1(\nu)$',/OVERPLOT)
p6c = plot(f/1e6, abs(cs2_fft), color='#0000ff', name='$cs_2(\nu)$',/OVERPLOT)
p6a.xtitle='Frequency !9n !3 (MHz)'
p6a.ytitle='Amplitude Spectrum !C of Compressed Signal (a.u.)'
p6a.XTICKVALUES = [0, 20, 40]
p6a.YTICKVALUES = [0, 5e5, 1e6] 
p6a.xminor = 1
p6a.yminor = 0
p6a.xticklen=0.025
p6a.yticklen=0.025
legend6 = legend(target=[p6a, p6b,p6c],position=[lxp[2],lyp[1]-0.12], font_size=8, sample_width=0.03, HORIZONTAL_ALIGNMENT='right', vertical_spacing=0.005, shadow=0)
p6a.xrange=[-20, 40]