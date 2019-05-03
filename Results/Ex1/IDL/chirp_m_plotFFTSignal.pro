p4a = plot(f/1e6, abs(s1_fft), color='#00ff00', /CURRENT, position=[x0g[0], y0g[1],x1g[0],y1g[1]], name='$|s1(\nu)|$')
p4b = plot(f/1e6, abs(s2_fft), color='#0000ff',/OVERPLOT, name='$|s2(\nu)|$')
p4c = plot(f/1e6, abs(shift(s_fft,N2)),  color='#ff0000',/OVERPLOT, name='$|s(\nu)|$')

p4a.xtitle='Frequency !9n !3 (MHz)'
p4a.ytitle='Amplitude Spectrum !C of Signal (a.u.)'
p4a.XTICKVALUES = [0, 20, 40]
p4a.YTICKVALUES = [0, 500, 1000] 
p4a.xminor = 1
p4a.yminor = 0
p4a.xticklen=0.025
p4a.yticklen=0.025

legend4 = legend(target=[p4a, p4b, p4c],position=[lxp[0],lyp[1]-0.12], font_size=8, sample_width=0.03, HORIZONTAL_ALIGNMENT='right', vertical_spacing=0.005, shadow=0)
p4a.xrange=[-20, 40]