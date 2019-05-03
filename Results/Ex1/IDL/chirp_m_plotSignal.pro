p1a = plot(t*1e6, real(s1), color='#00ff00', location=[0,0], $
      dimensions=[1000, 700], position=[x0g[0], y0g[0],x1g[0],y1g[0]],name='real $s_1(t)$')
p1b = plot(t*1e6, real(s2), color='#0000ff', linestyle=0,/OVERPLOT,name='real $s_2(t)$ ')
p1a.xtitle = 'time !5t!3 (!9m!3s)'
p1a.ytitle = 'Amplitude of Signal !5s!3(!5t!3) !C (a.u.)'
p1a.XTICKVALUES = [-1,0, 1]
p1a.YTICKVALUES = [-1,0,   1]
p1a.xticklen=0.025
p1a.yticklen=0.025
p1a.xminor = 1
p1a.yminor = 0
legend1 = legend(target=[p1a, p1b],position=[lxp[0],lyp[0]-0.19], font_size=8, sample_width=0.03, HORIZONTAL_ALIGNMENT='right', vertical_spacing=0.005, shadow=0)
p1a.yrange = [-1.2,1.5]