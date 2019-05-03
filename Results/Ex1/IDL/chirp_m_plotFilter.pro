p2a = plot(t*1e6 , real(h1), color='#00ff00', thick=1,/CURRENT, position=[x0g[1], y0g[0],x1g[1],y1g[0]], name='real $h_1(t)$')
p2b = plot(t*1e6 , real(h2), color='#0000ff', linestyle=0,/OVERPLOT, name='real $h_2(t)$')
p2a.xtitle='time !5 t !3(!9m!3s)'
p2a.ytitle='Amplitude !C of Filter Function !5h!3(!5t!4) (a.u.)' 
p2a.XTICKVALUES = [-1,0, 1]
p2a.YTICKVALUES = [-1, 0, 1] 
p2a.yminor = 0
p2a.xminor = 4
p2a.yrange = [-1.2, 1.5]
p2a.xticklen=0.025
p2a.yticklen=0.025
legend2 = legend(target=[p2a, p2b],position=[lxp[1],lyp[0]-0.19], font_size=8, sample_width=0.03, HORIZONTAL_ALIGNMENT='right', vertical_spacing=0.005, shadow=0)