yran = [min(phase_filter[frange]), max(phase_filter[frange])]-phase_filter[f0I]
 p8a = plot(f/1e6, phase_filter-phase_filter[f0I], color='#0000ff', Linestyle=0, xrange=[0, 50], yrange = yran, $
           /CURRENT, position=[x0g[1], y0g[2],x1g[1],y1g[2]],name='$\phi(\nu)$')
p8a.xtitle='Frequency !9n !3 (MHz)'
p8a.ytitle='Phase !9f!C !3 of Filter Function (rad)'           
p8a.XTICKVALUES = [0, 20, 40]
p8a.xminor = 1
p8a.yminor = 0
p8a.xticklen=0.025
p8a.yticklen=0.025
