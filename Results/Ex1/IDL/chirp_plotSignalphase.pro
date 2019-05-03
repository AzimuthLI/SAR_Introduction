yran = [min(phase_signal[where(f ge 0 and f le 50e6)]), max(phase_signal[where(f ge 0 and f le 50e6)])]-phase_signal[f0I]
p7a = plot(f/1e6,phase_signal-phase_signal[f0I], color='#ff0000', Linestyle=0, xrange=[0, 50], yrange = yran, $
           /CURRENT, position=[x0g[0], y0g[2],x1g[0],y1g[2]])
p7a.xtitle='Frequency !9n !3 (MHz)'
p7a.ytitle='Phase !9f!C !3 of Signal (rad)'
p7a.XTICKVALUES = [0, 20, 40]
p7a.xminor = 1
p7a.yminor = 0
p7a.xticklen=0.025
p7a.yticklen=0.025
