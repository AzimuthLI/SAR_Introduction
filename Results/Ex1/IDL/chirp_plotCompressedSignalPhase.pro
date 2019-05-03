yran = [min(phase_compressed[where(f ge 0 and f le 50e6)]), max(phase_compressed[where(f ge 0 and f le 50e6)])]-phase_compressed[f0I]
p9a = plot(f/1e6, phase_compressed-phase_compressed[f0I], color='#bb00bb', Linestyle=0, xrange=[0, 50], yrange = yran, $
           /CURRENT, position=[x0g[2], y0g[2],x1g[2],y1g[2]], name='$\phi_c(\nu)$')           
p9a.xtitle = 'Frequency !9n !3 (MHz)'
p9a.ytitle = 'Phase !9f!C !3 of Compressed Signal !C !5sc!3(!5t!3) (rad)'
p9a.XTICKVALUES = [0, 20, 40]
;p9a.YTICKVALUES = [-!Pi/4, 0, !Pi/4] 
;p9a.YTICKname = ['-1/4!9p', '0', '1/4!9p'] 
p9a.xminor = 1
p9a.yminor = 0
p9a.xticklen=0.025
p9a.yticklen=0.025
