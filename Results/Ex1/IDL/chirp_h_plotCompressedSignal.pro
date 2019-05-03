; plot compressed signal
p3a.select

; make old plot 90% transparent
p3a.transparency = 90
p3b.transparency = 90
p3c.transparency = 90

; add compressed signal with hamming window filter
p3e = plot(t*1e6, abs(shift(csh,round(N/2))), color='#ff8800', thick=1,/OVERPLOT, name='|csh(t)|')
p3a.xrange = [-0.25,0.25]
legend3.add, p3e
