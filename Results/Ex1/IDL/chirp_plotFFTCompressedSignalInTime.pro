p3a.select
p3d = plot(t*1e6, abs(shift(cs_ifft,round(N/2))), color='#00bb00', yrange=[min(real(cs)),max(real(cs))],linestyle=0, thick=0.33,/OVERPLOT)
p3d.name = '|cs(t), fft|'
p3a.Yrange = [-2000, 5000]
p3a.xrange = [-0.25, 0.25]
legend3.add, p3d
