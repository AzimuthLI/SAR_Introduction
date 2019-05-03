pro example_histogram

x=[2,5,7,3,4,11,14,2,2,6,7,7,8,9,11,4,4,4,4,4,5,5,5]

hist=histogram(x,MAX=10,MIN=3)

window,1
plot, findgen(10-3+1)+3,hist

y=[2,5.2,7.1,3.6,4.4,11,14.1,2,2,6.2,7,7,8,9,11,4,4,4,4,4,5,5,5]

hist2=histogram(y,BINSIZE=0.1,MAX=10,MIN=3)

window,2
plot, findgen(10-3+1)+3,hist2

stop

end