@colorbar_define

pro pol3
;COMPUTE COVARIANCE MATRIX(USE E.G. WINSZ=7X7)
;Plot:
;             phase(HH VV*),phase(HH HV*), phase(VV HV*)
;                     (using TV,color plot, WITH colorbar)
;             coherences (HH VV, HH HV, VV HV)            
;                   (using TV,black and white scaled 0 to 1 WITH colorbar)
;             histograms: phases and coherences (2 plots)

device,retain=2,decomposed=0

OPENR, uh, '/home/pisc_io/Documents/i_al_af_1206_hh_corr.dat', /xdr, /get_lun; 
dim = lonarr(2);  
READU, uh, dim;  
hh = complexarr(dim);  
READU, uh, hh;  
OPENR, uv, '/home/pisc_io/Documents/i_al_af_1206_vv_corr.dat', /xdr, /get_lun; 
vv = complexarr(dim);  
READU, uv, vv;  
OPENR, ux, '/home/pisc_io/Documents/i_al_af_1206_xx_corr.dat', /xdr, /get_lun; 
xx = complexarr(dim);  
READU, ux, xx;  

;LEXICOGRAPHIC SCATTERING VECTOR
im=complexarr(dim(0),dim(1),3)
im[*,*,0]=hh
im[*,*,1]=sqrt(2)*xx
im[*,*,2]=vv

;COVARIANCE MATRIX
cov=complexarr(dim(0),dim(1),3,3)
v=complexarr(3)
for i=0,dim(0)-1 do begin
    for j=0,dim(1)-1 do begin
    	v[0]=im[i,j,0]
	v[1]=im[i,j,1]
	v[2]=im[i,j,2]
	cov(i,j,*,*)=(v)#transpose(conj(v))
    end
end
cov_s=smooth(cov,[7,7,1,1]);expected value

;phase(HH VV*),phase(HH HV*), phase(VV HV*)
hh_hv=atan(cov_s(*,*,0,1),/phase)
hh_vv=atan(cov_s(*,*,0,2),/phase)
hv_vv=atan(cov_s(*,*,1,2),/phase)
;window,1,title= "phases HH HV*" 
;tv,bytscl(hh_hv,-!pi,!pi)
;window,2,title= "phases HH VV*" 
;tv,bytscl(hh_vv,-!pi,!pi)
;window,3,title= "phases HV VV*" 
;tv,bytscl(hv_vv,-!pi,!pi)

loadct,39
device, decompose=0
WINDOW,4,xsize=200+320,ysize=500,title='phases HH HV*'
tv,congrid(bytscl(hh_hv,-!pi,!pi),200,500)
colorbar = Obj_New("COLORBAR", Format='(F8.1)',/vertical)
colorbar->Draw
colorbar->SetProperty, Range=[-!pi,!pi], /Erase, /Draw
Obj_Destroy, colorbar

WINDOW,5,xsize=200+320,ysize=500,title='phases HH VV*'
tv,congrid(bytscl(hh_vv,-!pi,!pi),200,500)
colorbar = Obj_New("COLORBAR", Format='(F8.1)',/vertical)
colorbar->Draw
colorbar->SetProperty, Range=[-!pi,!pi], /Erase, /Draw
Obj_Destroy, colorbar

WINDOW,6,xsize=200+320,ysize=500,title='phases HV VV*'
tv,congrid(bytscl(hv_vv,-!pi,!pi),200,500)
colorbar = Obj_New("COLORBAR", Format='(F8.1)',/vertical)
colorbar->Draw
colorbar->SetProperty, Range=[-!pi,!pi], /Erase, /Draw
Obj_Destroy, colorbar
;!p.background=255
;device, decompose=1
;loadct,0

;Plot coherences (HH VV, HH HV, VV HV)            
;  (using TV,black and white scaled 0 to 1 WITH colorbar)

a1=complexarr(dim(0),dim(1))
for i=0,dim(0)-1 do begin
    for j=0,dim(1)-1 do begin
    a1(i,j)=(hh(i,j))#conj(xx(i,j))
    end
end
a1=smooth(a1,[7,7])
b1=smooth((abs(hh)^2),[7,7])
c1=smooth((abs(xx)^2),[7,7])
g1=a1/sqrt(b1*c1)

a2=complexarr(dim(0),dim(1))
for i=0,dim(0)-1 do begin
    for j=0,dim(1)-1 do begin
 a2(i,j)=(hh(i,j))#conj(vv(i,j))
    end
end
a2=smooth(a2,[7,7])
c2=smooth((abs(vv)^2),[7,7])
g2=a2/sqrt(b1*c2)

a3=complexarr(dim(0),dim(1))
for i=0,dim(0)-1 do begin
    for j=0,dim(1)-1 do begin
    a3(i,j)=(vv(i,j))#conj(xx(i,j))
    end
end
a3=smooth(a3,[7,7])
g3=a3/sqrt(c2*c1)

loadct,0
;device, decompose=0
WINDOW,7,xsize=dim(0)/3.+320,ysize=dim(1)/3.,title='coherence HH HV'
tv,congrid(bytscl(abs(g1),0,1),dim(0)/3.,dim(1)/3.)
colorbar = Obj_New("COLORBAR", Format='(F8.1)',/vertical)
colorbar->Draw
colorbar->SetProperty, Range=[0,1], /Erase, /Draw
Obj_Destroy, colorbar

WINDOW,8,xsize=dim(0)/3.+320,ysize=dim(1)/3.,title='coherence HH VV'
tv,congrid(bytscl(abs(g2),0,1),dim(0)/3.,dim(1)/3.)
colorbar = Obj_New("COLORBAR", Format='(F8.1)',/vertical)
colorbar->Draw
colorbar->SetProperty, Range=[0,1], /Erase, /Draw
Obj_Destroy, colorbar

WINDOW,9,xsize=dim(0)/3.+320,ysize=dim(1)/3.,title='coherence VV HV'
tv,congrid(bytscl(abs(g3),0,1),dim(0)/3.,dim(1)/3.)
colorbar = Obj_New("COLORBAR", Format='(F8.1)',/vertical)
colorbar->Draw
colorbar->SetProperty, Range=[0,1], /Erase, /Draw
Obj_Destroy, colorbar

!p.background=255
device, decompose=1

;histograms: phases and coherences (2 plots);
loadct,0
window,10,title= "histogram of phases(HH HV*)(HH VV*)(HV VV*)"
plot,findgen((2*!pi+1)*10)/10.-!pi,histogram(hh_hv,binsize=0.1,max=!pi,min=-!pi),yrange=[0, 500000],LINESTYLE=0
oplot,findgen((2*!pi+1)*10)/10.-!pi,histogram(hh_vv,binsize=0.1,max=!pi,min=-!pi),LINESTYLE=1
oplot,findgen((2*!pi+1)*10)/10.-!pi,histogram(hv_vv,binsize=0.1,max=!pi,min=-!pi),LINESTYLE=2


window,11,title= "histogram of coherence(HH VV)(HH HV)(VV HV)"
plot,findgen(1*100)/100.,histogram(abs(g1),binsize=0.01,max=1,min=-0),LINESTYLE=0
oplot,findgen(1*100)/100.,histogram(abs(g2),binsize=0.01,max=1,min=-0),LINESTYLE=1
oplot,findgen(1*100)/100.,histogram(abs(g3),binsize=0.01,max=1,min=-0),LINESTYLE=5

stop
CLOSE, uh
 free_lun, uh
CLOSE, uv
 free_lun, uv
CLOSE, ux
 free_lun, ux


end

