@colorbar_define
pro pol6
; Compute H, A and alpha (total) and alpha1 and visualize (use logical min/max values, e.g. 0 < H < 1). 
;         Plot histograms of each. 
;         Plot 2-D histogram (alpha vs. H as in Cloude97 paper but with
;              histograms, use IDL function hist2d).    
;         Overplot theoretical bounds (see Cloude97) on H and alpha and add 
;                  lines to show classification

OPENR, uh, '../i_al_af_1206_hh_corr.dat', /xdr, /get_lun
dim = lonarr(2)
READU, uh, dim
hh = complexarr(dim)
READU, uh, hh
OPENR, uv, '../i_al_af_1206_vv_corr.dat', /xdr, /get_lun
vv = complexarr(dim)
READU, uv, vv
OPENR, ux, '../i_al_af_1206_xx_corr.dat', /xdr, /get_lun
xx = complexarr(dim)
READU, ux, xx

;PAULI SCATTERING VECTOR
im=complexarr(3,dim(0),dim(1))
im[0,*,*]=(hh+vv)
im[1,*,*]=(hh-vv)
im[2,*,*]=(2.*xx)
im=im/(sqrt(2))
;coherency matrix.
t=complexarr(3,3,dim(0),dim(1))
v=complexarr(3)
for i=0,dim(0)-1 do begin
    for j=0,dim(1)-1 do begin
    	v[0]=im[0,i,j]
	v[1]=im[1,i,j]
	v[2]=im[2,i,j]
	t(*,*,i,j)=transpose(v)##(conj(v))
    end
end
t=smooth(t,[1,1,7,7])

;Compute eigenvalues and eigenvectors using
;built-in IDL function (LA_EIGENQL). 
eigen_value=complexarr(3,dim(0),dim(1))
eigen_vectors=complexarr(3,3,dim(0),dim(1))
for i=0,dim(0)-1 do begin;10 do begin;
    for j=0,dim(1)-1 do begin;10 do begin;
;    print,i,j,' eigenvalues and eigenvectors'
       eigen_value[*,i,j]=la_eigenql(t[*,*,i,j],eigenvectors=vectors)
       eigen_vectors[*,*,i,j]=vectors
endfor
endfor

; Compute H, A and alpha (total) and alpha1 and visualize (use logical min/max values, e.g. 0 < H < 1). ]
p=fltarr(3,dim(0),dim(1))
for i=0,dim(0)-1 do begin;10 do begin;
    for j=0,dim(1)-1 do begin;10 do begin;
;     print,i,j," p"
    p[*,i,j]=eigen_value[*,i,j]/total(eigen_value[*,i,j]); N.B. e_v[0]<e_v[1]<e_v[2]
     endfor
endfor
;Entropy H
;h=fltarr(1,dim(0),dim(1))
h=fltarr(dim(0),dim(1))
for i=0,dim(0)-1 do begin;10 do begin;
    for j=0,dim(1)-1 do begin;10 do begin;
       for k=0,2 do begin
;       print,i,j," h"
          ;h[*,i,j]=h[*,i,j]-p[k,i,j]*(alog10(p[k,i,j])/alog10(3))
          h[i,j]=h[i,j]-p[k,i,j]*(alog10(p[k,i,j])/alog10(3))
        endfor
     endfor
endfor
;Anisotropy a
a=fltarr(dim(0),dim(1))
for i=0,dim(0)-1 do begin;10 do begin;
    for j=0,dim(1)-1 do begin;10 do begin;
;     print,i,j," a"
       a[i,j]=(eigen_value[1,i,j]-eigen_value[0,i,j])/(eigen_value[1,i,j]+eigen_value[0,i,j]); N.B. e_v[0]<e_v[1]<e_v[2]
        ;a=(p[1,i,j]-p[0,i,j])/(p[1,i,j]+p[0,i,j])
    endfor
endfor

;alpha1
alpha_00=fltarr(dim(0),dim(1))
alpha_10=fltarr(dim(0),dim(1))
alpha_20=fltarr(dim(0),dim(1))
for i=0,dim(0)-1 do begin;10 do begin;
    for j=0,dim(1)-1 do begin;10 do begin;
;     print,i,j," alpha123"
alpha_00[i,j]=acos(abs(eigen_vectors[0,0,i,j])); N.B. e_v[0]<e_v[1]<e_v[2]
alpha_10[i,j]=acos(abs(eigen_vectors[1,0,i,j]))
alpha_20[i,j]=acos(abs(eigen_vectors[2,0,i,j]))
    endfor
endfor

;alpha total
alpha=fltarr(dim(0),dim(1))
for i=0,dim(0)-1 do begin;10 do begin;
    for j=0,dim(1)-1 do begin;10 do begin;
;     print,i,j," alpha"
alpha[i,j]=p[0,i,j]*alpha_00[i,j]+p[1,i,j]*alpha_10[i,j]+p[2,i,j]*alpha_20[i,j]
    endfor
endfor
alpha=alpha*!RADEG
alpha1=alpha_20*!RADEG; N.B. e_v[0]<e_v[1]<e_v[2]


loadct,0
device, decompose=0
WINDOW,1,xsize=dim(0)/3.+320,ysize=dim(1)/3.,title='entropy H'
tv,congrid(bytscl(abs(h),0,1),dim(0)/3.,dim(1)/3.)
colorbar = Obj_New("COLORBAR", Format='(F8.1)',/vertical)
colorbar->Draw
colorbar->SetProperty, Range=[0,1], /Erase, /Draw
Obj_Destroy, colorbar

WINDOW,2,xsize=dim(0)/3.+320,ysize=dim(1)/3.,title='anisotropy A'
tv,congrid(bytscl(abs(a),0,1),dim(0)/3.,dim(1)/3.)
colorbar = Obj_New("COLORBAR", Format='(F8.1)',/vertical)
colorbar->Draw
colorbar->SetProperty, Range=[0,1], /Erase, /Draw
Obj_Destroy, colorbar

WINDOW,3,xsize=dim(0)/3.+320,ysize=dim(1)/3.,title='alpha1 angle'
tv,congrid(bytscl(alpha1,0,90),dim(0)/3.,dim(1)/3.)
colorbar = Obj_New("COLORBAR", Format='(F8.1)',/vertical)
colorbar->Draw
colorbar->SetProperty, Range=[0,90], /Erase, /Draw
Obj_Destroy, colorbar

WINDOW,4,xsize=dim(0)/3.+320,ysize=dim(1)/3.,title='alpha angle'
tv,congrid(bytscl(alpha,0,90),dim(0)/3.,dim(1)/3.)
colorbar = Obj_New("COLORBAR", Format='(F8.1)',/vertical)
colorbar->Draw
colorbar->SetProperty, Range=[0,90], /Erase, /Draw
Obj_Destroy, colorbar


window,11,title= "histogram of entropy"
plot,findgen((1)*100)/100.,histogram(h,binsize=0.01,max=1,min=0)
window,12,title= "histogram of anysotropy"
plot,findgen((1)*100)/100.,histogram(a,binsize=0.01,max=1,min=0)
window,13,title= "histogram of alpha1"
plot,findgen(91),histogram(alpha1,binsize=1,max=90,min=0)
window,14,title= "histogram of alpha"
plot,findgen(91),histogram(alpha,binsize=1,max=90,min=0)

;         Plot 2-D histogram (alpha vs. H as in Cloude97 paper but with
;              histograms, use IDL function hist2d).    
;         Overplot theoretical bounds (see Cloude97) on H and alpha and add 
;                  lines to show classification


 Result = HIST_2D(H, alpha , BIN1=0.001, BIN2=0.001, MIN1=0, MAX1=1,MIN2=0, MAX2=90)
 window,15,xsize=1000,ysize=1000,title='alpha-H distribution'
 tvscl,congrid(result,1000,1000)
  axis,xrange=[0,1],yrange=[0,90],XTITLE = 'Entropy',  YTITLE = 'Alpha'
stop 

  e=1.e-5
  m=findgen(101)/100+e ;0<=m<=1
  H_1=-(alog10(1/(1+2*m))+2*m*alog10(m(1+2*m)))/(alog10(3)+m*alog10(9))
  alpha_1=m*!pi/(1+2*m)
  m=findgen(51)/100+.5+e ;0<=m<=0.5 & 0.5<=m<=1
  H_2=(-2*alog10(1/(1+2*m))+(1-2*m)*alog10((2*m-1)/(1+2*m)))/(alog10(3)+m*alog10(9))
  alpha_2=!pi/(1+2*m)
  oplot,H_1,alpha_1/!PI*180,/t3d,/noclip,color=color,thick=thick
  oplot,H_2,alpha_2/!PI*180,/t3d,/noclip, color=color,thick=thick  
  oplot,findgen(64)/100,make_array(1,64,value=90),/noclip,color=color,thick=thick


!p.background=255
device, decompose=1

CLOSE, uh
 free_lun, uh
CLOSE, uv
 free_lun, uv
CLOSE, ux
 free_lun, ux
end