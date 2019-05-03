@colorbar_define

PRO q06_HAalpha

; DATA FILE READING

; Opening files

OPENR, u_, '../i_al_af_1206_hh_corr.dat', /xdr, /get_lun

; Initialization of a two-element array

dim = lonarr(2)

; Reading the first two values of the file (image size)

READU, u_, dim

; Initialization of a complex-element matrix

hh = complexarr(dim)

; Reading the rest of the values of the file

READU, u_, hh

; Closing file

CLOSE, u_

free_lun, u_

; Opening files

OPENR, u_, '../i_al_af_1206_vv_corr.dat', /xdr, /get_lun

; Initialization of a two-element array

dim = lonarr(2)

; Reading the first two values of the file (image size)

READU, u_, dim

; Initialization of a complex-element matrix

vv = complexarr(dim)

; Reading the rest of the values of the file

READU, u_, vv

; Closing file

CLOSE, u_

free_lun, u_

; Opening files

OPENR, u_, '../i_al_af_1206_vv_corr.dat', /xdr, /get_lun

; Initialization of a two-element array

dim = lonarr(2)

; Reading the first two values of the file (image size)

READU, u_, dim

; Initialization of a complex-element matrix

xx = complexarr(dim)

; Reading the rest of the values of the file

READU, u_, xx

; Closing file

CLOSE, u_

free_lun, u_

; COHERENCY MATRIX

quantities= complexarr([dim,9])

quantities[*,*,0]=(abs(hh+vv))^2

quantities[*,*,1]=(hh+vv)*conj(hh-vv)

quantities[*,*,2]=2*(hh+vv)*conj(xx)

quantities[*,*,3]=(hh-vv)*conj(hh+vv)

quantities[*,*,4]=(abs(hh-vv))^2

quantities[*,*,5]=2*(hh-vv)*conj(xx)

quantities[*,*,6]=2*(xx)*conj(hh+vv)

quantities[*,*,7]=2*(xx)*conj(hh-vv)

quantities[*,*,8]=4*(abs(xx))^2

quantities=quantities/2.;

coherency= complexarr([dim,9])

win_sz=7

for i=0,8 do coherency[*,*,i]=smooth(quantities[*,*,i],win_sz,/EDGE_TRUNCATE, /NAN)

m_SEL=100

n_SEL=100

; COMPUTATION OF EIGENVALUES with LA_EIGENQL

eigenvalues=complexarr([dim,3])

eigenvec=complexarr([dim,3,3])

alpha1=fltarr([dim])

alpha2=fltarr([dim])

alpha3=fltarr([dim])

coh=complexarr(3,3)

for m=0,dim[0]-1 do begin

        for n=0,dim[1]-1 do begin

                coh[*,0]=coherency[m,n,0:2]

                coh[*,1]=coherency[m,n,3:5]

                coh[*,2]=coherency[m,n,6:8]

                eigenvalues[m,n,*] = LA_EIGENQL(coh,eigenvectors=eigenvec)

                alpha1[m,n]=acos(abs(eigenvec[0,2]))

                alpha2[m,n]=acos(abs(eigenvec[0,1]))

                alpha3[m,n]=acos(abs(eigenvec[0,0]))

        endfor

endfor


; APPEARANCE PROBABILITIES, ENTROPHY AND ANISOTROPHY

prob=complexarr(3)

H=fltarr([dim])

A=fltarr([dim])

alpha=fltarr([dim])

tot=total(eigenvalues,3)

for m=0,dim[0]-1 do begin

        for n=0,dim[1]-1 do begin

                prob=eigenvalues[m,n,*]/tot[m,n]

                H[m,n]=-total(prob*alog10(prob)/alog10(3))

                A[m,n]=(eigenvalues[m,n,1]-eigenvalues[m,n,0])/(eigenvalues[m,n,1]+eigenvalues[m,n,0])

                alpha_[m,n]=alpha1[m,n]*prob[2]+alpha2[m,n]*prob[1]+alpha3[m,n]*prob[0]

        endfor

endfor


; COMPUTATION OF EIGENVALUES analitically

eigenvalues_AC=complexarr([dim,3])

eigenvec_AC=complexarr([dim,3,3])

alpha1_AC=fltarr([dim])

alpha2_AC=fltarr([dim])

alpha3_AC=fltarr([dim])

coh=complexarr(3,3)

for m=0,dim[0]-1 do begin

        for n=0,dim[1]-1 do begin

                coh[*,0]=coherency[m,n,0:2]

                coh[*,1]=coherency[m,n,3:5]

                coh[*,2]=coherency[m,n,6:8]

a=-(coherency[m,n,0]+coherency[m,n,4]+coherency[m,n,8])

b=-((abs(coherency[m,n,1]))^2+(abs(coherency[m,n,2]))^2+(abs(coherency[m,n,5]))^2-coherency[m,n,0]*coherency[m,n,4]-coherency[m,n,0]*coherency[m,n,8]-coherency[m,n,4]*coherency[m,n,8])

c=-(coherency[m,n,0]*coherency[m,n,4]*coherency[m,n,8]-coherency[m,n,0]*(abs(coherency[m,n,5]))^2-coherency[m,n,8]*(abs(coherency[m,n,1]))^2-coherency[m,n,4]*(abs(coherency[m,n,2]))^2$

+coherency[m,n,1]*conj(coherency[m,n,2])*coherency[m,n,5]+conj(coherency[m,n,1])*coherency[m,n,2]*conj(coherency[m,n,5]))

eigenvalues_AC[0]=-a/3.-1/3.*((2*a^3-9*a*b+27*c+sqrt(complex((2*a^3-9*a*b+27*c)^2-4*(a^2-3*b)^3,0)))/2.)^(1/3.)-1/3.*((2*a^3-9*a*b+27*c-sqrt(complex((2*a^3-9*a*b+27*c)^2-4*(a^2-3*b)^3,0)))/2.)^(1/3.)

eigenvalues_AC[1]=-a/3.+complex(1,sqrt(3))/6.*((2*a^3-9*a*b+27*c+sqrt(complex((2*a^3-9*a*b+27*c)^2-4*(a^2-3*b)^3,0)))/2.)^(1/3.)+complex(1,-sqrt(3))/6.*((2*a^3-9*a*b+27*c-sqrt(complex((2*a^3-9*a*b+27*c)^2-$

4*(a^2-3*b)^3,0)))/2.)^(1/3.)

eigenvalues_AC[2]=-a/3.+complex(1,-sqrt(3))/6.*((2*a^3-9*a*b+27*c+sqrt(complex((2*a^3-9*a*b+27*c)^2-4*(a^2-3*b)^3,0)))/2.)^(1/3.)+complex(1,sqrt(3))/6.*((2*a^3-9*a*b+27*c-sqrt(complex((2*a^3-9*a*b+27*c)^2-$

4*(a^2-3*b)^3,0)))/2.)^(1/3.)

for i=0,2 do begin

        Qn=complex((coherency[m,n,8]-eigenvalues_AC[i])*conj(coherency[m,n,1])-conj(coherency[m,n,2])*coherency[m,n,5])

        Qd=complex((coherency[m,n,4]-eigenvalues_AC[i])*conj(coherency[m,n,2])-conj(coherency[m,n,1])*conj(coherency[m,n,5]))

        Q=Qn/Qd

        eigenvec_AC[*,i]=[-(coherency[m,n,5]/conj(coherency[m,n,1]))-(coherency[m,n,4]-eigenvalues_AC[i])*Q/conj(coherency[m,n,1]),Q,1]

        eig_norm=sqrt(total(abs(eigenvec_AC[*,i])^2))

        eigenvec_AC[*,i]=eigenvec_AC[*,i]/eig_norm

endfor

          

                alpha1_AC[m,n]=acos(abs(eigenvec_AC[0,2]))

                alpha2_AC[m,n]=acos(abs(eigenvec_AC[0,1]))

                alpha3_AC[m,n]=acos(abs(eigenvec_AC[0,0]))

        endfor

endfor

; APPEARANCE PROBABILITIES, ENTROPHY AND ANISOTROPHY

prob_AC=complexarr(3)

H_AC=fltarr([dim])

A_AC=fltarr([dim])

alpha_AC=fltarr([dim])

tot_AC=total(eigenvalues_AC,3)

for m=0,dim[0]-1 do begin

        for n=0,dim[1]-1 do begin

                prob_AC=eigenvalues_AC[m,n,*]/tot_AC[m,n]

                H_AC[m,n]=-total(prob_AC*alog10(prob_AC)/alog10(3))

                A_AC[m,n]=(eigenvalues_AC[m,n,1]-eigenvalues_AC[m,n,0])/(eigenvalues_AC[m,n,1]+eigenvalues_AC[m,n,0])

                alpha_AC[m,n]=alpha1_AC[m,n]*prob_AC[2]+alpha2_AC[m,n]*prob_AC[1]+alpha3_AC[m,n]*prob_AC[0]

        endfor

endfor
stop
; ; Visualisation
; 
; loadct,39
; 
; device, decompose=0
; 
; window, 1, xsize=1000+320, ysize=600, TITLE ="Entrophy"
; 
; TV,congrid(bytscl(H,0,1),1000,600)
; 
; colorbar = Obj_New("COLORBAR", Format='(F8.2)',/vertical)
; 
; colorbar->Draw
; 
; colorbar->SetProperty, Range=[0,1], /Erase, /Draw
; 
; Obj_Destroy, colorbar
; 
; !p.background=255
; 
; device, decompose=1
; 
; loadct,39
; 
; device, decompose=0
; 
; window, 2, xsize=1000+320, ysize=600, TITLE ="Anisotrophy"
; 
; TV,congrid(bytscl(A,0,1),1000,600)
; 
; colorbar = Obj_New("COLORBAR", Format='(F8.2)',/vertical)
; 
; colorbar->Draw
; 
; colorbar->SetProperty, Range=[0,1], /Erase, /Draw
; 
; Obj_Destroy, colorbar
; 
; !p.background=255
; 
; device, decompose=1
; 
; loadct,39
; 
; device, decompose=0
; 
; window, 3, xsize=1000+320, ysize=600, TITLE ="Alpha (total)"
; 
; TV,congrid(bytscl(alpha*180/!pi,0,90),1000,600)
; 
; colorbar = Obj_New("COLORBAR", Format='(F8.2)',/vertical)
; 
; colorbar->Draw
; 
; colorbar->SetProperty, Range=[0,90], /Erase, /Draw
; 
; Obj_Destroy, colorbar
; 
; !p.background=255
; 
; device, decompose=1
; 
; loadct,39
; 
; device, decompose=0
; 
; window, 4, xsize=1000+320, ysize=600, TITLE ="Alpha 1"
; 
; TV,congrid(bytscl(alpha1*180/!pi,0,90),1000,600)
; 
; colorbar = Obj_New("COLORBAR", Format='(F8.2)',/vertical)
; 
; colorbar->Draw
; 
; colorbar->SetProperty, Range=[0,90], /Erase, /Draw
; 
; Obj_Destroy, colorbar
; 
; !p.background=255
; 
; device, decompose=1
; 
; ; Histograms
; 
; loadct,0
; 
; device, decompose=0
; 
; window,5,xsize=1000, ysize=600, TITLE ="Histogram - Entrophy"
; 
; x_ref=findgen(101)/100
; 
; plot, x_ref,100*HISTOGRAM(H*100)/total(histogram(H)), XTITLE = "H", YTITLE = "pdf(H)"
; 
; loadct,0
; 
; device, decompose=0
; 
; window,6,xsize=1000, ysize=600, TITLE ="Histogram - Anisotrophy"
; 
; x_ref=findgen(101)/100
; 
; plot, x_ref,100*HISTOGRAM(A*100)/total(histogram(A)), XTITLE = "A", YTITLE = "pdf(A)"
; 
; loadct,0
; 
; device, decompose=0
; 
; window,7,xsize=1000, ysize=600, TITLE ="Histogram - Alpha (total)"
; 
; x_ref=findgen(101)/100*90
; 
; plot, x_ref,HISTOGRAM(alpha*180/!pi)/total(histogram(alpha*180/!pi)), XTITLE = "alpha", YTITLE = "pdf(alpha)"
; 
; loadct,0
; 
; device, decompose=0
; 
; window,8,xsize=1000, ysize=600, TITLE ="Histogram - Alpha 1"
; 
; x_ref=findgen(91)
; 
; plot, x_ref,HISTOGRAM(alpha1*180/!pi)/total(histogram(alpha1*180/!pi)), XTITLE = "alpha 1", YTITLE = "pdf(alpha 1)"
; 
; device, decompose=0
; 
; window,9,xsize=1000, ysize=600, TITLE ="Two-dimensional histogram - alpha vs. H"
; 
; shade_surf, hist_2D(alpha*180/!pi,H*100)/total(hist_2D(alpha*180/!pi,H*100)),findgen(88),findgen(100)/100,XTITLE = "alpha", YTITLE = "H",ZTITLE="pdf(alpha, H)"

END