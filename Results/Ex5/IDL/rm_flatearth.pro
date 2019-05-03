function rm_flatearth, img, wavevector, sensorheight, baseline, theta, gnd_rng_res, deltar

; this function calculates phase shift on a line on the ground, parallel to the baseline
; ONLY IF the BASELINE IS PERPENDICULAR TO THE ORBIT, this function is VALID and the fringes are in the direction of range.
; elseif the baseline is parallel to the orbit, THIS FUNCTION HAS TO BE MODIFIED, such that the fringes are in the direction of azimut.
; else if the baseline is neither parallel nor perpenticular to the orbit, 
; the angle theta is a COMBINATION OF AZIMUT AND RANGE. 
;
; this function is also ONLY VALID, if the baseline is much smaller that the height of the sensor.     


dim = size(img,/dim)
R0  = sensorheight/cos(theta) ; total sensor - target distance 
img_rglength = (dim[0]-1)*gnd_rng_res

y0   = R0*sin(theta) ; ground distance nadir - image center
ymin = y0 - 0.5*img_rglength
ymax = y0 + 0.5*img_rglength
theta_min = atan(ymin/sensorheight)
theta_max = atan(ymax/sensorheight)

theta_arr = findgen(dim[0])/dim[0]*(theta_max-theta_min) + theta_min ; change of angle across image
d_phi     = wavevector*abs(baseline)*sin(theta_arr) ; change in phase across image: wavevector * change of distance across image
delta_phi = fltarr(dim)
for i=0,dim(1)-1 do delta_phi[*,i] = d_phi + deltar*wavevector

return, img*exp(-complex(0.0,1.0)*delta_phi)

end