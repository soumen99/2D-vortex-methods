function [Cn,Ct]= wake_influence(vortx,vortz,xcoll,zcoll,theta)
r=sqrt((vortx-xcoll)^2+(zcoll-vortz)^2);
thetai=atan((zcoll-vortz)/(xcoll-vortx));
Cn=-cos(theta-thetai)/(2*pi*r);
Ct=-sin(theta-thetai)/(2*pi*r);
end
