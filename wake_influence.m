function [Cn,Ct]= wake_influence(vortx,vortz,x,z,theta)
pi=3.14;
r=sqrt((vortx-x)^2+(z-vortz)^2);
thetai=atan((z-vortz)/(x-vortx));
Cn=-cos(theta-thetai)/(2*pi*r);
Ct=-sin(theta-thetai)/(2*pi*r);
end
