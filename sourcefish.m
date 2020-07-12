function[An,At]=sourcefish(xm,ym,x1,y1,x2,y2,theta1,theta2)
pi=3.14;
r1=sqrt((x1-xm)^2+(y1-ym)^2);
r2=sqrt((x2-xm)^2+(y2-ym)^2);
beta=atan((ym-y2)/(xm-x2))-atan((ym-y1)/(xm-x1));
An=(1/2*pi)*(sin(theta1-theta2)*log(r2/r1)+cos(theta1-theta2)*beta);
At=(1/2*pi)*(-cos(theta1-theta2)*log(r2/r1)+sin(theta1-theta2)*beta);
end


