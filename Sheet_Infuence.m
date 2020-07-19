function [B_sheetn,B_sheett]=Sheet_Infuence(xm,zm,x1,z1,x2,z2,theta,theta_trail,n1,it)
B_sheetn=zeros(n1,1);
B_sheett=B_sheetn;

switch panel_type
    
    case'airfoil'

for i=1:n1
    [a,b]=sourcefish(xm(i),zm(it,i),x1,z1,x2,z2,theta(it,i),theta_trail);
    B_sheetn(i)=-b;
    B_sheett(i)=a;
end
    case'wake'
for i=1:n1
    [a,b]=sourcefish(xm(i),zm(it,i),x1,z1,x2,z2,0,theta_trail);
    B_sheetn(i)=-b;
    B_sheett(i)=a;
end
end
end

