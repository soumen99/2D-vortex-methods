function [xmain,zmain,vmain,xcoll,zcoll,vcoll,enx,enz,theta,p] = import_paneldata(l,omega,k,t,m,n,h0)

path(path,'C:\Users\ACER PC\Desktop\Projects\Flexible Bodies\Airfoil data');
[data]=xlsread('naca 0010.xlsx');
num=length(data);
xmain=zeros(1,num);
z=xmain;
zmain=zeros(n,num);
vmain=zmain;
p=zeros(1,n);
for i=1:num
    xmain(i)=l*data(num+1-i,1);
    z(i)=l*data(num+1-i,2);
    zmain(:,i)= z(i)-h0*sin(k*xmain(i)-omega*t)*exp(m*(xmain(i)/l-1));
    vmain(:,i)= h0*omega*cos(k*xmain(i)-omega*t)*exp(m*(xmain(i)/l-1));
end

xcoll=zeros(1,num-1);
zcoll=zeros(n,num-1);
vcoll=zcoll;
theta=zcoll;
enx=zcoll;
enz=zcoll;

for i=1:num-1
xcoll(i)=(xmain(i)+xmain(i+1))/2;
zcoll(:,i)=(zmain(:,i)+zmain(:,i+1))/2;
vcoll(:,1)=(vmain(:,1)+vmain(:,i+1))/2;
theta(:,i)=atan((zmain(:,i+1)-zmain(:,i))/(xmain(i+1)-xmain(i)));
enx(:,i)=-sin(theta(i));
enz(:,i)=cos(theta(i));
end
i=1:num-1;
for j=1:n
    p(j)=sum(((xmain(i+1)-xmain(i)).^2+(zmain(j,i+1)-zmain(j,i)).^2).^(0.5));
end
end



