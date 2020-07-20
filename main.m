clear
clc
%no of cycles
n=input('no of steps>');
%chord length
l=input('chord length>');
%freestream velocity
u_inf=input('freestream velocity>');
%frequency
f=input('frequency>');
%model constants
pi=3.14;
rho=1000;
%angular frequency
omega=2*pi*f;
%max amplitude
h0=0.05*l;
%strouhal number
st=f*h0/u_inf;
%wavelength
lambda=0.6*l;
%wave number
k=2*pi/lambda;
%stepsize
dt=0.001;
%total timesteps
%n=round(N/(f*dt));
%time array
t=0:dt:(n-1)*dt;
%amplitude growth rate
m=2;
%wake array
vortic=zeros(n,2);
%import kinematics
[xmain,zmain,vmain,xcoll,zcoll,vcoll,enx,eny,theta,p] = import_paneldata(l,omega,k,t,m,n,h0);
%number of collocation pts
n_coll=size(xcoll,2);
%locomotion in x
sx=-u_inf*t;
%vortex location wrt inertial coordinates
vortic(:,1)=l+sx;
vortic(:,2)=zmain(:,1);
tau=zeros(n,1);
%iterations to calculate sheet parameters
num_iters=10;
%initial parameters of attached sheet
u_wake=u_inf;
v_wake=vmain(1,1);
tan_trail=v_wake/u_wake;
theta_trail=atan(tan_trail);
delta=dt*sqrt(u_wake^2+v_wake^2);
xpanel1=xmain(1);
xpanel2=xmain(1)+delta/(1+tan_trail^2);
xmid=(xpanel1+xpanel2)/2;
zpanel1=zmain(1,1);
zpanel2=zmain(1,1)+delta*tan_trail/(1+tan_trail^2);
zmid=(zpanel1+zpanel2)/2;

counter=zeros(10,1);
%start time iteration
for it=1:n

if it==1
            circ=0;
        else
            circ=tau(it-1);
end
        
for iter=1:num_iters  
C_n=zeros(n_coll,1);
C_t=C_n;

%B_bound=zeros(n_coll);
%B_sheet=zeros(n_coll,1);
%Wake=zeros(n_coll,1);

%setup attached wake panel parameters for first timestep
if it==1                                                                         
    %for 1st iteration  
    
else
                                                            
it1=it-1;                                                     
%wake influence on collocation pts
[C_n,C_t] = C_matrix(n_coll,it1,vortic,xcoll+sx(it),zcoll,theta,tau,it,'airfoil');
end
%attached wake panel influence on collocation pts
[B_sheetn,B_sheett] = Sheet_Infuence(xcoll+sx(it),zcoll,xpanel1,zpanel1,xpanel2,zpanel2,theta,theta_trail,n_coll,it,'airfoil');
%freestream influence on collocation pts
V_stream_n = u_inf*enx(it,:)+vcoll(it,:).*eny(it,:);   
V_stream_t = u_inf*eny(it,:)-vcoll(it,:).*enx(it,:);
%gather terms to form 2nd right side matrix
B3=circ*B_sheetn/delta;
rhs2=-(B3+V_stream_n.'+C_n);



%bound influence coefficients
[An,At,Bn,Bt]=Influence_coeff(n_coll,n_coll,it,xcoll,zcoll,xmain,zmain,theta,'airfoil');
        B=sum(Bn,2);
        %gather terms to form 1st rhs matrix
        rhs1=(p(it)*B_sheetn/delta)-B;
        %solve system of eqn in terms of TAUk
        [sol1,sol2]=solve(An,rhs1,rhs2,n_coll);
        %setup unsteady kutta condition
        T11=At(1,:)*sol1+sum(Bt(1,:),2)-p(it)*B_sheett(1)/delta;
        Tn1=At(n_coll,:)*sol1+sum(Bt(n_coll,:),2)-p(it)*B_sheett(n_coll)/delta;
        T12=V_stream_t(1)+(circ*B_sheett(1)/delta)+C_t(1)+At(1,:)*sol2;
        Tn2=V_stream_t(n_coll)+(circ*B_sheett(n_coll)/delta)+C_t(n_coll)+At(n_coll,:)*sol2;
        %setup quadratic coefficients
        f1 = T11^2-Tn1^2;
        f2 = 2*(T11*T12-Tn1*Tn2-2/dt);
        f3 = T12^2-Tn2^2+2*circ/dt;
        %solve non linear equation
        s1 = (-f2-sqrt(f2^2-4*f1*f3))/(2*f1);
        s2 = (f2-sqrt(f2^2-4*f1*f3))/(2*f1);
        if abs(s1)>abs(s2)
            tau(it)=s2;
        else
            tau(it)=s1;
        end
        counter(iter)=tau(it);
        q=(tau(it)*sol1/p(it))+sol2;
        
        %calculate sheet velocity for next iteration
        [C_x,C_y] = C_matrix(1,it-1,vortic,xmid,zmid,0,tau,it,'sheet');
        [A_x,A_y,B_x,B_y]=Influence_coeff(n_coll,1,it,xmid,zmid,xmain+sx(it),zmain,theta,'sheet');
        
        %update parameters for next iteration
        u_wake = u_inf+A_x*q+(tau(it)/p(it))*sum(B_x,2)+C_x;
        v_wake = A_y*q+(tau(it)/p(it))*sum(B_y,2)+C_y;
        u_wake(isnan(u_wake))=0;
        v_wake(isnan(v_wake))=0;
        tan_trail=v_wake/u_wake;
        theta_trail=atan(tan_trail);
        delta=dt*sqrt(u_wake^2+v_wake^2);
        xpanel1=xmain(1)+sx(it);
        xpanel2=xpanel1+delta/(1+tan_trail^2);
        xmid=(xpanel1+xpanel2)/2;
        zpanel1=zmain(it,1);
        zpanel2=zpanel1+delta*tan_trail/(1+tan_trail^2);
        zmid=(zpanel1+zpanel2)/2;
        vortic(it,1)=xmid;
        vortic(it,2)=zmid;
        
end
        tau_sheet=(circ-tau(it))/delta;
       %calculate airfoil influence on wake
        [P_x,P_z,Q_x,Q_z]=Influence_coeff(n_coll,it-1,it,vortic(:,1),vortic(:,2),xmain,zmain,theta,'bound');
        %calculate sheet influence on wake
        [S_x,S_z]=Sheet_Infuence(vortic(:,1),vortic(:,2),xpanel1,zpanel1,xpanel2,zpanel2,0,theta_trail,it-1,it,'wake');
        ucore=P_x*q+tau(it)*sum(Q_x,2)/p(it)+tau_sheet*S_x;
        vcore=P_z*q+tau(it)*sum(Q_z,2)/p(it)+tau_sheet*S_z;
        ucore(isnan(ucore)) = 0;
        vcore(isnan(vcore)) = 0;
        %calculate wake-wake influence
        [W_x,W_y]=C_matrix(it-1,it-1,vortic,vortic(:,1),vortic(:,2),0,tau,it,'wake');
        %update total wake speed
        ucore=ucore+W_x+u_inf;
        vcore=vcore+W_y;
        %update position
        vortic(1:it-1,1)=vortic(1:it-1,1)+ucore*dt;
        vortic(1:it-1,2)=vortic(1:it-1,2)+vcore*dt;
end 



        
    










