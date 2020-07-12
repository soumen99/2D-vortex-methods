clear
clc
%no of cycles
N=input('no of cycles>');
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
n=round(N/(f*dt));
%time array
t=0:dt:(n-1)*dt;
%amplitude growth rate
m=2;
%import kinematics
[xmain,zmain,vmain,xcoll,zcoll,vcoll,enx,eny,theta,p] = import_paneldata(l,omega,k,t,m,n,h0);
%number of collocation pts
n_coll=size(x_main,1)-1;
%vortex location wrt inertial coordinates
vortic(:,1)=l+sx;
vortic(:,2)=zmain(:,1);
%locomotion in x
sx=-u_inf*t;
%start time iteration
for it=1:n
    
%A=zeros(n_coll);
%B_bound=zeros(n_coll);
%B_sheet=zeros(n_coll,1);
%Wake=zeros(n_coll,1);

%setup attached wake panel parameters for first timestep
if it==1
    %for 1st iteration
u_wake=u0;
v_wake=vmain(1,1);
tan_trail=u_wake/v_wake;
theta_trail=atan(tan_trail);
delta=dt*sqrt(u_wake^2+v_wake^2);
xpanel1=xmain(1);
xpanel2=xmain(1)+delta/(1+tan_trail^2);
zpanel1=zmain(it,1);
zpanel2=zmain(it,1)+delta*tan_trail/(1+tan_trail^2);
    
else
    it1=it-1;
         %to find downwash and total circulation due to n-1 vortexes                                       
%              %for i_coll=1:n_coll
%                  for j = 1:it1
%                  [Cn,Ct]=wake_influence(vortic(j,1),vortic(j,2),xcoll(i_coll),zcoll(it,i_coll),theta(it,i_coll));
%                  Wake(i_coll)=Wake(i_coll)+Cn*(tau(j-1)-tau(j));
%                  if i_coll==1
%                      C_1=C_1+Ct*(tau(j-1)-tau(j));
%                  end
%                  if i_coll==n_coll
%                      C_n=C_n+Ct*(tau(j-1)-tau(j));
%                  end
%                  [B_sheet(i_coll), ]=sourcefish(xcoll(i_coll),zcoll(it,i_coll),xpanel1,zpanel1,xpanel2,zpanel2,theta(it,i_coll),theta_trail);
%                  end
%              
%              end
%wake influence on collocation pts
[C_n,C_t] = C_matrix(n_coll,it1,vortic,xcoll,zcoll,theta,tau,it);
%attached wake panel influence on collocation pts
[B_sheetn,B_sheett] = Sheet_Infuence(xcoll,zcoll,xpanel1,zpanel1,xpanel2,zpanel2,theta,theta_trail,n_coll,it);
%freestream influence on collocation pts
V_stream_n = u_inf*enx(it,:)+vcoll(it,:).*eny(it,:);   
V_stream_t = u_inf*eny(it,:)-vcoll(it,:).*enx(it,:);
%gather terms to form 2nd right side matrix
B3=tau(it-1)*B_sheetn/delta;
rhs2=-(B3+V_stream_n+C_n);
end

% for i_num=1:n_coll
%         for k_num=1:n_coll
%             if i_num==k_num
%                 A(i_num,k_num)=0.5;
%                 B_bound(i_num,k_num)=0;
%             else
%             [a,b]=sourcefish(xcoll(i_num),zcoll(it,i_num),xmain(k_num),zmain(it,k_num),xmain(k_num+1),zmain(it,k_num+1),theta(it,i_num),theta(it,k_num));
%             A(i_num,k_num)=a;             
%             B_bound(i_num,k_num)=-b;
%             if i_num==1
%                 A_1=A_1-B_bound(i_num,k_num);
%             else if i_num==n_coll
%                     A_n=A_n-B_bound(i_num,k_num);
%                 end
%             end
%             end
%         end
%bound influence coefficients
[An,At,Bn,Bt]=Influence_coeff(n_coll,it,xcoll,zcoll,xmain,zmain,theta);
        B=sum(Bn,2);
        %gather terms to form 1st rhs matrix
        rhs1=(p(it)*B_sheetn/delta)-B;
        %solve system of eqn in terms of TAUk
        [sol1,sol2]=solve(An,rhs1,rhs2,n_coll);
        %setup unsteady kutta condition
        T11=At(1,:)*sol1+sum(Bt(1,:),2)-p(it)*Bsheett(1)/delta;
        Tn1=At(n_coll,:)*sol1+sum(Bt(n_coll,:),2)-p(it)*Bsheett(n_coll)/delta;
        T12=V_stream_t(1)+(tau(it-1)*Bsheett(1)/delta)+C_t(1);
        Tn2=V_stream_t(n_coll)+(tau(it-1)*Bsheett(n_coll)/delta)+C_t(n_coll);
        
end 


        
    










