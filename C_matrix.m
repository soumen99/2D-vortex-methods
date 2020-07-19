function [C_n,C_t] = C_matrix(n1,it1,vortic,x,z,theta,tau,it)
C_n=zeros(n1,1);
C_t=C_n;
%to find downwash and total circulation due to n-1 vortexes 

switch panel_type
    
    case'airfoil'
        
             for i_coll=1:n1
                 for j = 1:it1
                     
                 [Normal,Tang]=wake_influence(vortic(j,1),vortic(j,2),x(i_coll),z(it,i_coll),theta(it,i_coll));
                 C_n(i_coll)=C_n(i_coll)+Normal*(tau(j-1)-tau(j));
                 C_t(i_coll)=C_t(i_coll)+Tang*(tau(j-1)-tau(j));
                     
                 end
             end
             
    case'wake'
        
        for i_coll=1:n1
                 for j = 1:it1
                     if i_coll~=j
                 [Normal,Tang]=wake_influence(vortic(j,1),vortic(j,2),x(i_coll),z(it,i_coll),0);
                 C_n(i_coll)=C_n(i_coll)+Normal*(tau(j-1)-tau(j));
                 C_t(i_coll)=C_t(i_coll)+Tang*(tau(j-1)-tau(j));
                     end 
                 end
        end
end
        
end