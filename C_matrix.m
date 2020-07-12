function [C_n,C_t] = C_matrix(n_coll,it1,vortic,xcoll,zcoll,theta,tau,it)
C_n=zeros(n_coll,1);
C_t=C_n;
%to find downwash and total circulation due to n-1 vortexes                                       
             for i_coll=1:n_coll
                 for j = 1:it1
                 [Normal,Tang]=wake_influence(vortic(j,1),vortic(j,2),xcoll(i_coll),zcoll(it,i_coll),theta(it,i_coll));
                 C_n(i_coll)=C_n(i_coll)+Normal*(tau(j-1)-tau(j));
                 C_t(i_coll)=C_t(i_coll)+Tang*(tau(j-1)-tau(j));
                 end
             end
end