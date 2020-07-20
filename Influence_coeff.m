function[An,At,Bn,Bt]=Influence_coeff(n_coll,n1,it,xm,zm,x,z,theta,panel_type)

An=zeros(n1,n_coll);
At=An;
Bn=An;
Bt=An;


switch panel_type

    case 'airfoil'



for i=1:n1
        for k=1:n_coll
            if i==k
                An(i,k)=0.5;
                Bn(i,k)=0;
            else
            [a,b]=sourcefish(xm(i),zm(it,i),x(k),z(it,k),x(k+1),z(it,k+1),theta(it,i),theta(it,k));
            An(i,k)=a;             
            Bn(i,k)=-b;
            end
        end
        At=-Bn;
        Bt=An;
end
        
     case 'sheet'
    
    for i=1:n1
        for k=1:n_coll
            
            [a,b]=sourcefish(xm,zm,x(k),z(it,k),x(k+1),z(it,k+1),0,theta(it,k));
            An(i,k)=a;             
            Bn(i,k)=-b;
            
        end
        At=-Bn;
        Bt=An;
    end


case 'wake'
    
    for i=1:n1
        for k=1:n_coll
            
            [a,b]=sourcefish(xm(i,1),zm(i,2),x(k),z(it,k),x(k+1),z(it,k+1),0,theta(it,k));
            An(i,k)=a;             
            Bn(i,k)=-b;
            
        end
        At=-Bn;
        Bt=An;
    end
end

    
        
end
