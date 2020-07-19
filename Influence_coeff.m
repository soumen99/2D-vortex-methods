function[An,At,Bn,Bt]=Influence_coeff(n_coll,n1,it,xcoll,zcoll,xmain,zmain,theta)

An=zeros(n1,n_coll);
At=An;
Bn=An;
Bt=An;


switch panel_type

    case 'bound'



for i=1:n1
        for k=1:n_coll
            if i==k
                An(i,k)=0.5;
                Bn(i,k)=0;
            else
            [a,b]=sourcefish(xcoll(i),zcoll(it,i),xmain(k),zmain(it,k),xmain(k+1),zmain(it,k+1),theta(it,i),theta(it,k));
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
            
            [a,b]=sourcefish(xmid,zmid,xmain(k),zmain(it,k),xmain(k+1),zmain(it,k+1),0,theta(it,k));
            An(i,k)=a;             
            Bn(i,k)=-b;
            
        end
        At=-Bn;
        Bt=An;
    end


case 'wake'
    
    for i=1:n1
        for k=1:n_coll
            
            [a,b]=sourcefish(vortic(i,1),vortic(i,2),xmain(k),zmain(it,k),xmain(k+1),zmain(it,k+1),0,theta(it,k));
            An(i,k)=a;             
            Bn(i,k)=-b;
            
        end
        At=-Bn;
        Bt=An;
    end
end

    
        
end
