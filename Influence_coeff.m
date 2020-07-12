function[An,At,Bn,Bt]=Influence_coeff(n_coll,it,xcoll,zcoll,xmain,zmain,theta)

An=zeros(n_coll);
At=An;
Bn=An;
Bt=An;
for i=1:n_coll
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
