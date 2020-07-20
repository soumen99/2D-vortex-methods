function[sol1,sol2]=solve(An,rhs1,rhs2,n_coll)
Ab1=[An rhs1];
Ab2=[An rhs2];
         for i=1:n_coll-1
             for j=i+1:n_coll
                 alpha1=Ab1(j,i)/Ab1(i,i);
                 alpha2=Ab2(j,i)/Ab2(i,i);
                 Ab1(j,:)=Ab1(j,:)-alpha1*Ab1(i,:);
                 Ab2(j,:)=Ab2(j,:)-alpha2*Ab2(i,:);
             end
         end
         sol1=zeros(n_coll,1);
         sol2=sol1;
         for i=n_coll:-1:1
             sol1(i)= (Ab1(i,end)-Ab1(i,i+1:n_coll)*sol1(i+1:n_coll))/Ab1(i,i);
             sol2(i)= (Ab2(i,end)-Ab2(i,i+1:n_coll)*sol2(i+1:n_coll))/Ab2(i,i);
         end