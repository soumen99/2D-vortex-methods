function [B_sheetn,B_sheett]=Sheet_Infuence(xcoll,zcoll,xpanel1,zpanel1,xpanel2,zpanel2,theta,theta_trail,n_coll,it)
B_sheetn=zeros(n_coll,1);
B_sheett=B_sheetn;
for i=1:n_coll
    [a,b]=sourcefish(xcoll(i_coll),zcoll(it,i_coll),xpanel1,zpanel1,xpanel2,zpanel2,theta(it,i_coll),theta_trail);
    B_sheetn(i)=-b;
    B_sheett(i)=a;
end
end
