function [ZH Zh0]=k_calcZH(flier, ra)
% racuna matrice za izracunavanje ukupnog momenta sistema za tacku ra
ZH=zeros(3,flier.N);
Zh0=zeros(3,1);
al=zeros(3,flier.N);
bl=zeros(3,flier.N);
for ln=1:flier.N
    al(:,:)=flier.alf(:,ln,:);
    bl(:,:)=flier.bet(:,ln,:);
    gl=flier.gam(:,ln);
    dl=flier.del(:,ln);
    
    AA=flier.A(:,:,ln);
    alfmod=AA*flier.J(:,:,ln)*(AA')*al;
    gammod=AA*flier.J(:,:,ln)*(AA')*gl+AA*cross(AA'*flier.omg(:,ln), flier.J(:,:,ln)*AA'*flier.omg(:,ln));
    
    rcnov=(flier.rc(:,ln)-ra);
    ZH=ZH-cross(rcnov*ones(1,flier.N),bl)*flier.ms(ln)-alfmod;
    Zh0=Zh0+cross(rcnov,[0;0;-9.81]-dl)*flier.ms(ln)-gammod; %%%%% zasto jedan plus a drugi minus !!!!!!!!!!!
    
end;