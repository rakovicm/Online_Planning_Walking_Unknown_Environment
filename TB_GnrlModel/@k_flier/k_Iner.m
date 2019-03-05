function [Fi Mi]=k_Iner(a,qa)
Fi=zeros(3,a.N);
Mi=zeros(3,a.N);
al=zeros(3,a.N);
bl=zeros(3,a.N);
for ln=1:a.N
    al(:,:)=a.alf(:,ln,:);
    bl(:,:)=a.bet(:,ln,:);
    gl=a.gam(:,ln);
    dl=a.del(:,ln);
    Fi(:,ln)=-(bl(:,:)*qa+ dl)*a.ms(ln);
    
    alfmod=-a.A(:,:,ln)*a.J(:,:,ln)*(a.A(:,:,ln)')*al;
    gammod=-a.A(:,:,ln)*a.J(:,:,ln)*(a.A(:,:,ln)')*gl-a.A(:,:,ln)*cross((a.A(:,:,ln)')*a.omg(:,ln), a.J(:,:,ln)*((a.A(:,:,ln)')*a.omg(:,ln)));
    Mi(:,ln)=alfmod*qa+gammod;
end;