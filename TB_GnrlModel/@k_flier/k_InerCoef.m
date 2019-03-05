function [am a0m bm b0m]=k_InerCoef(a)
%%%% proveriti
am=a.bet;
a0m=a.del;
bm=a.alf;
b0m=a.gam;

for ln=1:a.N    
    am(:,ln,:)=-a.ms(ln)*am(:,ln,:);
    a0m(:,ln) =-a.ms(ln)*a0m(:,ln);
    bbm=zeros(3,a.N);
    bbm(:,:)=bm(:,ln,:);
    bm(:,ln,:)=-a.A(:,:,ln)*a.J(:,:,ln)*(a.A(:,:,ln)')*bbm;
    b0m(:,ln) =-a.A(:,:,ln)*a.J(:,:,ln)*(a.A(:,:,ln)')*b0m(:,ln)-a.A(:,:,ln)*cross((a.A(:,:,ln)')*a.omg(:,ln), a.J(:,:,ln)*((a.A(:,:,ln)')*a.omg(:,ln)));      
end;
