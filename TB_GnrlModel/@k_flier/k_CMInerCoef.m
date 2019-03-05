function [Fk Fk0 Mk Mk0]=k_CMInerCoef(a)
%%%% proveriti
am=a.bet;
a0m=a.del;
bm=a.alf;
b0m=a.gam;
[rcm ]=k_CentarMase(a);
Fk=zeros(3,a.N);
Fk0=zeros(3,1);
Mk=Fk;
Mk0=Fk0;
aam=zeros(3,a.N);
bbm=aam;
for ln=1:a.N    
    aam(:,:)=a.ms(ln)*am(:,ln,:);
    Fk =Fk-aam;
    Fk0=Fk0-a.ms(ln)*a0m(:,ln);  
    bbm(:,:)=bm(:,ln,:);
    Mk=Mk-a.A(:,:,ln)*a.J(:,:,ln)*(a.A(:,:,ln)')*bbm-cross((a.rc(:,ln,1)-rcm)*ones(1,a.N),aam);
    Mk0=Mk0-a.A(:,:,ln)*a.J(:,:,ln)*(a.A(:,:,ln)')*b0m(:,ln)-a.A(:,:,ln)*cross((a.A(:,:,ln)')*a.omg(:,ln), a.J(:,:,ln)*((a.A(:,:,ln)')*a.omg(:,ln)));      
    Mk0=Mk0-cross(a.rc(:,ln,1)-rcm,a.ms(ln)*a0m(:,ln));
end;
