function [I]=k_InerSys(a)
I=eye(3);
[rcm,vcm]=k_CentarMase(a);
for ln=1:a.N
    JG=a.A(:,:,ln)*a.J(:,:,ln)*a.A(:,:,ln)';
    r=a.rc(:,ln,1)-rcm;
    I(1,1)=I(1,1)+JG(1,1)+(r(2)^2+r(3)^2)*a.ms(ln);
    I(2,2)=I(2,2)+JG(2,2)+(r(1)^2+r(3)^2)*a.ms(ln);
    I(3,3)=I(3,3)+JG(3,3)+(r(2)^2+r(1)^2)*a.ms(ln);
end;
