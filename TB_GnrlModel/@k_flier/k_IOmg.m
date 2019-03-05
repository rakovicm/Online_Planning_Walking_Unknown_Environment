function [L]=k_IOmg(a)
%daje i*omega u globalnim koordinatama
%sracunati za centar mase
[rcm Vcm]=k_CentarMase(a);
L=0;
for i=1:a.N
    Q=a.A(:,:,i);
    L=L + cross(a.rc(:,i,1)-rcm,a.v(:,i)-Vcm)*a.ms(i) + Q*((a.J(:,:,i)*inv(Q)*a.omg(:,i)));
end;