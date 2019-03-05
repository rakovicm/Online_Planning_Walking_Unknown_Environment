function E=k_Energija(a)
E=0;
g=9.81;
for i=1:a.N
    V=a.v(:,i);
    Omg=a.omg(:,i);
    Q=a.A(:,:,i);
    ms=a.ms(i);
    J=a.J(:,:,i);
    h=a.rc(3,i,1);
    E=E+ms*V'*V/2+Omg'*Q*J*Q'*Omg/2+ms*g*h;
end;