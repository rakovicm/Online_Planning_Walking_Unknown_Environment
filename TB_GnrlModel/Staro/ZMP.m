function rzmp=ZMP(F,rk)
rk(3,:)=0;
F([1,2,6],:)=0;
M=F(4:6,:);
F=F(1:3,:);
L=cross(M,F);
L=L/VecMod(L);
rzmp=rk-L*VecMod(M)./VecMod(F);