function xp=slackras(s,x)
global flier TAC Kn Cn

xx=x(1:end-1);
t=x(end);
f=fliergeneral12(t,xx);

N=flier.N;
qv=xx(N+1:2*N);
qa=f(N+1:2*N);
[Jc Ac]=km_jakP(flier,TAC); 


sum=-Kn*(Jc(3,:)*qv)-Cn* (Jc(3,:)*qa+Ac(3));
xp=[f;1]*inv(sum);