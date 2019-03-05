function xp=slackkon(s,x)
global flier TAC 

xx=x(1:end-1);
t=x(end);
f=fliergeneral12(t,xx);
N=flier.N;
qv=x(N+1:2*N);
[Jc Ac]=km_jakP(flier,TAC); %#ok<NASGU>   
dd=zeros(3,12);
dd(:)      =f(3*N-5:3*N+30); %%%%%%%%delta
del        =dd(3,TAC);
%%%%%%%%%%%%%%%%%%%5 suma je izvod uslova po vremenu u nasem slucaj
%%%%%%%%%%%%% brzina kontaktne tacke
sum=Jc(3,:)*qv+del;
xp=[f;1]*inv(sum);