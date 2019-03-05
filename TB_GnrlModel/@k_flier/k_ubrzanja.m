function [qanew]=k_ubrzanja(a, azel, epszel, cn, qa)
% funkcija vraca zeljena ubrzanja centra mase baznog segmenta tako da 
% ubrzanje tacke na segmentu ln odredjena vektorom cvec (u lokalnom
% koordinatnom sistemu) bude azel i ugaono ubrzanje  bude epszel
% pri tome su generalisana ubrzanja qa

ln=a.Con.lnr(cn);
Q=a.A(:,:,ln);

al=zeros(3,a.N);
bt=zeros(3,a.N);
al(:,:)=a.alf(:,ln,:);
bt(:,:)=a.bet(:,ln,:);
g=a.gam(:,ln);
d=a.del(:,ln);
alfbet=[al;bt];
epscmzel=epszel;  %zelenja ugaona brzina centra mase
cvecGl=Q*a.Con.cvec(:,cn);

acmzel=azel-cross(epszel,cvecGl)-cross(a.omg(:,ln),cross(a.omg(:,ln),cvecGl)); 
% zeljeno linearno ubrzanje centra mase

qanew=alfbet(1:6,1:6)\([epscmzel-g;acmzel-d]-alfbet(1:6,7:a.N)*qa(7:a.N));


