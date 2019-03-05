function [qvnew]=k_brzine(a, vzel, wzel, cn)
% funkcija vraca zeljene brzine centra mase baznog segmenta tako da 
% brzina tacke na segment ln odredjena vektorom cvec (u lokalnom
% koordinatnom sistemu) bude vzel i ugaona brzine bude wzel
ln=a.Con.lnr(cn);
Q=a.A(:,:,ln);


al=zeros(3,a.N);
bt=zeros(3,a.N);
al(:,:)=a.alf(:,ln,:);
bt(:,:)=a.bet(:,ln,:);
alfbet=[al;bt];
wcmzel=wzel;  %zelenja ugaona brzina centra mase
vcmzel=vzel-cross(wzel,Q*a.Con.cvec(:,cn)); % zeljena linearna brzina centra mase

qvnew=alfbet(1:6,1:6)\([wcmzel;vcmzel]-alfbet(1:6,7:a.N)*a.dq(7:a.N));


