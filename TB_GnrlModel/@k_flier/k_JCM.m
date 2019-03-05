function [Jcm Acm] =k_JCM(a)
%[Jcm Acm]=k_CentarMase(flier)
% Ulazi:
%   flier - objekat tipa flier
%   
% Rezultati
%   rcm   - Polozaj centra mase sistema u globalnom koordinatnom sistemu
%   vcm   - Brzina centra mase sistema u globalnom koordinatnom sistemu
%   qnew  - Vrednost prvih 6 koordinata flier-a tako da pozicija kontaktne
% Napomena
%   Za izracunavanje rcm potrebno je prethodno pozvati funkciju k_geo,
%       dok je za izracunavanje vcm potrebno prethodno pozvati i k_kine i
%       k_dyn
Jcm=zeros(3,a.N);
Acm=zeros(3,1);

bl=zeros(3,a.N);
for i=1:a.N
    bl(:,:)=a.bet(:,i,:);
    Jcm=Jcm+a.ms(i)*bl;
    Acm=Acm+a.ms(i)*a.del(:,i);
end
M=sum(a.ms(:));
Jcm=Jcm/M;
Acm=Acm/M;
