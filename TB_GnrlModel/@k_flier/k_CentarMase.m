function [rcm Vcm] =k_CentarMase(a)
%[rcm vcm]=k_CentarMase(flier)
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
rcm=[0;0;0];
Vcm=[0;0;0];
for i=1:a.N
    rcm=rcm+a.rc(:,i,1)*a.ms(i);
    Vcm=Vcm+a.v(:,i,1)*a.ms(i);
end
M=sum(a.ms(:));
rcm=rcm/M;
Vcm=Vcm/M;

