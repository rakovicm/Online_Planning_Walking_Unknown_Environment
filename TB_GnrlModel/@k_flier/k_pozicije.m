function [qnew]=k_pozicije(a, Rzel, Azel, cn)
%[qnew]=k_pozicije(flier, Rzel, Qzel, Cn)
% Ulazi:
%   flier - objekat tipa flier
%   Rzel  - Zeljena pozicija kontaktne tacke
%   Qzel  - Zeljena matrica transformacije izmedju koordinatnog sistema
%               kontaktne tacke i globalnog koordinatnog sistema
%   Cn    - Indeks zeljene kontaktne tacke
%
% Rezultat
%   qnew  - Vrednost prvih 6 koordinata flier-a tako da pozicija kontaktne
%               tacke Cn bude Rzel i da matrica rotacija izmedju koordinatnog 
%               sistema kontaktne tacke i globalnog bude Qzel
%
% Napomena
%   prethodno je potrebno pozvati funkciju k_geo
ln=a.Con.lnr(cn);
A1ln=a.A(:,:,ln);
A16 =a.A(:,:,6);
A6ln=A16'*A1ln;

R6 = a.rc(:,6);
Rln= a.rc(:,ln);
R6ln=Rln-R6;
R6lnloc=A16'*R6ln;

A16=Azel/(A6ln*a.Con.cQ(:,:,cn));
A1ln=Azel/a.Con.cQ(:,:,cn);
R1ln=Rzel-A1ln*a.Con.cvec(:,cn);

switch a.ors
case 'xyz'
	qnew(4:6) = TRotXYZ(A16);
case 'xzy'
	qnew(4:6) = TRotXZY(A16);
case 'yxz'
	qnew(4:6) = TRotYXZ(A16);
case 'yzx'
	qnew(4:6) = TRotYZX(A16);
case 'zxy'
	qnew(4:6) = TRotZXY(A16);
case 'zyx'
	qnew(4:6) = TRotZYX(A16);
otherwise
	error('The requested order ''%s'' is not supported!',a.ors);
end

qnew(1:3)=R1ln-A16*R6lnloc;

