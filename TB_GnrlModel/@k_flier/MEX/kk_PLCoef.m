function [Pc Lc]=kk_PLCoef(a,r0)
% [Pc Lc]=kk_PLCoef(a,r)
% Racuna koeficiente Pc Lc koji sluze za racunanje kolicine kretanja
% i momenta kolicine kretanja za tacku r
% 
%    Kolicina kretanja:
%         P=Pc*q
%         
%    Moment kolicine kretanja:
%         L=Lc*q

[Pc Lc]=kPLCoef(a.N,a.ch,a.M,a.bs,a.cl,a.rc(:,:,1),a.ac,a.bc, r0 );
