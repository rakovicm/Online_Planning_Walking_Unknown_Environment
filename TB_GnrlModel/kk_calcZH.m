function [ZH Zh0]=kk_calcZH(a, ra)
% racuna matrice za izracunavanje ukupnog momenta sistema za tacku ra
[ZH Zh0]=kcalcZH(a.N,a.ch,a.M, a.bs, a.cl,a.rc,a.A, a.alf,a.bet,a.gam,a.del, a.omg,a.ms, a.J, ra);
