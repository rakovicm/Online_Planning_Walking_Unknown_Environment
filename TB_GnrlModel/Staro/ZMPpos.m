function [Zx, Zy] = ZMPpos(F)
% Finds the position of the ZMP in the contact plane
% Contact plane : the plane formed by the axes Sx_Sy of the
% functional coordinate system. The Fx, Fy and M3 components
% of the reaction generalized forces are assumed to be
% compenzated by friction.
% Usage : [Zx, Zy] = ZMPpos(F)
%

% vector F : Fx Fy Fz M1 M2 M3. Which M is which is defined
% by the contact object
Fz = [0 0 F(3)]';
Mr = [F(4:5); 0];
ri = cross(Mr,Fz);
ri = -ri / vmod(ri);
R = ri * ( vmod(Mr) / abs(F(3)) );
Zx = R(1);
Zy = R(2);
