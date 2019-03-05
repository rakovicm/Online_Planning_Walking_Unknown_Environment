function [J,JA] = k_jakL(a, lnk)
% Calculates the Jacobi matrix of one of the links
% Usage: J = k_jakL(flier_object, lnk)
% dX = J*dq : where X is the orientation of the lnk-th link of the system
% this Jacobian gives the motion of the link's COM
%

J = zeros(6,a.N);

J(1:3,:) = a.bet(:,lnk,:);
J(4:6,:) = a.alf(:,lnk,:);

JA = [ a.del(:,lnk); a.gam(:,lnk) ];
