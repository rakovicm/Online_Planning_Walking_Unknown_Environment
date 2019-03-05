function [lnr, vec, Qtr] = k_getcontact(a, cno)
% Gets the parameters of one of the contact points
% Usage : [lnr, vec, Qtr] = k_getcontact(k_flier, cno)
% k_flier - a k_flier object
% cno - index number of the contact
% lnr - number of the link that makes the contact
% vec - vector in the local frame of the lnr-th link
%		that points to the contact point
% Qtr - homognous transformation to the contact point
%		local frame
%
vec = a.Con.cvec(:,cno);
Qtr = a.Con.cQ(:,:,cno);
lnr = a.Con.lnr(cno);
