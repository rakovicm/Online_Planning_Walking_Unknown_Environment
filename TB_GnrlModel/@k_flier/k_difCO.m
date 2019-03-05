function [dif] = k_difCO(obf, cif, obo, cio)
% 
% Usage: dif = k_difCO(obf, cif, obo, cio)
%
% returns the spatial vector between the two contact points
% One of them is on the flier (obf/cif),
% the other is on a contact object (obo/cio)
% IMPORANT !!! the resulting vector is in s-coordinates (of 'obo')
%
% obs - the object carrying the s coordinate frame
% obx - the other object
% cis, cix - the corresponding contact indices
%

ocontacts = obo.contact;
ocoords = obo.q;
ocQ = obo.cQ;
Xs = obo.Q * ocontacts(:,cio) + ocoords(1:3,1);
linn = obf.Con.lnr(cif);
Xx = obf.rc(:,linn,1) + obf.A(:,:,linn) * obf.Con.cvec(:,cif);
dif = Xx - Xs;
dif = ( obo.Q * ocQ(:,:,cio) ) \ dif;
