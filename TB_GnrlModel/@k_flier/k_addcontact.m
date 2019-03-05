function a = k_addcontact(a, cno, lnr, vec, Qtr)
% Adds a contact point to a selected link of the flier
% Usage : k_flier = k_addcontact(k_flier, cno, lnr, vec, Qtr)
% k_flier - a k_flier object
% cno - index number of the contact
% lnr - number of the link that makes the contact
% vec - vector in the local frame of the lnr-th link
%		that points to the contact point
% Qtr - homognous transformation to the contact point
%		local frame
%
a.Con.cvec(:,cno) = vec;
a.Con.cQ(:,:,cno) = Qtr;
a.Con.lnr(cno) = lnr;
