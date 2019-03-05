function X = k_linkX(a, cn, rorder)
% returns the external coordinates of one point of the link
% Usage : X = k_linkX(k_flier, cn, rorder)
% X - position of the contact point of the link
% k_flier - an instance of the k_flier object
% cn - inedex number of the flier contact
% rorder - a string describing the rotation order - p.e. 'xyz'
%

% Find the corresponding link number
ln = a.Con.lnr(cn);

% Real transformation
RT = a.A(:,:,ln) * a.Con.cQ(:,:,cn);
% Real vector
RV = a.rc(:,ln,1) + a.A(:,:,ln)*a.Con.cvec(:,cn);

switch rorder
case 'xyz'
	As = TRotXYZ(RT);
case 'xzy'
	As = TRotXZY(RT);
case 'yxz'
	As = TRotYXZ(RT);
case 'yzx'
	As = TRotYZX(RT);
case 'zxy'
	As = TRotZXY(RT);
case 'zyx'
	As = TRotZYX(RT);
otherwise
	error(sprintf('The requested order ''%s'' is not supported!',rorder));
end

X = [RV; As'];
