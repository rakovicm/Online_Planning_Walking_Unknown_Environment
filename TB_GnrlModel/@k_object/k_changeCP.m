function a = k_changeCP(a,cpoint,varargin)
% Usage : k_object = k_changeCP(k_object,cpoint[,contact_index])
% sets up the CoM to contact point vector
% defined in the local coordinate frame
% 	k_object	- 	the object the contact point is defined for
%	cpoint		- 	a vector in the local coordinate frame that
%					defines the contact point
%	contact_index	-	optional, if not given, considered 1
%

switch nargin
case 2
	i = 1;
case 3
	i = varargin{1};
otherwise
	error('Invalid number of aguments!');
end

a.contact(:,i) = cpoint;
[a.surf.X,a.surf.Y,a.surf.Z] = GenSphr(vmod(cpoint),10);
