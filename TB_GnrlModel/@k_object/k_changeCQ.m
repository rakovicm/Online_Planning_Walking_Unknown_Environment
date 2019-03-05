function a = k_changeCQ(a,trf,varargin)
% Usage : k_object = k_changeCQ(k_object,CQ[,contact_index])
% sets up the coordinate transformation matrix
% form the local coordinate frame to the
% contact potint i.e. s coordinate frame
% 	k_object	- 	the object the contact point is defined for
%	CQ			- 	a transformation matrix form the object local
%					coordinate frame to the s-cordinate frame
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

a.cQ(:,:,i) = trf;
