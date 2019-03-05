function a = k_changeGrav(a,gravacc)
% Usage : k_object = k_changeCP(k_object,gravacc)
% sets up the CoM to contact point vector
% defined in the local coordinate frame
% 	k_object	- 	the object the contact point is defined for
%	gravacc		-	the new value of the gravitational acceleration

if nargin == 2
	a.g = gravacc;
else
	error('Invalid number of aguments!');
end
