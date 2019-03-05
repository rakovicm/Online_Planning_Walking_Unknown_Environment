function [a, X] = k_objectX(a, q, varargin)
% returns the object's contact point -
% (the s coordinate frame) in the global
% coordinate frame. Negation is implicitly
% included so:
% s_dif = k_objectX(k_object, q[, contact_index]) + k_linkX(k_flier, contact_ndx)
% 
% Usage:	[k_object, X] = k_objectX(k_object, q[, contact_index])
% q - the coordinate vector - will be stored in the object
% contact_index - optional. If not specified 1 is assumed
% the object coordinates are updated only if contact_index equals 1 or isn't specfied
%

switch nargin
case 2
	ci = 1;
case 3
	ci = varargin{1};
otherwise
	error('Invalid number of aguments!');
end

if ci==1
	% the object stores the coordinates
	a.q = q;
end

% helper vectors
rs = eye(3);
lt = a.cQ(:,:,ci);
for i=1:3
	switch(a.rots(i))	% calculate rotation axes
		case 'x'
			e(:,i) = rs * lt * [1 0 0]';
			rs = rs * rotx(a.q(i+3));
		case 'y'
			e(:,i) = rs * lt * [0 1 0]';
			rs = rs * roty(a.q(i+3));
		case 'z'
			e(:,i) = rs * lt * [0 0 1]';
			rs = rs * rotz(a.q(i+3));
	end
end

a.Q = rs;
a.e = e;

a.rc = zeros(3,6);
a.rc(:,1:3) = [a.q(1) a.q(2) a.q(3); 0 a.q(2) a.q(3); 0 0 a.q(3)]';

switch a.rots
case 'xyz'
	As = TRotXYZ(rs*lt);
case 'xzy'
	As = TRotXZY(rs*lt);
case 'yxz'
	As = TRotYXZ(rs*lt);
case 'yzx'
	As = TRotYZX(rs*lt);
case 'zxy'
	As = TRotZXY(rs*lt);
case 'zyx'
	As = TRotZYX(rs*lt);
otherwise
	error(sprintf('The requested order ''%s'' is not supported!',rorder));
end

% the sign of the result has to be changed in odrer
% to be   s = k)objectX + k_linkX (instead of - inbetween)
X = [ rs*a.contact(:,ci)+a.q(1:3,1); As' ];
