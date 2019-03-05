function [J, A] = k_Jsl(a,varargin)
% returns J [6x6] and A [6x1]
% ddot_s = J * ddot_X + A  -  Jacobian form
% connects the global coordinates of one of the Mechanisms' LINKS - X
% in the global frame with the local s-frame of the object
% the local s-frame is defined by a vector and a homogenous transform
% with respect to the object's COM.
% Usage : [J, A] = k_Jsl(a[,contact_index])
% a - 'k_object' object instance
% object_index - optional, if not given considered 1
%
% WARNING - k_Jsb should be called previously to set up the object motion
%

switch nargin
case 1
	ci = 1;
case 2
	ci = varargin{1};
otherwise
	error('Invalid number of aguments!');
end

% get the already calculated values from the object
Q = a.Q;

% the sole purpose of this Jacobian is to 
tQ = a.cQ(:,:,ci);
sQ = Q * tQ;
Qx = [sQ' zeros(3); zeros(3) sQ'];
J = Qx;
A = zeros(6,1);
