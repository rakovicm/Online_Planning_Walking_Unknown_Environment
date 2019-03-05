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

% helper vectors
rs = eye(3);
lt = a.cQ;

% THEORY : ddot_X = Ja*ddot_s + Aa =>
%		=> ddot_s = inv(Ja)*ddot_X - inv(Ja)*Aa
% therefore the returned result should be
% J = inv(Ja) since we calculate Ja (even if the variable is called J)
% A = inv(Ja)*Aa since we caclulate Aa (even if we call it A)

% get the already calculated values from the object
rs = a.Q;
omg = a.omg;
e = a.e;

% Jacobian form 6x6 matrix
rp = rs * a.contact(:,ci);
J = [ rs*[1 0 0]' rs*[0 1 0]' rs*[0 0 1]' cross(e(:,1),rp) cross(e(:,2),rp) cross(e(:,3),rp);
		zeros(3) [e(:,1) e(:,2) e(:,3)] ];
J = inv(J);

% Adjoint vector
A = zeros(6,1);
A(1:3,1) = a.dq(6)*cross( cross( omg(:,2),e(:,2) ),rp ) + ...
	cross( omg(:,3),cross( omg(:,3),rp ) );
A(4:6,1) = a.dq(5)*cross(omg(:,1),e(:,2)) + a.dq(6)*cross(omg(:,2),e(:,3));
A = -J*A;

% correction because of the different orientation of the s-frame
tQ = rs * a.cQ(:,:,ci);
Qx = [tQ' zeros(3); zeros(3) tQ'];
J = Qx * J;
A = Qx * A;
