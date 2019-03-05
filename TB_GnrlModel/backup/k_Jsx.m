function [a, J, A] = k_Jsx(a,dq,varargin)
% returns J [6x6] and A [6x1]
% ddot_s = J * ddot_X + A  -  Jacobian form
% connects the global coordinates of the object - X in global frame
% with the local s-frame of the object
% the local s-frame is defined by a vector and a homogenous transform
% with respect to the object's COM.
% Usage : [a, J, A] = k_Jsb(a, dq)
% a - 'k_object' object instance
% dq - first derivatives of the coordinates
% object_index - optional, if not given considered 1
% !!! It sets up the first derivatives of the object coordinates only if
%	the contact_index equals 1 !!!!
% WARNING - k_objectX should be called first
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
	% fill up the object structures
	a.dq = dq;
end

% THEORY : ddot_X = Ja*ddot_s + Aa =>
%		=> ddot_s = inv(Ja)*ddot_X - inv(Ja)*Aa
% therefore the returned result should be
% J = inv(Ja) since we calculate Ja (even if the variable is called J)
% A = -inv(Ja)*Aa since we caclulate Aa (even if we call it A)

% get the already calculated values from the object
rs = a.Q;
e = a.e;

% calculate the kinematics koefficients
a.omg(:,1:3) = zeros(3,3);
a.alf(:,1:3) = zeros(3,3);	% no angular motion depending on dq(1:3)
a.gam = zeros(3,1);	% no angular motion depending on dq(1:3)
a.v(:,1:3) = [ a.dq(1) 0 0; a.dq(1) a.dq(2) 0; a.dq(1) a.dq(2) a.dq(3) ]';	% linear velocity of the COM
a.bet(:,1:3) = [1 0 0; 0 1 0; 0 0 1]';
a.del = zeros(3,1);	% these linear components would depend on rotational compopnents
							% that are equal to zero
for i = 4:6
	a.alf(:,i) = a.e(:,i-3);
	a.omg(:,i) = a.omg(:,i-1) + a.dq(i)*a.e(:,i-3);
	a.gam = a.gam + a.dq(i)*cross(a.omg(:,i-1),a.e(:,i-3));
	a.v(:,i) = a.v(i-1);
	a.bet(:,i) = [0 0 0]';
end

% Jacobian form 6x6 matrix
rp = rs * a.contact(:,ci);
J = [ rs*[1 0 0]' rs*[0 1 0]' rs*[0 0 1]' cross(e(:,1),rp) cross(e(:,2),rp) cross(e(:,3),rp);
		zeros(3) [e(:,1) e(:,2) e(:,3)] ];

% Adjoint vector
A = zeros(6,1);
A(1:3,1) = a.dq(6)*cross( cross( a.omg(:,5),e(:,2) ),rp ) + ...
	cross( a.omg(:,6),cross( a.omg(:,6),rp ) );
A(4:6,1) = a.dq(5)*cross(a.omg(:,4),e(:,2)) + a.dq(6)*cross(a.omg(:,5),e(:,3));
