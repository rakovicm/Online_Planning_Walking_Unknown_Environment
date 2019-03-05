function [a, J, A] = k_Jso(a,varargin)
% returns J [6x6] and A [6x1]
% ddot_s = J * ddot_X + A  -  Jacobian form
% connects the global coordinates of the object - X in global frame
% with the local s-frame of the object
% THE "CONTACT POINT" IS ALWAYS IN THE COM OF THE OBJECT
% THE CONTACT COORDINATE FRAME COINCIDES WITH THE OBJECT'S CF
% Usage : [a, J, A] = k_Jsb(a, [dq])
% a - 'k_object' object instance
% dq - first derivatives of the coordinates
% !!! It sets up the first derivatives of the object coordinates only if
%	the argument  dq  is given
% WARNING - k_objectX should be called first
%

if nargin > 2
	error('Invalid number of aguments!');
end

% get the already calculated values from the object
Q = a.Q;
e = a.e;

if nargin == 2
		a.dq = varargin{1};
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
end

% Jacobian form 6x6 matrix
J = [	[1 0 0]' [0 1 0]' [0 0 1]' zeros(3,3)
		zeros(3,3) e(:,1) e(:,2) e(:,3) ];

% Adjoint vector
A = zeros(6,1);
A(1:3,1) = zeros(3,1);
A(4:6,1) = a.dq(5)*cross( a.omg(:,4),e(:,2) ) + ...
	a.dq(5)*cross( a.omg(:,4),e(:,2) );
