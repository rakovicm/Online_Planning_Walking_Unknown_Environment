function a = k_object(varargin)
% CONSTRUCTOR of the k_object object
% Usage :
% 1. a = k_object(k_object)	- copy constructor
% 2. a = k_object(string)
% 3. a = k_object(string, mass, inertia_matrix)
%	string - e.g. 'xyz' - defines the rotation order
%	creates a new default object
%

switch(nargin)
	case 1
		if( isa(varargin{1},'k_object') )
%		copy constructor
			a = varargin{1};
		elseif( isa(varargin{1},'char') )
%		create default object with defined rot. order
			a = CreateDefault(varargin{1}, 0, zeros(3));
			a = class(a, 'k_object');
		else
			error('The only argument has to be either k_object or 3-character string');
		end
	case 3
		if( isa(varargin{1}, 'char') & isa(varargin{2}, 'double') & ...
				isa(varargin{3}, 'double') )
			a = CreateDefault(varargin{1}, varargin{2}, varargin{3});
			a = class(a, 'k_object');
		else
			error('The 3 arguments has to be: rots string, mass, inertia matrix');
		end
	otherwise
		error('Wring arguments');
end;

function a = CreateDefault(rots, m, J, idstr)
a.q = zeros(6,1);	% coordinates
a.dq = zeros(6,1);	% first derivatives
a.rots = rots;	% rotation order
a.m = m;
a.J = J;
a.contact = [0 0 0]';	% COM to contatct point vector
a.cQ = eye(3);	% coord. transformation to the functional frame S
a.Q = eye(3);	% homogenous tansform to the link CoM
a.omg = zeros(3,6);	% angular velocities of the object
a.v = zeros(3,6);	% linear velocity params of the object CoM
a.e = zeros(3,3);	% object rotation axes
a.alf = zeros(3, 6);
a.bet = zeros(3, 6);
a.gam = zeros(3,1);
a.del = zeros(3,1);
a.ac = zeros(3,6);
a.bc = zeros(3,6);
a.ac0 = zeros(3,1);
a.bc0 = zeros(3,1);
a.rc = zeros(3,6);
a.surf = struct('X',[],'Y',[],'Z',[]);		% for object drawing
a.g = 9.81;
