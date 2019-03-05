function k_drawobject(a, varargin)
%	Function k_draw object: draws the object as a sphere which radius is such
%							that the contact point is on the surface of the sphere
%	Usage : k_drawobject(k_object_object, optional_arguments);
%	optional 1 : 0 (default) - hold off, 1 - hold on
%	optional 2 : [az, elev], don't change if missing
%	optional 3 : [xmin xmax ymin ymax zmin zmax], don't change if missing
%	the optional arguments must be in the given order to maintain their	meaning
%	k_object_object is an object of type 'k_object'
%

if nargin == 1
		hold off;
else
	if varargin{1}
		hold on;
	else
		hold off;
	end
end

% make the coordinate transformation and relocation
X = a.surf.X;
Y = a.surf.Y;
Z = a.surf.Z;
for i = 1:size(X,1)
	for j = 1:size(X,2)
		V = a.Q * [X(i,j) Y(i,j) Z(i,j)]' + a.q(1:3);
		X(i,j) = V(1);
		Y(i,j) = V(2);
		Z(i,j) = V(3);
	end
end

% DRAW IT OUT
surf(X,Y,Z);

switch nargin
	case 1
		axis equal;
		axis([-2 2 -2 2 -1 3]);
	case 2
		axis equal;
		axis([-2 2 -2 2 -1 3]);
	case 3
		view(varargin{2});
		axis equal;
		axis([-2 2 -2 2 -1 3]);
	case 4
		view(varargin{2});
		axis(varargin{3})
	otherwise
		error('Illegal arguments!');
end
