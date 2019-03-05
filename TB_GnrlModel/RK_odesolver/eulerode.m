function yn = eulerode(t,y,dy1,dt)
% ORDINARY DIFFERENTIAL EQUATION SOLVER
% makes one integration step using the basic Euler formula
% Usage:	yn = eulerode(t,y,dy1,dt,eqfun)
%	yn	-	result, solution in the (n+1)-th point
%	t	-	independent variable in the n-th point
%	y	-	value of the function in the n-th point
%	dy1	-	derivative of the function in the n-th point
%	dt	-	integration step size
% Note: y, dy1 and yn can be N-element vectors
%

yn = y + dt * dy1;
