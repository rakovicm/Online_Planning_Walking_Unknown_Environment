function yn = rkode2(t,y,dy1,dt,eqfun)
% ORDINARY DIFFERENTIAL EQUATION SOLVER
% makes one integration step using second-order Runge-Kutta formulae
% Usage:	yn = rkode2(t,y,dy1,dt,eqfun)
%	yn	-	result, solution in the (n+1)-th point
%	t	-	independent variable in the n-th point
%	y	-	value of the function in the n-th point
%	dy1	-	derivative of the function in the n-th point
%	dt	-	integration step size
%	eqfun	-	handle to the equation function
% eqfun form:
% dy = eqfun(t,y)
% Note: y, dy1 and yn can be N-element vectors
%

k1 = dt * dy1;
k2 = dt * feval(eqfun,t+dt/2,y+k1/2);
yn = y + k2;
