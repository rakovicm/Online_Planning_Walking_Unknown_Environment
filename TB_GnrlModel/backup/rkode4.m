function yn = rkode4(t,y,dy1,dt,eqfun)
% ORDINARY DIFFERENTIAL EQUATION SOLVER
% makes one integration step using fourth-order Runge-Kutta formulae
% Usage:	yn = rkode4(t,y,dy1,dt,eqfun)
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
k3 = dt * feval(eqfun,t+dt/2,y+k2/2);
k4 = dt * feval(eqfun,t+dt,y+k3);
yn = y + (k1+k4)/6 + (k2+k3)/3;
