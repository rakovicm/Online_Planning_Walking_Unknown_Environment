function [x, v, a, t] = ADprofile(tm, xb, xe, ak, dk, n)
% creates a smooth position - velocity - acceleration profile
% zero start/ zero stop velocity
% Usage:
% 	[x, v, a, t] = ADprofile(tm, xb, xe, ak, dk, n)
% tm  --  duration of the process
% xb  --  start position
% xe  --  end position
% ak  --  acceleration time partition
% dk  --  deceleration time partition
%      !!! (ak + dk) < 1
% n   --  time points, (n-1) intervals
%
% xb, xe can be vector colums of same size
% th length of the generated x, v and a will be the same as
% the length of the vectors xb and xe

ne = size(xe,1);	% vector length
dx = xe - xb;
t = linspace(0, tm, n);
ta = tm * ak;
td = tm * dk;
tdl = tm - td;
omgmax = (2*dx) / (tm * (2 - ak - dk) );
alfa = omgmax / (ta);
alfd = omgmax / (td);
omgka = omgmax / ta;
omgkd = omgmax / td;
linb = 0.5 * alfa * ta^2;
line = dx - 0.5 * alfd * td^2;
link = (line - linb)/(tm - ta - td);

x = zeros(ne,n);
v = x;
a = x;

for i = 1:n
	tt = t(i);
	if tt < ta
		x(:,i) = 0.5 * alfa * tt^2 + xb;
		v(:,i) = tt * omgka;
		a(:,i) = alfa;
	elseif tt < tdl
		x(:,i) = linb + link * (tt - ta) + xb;
		v(:,i) = omgmax;
%		a(i) = 0	% zero anyway
	else
		x(:,i) = dx - 0.5 * alfd * (tm - tt)^2 + xb;
		v(:,i) = omgmax - (tt - tdl)*omgkd;
		a(:,i) = -alfd;
	end
end
