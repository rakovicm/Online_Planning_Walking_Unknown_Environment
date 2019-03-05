function [nq, ndq, nddq] = AccelProfile(q0,dq0,dqe,dt,tn,n)
% makes a speed transition

a = (dqe - dq0)/((tn-1)*dt);
r = size(a,1);
nq = zeros(r,n);
ndq = zeros(r,n);
nddq = zeros(r,n);

for i=1:n
	if i<=tn
		nddq(:,i) = a;
		ndq(:,i) = dq0 + a*(i-1)*dt;
		nq(:,i) = q0 + dq0*(i-1)*dt + 0.5*a*((i-1)*dt)^2;
	else
		ndq(:,i) = dqe;
		nq(:,i) = nq(:,tn) + dqe*(i-tn)*dt;
	end
end
