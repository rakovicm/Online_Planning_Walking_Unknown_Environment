function a = k_dyn(a)
% calculates the dynamics coeffcients of the system
%  Usage : k_flier_object = k_dyn(k_flier_object)
% thos function calculates a bunch of coefficients needed for later calculation
% of different dynamics-related quantities:
%	* inertia matrices H and h0 <-> k_inemat
%	* inertial forces and torques <-> k_infor
% according to this, MUST be called before these functions!
%

% *** intermediate coefficients ***
T = zeros(3);
lamb = zeros(3,1);

% *** calculation of the a, b, a0 and b0 coefficients
for c_chain = 1:a.ch	% joints are in different chains
	chain = a.M(c_chain,:);		% current chain
	for c_si = a.bs(c_chain):a.cl(c_chain)		% cycle though the current chain
		i = chain(c_si);
		% calculate the ac parameters for the dynamics (linear motion component)
		a.ac(:,i,chain(1:c_si)) =  -a.ms(i) * a.bet(:,i,chain(1:c_si));
		a.ac0(:,i) = -a.ms(i) * a.del(:,i);
		% calculate the bc parameters for the dynamics (rotational motion component)
		T = a.A(:,:,i) * a.J(:,:,i) * a.A(:,:,i)';
		lamb = cross( a.omg(:,i),(a.A(:,:,i)*a.J(:,:,i)*a.A(:,:,i)'*a.omg(:,i)) );
		for j = chain(1:c_si)
			a.bc(:,i,j) = -T*a.alf(:,i,j);
		end
		a.bc0(:,i) = -T*a.gam(:,i) - lamb;
	end
end
