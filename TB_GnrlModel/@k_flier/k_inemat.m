function [H, h0] = k_inemat(a);
% Calculates the inertia matrix and adjoint vector of the system
%  Usage :  [H, h0] = k_inemat(k_flier_object)
%

% *** result space ***
H = zeros(a.N);
h0 = zeros(a.N,1);		% make room for the result
G = [0 0 -9.81]';		% gravity, fixed (for now)

dc = 0;		% correction
for c_chain = 1:a.ch	% joints are in different chains
	chain = a.M(c_chain,:);		% current chain
	if (c_chain==2)		% only the first chain should include the
		dc = 1;			% starting segment of the particular chain
	end					% the others should jump over it .... dc
	for c_si = (a.bs(c_chain)+dc):a.cl(c_chain)		% cycle though the current chain
		i = chain(c_si);
		for t = chain(1:c_si)
			if(t>3)		% rotational DOF, torque calculation
				h0(t,1) =  h0(t,1) - dot( a.ge(:,t),(cross( a.rc(:,i,t),(a.ac0(:,i)+G*a.ms(i)) )+a.bc0(:,i)) );
				for s = chain(1:c_si)
					H(t,s) = H(t,s) - dot( a.ge(:,t),(a.bc(:,i,s)+cross( a.rc(:,i,t),a.ac(:,i,s) )) );
				end
			else		% linear DOF, force calculation
				h0(t,1) = h0(t,1) - dot( a.ge(:,t),( a.ac0(:,i) + G*a.ms(i) ) );
				for s = chain(1:c_si)
					H(t,s) = H(t,s) - dot( a.ge(:,t),a.ac(:,i,s) );
				end
			end
		end
	end
end
