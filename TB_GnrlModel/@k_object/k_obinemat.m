function [OH,Oh0] = k_obinemat(a)

% CALCULATION OF THE DYNAMICS PARAMETERS
% because only the sixth segment is not fixious
a.ac = -a.m * a.bet;
a.ac0 = -a.m * a.del;
a.bc = zeros(3,6);
T = a.Q * a.J * a.Q';
lamb = cross( a.omg(:,6),T*a.omg(:,6) );
for j = 1:6
	a.bc(:,j) = -T*a.alf(:,j);
end
a.bc0 = -T*a.gam - lamb;

% CALCULATION OF THE INERTIA MATRICES AND ADJOINT VECTOR - OH, Oh0
le = [1 0 0; 0 1 0; 0 0 1]';
OH = zeros(6,6);
Oh0 = zeros(6,1);
G = [0 0 -a.g]';	% fixed gravitational acceleration
for j = 1:6
	if(j>3)	% translational motion
		Oh0(j,1) =  Oh0(j,1) - dot( a.e(:,j-3),(cross( a.rc(:,j),(a.ac0+G*a.m) )+a.bc0) );
		for s = 1:6
			OH(j,s) = OH(j,s) - dot( a.e(:,j-3),(a.bc(:,s)+cross( a.rc(:,j),a.ac(:,s) )) );
		end
	else	% rotational motion
		Oh0(j,1) = Oh0(j,1) - dot( le(:,j),( a.ac0 + G*a.m ) );
		for s = 1:6
			OH(j,s) = OH(j,s) - dot( le(:,j),a.ac(:,s) );
		end
	end
end
