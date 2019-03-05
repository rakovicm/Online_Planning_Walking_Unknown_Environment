function a = k_geo(a, q)
% calculates the geometry and configuration of the mechanism
% Usage :  k_flier_object = k_geo( k_flier_object, coordinate vector q )
%
if( length(q) ~= (a.N) )
	error('Coordinate vector length illegal error!');
end

a.q = q;	% internal angles (1:3 base translation, 4:6 - base orient., 7.. internal ang.)

a.A(:,:,1) = eye(3);	% transformation to the first (dummy) segment must be unity
a.rc = zeros(3, a.N, a.N);	% clear the rc matrices
% cycle throug the chains and apply Rodrigue's famous formula
for c_chain = 1:a.ch	% joints are in different chains
	chain = a.M(c_chain,:);
	TT = zeros(3,3);
	i = chain(a.bs(c_chain));
	a.gpr(:,i,c_chain) = a.A(:,:,i) * a.lpr(:,i,c_chain);	
%	disp(sprintf('%d,%d',i,c_chain));	% debug
	for c_si = (a.bs(c_chain)+1):a.cl(c_chain)		% cycle though the current chain
		i = chain(c_si);
		ip = chain(c_si-1);
		if( i>6 )
			for k=1:3		% Rodrigue's formula
				TT(:,k) = a.As0(:,k,i)*cos(a.q(i)) +...
					(1-cos(a.q(i))) * dot( a.lne(:,ip,c_chain),a.As0(:,k,i) ) * a.lne(:,ip,c_chain) +...
					cross( a.lne(:,ip,c_chain),a.As0(:,k,i) )*sin(a.q(i));
			end
			a.A(:,:,i) = a.A(:,:,ip) * TT;
			a.ge(:,i) = a.A(:,:,i) * a.le(:,i);
			a.gdr(:,i) = a.A(:,:,i) * a.ldr(:,i);
		else
			switch i
			case 2
				a.gdr(:,1) = [a.q(1) 0 0]';
				a.rc(:,1,1) = a.gdr(:,1);
				a.A(:,:,1) = eye(3);
				a.ge(:,1) = a.A(:,:,1) * a.le(:,1);
				a.gdr(:,2) = [0 a.q(2) 0]';
				a.A(:,:,2) = eye(3);
            case 3
				a.gdr(:,3) = [0 0 a.q(3)]';
				a.A(:,:,3) = eye(3);
			otherwise
				switch a.ors(i-3)
				case 'x'
					TT = rotx(a.q(i));
                case 'y'
					TT = roty(a.q(i));
                case 'z'
					TT = rotz(a.q(i));
                end
				a.A(:,:,i) = a.A(:,:,ip) * TT;
			end
			a.ge(:,i) = a.A(:,:,i) * a.le(:,i);
		end
% --- calculate the gpr vector for the current link, current chain ---
		a.gpr(:,i,c_chain) = a.A(:,:,i) * a.lpr(:,i,c_chain);
%		disp(sprintf('%d,%d',i,c_chain));	% debug
%		disp(a.gpr(:,i,c_chain));			% debug
% --- calculate the JOINT to COM vectors rc (related to i-th link) ---
		a.rc(:,i,i) = a.gdr(:,i);
		AD = - a.gpr(:,ip,c_chain) + a.gdr(:,i);
		for k = chain( 1:(c_si-1) )
			a.rc(:,i,k) = a.rc(:,ip,k) + AD;
		end
	end
	for i = chain(6:(a.bs(c_chain)-1))		% fix gpr vectors (links preceding the fork)
		a.gpr(:,i,c_chain) = a.A(:,:,i) * a.lpr(:,i,c_chain);
%		disp(sprintf('%d,%d',i,c_chain));	% debug
	end
end
