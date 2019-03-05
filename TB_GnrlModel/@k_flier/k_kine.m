function a = k_kine(a, dq)
% calculates the kinematics of the system ...
%
% Usage: k_kine_object = k_kine(k_kine_object, dq)
%
% operates on a k_flier object
% WARNING - for this method to operate properly, the 'k_geo' method MUST be
% WARNING - called previously the k_geo method inserts the adequate internal
% WARNING - coordinates 'q'
% WARNING - through k_kine, the first derivatives only will be inserted 'dq'
%
if( length(dq) ~= (a.N) )
	error('Coordinate derivative vector length illegal error!');
end

a.dq = dq;	% update internal data

% first link parameters, set separately in advance
a.omg(:,1:3) = zeros(3,3);	% set the rotational motion koefficients of links 1 - 3 in advance
a.alf(:,1,1) = [ 0 0 0 ]';
a.alf(:,2,1:2) = [ 0 0 0; 0 0 0 ]';
a.alf(:,3,1:3) = [ 0 0 0; 0 0 0; 0 0 0 ]';
a.gam(:,1:3) = [ 0 0 0; 0 0 0; 0 0 0 ]';

a.v(:,1) = [dq(1) 0 0]';	% set the linear motion koefficients of links 1 - 3 in advance
a.v(:,2) = a.v(:,1) + [0 dq(2) 0]';
a.v(:,3) = a.v(:,2) + [0 0 dq(3)]';
a.bet(:,1,1) = [ 1 0 0 ]';
a.bet(:,2,1:2) = [ 1 0 0; 0 1 0 ]';
a.bet(:,3,1:3) = [ 1 0 0; 0 1 0; 0 0 1 ]';
a.del(:,1:3) = [ 0 0 0; 0 0 0; 0 0 0 ]';

for c_chain = 1:a.ch	% joints are in different chains
%(params for 1)	% base segment should be processed in advance
	chain = a.M(c_chain,:);		% current chain
	for c_si = (a.bs(c_chain)+1):a.cl(c_chain)		% cycle though the current chain
		i = chain(c_si);
		ip = chain(c_si-1);
		if(i>3)
% --- link angular velocities and acceleration ---
			a.alf(:,i,chain(1:(c_si-1))) = a.alf(:,ip,chain(1:(c_si-1)));	% i-th link pre i-th alfa coeffs
			a.alf(:,i,i) = a.ge(:,i);	% i-th link last alfa coeffcient
			a.omg(:,i) = a.omg(:,ip) + dq(i)*a.ge(:,i);	% i-th link angular velocity
			a.gam(:,i) = a.gam(:,ip) + dq(i)*cross(a.omg(:,ip),a.ge(:,i));	% i-th link gamma (alfa0) coefficient
% --- link linear velocities and accelerations ---
			a.v(:,i) = a.v(:,ip) - cross( a.omg(:,ip),a.gpr(:,ip,c_chain) ) + cross( a.omg(:,i),a.gdr(:,i) );
			difR = a.gdr(:,i)-a.gpr(:,ip,c_chain);
			for j=chain(1:(c_si-1))
				a.bet(:,i,j) = a.bet(:,ip,j) + cross( a.alf(:,ip,j),difR );
			end
			a.bet(:,i,i) = cross( a.ge(:,i),a.gdr(:,i) );
				a.del(:,i) = a.del(:,ip) + cross( a.gam(:,ip),difR ) + ...
					dq(i)*cross( cross( a.omg(:,ip),a.ge(:,i) ),a.gdr(:,i) ) - ...
					cross( a.omg(:,ip),cross( a.omg(:,ip),a.gpr(:,ip,c_chain) ) ) + ...
					cross( a.omg(:,i),cross( a.omg(:,i),a.gdr(:,i) ) );
		end
	end
end
