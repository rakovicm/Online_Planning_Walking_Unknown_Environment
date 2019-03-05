function a = kk_kine(a, dq)
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

%%%
a.dq = dq;	% update internal data

[a.alf a.bet a.gam a.del a.omg a.v]=kkine(a.N,a.ch,a.M,a.bs,a.cl,a.ge,a.dq,a.gpr,a.gdr);
%%%%%%%%%%%% treba odrediti 
%%% omega, ve  alfa beta gama delta
% first link parameters, set separately in advance
%%a.alf=alf;

