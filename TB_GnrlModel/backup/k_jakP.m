function [J,JA] = k_jakP(a, cn)
% Calculates the Jacobi matrix of a specific point one of the links
% Usage: J = k_jakP(flier_object, cn)
% dX = J*dq : where X is the orientation of a specific point
% in the lnk-th link of the system.
% This Jacobian gives the motion of a specific point located in the
% coordinate frame of the lnk-th link. The point is defined
% by the cn index. This index refers to a specific contact point
% defined in advance using function k_addcontact.
%

J = zeros(6,a.N);
JA = zeros(6,1);

% retrieve the actual vector and link number from the object structures
lnk = a.Con.lnr(cn);
svec = a.Con.cvec(:,cn);

% automaticly find the chain that contains the link
% ip - the number of the preceeding link
% ch - the number of the chain containing the link
if(lnk < 4)
	error('''lnk'' is not allowed to be 1,2 or 3.');
end;

for ch = 1:a.ch
	j = find(a.M(ch,:) == lnk);
	if(j)
		ip = a.M(ch,j-1);
		break;
	end
end

gsv = a.A(:,:,lnk)*svec;
J(1:3,:) = a.bet(:,lnk,:) ;
J(4:6,:) = a.alf(:,lnk,:);
J(1:3,lnk) = J(1:3,lnk) + cross( a.ge(:,lnk),gsv );

JA(1:3,1) = a.del(:,lnk) + cross( a.gam(:,lnk-1),cross( gsv,a.gpr(:,ip,ch) ) ) +...
				a.dq(lnk) * cross( cross( a.omg(:,lnk),a.ge(:,lnk) ),gsv ) +...
				cross( a.omg(:,lnk),cross( a.omg(:,lnk),gsv ) );
JA(4:6,1) = a.gam(:,lnk);
