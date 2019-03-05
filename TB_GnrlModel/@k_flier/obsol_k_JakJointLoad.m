function J = k_JakJointLoad(a, cp, jnt)
% Usage: J = k_JakJointLoad(a, cp, jnt)
% returns a matrix that can be used to calcluate
% load exerted to a joint, when contact force in a
% in a contact point is known. 6 dimensional force
% vector is needed
% * cp is the contact point
% * jnt is the slected joint number
% * a is a flier object
% the result is used as follows
% Fj = J' * Fc
% where Fc is the contact force vector and Fj,
% is the joint load vector, both are 6-dimensional

% contact point link
i = a.Con.lnr(cp);

ACP = a.A(:,:,i)*a.Con.cQ(:,:,cp);

eg = a.A(:,:,jnt)*[1 0 0]';
ag = a.A(:,:,jnt)*[0 1 0]';
bg = a.A(:,:,jnt)*[0 0 1]';

J = zeros(6);

rcp = a.rc(:,i,jnt) + a.A(:,:,i)*a.Con.cvec(:,cp);
J(1:3,1:3) = ACP;
J(4:6,4:6) = ACP;
J(1:3,4) = ACP*cross(eg,rcp);
J(1:3,5) = ACP*cross(ag,rcp);
J(1:3,6) = ACP*cross(bg,rcp);
