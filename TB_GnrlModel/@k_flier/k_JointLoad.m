function LF = k_JointLoad(a, ddq, jnt, cps, FS)
% Usage: LF = k_JointLoad(flier, ddq, jnt, cps, FS)
% returns joint load forces and moments relative
% to segment's coordinate frame

% * flier - the flier object
% * ddq - acceleration vector
% * jnt - the slected joint number
% * cp - vecotr of selected (usually all) contact point indices [1 x ncp]
% * FS - the calculated reaction forces matrix [6 x ncp], relative to
%			world cooridinate frame

[k1, k2] = find((a.M==jnt));
j = [];
for k=1:length(k1)
	i = k1(k);
	if(k2(k)>a.bs(i))
		chunk = a.M(i,k2(k)+1:end);
	else
		chunk = a.M(i,a.bs(i)+1:end);
	end
	j = cat(2,j,(chunk(1:length(find(chunk)))));
end
j = cat(2,jnt,j);

lj = [0 0 0]';
rj = lj;
for s=j
	la = [0 0 0]';
	ra = la;
	for b=1:a.N
		la = la + a.ac(:,s,b) * ddq(b);
		ra = ra + a.bc(:,s,b) * ddq(b);
	end
	la = la + a.ac0(:,s) + a.ms(s)*[0 0 -9.81]';
	ra = ra + a.bc0(:,s);
	lj = lj + la;
	rj = rj + ra + cross( a.rc(:,s,jnt) , la );
end

LF1 = [ a.A(:,:,jnt)'*lj; a.A(:,:,jnt)'*rj ];

% addition of the foreign forces
if length(cps) ~= size(FS,2)
	error('The number of forces must match the number of contact points');
end

JF = a.A(:,:,jnt);
LF2 = zeros(6,1);
for k=1:length(cps)
	cs = a.Con.lnr(cps(k)); %%%%% segment na kome se nalazi kontaktna tacka
	if ~isempty(find(j == cs))
		CF = a.A(:,:,cs);
		CsF = a.A(:,:,cs)*a.Con.cQ(:,:,cps(k));
		cpF = zeros(6,1);
		exF = cpF;
		exF(1:3) = FS(1:3,k); %%%% deluje da je u eksternim
		exF(4:6) = FS(4:6,k); %%%% deluje da je u eksternim
		gV = a.rc(:,cs,jnt) + CF*a.Con.cvec(:,cps(k));
		cpF(1:3,1) = JF \ exF(1:3,1);
		cpF(4:6,1) = JF \ ( cross(gV,exF(1:3,1)) + exF(4:6) );
		LF2 = LF2 + cpF;
	end
end

LF = LF1 + LF2;
