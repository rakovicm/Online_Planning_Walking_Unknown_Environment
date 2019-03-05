function a = k_flier(varargin)
%
% CONSTRUCTS THE K_flier OBJECT
%
% a = k_flier( k_flier_object )		-	copy constructor
%
% a = k_flier( links structure )	-	constructs a new object from a
%												'links' structure data
% links structure:
%	.M - [M x k] matrix, defining different chains (M - number of chains, k - 
%												length of the lingest chain)
%	.ors - char([3 x 1]) x,y,z - defines the 3 base rotation angle axes sequence
%	.linkgeo[N]		- geometry structure array, N - number of links
%		.le		- [3 x 1] joint axis vector ( j-th axis in j-th link frame )
%		.lne	- [3 x M] next joint axis vector ( (j+1)-th axis in j-th link frame )
%							M vectors,one for each chain
%		.ldr	- [3 x 1] r_jj vector, local frame
%		.lpr	- [3 x M] r_j_j+1 vectors, local frame, M vectors, one for each chain
%	.linkdyn[N]	- dynamic parameters of the N links, structure array
%		.m		- segment masses
%		.J		- [3 x 3] - mimnets of inertia matrices
%
if( nargin == 1)
	if( isa(varargin{1}, 'k_flier'))
		a = varargin{1};
	elseif( isa(varargin{1}, 'struct') )
		a = Decompose( varargin{1} );
		a = AnalyseAssemble( a );
% - coordinates
		a.q = zeros(a.N, 1);		% internal coordinates (angles)
		a.dq = zeros(a.N, 1);		% first derivatives of internal coord.
		a.ddq = zeros(a.N, 1);	% second derivatives of internal coord.
% - transformation matrices
		a.A = zeros(3,3,a.N);	% transf. matrix from global frame to link j
% - global link geometry 
		a.ge = zeros(3, a.N);	% joint axes vectors
		a.gdr = zeros(3, a.N);	% global dr vectors of the links
		a.gpr = zeros(3, a.N, a.ch);	% global pr vectors of links for each chain
		a.rc = zeros(3, a.N, a.N);		% global CM position vectors relative to joints
%							rc(:,i,j) : vector form j-th joint to i-th link mass center
		a.Con = struct('cvec',[0 0 0]','cQ',zeros(3),'lnr',0);
% - kinematics, koefficients and angular/linear velocities/accelerations
		a.alf = zeros(3, a.N, a.N);	% alfa coefficients
		a.bet = zeros(3, a.N, a.N);	% beta coefficients
		a.gam = zeros(3, a.N);	% gamma koefficients (alfa null)
		a.del = zeros(3, a.N);	% delta koefficients (beta null)
		a.omg = zeros(3, a.N);	% angular velocities of the links
		a.v = zeros(3, a.N);	% linear velocities of the links
% - dynamics, for now, no parameters are stored
		a.ac = zeros(3,a.N,a.N);
		a.bc = zeros(3,a.N,a.N);
		a.ac0 = zeros(3,a.N);
		a.bc0 = zeros(3,a.N);
% - finally, create the k_flier object
		a = class(a, 'k_flier');
        
	else
		error('Input argument must be of class k_flier or a structure ''link''!\nNO DEFAULT OBJECT SUPPORTED');
    end
 elseif (nargin ==2)
    if ( isa(varargin{1}, 'struct') ) && (varargin{2}==1)
        a = Decompose( varargin{1} );
		a = AnalyseAssembleNK( a );
% - coordinates
		a.q = zeros(a.N, 1);		% internal coordinates (angles)
		a.dq = zeros(a.N, 1);		% first derivatives of internal coord.
		a.ddq = zeros(a.N, 1);	% second derivatives of internal coord.
% - transformation matrices
		a.A = zeros(3,3,a.N);	% transf. matrix from global frame to link j
% - global link geometry 
		a.ge = zeros(3, a.N);	% joint axes vectors
		a.gdr = zeros(3, a.N);	% global dr vectors of the links
		a.gpr = zeros(3, a.N, a.ch);	% global pr vectors of links for each chain
		a.rc = zeros(3, a.N, a.N);		% global CM position vectors relative to joints
%							rc(:,i,j) : vector form j-th joint to i-th link mass center
		a.Con = struct('cvec',[0 0 0]','cQ',zeros(3),'lnr',0);
% - kinematics, koefficients and angular/linear velocities/accelerations
		a.alf = zeros(3, a.N, a.N);	% alfa coefficients
		a.bet = zeros(3, a.N, a.N);	% beta coefficients
		a.gam = zeros(3, a.N);	% gamma koefficients (alfa null)
		a.del = zeros(3, a.N);	% delta koefficients (beta null)
		a.omg = zeros(3, a.N);	% angular velocities of the links
		a.v = zeros(3, a.N);	% linear velocities of the links
% - dynamics, for now, no parameters are stored
		a.ac = zeros(3,a.N,a.N);
		a.bc = zeros(3,a.N,a.N);
		a.ac0 = zeros(3,a.N);
		a.bc0 = zeros(3,a.N);
% - finally, create the k_flier object
		a = class(a, 'k_flier');
    end
else
	error('Constructor has to be called with one structure array or k_flier object argument!');
end

%-----------------------------------------------------------------------------------
function a = Decompose(links)
if( all( ismember( { 'M' 'angs' 'linkgeo' 'linkdyn' }, fieldnames(links) ) ) )
% -------- setting general parameters -----------
	a.N = length(links.linkgeo) + 5;	% No. of links (+5 because of base positioning)
	a.ch = size(links.M,1); % No. of chains
	a.M = [ repmat(1:5, a.ch, 1) links.M+(5 * (links.M ~= 0)) ];	% chain configuration
	a.ors = links.angs;		% orientation angles sequence
% -------- decomposing dynamic parameters -----------
	if((a.N-5) ~= max(size(links.linkdyn)))
		error('linkgeo and linkdyn has to be of same dimension');
	end
	a.J = zeros(3,3,a.N);
	a.ms = zeros(1,a.N);
	for i=6:a.N		% links until the 6-th are just dummy links, dynamics remain ZERO
		a.J(:,:,i) = links.linkdyn(i-5).J;
		a.ms(1,i) = links.linkdyn(i-5).m;
	end
% -------- decomposing geometry parameters -----------
	a.le = [ 1 0 0; 0 1 0; 0 0 1 ]';
	a.lne = repmat([ 0 1 0; 0 0 1 ]',[1 1 a.ch]);
	a.ldr = zeros(3,a.N);
	a.lpr = zeros(3,a.N,a.ch);
	a.su = zeros(a.ch,a.N);
	a.sl = zeros(1,a.N);
% le (local joint axes) vectros for 1 - 3 remain 0 since there is no rotation around them
	for i=4:6
        switch a.ors(i-3)
		case 'x'
			a.le(:,i) = [1 0 0]';
			a.lne(:,i-1,:) = repmat([1 0 0]',[1 1 a.ch]);
		case 'y'
			a.le(:,i) = [0 1 0]';
			a.lne(:,i-1,:) = repmat([0 1 0]',[1 1 a.ch]);
		case 'z'
			a.le(:,i) = [0 0 1]';
			a.lne(:,i-1,:) = repmat([0 0 1]',[1 1 a.ch]);
		end
	end
	ju = 1;
	jl = 1;
	for i=6:a.N		% links until the 6-th are just dummy links, vectors remain ZERO
		if (i>6)
			a.le(:,i) = links.linkgeo(i-5).le;
		end
		a.lne(:,i,:) = links.linkgeo(i-5).lne;
		a.ldr(:,i) = links.linkgeo(i-5).ldr;
		a.lpr(:,i,:) = links.linkgeo(i-5).lpr;
		% check if the i-th joint has a specificity on its LOWER side
		if(vmod(links.linkgeo(i-5).aldr)>0)
			a.sl(i) = jl;
			a.aldr(:,jl) = links.linkgeo(i-5).aldr;
			jl = jl + 1;
		end
		% check if the i-th joint has a specificity on its UPPER side
		ac = vmod(links.linkgeo(i-5).alpr);
		while(~isempty(find(ac>0,1)))
			a.su(find(ac>0,1),i) = ju;
			a.alpr(:,ju) = links.linkgeo(i-5).alpr(:,find(ac>0,1));
			ac(find(ac>0,1)) = 0;
			ju = ju + 1;
		end
	end
else
	error('Invalid ''links'' structure fields');
end

%-----------------------------------------------------------------------------------
function a = AnalyseAssemble(a)
% -------- find chain forks (fork segment index) & chain lengths --------
a.bs = [1];		% chain fork indices
a.cl = [];		% chain lengths
for i=1:a.ch
	if( i>1 )
		tf = find( ~( a.M(i-1,:) == a.M(i,:) ) );
		a.bs = [ a.bs tf(1)-1 ];
	end
	tf = find( ~(a.M(i,:) ) );
	if( length(tf) )
		a.cl = [ a.cl tf(1)-1 ];
	else
		a.cl = [ a.cl size(a.M,2) ];
	end
end
% ---------      Assemble joints      -----------
% find transfrom matrices in their zero positions 
% c_chain - current chain that is being analysed
% c_si - current segment index that is being analysed
a.As0(:,:,1) = eye(3);	% sequential transf. matrix to the first segment always unity
a.Ad0(:,:,1) = eye(3);	% direct transf. matrix to the first segment is unity
for c_chain = 1:a.ch	% joints are in different chains
	chain = a.M(c_chain,:);
	for c_si = (a.bs(c_chain)+1):a.cl(c_chain)		% cycle though the current chain
		i = chain(c_si);
		ip = chain(c_si-1);
		if( i>6 )
			if(a.su(c_chain,ip))
				tap = cross3( -a.lne(:,ip,c_chain) , cross3( a.alpr(:,a.su(c_chain,ip)) , a.lne(:,ip,c_chain) ) );
			else
				tap = cross3( -a.lne(:,ip,c_chain) , cross3( a.lpr(:,ip,c_chain) , a.lne(:,ip,c_chain) ) );
			end
			if(a.sl(i))
				tan = cross3( a.le(:, i) , cross3( a.aldr(:,a.sl(i)), a.le(:, i) ) );
			else
				tan = cross3( a.le(:, i) , cross3( a.ldr(:,i), a.le(:, i) ) );
			end
			tap = tap / vmod(tap);
			tan = tan / vmod(tan);
			tbp = cross3(a.lne(:,ip,c_chain), tap);
			tbn = cross3(a.le(:,i), tan);
			a.As0(:,:,i) = [a.lne(:,ip,c_chain) tap tbp]*[a.le(:,i) tan tbn]';
		else
			a.As0(:,:,i) = eye(3);
		end
		a.Ad0(:,:,i) = a.Ad0(:,:,ip) * a.As0(:,:,i);
	end
end

function a = AnalyseAssembleNK(a)
% -------- find chain forks (fork segment index) & chain lengths --------
a.bs = [1];		% chain fork indices
a.cl = [];		% chain lengths
for i=1:a.ch
	if( i>1 )
		tf = find( ~( a.M(i-1,:) == a.M(i,:) ) );
		a.bs = [ a.bs tf(1)-1 ];
	end
	tf = find( ~(a.M(i,:) ) );
	if( length(tf) )
		a.cl = [ a.cl tf(1)-1 ];
	else
		a.cl = [ a.cl size(a.M,2) ];
	end
end
% ---------      Assemble joints      -----------
% find transfrom matrices in their zero positions 
% c_chain - current chain that is being analysed
% c_si - current segment index that is being analysed
a.As0(:,:,1) = eye(3);	% sequential transf. matrix to the first segment always unity
a.Ad0(:,:,1) = eye(3);	% direct transf. matrix to the first segment is unity
for c_chain = 1:a.ch	% joints are in different chains
	chain = a.M(c_chain,:);
	for c_si = (a.bs(c_chain)+1):a.cl(c_chain)		% cycle though the current chain
		i = chain(c_si);
		ip = chain(c_si-1);
		if( i>6 )
			if(a.su(c_chain,ip))
				tap = cross3( -a.lne(:,ip,c_chain) , cross3( a.alpr(:,a.su(c_chain,ip)) , a.lne(:,ip,c_chain) ) );
			else
				tap = cross3( -a.lne(:,ip,c_chain) , cross3( a.lpr(:,ip,c_chain) , a.lne(:,ip,c_chain) ) );
			end
			if(a.sl(i))
				tan = cross3( a.le(:, i) , cross3( a.aldr(:,a.sl(i)), a.le(:, i) ) );
			else
				tan = cross3( a.le(:, i) , cross3( a.ldr(:,i), a.le(:, i) ) );
			end
			tap = tap / vmod(tap);
			tan = tan / vmod(tan);
			tbp = cross3(a.lne(:,ip,c_chain), tap);
			tbn = cross3(a.le(:,i), tan);
			a.As0(:,:,i) = [a.lne(:,ip,c_chain) tap tbp]*[a.le(:,i) tan tbn]';
		else
			a.As0(:,:,i) = eye(3);
        end
        a.As0(:,:,i) = eye(3);
		a.Ad0(:,:,i) = a.Ad0(:,:,ip) * a.As0(:,:,i);
	end
end
