function rv = subsref(a, S)
switch S(1).type
case '.'
	switch S(1).subs
	case 'N'
		rv = a.N;
	case 'omg'
		rv = a.omg;
	case 'v'
		rv = a.v;
	case 'q'
		rv = a.q;
	case 'dq'
		rv = a.dq;
	case 'ddq'
		rv = a.ddq;
	case 'alf'
		rv = a.alf;
	case 'bet'
		rv = a.bet;
	case 'gam'
		rv = a.gam;
	case 'del'
		rv = a.del;
	case 'ac'
		rv = a.ac;
	case 'bc'
		rv = a.bc;
	case 'ac0'
		rv = a.ac0;
	case 'bc0'
		rv = a.bc0;
	case 'ms'
		rv = a.ms;
	case 'J'
		rv = a.J;
	case 'ge'
		rv = a.ge;
	case 'bs'
		rv = a.bs;
	case 'cl'
		rv = a.cl;
	case 'ch'
		rv = a.ch;
	case 'M'
		rv = a.M;
	case 'rc'
		rv = a.rc;
	case 'A'
		rv = a.A;
    case 'gdr'
        rv=a.gdr;
    case 'gpr'
        rv=a.gpr;
    case 'le'
        rv=a.le;
otherwise
		error('Getting that data is not inmplemented YET!');
	end
otherwise
	error('That can''t be displayed!');
end
