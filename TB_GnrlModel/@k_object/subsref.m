function rv = subsref(a, S)
switch S(1).type
case '.'
	switch S(1).subs
	case 'e'
		rv = a.e;
	case 'q'
		rv = a.q;
	case 'dq'
		rv = a.dq;
	case 'rots'
		rv = a.rots;
	case 'contact'
		rv = a.contact;
	case 'cQ'
		rv = a.cQ;
	case 'Q'
		rv = a.Q;
	case 'omg'
		rv = a.omg;
	case 'v'
		rv = a.v;
	case 'alf'
		rv = a.alf;
	case 'bet'
		rv = a.bet;
	case 'gam'
		rv = a.gam;
	case 'del'
		rv = a.del;
	case 'm'
		rv = a.m;
	case 'J'
		rv = a.J;
	case 'rc'
		rv = a.rc;
	case 'ac'
		rv = a.ac;
	case 'ac0'
		rv = a.ac0;
	case 'bc'
		rv = a.bc;
	case 'bc0'
		rv = a.bc0;
	otherwise
		error('Getting that data is not inmplemented YET!');
	end
otherwise
	error('That can''t be displayed!');
end
