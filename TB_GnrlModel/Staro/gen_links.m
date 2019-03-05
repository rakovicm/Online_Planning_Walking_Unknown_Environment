function lk = gen_links(fn)
% Opens a text file, reads the information in it and generates
% the 'links' structure of the flier
% Usage : 	links = gen_links(text_file_name)
%

% ----- Open the input file ----
fi = fopen(fn, 'r');
if( fi == -1 )
	error('Couldn''t open input file!');
end

% ----- Check the file signature : 'K_FLIER' at the top ----
inl = fgetl(fi);
vi = sscanf(inl, '%s', Inf);
if( ~strcmp( vi,'K_FLIER' ) )
	fclose(fi);
	error('K_FLIER signature must be present at the top of the file!');
end

% ----- Input parameters for each link one by one -----
i = 0;		% current link index
stp = 0;	% regular stop flag
alf = 1;	% all fields of the link completed
chs = 0;	% chain structure description in progress
ors = 0;	% base orientation angles sequence
orss = 0;	% base orientation angels sequence defined
lk.M = [];
while( (~feof(fi)) & (~stp) )
	vi = fgetl(fi);
	[pi,n] = sscanf(vi,'%s',1);
	switch pi
	case '$'
		ors = 0;
		if( i == 0 )
			chs = 1;
		else
			fclose(fi);
			error( 'Chain description has to be before description of the links!' );
		end
	case '&'
		ors = 1;
	case '#'
		ors = 0;
		if(alf)
			% fill in the optional field - joint specificity
			if(exist('flds'))
				if(~flds(7))
					lk.linkgeo(i).aldr = zeros(3,1);
				end
				if(~flds(8))
					lk.linkgeo(i).alpr = zeros(3,1);
				end
			end
			i = i + 1;
			alf = 0;
			flds = zeros(1,8);
			chs = 0;
		else
			fclose(fi);
			error( strcat('Not all fields filled for link # ',num2str(i)) );
		end
		if( max(size(lk.M)) == 0 )
			fclose(fi);
			error( 'Link definition can''t commence before chain definition!' );
		end
	case '@'
		stp = 1;
		if( ~alf)
			fclose(fi);
			error( strcat('Data end reached, but not all fields filled for link # ',num2str(i)) );
		elseif( ~orss )
			fclose(fi);
			error( 'Base orientation angels sequence not set!' );
		end
	case ''
		% do nothing, it's not an error
	otherwise
		if( ~chs )
			switch pi
			case 'le'
				try
					flds(1) = 1;
					[lk.linkgeo(i).le(1,1), lk.linkgeo(i).le(2,1), lk.linkgeo(i).le(3,1)] = ...
						strread(vi,'le %f %f %f','commentstyle','c++');
				catch
					fclose(fi);
					error( strcat('Error processing ''le'' field of link # ',num2str(i)) );
				end
			case 'lne'
				try
					flds(2) = 1;
					vi = strread(vi,'lne %s','delimiter','|');
					[lne1, lne2, lne3] = strread(vi{1}, '%f %f %f;', 'commentstyle', 'c++');
					lk.linkgeo(i).lne = [lne1 lne2 lne3]';
				catch
					fclose(fi);
					error( strcat('Error processing ''lne'' field of link # ',num2str(i)) );
				end
			case 'ldr'
				try
					flds(3) = 1;
					[lk.linkgeo(i).ldr(1,1), lk.linkgeo(i).ldr(2,1), lk.linkgeo(i).ldr(3,1)] = ...
						strread(vi,'ldr %f %f %f','commentstyle','c++');
				catch
					fclose(fi);
					error( strcat('Error processing ''ldr'' field of link # ',num2str(i)) );
				end
			case 'lpr'
				try
					flds(4) = 1;
					vi = strread(vi,'lpr %s','delimiter','|');
					[lpr1, lpr2, lpr3] = strread(vi{1}, '%f %f %f;', 'commentstyle', 'c++');
					lk.linkgeo(i).lpr = [lpr1 lpr2 lpr3]';
				catch
					fclose(fi);
					error( strcat('Error processing ''lpr'' field of link # ',num2str(i)) );
				end
			case 'm'
				try
					flds(5) = 1;
					lk.linkdyn(i).m = ...
						strread(vi,'m %f','commentstyle','c++');
				catch
					fclose(fi);
					error( strcat('Error processing ''m'' field of link # ',num2str(i)) );
				end
			case 'J'
				try
					flds(6) = 1;
					[lk.linkdyn(i).J(1,1), lk.linkdyn(i).J(1,2), lk.linkdyn(i).J(1,3),...
						lk.linkdyn(i).J(2,1), lk.linkdyn(i).J(2,2), lk.linkdyn(i).J(2,3),...
						lk.linkdyn(i).J(3,1), lk.linkdyn(i).J(3,2), lk.linkdyn(i).J(3,3) ] = ...
							strread(vi,'J %f %f %f; %f %f %f; %f %f %f;','commentstyle','c++');
				catch
					fclose(fi);
					error( strcat('Error processing ''J'' field of link # ',num2str(i)) );
				end
			case 'aldr'
				try
					flds(7) = 1;
					[lk.linkgeo(i).aldr(1,1), lk.linkgeo(i).aldr(2,1), lk.linkgeo(i).aldr(3,1)] = ...
						strread(vi,'aldr %f %f %f','commentstyle','c++');
				catch
					fclose(fi);
					error( strcat('Error processing ''aldr'' field of link # ',num2str(i)) );
				end
			case 'alpr'
				try
					flds(8) = 1;
					vi = strread(vi,'alpr %s','delimiter','|');
					[alpr1, alpr2, alpr3] = strread(vi{1}, '%f %f %f;', 'commentstyle', 'c++');
					lk.linkgeo(i).alpr = [alpr1 alpr2 alpr3]';
				catch
					fclose(fi);
					error( strcat('Error processing ''alpr'' field of link # ',num2str(i)) );
				end
			end
			if( all(flds(1:6)) )
				if(~flds(7))	% if the alternative vector is not set, just set it to zero
					lk.linkgeo(i).aldr = zeros(3,1);
				end
				if(~flds(8))
					lk.linkgeo(i).alpr = zeros(3,1);
				end
				alf =1;
			end
		elseif( ors )
			ors = 0;
			try
				[lk.angs(1), lk.angs(2), lk.angs(3)] = ...
					strread(vi,'%c %c %c','commentstyle','c++');
			catch
				fclose(fi);
				error( 'Error processing base orientation angle sequence!' );
			end
			orss = 1;
		else
			imed = strread(vi,'%f','commentstyle','c++')';
			dif = size(imed,2) - size(lk.M,2);
			if( dif >= 0 )
				lk.M = [lk.M zeros(size(lk.M,1),dif)];
			else
				imed = [imed zeros(1,-dif)];
			end
			lk.M = [lk.M; imed];
		end
	end
end
if( ~stp )
	fclose(fi);
	error('Unexpected end of file! End data with ''@''!');
end

fclose(fi);
