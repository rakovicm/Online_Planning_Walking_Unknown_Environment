function [J,JA] = kk_jakP(a, cn)
% calculates the geometry and configuration of the mechanism
% Usage :  k_flier_object = k_geo( k_flier_object, coordinate vector q )
%

	% internal angles (1:3 base translation, 4:6 - base orient., 7.. internal ang.)
lnk = a.Con.lnr(cn);
svec = a.Con.cvec(:,cn);
gsv = a.A(:,:,lnk)*svec;
if(lnk < 4)
	error('''lnk'' is not allowed to be 1,2 or 3.');
end;


[J JA]=kjakp(a.N,a.ch,a.M,a.bs,a.cl,a.alf,a.bet,a.gam,a.del,a.ge,a.omg,a.dq,lnk,gsv);
    
