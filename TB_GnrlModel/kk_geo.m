function a = kk_geo(a, q)
% calculates the geometry and configuration of the mechanism
% Usage :  k_flier_object = k_geo( k_flier_object, coordinate vector q )
%
if( length(q) ~= (a.N) )
	error('Coordinate vector length illegal error!');
end

a.q = q;	% internal angles (1:3 base translation, 4:6 - base orient., 7.. internal ang.)
[a.A a.rc a.ge a.gdr a.gpr]=kgeo(a.N,a.ch,a.M,a.bs,a.cl,a.ors,a.q,a.lne,...
    a.As0,a.le,a.ldr, a.lpr);
