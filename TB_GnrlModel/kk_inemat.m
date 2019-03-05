function [H, h0] = kk_inemat(a)
% Calculates the inertia matrix and adjoint vector of the system
%  Usage :  [H, h0] = k_inemat(k_flier_object)
		% gravity, fixed (for now)


[H,h0]=kinemat(a.N,a.ch,a.M,a.bs,a.cl,a.ms,a.ge,a.rc,a.ac,a.ac0,a.bc,a.bc0);
