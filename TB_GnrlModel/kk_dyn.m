function a = kk_dyn(a)
% calculates the dynamics coeffcients of the system
%  Usage : k_flier_object = k_dyn(k_flier_object)
% thos function calculates a bunch of coefficients needed for later calculation
% of different dynamics-related quantities:
%	* inertia matrices H and h0 <-> k_inemat
%	* inertial forces and torques <-> k_infor
% according to this, MUST be called before these functions!
%

%Prilagodi dimenzije ulaznim podacima ko zna kako on matricu koja ima tri
%dimenzije spakuje u dvodimenzionalnu matricu

[a.ac, a.ac0, a.bc, a.bc0] = kdyn(a.N,a.ch,a.M,a.bs,   a.cl,a.ms,a.J, a.A,...
                           a.alf,  a.bet,a.gam,a.del ,a.omg,a.v);
  
    
  
    
  

   
	
  
   
  
	
 