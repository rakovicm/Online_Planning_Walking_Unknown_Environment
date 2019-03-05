function [dif] = k_dif(obs, cis, obx, cix)
% 
% Usage: dif = k_dif(obs, cis, obx, cix)
%
% returns the vector between the two contact points
% 
%
% obs - the object carrying the s coordinate frame
% obx - the other object
% cis, cix - the corresponding contact indices
%

Xs = [ obs.Q * obs.contact(:,cis) + obs.q(1:3,1) ];
Xx = [ obx.Q * obx.contact(:,cix) + obx.q(1:3,1) ];
dif = Xx - Xs;
dif = ( obs.Q * obs.cQ(:,:,cis) ) \ dif;
