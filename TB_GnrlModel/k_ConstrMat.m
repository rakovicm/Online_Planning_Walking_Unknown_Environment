function [Jc, Jf, Ac, Af K m] = k_ConstrMat(flier, cn, object, ovel, oacc, ccs)
% forms the constrained motion (also impact analysis) matrices
% Usage : [Jc, Jf, Ac, Af] = k_ConstrMat(flier, cn, object, ovel, oacc, cs, fs)
% flier - the flier object
% cn - flier object contact index
% object - the contacted object
% ovel - 6 element vector - object coordinates first derivative
% oacc - 6 element vector - object coordinates second derivative
% cs - the constrained components of the contact
% fs - complement of the previous, the not constrained, free components
%
m=sum(ccs);
cs=[];
fs=[];
for i=1:6
    if ccs(i)
      cs=[cs i];
    else
        fs=[fs i];
    end;
end;
        

[Jl, Al] = k_jakP(flier, cn);
[Jsl, Asl] = k_Jsl(object);
[object, Jsb, Asb] = k_Jsb(object, ovel);		% + first coordinate derivatives
Jslc = Jsl(cs,:);
Jslf = Jsl(fs,:);
Jsbc = Jsb(cs,:);
Jsbf = Jsb(fs,:);
Asc = Asl(cs,:) + Asb(cs,:);
Asf = Asl(fs,:) + Asb(fs,:);
Jc = Jslc * Jl;
Jf = Jslf * Jl;
Ac = Jslc * Al + Jsbc * oacc + Asc;
Af = Jslf * Al + Jsbf * oacc + Asf;
K=zeros(6);

for i=1:m
    K(i,cs(i))=1;
end;
for i=m+1:6
    K(i,fs(i-m))=1;
end;