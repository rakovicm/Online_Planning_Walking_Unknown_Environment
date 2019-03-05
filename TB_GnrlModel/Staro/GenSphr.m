function [x, y, z] = GenSphr(radius,poly)
%    new_sphere=GenSphr(radius,poly)
%    This function generates a new sphere 
%    Input variables:
%    		radius: sphere radius in (m)
%			poly:	no. of segments forming the sphere
%    Output variables:    
%           new_sphere: Structure od sphere x,y,and t coordinates

[x,y,z]=sphere(poly);
for i=1:poly+1 
   for j=1:poly+1
       p_B(:,1)=[radius*x(i,j) radius*y(i,j) radius*z(i,j) 1]';
       x(i,j)=p_B(1);
       y(i,j)=p_B(2);
       z(i,j)=p_B(3);
   end
end
