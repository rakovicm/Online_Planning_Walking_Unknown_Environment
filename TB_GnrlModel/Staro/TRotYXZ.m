function RS = TRotYXZ(T);
% RS = TRotYXZ(T) returns a series of rotation around the Y, X and Z axes respectively
% that are equivalent to the given homogenous transformation

aY = atan(T(1,3)/T(3,3));
aX = atan(-T(2,3)/(T(3,3)*cos(aY) +  T(1,3)*sin(aY)));
aZ = atan(T(2,1)/T(2,2));
RS = [aY aX aZ];
