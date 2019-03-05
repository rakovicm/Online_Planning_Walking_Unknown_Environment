function RS = TRotZXY(T);
% RS = TRotZXY(T) returns a series of rotation around the Z, X and Y axes respectively
% that are equivalent to the given homogenous transformation

aZ = atan(-T(1,2)/T(2,2));
aX = atan(T(3,2)/(T(2,2)*cos(aZ) - T(1,2)*sin(aZ)));
aY = atan(-T(3,1)/T(3,3));
RS = [aZ aX aY];
