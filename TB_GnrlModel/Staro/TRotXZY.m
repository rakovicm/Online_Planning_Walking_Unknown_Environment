function RS = TRotXZY(T);
% RS = TRotYZX(T) returns a series of rotation around the X, Z and Y axes respectively
% that are equivalent to the given homogenous transformation

aX = atan(T(3,2)/T(2,2));
aZ = atan(-T(1,2)/(T(2,2)*cos(aX) + T(3,2)*sin(aX)));
aY = atan(T(1,3)/T(1,1));
RS = [aX aZ aY];
