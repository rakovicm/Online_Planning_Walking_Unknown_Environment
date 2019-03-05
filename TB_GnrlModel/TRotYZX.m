function RS = TRotYZX(T);
% RS = TRotYZX(T) returns a series of rotation around the Y, Z and X axes respectively
% that are equivalent to the given homogenous transformation

aY = atan(-T(3,1)/T(1,1));
aZ = atan(T(2,1)/(T(1,1)*cos(aY) - T(3,1)*sin(aY)));
aX = atan(-T(2,3)/T(2,2));
RS = [aY aZ aX];
