function RS = TRotZYX(T);
% RS = TRotZYX(T) returns a series of rotation around the Z, Y and X axes respectively
% that are equivalent to the given homogenous transformation

aZ = atan(T(2,1)/T(1,1));
aY = atan(-T(3,1)/(T(1,1)*cos(aZ) + T(2,1)*sin(aZ)));
aX = atan(T(3,2)/T(3,3));
RS = [aZ aY aX];
