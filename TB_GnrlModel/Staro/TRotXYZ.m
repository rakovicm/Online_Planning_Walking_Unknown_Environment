function RS = TRotXYZ(T);
% RS = TRotXYZ(T) returns a series of rotation around the X, Y and Z axes respectively
% that are equivalent to the given homogenous transformation

aX = atan(-T(2,3)/T(3,3));
aY = atan(T(1,3)/(T(3,3)*cos(aX) - T(2,3)*sin(aX)));
aZ = atan(-T(1,2)/T(1,1));
RS = [aX aY aZ];
