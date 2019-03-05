function y=func(x)
global Qtil S
A(:,1)=x(1:3);
p=x(4);
y(1)=A'*A-1;
y(2:4)=S-(p*eye(3)+Qtil)*A;
