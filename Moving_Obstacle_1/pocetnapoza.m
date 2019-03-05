%%%%%%%%%%load poza2

global Qtil S
letac_motor_kontakti;
links=gen_links('anthr46ss.m');
flier=k_flier(links);
N=flier.N;
qsave=zeros(N,1);
q=qsave;
q(1:6)=0;
seg=[ 6 1; 9 1; 10 1; 13 1; 14 1; 6 2; 17 2; 18 2; 21 2; 22 2];
seg=[seg; 6 3; 23 3; 25 3; 27 3; 29 3; 31 3; 33 3; 35 3; 37 3];
seg=[seg; 39 3; 41 3; 44 4; 44 3; 45 3; 46 3; 47 3;  48 4];
seg=[seg; 49 4; 50 4; 51 4];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%    KONTAKTNE TACKE    %%%%%%%%%%%%%%%%%%%%
%%%                                                         %%%        
%%%                                           ^ X           %%%
%%%                                           |             %%%
%%%                     5 _______ 6         11|_______ 12   %%%
%%%                      |       |            |       |     %%%
%%% Y <________________ 3|_______|4__________9|_______|10   %%%
%%%                      |       |            |       |     %%%
%%%                      |levo   |            |desno  |     %%%    
%%%                      |stopalo|            |stopalo|     %%%
%%%                      |       |            |       |     %%%
%%%                      |  x 14 |            |  x 15 |     %%%
%%%                     1|_______|2          7|_______|8    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%
flier=k_addcontact(flier,  1, 21, [-0.112885  0.05  -0.041441],eye(3)); 
flier=k_addcontact(flier,  2, 21, [-0.112885 -0.05  -0.041441],eye(3));
flier=k_addcontact(flier,  3, 21, [ 0.082885  0.05  -0.041441],eye(3));
flier=k_addcontact(flier,  4, 21, [ 0.082885 -0.05  -0.041441],eye(3));
flier=k_addcontact(flier,  5, 22, [ 0.06      0.05   0       ],eye(3));
flier=k_addcontact(flier,  6, 22, [ 0.06     -0.05   0       ],eye(3));
%%%%%%%%%%%%%%%%%%%%%%
flier=k_addcontact(flier,  7, 13, [-0.112885  0.05  -0.041441],eye(3)); 
flier=k_addcontact(flier,  8, 13, [-0.112885 -0.05  -0.041441],eye(3));
flier=k_addcontact(flier,  9, 13, [ 0.082885  0.05  -0.041441],eye(3));
flier=k_addcontact(flier, 10, 13, [ 0.082885 -0.05  -0.041441],eye(3));
flier=k_addcontact(flier, 11, 14, [ 0.06      0.05   0       ],eye(3));
flier=k_addcontact(flier, 12, 14, [ 0.06     -0.05   0       ],eye(3));

%bazni segment
flier=k_addcontact(flier, 13, 6, [0 0 0],eye(3));

%skocni zglob leve i desne noge respektivno
flier=k_addcontact(flier, 14, 18, [0 0 -0.20955],eye(3));
flier=k_addcontact(flier, 15, 10, [0 0 -0.20955],eye(3));

%kuk leve i desne noge respektivno
flier=k_addcontact(flier, 16, 15, [0  0  0],eye(3));
flier=k_addcontact(flier, 17,  7, [0  0  0],eye(3));


%peta
flier=k_addcontact(flier,  18, 21, [-0.112885  0  -0.041441],eye(3)); 
flier=k_addcontact(flier,  19, 13, [-0.112885  0  -0.041441],eye(3)); 

%rame
flier=k_addcontact(flier,  20, 43, [0 0.1786 0.0313],eye(3)); 
flier=k_addcontact(flier,  21, 43, [0 -0.1786 0.0313],eye(3)); 

%ruke
flier=k_addcontact(flier,  22, 47, [0 0 -0.132],eye(3)); 
flier=k_addcontact(flier,  23, 51, [0 0 -0.132],eye(3)); 

%trup
flier=k_addcontact(flier,  24, 43, [0 0 0],eye(3)); 

%prsti
flier=k_addcontact(flier,  25, 14, [0 0 0],eye(3)); 
flier=k_addcontact(flier,  26, 22, [0 0 0],eye(3)); 


brojkontakata=12;

mass=sum(flier.ms);
q=zeros(size(q));
% q(3)=0.9;
q(13)=-0.3;
q(21)=-0.3;

q(10)=0.4;
q(18)=0.4;

q( 9)=-0.1;
q(17)=-0.1;

q(46)=deg2rad(10);
q(50)=deg2rad(10);

q(47)=-deg2rad(20);
q(51)=-deg2rad(20);

% q( 8)=-0.01;
% q(16)=-0.01;
% 
% q(11)=-0.01;
% q(19)=-0.01;

% q(15)=0.1;
% q(23)=0.1;


qsave=q;
km_geo(flier,q);
km_kine(flier,zeros(N,1));
km_dyn(flier);
rcm=k_CentarMase(flier);
Q=zeros(3,3);
RR=zeros(3,1);
for i=1:12
    R(:,i)=km_linkX(flier,i,'xyz');    
end;

M=zeros(15,16);
sel=ones(12,1);
ss=12;
del=-1*ones(12,1);
K=Kn;
while any(del<0)
    RR=zeros(3,1);
    Q=zeros(3,3);
    for i=1:12
        if sel(i)
            RR=RR+R(1:3,i);
            Q =Q+R(1:3,i)*R(1:3,i)';
        end;
    end;
    S=(RR/ss-rcm)*mass*9.81/K;
    Qtil=Q-1/ss*RR*RR';
    RES=fsolve(@func,[0;0;1;0]);
    A=RES(1:3);
    p=RES(4);
    Zb=-1/ss*(RR'*A+mass*9.81/K);    
    for i=1:12;
        if sel(i)
            del(i)=-Zb-R(1:3,i)'*A;
        else
            del(i)=0;
        end;
    end;
    [mm pos]=min(del);
    sel(pos)=0;
    ss=ss-1;
    
    
end;
A1New=A(1)*cos(qsave(6))-A(2)*sin(qsave(6));
A2New=A(1)*sin(qsave(6))+A(2)*cos(qsave(6));
A3New=A(3);
qsave(4)=asin(A2New);
qsave(5)=atan(-A1New/A3New);
qsave(1:2)=0;
qsave(3)=Zb;
qvsave=0*qsave;
km_geo(flier,qsave);
km_kine(flier,qvsave);
km_dyn(flier);
[H,h0]=km_inemat(flier);
qa=qvsave;
taudod=zeros(N,1);

for i=1:12    
    [Jk Ak]=k_jakP(flier,i);    
    F(:,i)=diag([Kn,KK,KK])*[0;0;del(i)];
    taudod=taudod+Jk(1:3,:)'*F(:,i);
end;
tau=H*qa+h0-taudod;
Ir=(Jr*qa+Bc*qvsave+tau)/Cm;
Irp=0*Ir;
Up=Lr*Irp+Rr*Ir+Ce*qvsave;
Ir=SKR*Ir;
Up=SKR*Up;
save poza3 qsave qvsave del Ir Up  -V6
