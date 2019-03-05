function xp=fliergeneral12(ts,x)
global Up flier t Fuks itr
global PRO SKR
global Lr Bc Jr Cm Ce Rr 
global mi Kn Cn KK CC
global CT CTPROM
global XC YC XCL YCL XCR YCR XCPL YCPL XCPR YCPR

N=flier.N;

q=x(1:N);
qv=x(N+1:2*N);
Ir=x(2*N+1:3*N-6);
del=zeros(3,12);
del(:)=x(3*N-5:3*N+30);

delP=zeros(size(del));
F=delP;
lam=zeros(12,1);

flier=kk_geo(flier,q);
flier=kk_kine(flier,qv);
flier= kk_dyn(flier);
[H h0]=kk_inemat(flier);


Ct=[CC 0; 0 CC];
Kt=[KK 0; 0 KK];
Q=[1 0 0; 0 1 0; 0 0 -mi^2];
Cmat=diag([CC,CC,Cn]);
Kmat=diag([KK,KK,Kn]);
taudod=zeros(N,1);
Fuks=zeros(6,1);

XC=[];
YC=[];
XCL=[];
YCL=[];
XCR=[];
YCR=[];
XCPL=[];
YCPL=[];
XCPR=[];
YCPR=[];

for ii=1:12
    Rk=km_linkX(flier,ii,'xyz');
    [Jk Ak]=km_jakP(flier,ii);     %#ok<NASGU>
    Vk=Jk*qv;
    
%     if CT(ii)
        %%%%%%%%%%%%%%%% ako je kontakt aktivan onda je delta+Rk=0
        delP(3,ii)=-Vk(3);
        F(3,ii)=Kn*del(3,ii)+Cn*delP(3,ii);
        
        %%%%%% stick    
        F(1:2,ii)=Kt*del(1:2,ii)-Ct*Vk(1:2);  
        delP(1:2,ii)=lam(ii)*F(1:2,ii)-Vk(1:2);            
        if F(3,ii)<=0
            CTPROM(ii)=1;
            F(:,ii)=0;            
        elseif  (F(:,ii)'*Q*F(:,ii)>0) %%%%%%%%%%%%%% slip            
            M=inv(Ct)*Kt*del(1:2,ii)-Vk(1:2);
            lam(ii)=1/CC-sqrt(M'*M/((mi*F(3,ii))^2));
            F(1:2,ii)=inv(inv(Ct)-lam(ii)*eye(2))*M;
            delP(1:2,ii)=lam(ii)*F(1:2,ii)-Vk(1:2);            
        end;
        
        
        if(F(3,ii)~=0)
            XC=[XC Rk(1)];
            YC=[YC Rk(2)];
            if ii>0 && ii<7 
                XCL=[XCL Rk(1)];
                YCL=[YCL Rk(2)];
                if ii>2
                    XCPL=[XCPL Rk(1)];
                    YCPL=[YCPL Rk(2)];
                end            
            end 
            if ii>6 && ii<13
                XCR=[XCR Rk(1)];
                YCR=[YCR Rk(2)];
                if ii>8
                    XCPR=[XCPR Rk(1)];
                    YCPR=[YCPR Rk(2)];
                end
            end  
        end
%         
%     else 
%         delP(1:3,ii)=-Cmat\Kmat*del(1:3,ii);
%         if Rk(3)+del(3,ii)<0
%             CTPROM(ii)=1;
%         end;
%     end

    taudod=taudod+Jk(1:3,:)'*F(:,ii);
    Fuks(1:3)=Fuks(1:3)+F(:,ii);
    Fuks(4:6)=Fuks(4:6)+cross2(Rk(1:3),F(:,ii));
end;



if(itr<3)
    Upnew=interp1([t 100],[Up';zeros(1,45)],ts,'linear')';
else
    Upnew=interp1([t(itr-2:itr+2)],[Up(:,itr-2:itr+2)'],ts,'linear')';
% Upnew2-Upnew
end

if(any(isnan(Upnew)))
    Upnew=interp1([t 100],[Up';zeros(1,45)],ts,'linear')';
end

M1=[H                                  zeros(N,N-6)        -PRO       ];
M4=[Jr*SKR                             zeros(N-6)          eye(N-6)   ];
M5=[zeros(N-6,N)                       Lr*eye(N-6)         zeros(N-6) ];



   
D1=-h0+taudod;
D4=Cm*Ir-Bc*SKR*qv;
D5=Upnew-Rr*Ir-Ce*SKR*qv;
% % % if EP
% % %     D4((19:21)-6)=0;
% % %     M4((19:21)-6,:)=0;
% % %     M4((19:21)-6,N+N-6+(19:21)-6)=eye(3);
% % % end;

M=[M1; M4; M5];
D=[D1; D4; D5];    

res=M\D;
qa=res(1:N);
Irp=res(N+1:2*N-6);
tau=res(2*N-5:3*N-12);
xp=[qv; qa; Irp; delP(:); tau; F(:); lam(:)];

