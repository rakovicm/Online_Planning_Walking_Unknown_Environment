function [XCp YCp]=polyoffset(XC,YC,dist)
    if length(XC)>3  

        T=[XC(1); YC(1)];
        P=[XC(end-1); YC(end-1)];
        S=[XC(2); YC(2)];    
        p=length(XC);
        vp=(T-P);
        vp=vp/sqrt(vp'*vp);
        np=[-vp(2); vp(1)];
        vs=(S-T);
        vs=vs/sqrt(vs'*vs);
        ns=[-vs(2); vs(1)];
        np=np*sign(np'*vs);
        ns=ns*sign(-ns'*vp);
        MAT(1,:)=np';
        MAT(2,:)=ns';
        VEC=dist+[P'*np;T'*ns];
        res=inv(MAT)*VEC;
        XCp=zeros(p,1);
        YCp=XCp;
        XCp(1)=res(1);
        YCp(1)=res(2);
        for i=2:p-1
            vp=vs;        
            %%%%%%%% pomeramo tacku i
           T=[XC(i); YC(i)];       
           S=[XC(i+1); YC(i+1)];       
           vs=(S-T);
           vs=vs/sqrt(vs'*vs);
           ns=[-vs(2); vs(1)];       
           ns=ns*sign(-ns'*vp);
           koeff=dist+T'*ns-[XCp(i-1) YCp(i-1)]*ns;
           if ((vp'*ns)==0)           
               disp('Greska, kolinearne tacke')          
           end
           koeff=koeff/(vp'*ns);     
           if (koeff<0)
               disp('Linija Degradira');
           end       
           XCp(i)=XCp(i-1)+koeff*vp(1);
           YCp(i)=YCp(i-1)+koeff*vp(2);
        end;
        %%%ind=convhull(XCp,YCp);
        %%%XCp=XCp(ind);
        %%%YCp=YCp(ind);    
        XCp(end)=XCp(1);
        YCp(end)=YCp(1);
        if ([XC(end)-XC(end-1);YC(end)-YC(end-1)]'*[XCp(end)-XCp(end-1);YCp(end)-YCp(end-1)]<0)
              disp('Linija Degradira');
        end;
    else
        disp('Gari, daj vise tacaka!');
    end
end