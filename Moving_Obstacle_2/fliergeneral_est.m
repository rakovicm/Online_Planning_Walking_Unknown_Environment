function xp=fliergeneral_est(ts,x)

    global flier PRO SKR seg tausim itr actseg Nseg t d qnew l1 l2 l3 v

    N=flier.N;
    qmereno=qnew(:,itr);
    qest=x(1:N);%ucitavam vrednosti vektora stanja u pomocne promenljive
    qvest=x(N+1:2*N);
    dest=x(2*N+1:3*N);
    xp=x;

%     km_geo(flier,qest);%racuna geometriju
%     km_kine(flier,qvest);%racuna kinematiku
%     km_dyn(flier);%racuna dinamiku
%     [H h0]=km_inemat(flier);%racuna matrice H i h0

    if(itr>2)
        Upr = interp1([t(itr-2:itr+1)],[v(:,itr-2:itr+1)]',ts,'linear')';%interpolacija upravljanja jer RKode racuna u 4 tacke izmedju dve periode odabiranja
%         if(t(i-1)>ts)
%             disp('err1');
%             t(i-1)-ts
%             pause(0.5);
%         end
%         if(ts>t(i))
%             disp('err2');
%             ts-t(i)
%             pause(0.5);                        
%         end
%     elseif(i>1)
%         tauupr = interp1([t(i-1:i)],[tausim(:,i-1:i)]',ts,'linear')';
    else
        Upr = v(:,itr);
    end

    
    
    
%     n=7;
    
    xp(1:N,1)=qvest+l1*(qmereno-qest);%izvod estimirane pozicije
    xp(N+1:2*N,1)=dest+Upr+l2*(qmereno-qest);%izvod estimirane brzine
    xp(2*N+1:3*N,1)=l3*(qmereno-qest);%izvod estimiranog poremecaja
%     for n=7:29% anti windup
%         if (dest(n,1)>=200) && ((qmereno(n)-qest(n))>0)
%             xp(2*N+n,1)=0;
%         elseif (dest(n,1)<=-200) && ((qmereno(n)-qest(n))<0)
%             xp(2*N+n,1)=0;
%         else
%             xp(2*N+n,1)=l3*(qmereno(n,1)-qest(n,1));%izvod estimiranog poremecaja
%         end
%     end
    
end