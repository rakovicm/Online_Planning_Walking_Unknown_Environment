function simulacija
    clc    
    clear global;
    clear variables;

    global Up dt t Fuks ZMPtr CMtr CMzelj itr l1 l2 l3 v
    global XC YC XCL YCL XCR YCR XCPL YCPL XCPR YCPR TAC
    global segdkuk segdkol segdskz segdprst seglkuk seglkol seglskz seglprst 
    global segtrupx segtrupy segdrame segdlakat seglrame segllakat 
    global seglnoga segdnoga segtrup seglruka segdruka 

    global primosloninanoguakt brzosloninalnogul brzosloninadnogul brzosloninalnogur  brzspustilnogul brzspustilnogur brzosloninadnogur primgrcinoguakt brzgrcidnogul primopruzinoguakt brzopruzidnogu brzopruzidnogul brzopruzidnogur 
    global primlnl primlnd primlrl primlrd  primlt nlbsy ndbsy tbsy rbsy rdbsy
    global rbaz rskzl rskzd rkukl rkukd rpetal rpetad rramel rramed rrukal rrukad rtrup rprstl rprstd 
    global qnew qvnew qanew
    global statusnaginjanjenapred

    global primnaginjanjnapredakt brznaginjanjelnapred
    global brzgrcilnogul primspustinoguakt brzspustidnogul brzspustidnogur 
    global brzpocetnanagninapred itpocnagninapred
    global brzpocetnaopruzi itpocopruzi
    global brzpocetnagrci itpocgrci
    global brzpocetnapodigninogu brzpocetnaosloninalnogu itpocosloninalnogu brzpocetnaosloninadnogu itpocosloninadnogu
    global brzpocetnaspusti itpocspusti brzpocetnaosloninanogu brzspustidpetu brzpocetnaspustidpetu
    global primodgurninoguakt
    global brzodgurnidnogu brzodgurnilnogu brzpocetnaodgurni itpocodgurni brzspustilprst brzpocetnaspustilprst
    global brzspustilnogu brzspustidnogu brzgrcidnogur brzgrcilnogur 
    global brzspustilpetu brzpocetnaspustilpetu brzspustidprst brzpocetnaspustidprst
    global brzpocetnaosloninadnogu brzpocetnaspustidnogu brzpocetnaspustilprst brzpocetnaosloninalnogu brzpocetnaspustilnogu brzpocetnaspustidprst
    global zeljeni_pravac
    global zeljeni_pravac_L
    global zeljeni_pravac_D
    
    global flier flier_m ikraj
    global PRO SKR
    global Lr Bc Jr Cm Ce Rr 
    global mi Kn Cn KK CC
    global CT
    global CTPROM
    global r_saka_d_p r_saka_l_p
    
    load poza3
    load pocetniparametri
%     load podeseni_parametri9
%     
% qsave(1) = 5.9;
% qsave(2) = -2.1;
qsave(6) = 0;
    
    statusopruzi=0;
    statusgrci=0;    
    statusgrci1=0;
    statusnanog=0;
    statusispnap=0;
    
    %Inicijalizacija promenljivih
    letac_motor_kontakti;

    %Definisanje vremena simulacije
    i=1;%Trenutna iteracija
    itr=i;
    dt=0.8/1200;%period integracije
    tk=60;%ukupno vreme simulacije
    brit=round(tk/dt);%Ukupan broj iteracija
    t=(0:brit-1)*dt;%ukupno vreme simulacije

    %Pocetna vrednost vektora stanja
    x0=zeros(4*N-12+7*brojkontakata,1);
    q0=qsave;
    qv0=qvsave;
    Ir0=Ir;
    d0=[zeros(1,12);zeros(1,12);del'];
    x0(1:3*N+30,1)=[q0;qv0;Ir0;d0(:)];

    %Upravljacka promenljiva
    Up=Up*ones(1,brit);


    %Prvo diferenciranje 
    x0p=fliergeneral12(t(1),x0);
    X=zeros(length(x0),brit);
    Xp=X;
    X(:,i)=x0;
    Xp(:,i)=x0p;

    %Inicijalizacija rpomenljivih
    dd=zeros(3,brojkontakata);FF=dd; ddp=dd;
    qnew=zeros(N,brit); qvnew=qnew; qanew=qnew; qvv=qnew; Upr=qnew;
    Irnew=zeros(N-6,brit); Irpnew=Irnew; taunew=Irnew; 
    del=zeros(3,brojkontakata,brit); delp=del; Fnew=del; 
    Rk=zeros(6,brojkontakata,brit);Vk=Rk;
    Fuk=zeros(6,brit); rcm=zeros(3,brit); vcm=rcm; ZMPpos=rcm;
    ZMPtr=ZMPpos(:,i);
    CMtr=rcm(:,i);
    lam=zeros(brojkontakata,brit);
    CT(:)=1;
    CTPROM=0*CT;
    bk=4;


    %%Izvlacenje promenljivih iz pocetnog stanja
    qnew(:,i)  =X(1:N,i);
    qvnew(:,i) =X(N+1:2*N,i);
    Irnew(:,i) =X(2*N+1:3*N-6,i);
    dd(:)      =X(3*N-5:3*N+30,i);
    del(:,:,i) =dd;
    qanew(:,i) =Xp(N+1:2*N,i);
    Irpnew(:,i)=Xp(2*N+1:3*N-6,i);
    ddp(:)     =Xp(3*N-5:3*N+30,i); 
    delp(:,:,i)=ddp;
    taunew(:,i)=Xp(3*N+31:4*N+24,i);
    FF(:)      =Xp(4*N+25:4*N+60,i);
    Fnew(:,:,i)=FF;
    lam(:,i)   =Xp(4*N+61:4*N+72,i);               
    [rcm(:,i) vcm(:,i)]=k_CentarMase(flier);   
    nknew=sum(Fnew(3,:,i)>0);    
    nkold=nknew;
    primlnl=0;primlnd=0;primlrl=0;primlrd=0;primlt=0;nlbsy=0;ndbsy=0;tbsy=0;
    rbsy=0;rdbsy=0;nlbsy=0;ndbsy=0;tbsy=0;rbsy=0;rdbsy=0;

    faza=-1;
    sacuvaj=0;
%     %%%Na pocetku je za svaki slucaj ukljucen CM i ZMP regulator
    xa=5;
    xb=5;
   CMzeljt(:,i)=[rcm(1:2,i);0];

    brzinakoraka=1.6;
    visinakoraka=1;
    duzinakoraka=1;
    skretanje=0;
    zeljeni_pravac =  qsave(6);
    
  
    CMzelj(:,i)=[CMzeljt(:,i)];
    CMzeljs=CMzelj(:,i);
    promenacm=0;
    promenacms=0;
 
    dCMzelj=[0;0;0];
    icmzeljtr=i;
    ikraj = i;

    
    KpzmpcmynL = 0.4;%4.17;
    KizmpcmynL = 0.0000;
    KdzmpcmynL = 0.014;%20.02;    
    KpzmpcmxnL = 0.4;%4.17;
    KizmpcmxnL = 0.0000;
    KdzmpcmxnL = 0.014;%20.02;

    KpcmynL = 0.5;%4;
    KicmynL = 0.0000;
    KdcmynL = 0.014;%50.50;    
    KpcmxnL = 0.5;%4;
    KicmxnL = 0.0000;
    KdcmxnL = 0.014;%50.50;
    
    
    KpzmpcmynD = 0.4;%4.17;
    KizmpcmynD = 0.0000;
    KdzmpcmynD = 0.014;%20.02;    
    KpzmpcmxnD = 0.4;%4.17;
    KizmpcmxnD = 0.0000;
    KdzmpcmxnD = 0.014;%20.02;


    KpcmynD = 0.5;%4;
    KicmynD = 0.0000;
    KdcmynD = 0.014;%50.50;    
    KpcmxnD = 0.5;%4;
    KicmxnD = 0.0000;
    KdcmxnD = 0.014;%50.50;
    
        
    KpzmpcmynLmax = 0.4;%4.17;
    KizmpcmynLmax = 0.0000;
    KdzmpcmynLmax = 0.014;%20.02;    
    KpzmpcmxnLmax = 0.4;%4.17;
    KizmpcmxnLmax = 0.0001;
    KdzmpcmxnLmax = 0.014;%20.02;

    KpcmynLmax = 0.5;%4;
    KicmynLmax = 0.0000;
    KdcmynLmax = 0.014;%50.50;    
    KpcmxnLmax = 0.5;%4;
    KicmxnLmax = 0.0000;
    KdcmxnLmax = 0.014;%50.50;
    
    
    KpzmpcmynDmax = 0.4;%4.17;
    KizmpcmynDmax = 0.0000;
    KdzmpcmynDmax = 0.014;%20.02;    
    KpzmpcmxnDmax = 0.4;%4.17;
    KizmpcmxnDmax = 0.0000;
    KdzmpcmxnDmax = 0.014;%20.02;


    KpcmynDmax = 0.5;%4;
    KicmynDmax = 0.0000;
    KdcmynDmax = 0.014;%50.50;    
    KpcmxnDmax = 0.5;%4;
    KicmxnDmax = 0.0000;
    KdcmxnDmax = 0.014;%50.50;
    
%     
%     
%     Kpzmpcmyt = 0.7;%0.016;
%     Kizmpcmyt = 0.0;%0.00;
%     Kdzmpcmyt = 0.001;%0.02;    
%     Kpzmpcmxt = 0.7;%0.016;
%     Kizmpcmxt = 0.0;%0.00;
%     Kdzmpcmxt = 0.001;%0.02;
% 
%     Kpcmyt = 0.05;%0.022;
%     Kicmyt = 0.00;%0.00;
%     Kdcmyt = 0.002;%0.012;    
%     Kpcmxt = 0.05;%0.022;
%     Kicmxt = 0.00;%0.00;
%     Kdcmxt = 0.002;%0.012;

    Kpzmpcmyt = 3.7;%0.016;
    Kizmpcmyt = 0.0;%0.00;
    Kdzmpcmyt = 0.04;%0.02;    
    Kpzmpcmxt = 3.7;%0.016;
    Kizmpcmxt = 0.0;%0.00;
    Kdzmpcmxt = 0.04;%0.02;

    Kpcmyt = 0.05;%0.022;
    Kicmyt = 0.00;%0.00;
    Kdcmyt = 0.02;%0.012;    
    Kpcmxt = 0.05;%0.022;
    Kicmxt = 0.00;%0.00;
    Kdcmxt = 0.02;%0.012;
    
    
    LevaNogaZMPReg=0;
    DesnaNogaZMPReg=0;
    broj_koraka=1;
    skreni=1;
    
    podesi_koeficijente_za_sinhronizaciju=0;
    
%     brit=8000;


%%Nelinearno upravljanje

   NL_Control=1;
   v1=zeros(N,1);%I dejstvo, sliding mode drugog reda
   
  eta=4;
%     k=zeros(N,1);
%     k(actseg)=abs(dest(actseg,i))+eta(actseg)';
%     K=eye(N).*[k k k k k k k k k k k k k k k k k k k k k k k k k k k k k k];
%     kI=2*[75 75 75 75 75 75 175 175 25 175 25 175 25 75 75 75 75 75 75 75 75 75 75 75 75 75 75 75 75 75 ]';
%     KI=eye(N).*[kI kI kI kI kI kI kI kI kI kI kI kI kI kI kI kI kI kI kI kI kI kI kI kI kI kI kI kI kI kI];
    
    KI=1.50;
%     k0=2*[0 0 0 0 0 0 900 900 900 900 900 900 900 0 0 986.9333 986.9333 986.9333 986.9333 986.9333 986.9333 986.9333 986.9333 986.9333 986.9333 986.9333 986.9333 986.9333 986.9333 0]';%pojacanje P dejstva pole placement
%     K0=eye(N).*[k0 k0 k0 k0 k0 k0 k0 k0 k0 k0 k0 k0 k0 k0 k0 k0 k0 k0 k0 k0 k0 k0 k0 k0 k0 k0 k0 k0 k0 k0];
    
%     k1=3*[0 0 0 0 0 0 40 40 40 40 40 40 40 0 0 43.98 43.98 43.98 43.98 43.98 43.98 43.98 43.98 43.98 43.98 43.98 43.98 43.98 43.98 0]';%pojacanje D dejstva pole placement
%     K1=eye(N).*[k1 k1 k1 k1 k1 k1 k1 k1 k1 k1 k1 k1 k1 k1 k1 k1 k1 k1 k1 k1 k1 k1 k1 k1 k1 k1 k1 k1 k1 k1];
   
    fi=3;%*[0.5 0.5 0.5 0.5 0.5 0.5 1.8 1.8 2 1.8 5 1.8 1.8 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5]';
    sigmaI=0;
%     beta=zeros(N,1);
%     
%     s=zeros(N,brit);%povrsina klizanja
%     b=[0 0 0 0 0 0 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]';
%     
%     qe(:,1)=zeros(N,1);%greska pozicije na nivou unutrsnjih koordinata
%     pozLe(:,1)=pozL(:,1)-pozLd(:,1);%greska pozicije TCP
%     pozKe(:,1)=pozK(:,1)-pozKd(:,1);%greska pozicije TCP
%     pozVe(:,1)=pozV(:,1)-pozVd(:,1);%greska pozicije TCP
   
    v = zeros(N,brit);%novo upravljanje
    l1=230;%pojacanja estimatora poremecaja
    l2=2141.0;
    l3=84.5290;

    Xest(:,i)=zeros(3*N,1);%inicijalizacija vektora estimiranih stanja
    Xest(1:N,i)=qnew(:,1);%pocetna estimirana stanja jednaka stvarnim
    Xest(N+1:2*N,i)=zeros(N,1);%pocetne estimirane brzine su nula
    Xest(2*N+1:3*N,i)=zeros(N,1);%pocetne vrednosti estimiranih poremecaja

    Xpest=fliergeneral_est(t(i),Xest(:,i));%poziva se funkcija koja racuna izvod vektora estimiranih stanja u prvoj iteraciji
    qnewest(:,i)=Xest(1:N,i);
    qvnewest(:,i)=Xest(N+1:2*N,i);
    qanewest(:,i)=Xpest(N+1:2*N,i);
    dest(:,i)=Xest(2*N+1:3*N,i);
    destd(:,i)=Xpest(2*N+1:3*N,i);    
    
    
    
    
% % % % % %    %%SIMULACIJA - NASTAVAK
    load('sim1_faza1.5b1v1d1n.mat');    
    tk=20;%ukupno vreme simulacije
    brit=round(tk/dt);%Ukupan broj iteracija
    t=(0:brit-1)*dt;%ukupno vreme simulacije
    Up=[Up(:,1:ikraj+1) ones(size(Up,1),brit-ikraj-1)];
    resetprimitivcontrolvar;
% % % %  


%     cilj = [1.5+0.5,1.5,pi/2;3+0.5,1.5,0;6+0.5,0,-pi/2;];
    cilj = [1.5+0.5,1.5,pi/2;3+0.5,0,-pi/2;6+0.5,-2,0;8+0.5,0,0;];
    brojkoraka = 1;
    brzinakoraka_nom=1;
    visinakoraka=1;
    duzinakoraka_nom=1;
    skretanje=0;
    zeljeni_pravac =  qsave(6);
    zeljeni_pravac_L =  qsave(6);
    zeljeni_pravac_D =  qsave(6);
    
    k_brz = 1;
    k_duz = 1;
    brzinakoraka=k_brz*brzinakoraka_nom;
    duzinakoraka=k_duz*duzinakoraka_nom;
    
    last_min_distance_index=20000;
    last_min_distance_index_D=20000;
    last_min_distance_index_L=20000;
    
    igen = i;
    put_i=0;
    while i<brit-2
        %%%Sledeci trenutak
        i=i+1;
        itr=i;
        t(i)=t(i-1)+dt;
        t(i)
        if(t(i)>0.05 && faza == -1)
            faza=0;
        end
            
%         disp(t(i));

        %%%Integracija i diferencijacija
        X(:,i)=rkode4(t(i-1),X(:,i-1),Xp(:,i-1),dt,@fliergeneral12,bk);    
        Xp(:,i)=fliergeneral12(t(i),X(:,i));    

        %%Izvlacenje promenljivih iz vektora stanja
        dd(:)=X(3*N-5:3*N+30,i);
        del(:,:,i) =dd;
        FF(:)=Xp(4*N+25:4*N+60,i);
        Fnew(:,:,i)=FF;      

% % %         nkold=nknew;
% % %         nknew=sum(Fnew(3,:,i)>0);
% % %         if nknew~=nkold
% % %             i=i-1;
% % %             bk=16;
% % %             disp('nazad');
% % %             continue;
% % %         else
% % %             bk=1;
% % %         end
% % %         
        %%Provera da li je doslo do promene kontakta
        %%Odredjivanje trenutka kontakta pomocu slack variable
        %%Integracija do tog trenutka

        %%Nastavak izvlacenje promenljivih iz vektora stanja
        qnew(:,i)  =X(1:N,i);
        qvnew(:,i) =X(N+1:2*N,i);
        Irnew(:,i) =X(2*N+1:3*N-6,i);
        dd(:)      =X(3*N-5:3*N+30,i);
        del(:,:,i) =dd;
        qanew(:,i) =Xp(N+1:2*N,i);
        Irpnew(:,i)=Xp(2*N+1:3*N-6,i);
        ddp(:)     =Xp(3*N-5:3*N+30,i); 
        delp(:,:,i)=ddp;
        taunew(:,i)=Xp(3*N+31:4*N+24,i);
        FF(:)      =Xp(4*N+25:4*N+60,i);
        Fnew(:,:,i)=FF;
        lam(:,i)   =Xp(4*N+61:4*N+72,i);               
        [rcm(:,i) vcm(:,i)]=k_CentarMase(flier);   

        %%Izracunavanje ZMP-a
        ZMPpos(:,i)=ZMPb(Fuks);
        ZMPtr=ZMPpos(:,i);
        CMtr=rcm(:,i);

        %%Upravljanje
        Upr(:,i)=PRO*Up(:,i);
        dUpr(:,i)=0*Upr(:,i);

        qvzelj(:,i) = zeros(N,1);     
    %     if(i<1000)
    %         qvzelj(segtrupy,i)=qvzelj(segtrupy,i-1)+0.0001;
    %     elseif(i<2000)
    %         qvzelj(segtrupy,i)=qvzelj(segtrupy,i-1)-0.0001;
    %     end

        rbaz = km_linkX(flier,13,'xyz');
        rskzl = km_linkX(flier,14,'xyz');
        rskzd = km_linkX(flier,15,'xyz');
        rkukl = km_linkX(flier,16,'xyz');
        rkukd = km_linkX(flier,17,'xyz');
        rpetal = km_linkX(flier,18,'xyz');
        rpetad = km_linkX(flier,19,'xyz');
        rramel = km_linkX(flier,20,'xyz');
        rramed = km_linkX(flier,21,'xyz');
        rrukad = km_linkX(flier,22,'xyz');
        rrukal = km_linkX(flier,23,'xyz');
        rtrup = km_linkX(flier,24,'xyz');
        rprstl = km_linkX(flier,25,'xyz');
        rprstd = km_linkX(flier,26,'xyz');      
         
        [Jbaza Abaza] = km_jakP(flier,13);
        [Jlpeta Alpeta] = km_jakP(flier,18);
        [Jdpeta Adpeta] = km_jakP(flier,19);

        if(i>1200)
            aba=5;
        end
        
        

        if(last_min_distance_index>9500)
            put_i=put_i+1;
            c_pose(1)=rbaz(1);
            c_pose(2)=rbaz(2);
            A=flier.A;
            rbazor = TRotXYZ(A(:,:,6));
            c_pose(3) = rbazor(3);
            goal = [8 0 0];
            [xp yp putanja] = search_point_for_clatoid(c_pose, goal);
            last_min_distance_index = 1;
            last_min_distance_index_L = 1;
            last_min_distance_index_D = 1;
            igen = i+size(putanja,2);
            putanja_s(:,:,put_i) = putanja;
        end
        
        
        
%         if(i==2)            
%             [zeljeni_pravac duzinakoraka putanja] = mytest(flier, (rbaz+[0 0 0 0 0 0]'), [3 0 0]); 
%             last_min_distance_index = 1;
%             last_min_distance_index_L = 1;
%             last_min_distance_index_D = 1;
%             putanja_s(:,:,1) = putanja;
%         elseif(i==5000)            
%             [zeljeni_pravac duzinakoraka putanja] = mytest(flier, (rbaz+[0 0 0 0 0 0]'), cilj(1,:)); 
%             last_min_distance_index = 1;
%             last_min_distance_index_L = 1;
%             last_min_distance_index_D = 1;
%             putanja_s(:,:,2) = putanja;
%             visinakoraka  = 1;
%         elseif(i==15000)            
%             [zeljeni_pravac duzinakoraka putanja] = mytest(flier, [putanja(1:3,min_distance_index)], cilj(2,:));     
%             last_min_distance_index = 1;
%             last_min_distance_index_L = 1;
%             last_min_distance_index_D = 1;      
%             putanja_s(:,:,3) = putanja;
%             visinakoraka  = 1;
%         elseif(i==25000)            
%             [zeljeni_pravac duzinakoraka putanja] = mytest(flier, [putanja(1:3,min_distance_index)], cilj(3,:));     
%             last_min_distance_index = 1;
%             last_min_distance_index_L = 1;
%             last_min_distance_index_D = 1;   
%             putanja_s(:,:,4) = putanja;
%             visinakoraka  = 1;
%         elseif(i==45000)            
%             [zeljeni_pravac duzinakoraka putanja] = mytest(flier, [putanja(1:3,min_distance_index)], cilj(4,:));     
%             last_min_distance_index = 1;
%             last_min_distance_index_L = 1;
%             last_min_distance_index_D = 1;   
%             putanja_s(:,:,5) = putanja;
%             visinakoraka  = 1;
%         end        
        
        
        x_robot=rbaz(1);
        y_robot=rbaz(2);
        x_robot_L=rpetal(1);
        y_robot_L=rpetal(2);
        x_robot_D=rpetad(1);
        y_robot_D=rpetad(2);
        A=flier.A;
        rbazor = TRotXYZ(A(:,:,6));

        min_distance=sqrt((putanja(1,1)-x_robot)^2+(putanja(2,1)-y_robot)^2);
        min_distance_L=sqrt((putanja(1,1)-x_robot_L)^2+(putanja(2,1)-y_robot_L)^2);
        min_distance_D=sqrt((putanja(1,1)-x_robot_D)^2+(putanja(2,1)-y_robot_D)^2);
        min_distance_index=1;
        min_distance_index_L=1;
        min_distance_index_D=1;

        for iii=max(last_min_distance_index-1000,1):min(last_min_distance_index+2000,10000+1)
           distance=sqrt((putanja(1,iii)-x_robot)^2+(putanja(2,iii)-y_robot)^2);
           distance_L=sqrt((putanja(1,iii)-x_robot_L)^2+(putanja(2,iii)-y_robot_L)^2);
           distance_D=sqrt((putanja(1,iii)-x_robot_D)^2+(putanja(2,iii)-y_robot_D)^2);
%            l=[x_robot-putanja(1,iii),y_robot-putanja(2,iii)]';
%            vec_distance(:,iii)=[x_robot-putanja(1,iii),y_robot-putanja(2,iii)]'; %#ok<AGROW>
%            vec_distance_L(:,iii)=[x_robot_L-putanja(1,iii),y_robot_L-putanja(2,iii)]'; %#ok<AGROW>
%            vec_distance_D(:,iii)=[x_robot_D-putanja(1,iii),y_robot_D-putanja(2,iii)]'; %#ok<AGROW>
%            poz_ort(:,iii) = vec_distance(:,iii)/VecMod2(vec_distance(:,iii)); %#ok<AGROW>
%            poz_ort_L(:,iii) = vec_distance_L(:,iii)/VecMod2(vec_distance_L(:,iii)); %#ok<AGROW>
%            poz_ort_D(:,iii) = vec_distance_D(:,iii)/VecMod2(vec_distance_D(:,iii)); %#ok<AGROW>
           if (distance<min_distance)
               min_distance=distance;
               min_distance_index=iii;
               last_min_distance_index = min_distance_index;
           end
           if (distance_L<min_distance_L)
               min_distance_L=distance_L;
               min_distance_index_L=iii;
               last_min_distance_index_L = min_distance_index_L;
           end
           if (distance_D<min_distance_D)
               min_distance_D=distance_D;
               min_distance_index_D=iii;
               last_min_distance_index_D = min_distance_index_D;
           end
        end
        
        
        min_distance_index = min(min_distance_index,10000-1);
        min_distance_index_L = min(min_distance_index_L,10000-1);
        min_distance_index_D = min(min_distance_index_D,10000-1);
        
        ugao=putanja(3,min_distance_index);
        ugao_L=putanja(3,min_distance_index_L);
        ugao_D=putanja(3,min_distance_index_D);
        
        curv=putanja(4,min_distance_index);
        curv_L=putanja(4,min_distance_index_L);
        curv_D=putanja(4,min_distance_index_D);
        
        if(brojkoraka>1)            
            k_brz = ((1/abs(curv))/2)^(1/3);
            k_brz = saturacija2(k_brz,0.5,1.5);
        end
        brzinakoraka = k_brz*brzinakoraka_nom;        
        
        k_duz = 0.5*(k_brz-0.5)+0.75;
        
        duzinakoraka = k_duz*duzinakoraka_nom;
        
        if(brzinakoraka>1.2)
            brzinakoraka=1.2;
        end
        if(brzinakoraka<0.9)
            brzinakoraka=0.9;
        end
        if(duzinakoraka>1.1)
            duzinakoraka=1.1;
        end
        if(duzinakoraka>0.9)
            duzinakoraka=0.9;
        end
        
        
        vektor=[x_robot-putanja(1,min_distance_index) y_robot-putanja(2,min_distance_index) 0];
        vektor_L=[x_robot_L-putanja(1,min_distance_index_L) y_robot_L-putanja(2,min_distance_index_L) 0];
        vektor_D=[x_robot_D-putanja(1,min_distance_index_D) y_robot_D-putanja(2,min_distance_index_D) 0];
        tangenta=[cos(ugao) sin(ugao) 0]';
        tangenta_L=[cos(ugao_L) sin(ugao_L) 0]';
        tangenta_D=[cos(ugao_D) sin(ugao_D) 0]';
        vek_proizvod=cross(vektor,tangenta);
        vek_proizvod_L=cross(vektor_L,tangenta_L);
        vek_proizvod_D=cross(vektor_D,tangenta_D);
        
        smer=sign(vek_proizvod(1,3));                       
        smer_L=sign(vek_proizvod_L(1,3));                       
        smer_D=sign(vek_proizvod_D(1,3));                       
        
        Kp=1*smer;
        Kp_L=1*smer_L;
        Kp_D=1*smer_D;

        at = atan2(min_distance,duzinakoraka);
        at_L = atan2(min_distance_L,duzinakoraka);
        at_D = atan2(min_distance_D,duzinakoraka);
        
        zeljeni_pravac=ugao+Kp*at;
        zeljeni_pravac_L_new=ugao_L+Kp_L*at_L;
        zeljeni_pravac_D_new=ugao_D+Kp_D*at_D;
        
        ett_zp = exp(-dt/0.05);

        zeljeni_pravac_L=zeljeni_pravac_L_new*(1-ett_zp)+ett_zp*zeljeni_pravac_L;
        zeljeni_pravac_D=zeljeni_pravac_D_new*(1-ett_zp)+ett_zp*zeljeni_pravac_D;
 
        
%         zeljeni_pravac_L=ugao_L-0*smer_L*Kp*(0-atan2(min_distance_L-0.1103,duzinakoraka));
%         zeljeni_pravac_D=ugao_D-0*smer_D*Kp*(0-atan2(min_distance_D-0.1103,duzinakoraka));
        
        zp(:,i) = zeljeni_pravac;
        zp_L(:,i) = zeljeni_pravac_L;
        zp_D(:,i) = zeljeni_pravac_D;
        
        brz_rec(:,i) = brzinakoraka;
        duz_rec(:,i) = duzinakoraka;
        
        rbaz_zelj(1) = rbaz(1)+brzinakoraka*cos(zeljeni_pravac)*dt;% putanja(1,min_distance_index-100);
        rbaz_zelj(2) = rbaz(2)+brzinakoraka*sin(zeljeni_pravac)*dt;% putanja(2,min_distance_index-100);
                
                
        
        
        CMzelj(:,i)=CMzelj(:,i-1);
        
        if(faza==0)
            sacuvaj=1;
            faza=0.5;
            resetprimitivcontrolvar;
        end
        if(faza==0.5)   
            CMzelj(:,i)=[rpetal(1:2)+(rprstl(1:2)-rpetal(1:2))*0.5; 0];
            CMkr = rotx(rpetal(4))*roty(rpetal(5))*rotz(rpetal(6))*[0 -0.03*satlin(brzinakoraka-1) 0]';
            CMzelj(:,i) = CMzelj(:,i)+CMkr;
            promenacm=1;    
            [qvzelj statusnanog]= primosloninanogu2(qvzelj,flier,0.2,0.93,CMzelj(:,i),i,'L');
            if(statusnanog==3)
                qvzelj(:,i) = zeros(N,1); 
                sacuvaj=1;
                faza=1.5;
                ikraj=i;
                resetprimitivcontrolvar;
                statusopruzi=0;
                statusgrci=0;    
                statusnanog=0;
                statusispnap=0;
            end            
        end
        if(faza==1)   
            CMzelj(:,i)=[rpetal(1:2)+(rprstl(1:2)-rpetal(1:2))*0.25; 0];
            CMkr = rotx(rpetal(4))*roty(rpetal(5))*rotz(rpetal(6))*[0 -0.03*satlin(brzinakoraka-1) 0]';
            CMzelj(:,i) = CMzelj(:,i)+CMkr;
%             [qvzelj statusnanog]= primosloninanogu(qvzelj,flier,brzinakoraka*0.2,0.94,CMzelj(:,i),i,'D');
            [qvzelj statusnanog]= primosloninanogu2(qvzelj,flier,brzinakoraka*Kb_osl_L,Kv_osl_L,CMzelj(:,i),i,'L');

            if(statusnanog==3)
                qvzelj(:,i) = zeros(N,1); 
                sacuvaj=1;
                faza=1.5;
                ikraj=i;
                resetprimitivcontrolvar;
                statusopruzi=0;
                statusgrci=0;    
                statusnanog=0;
                statusispnap=0;
            end            
        end
        if(faza==1.5)
                CMzelj(:,i)=[rpetal(1:2)+(rprstl(1:2)-rpetal(1:2))*0.5; 0];
                CMkr = rotx(rpetal(4))*roty(rpetal(5))*rotz(rpetal(6))*[0 -0.03*satlin(brzinakoraka-1) 0]';
                CMzelj(:,i) = CMzelj(:,i)+CMkr;
                
                promenacm=2; 
                [qvzelj statusispnap]= primnaginjanjenapred(qvzelj,flier,brzinakoraka*Kb_nn_L,Kv_nn_L,CMzelj(:,i),i,'L'); 
                [qvzelj statusgrci]= primgrcinogu(qvzelj,flier,brzinakoraka*Kb_gn_D,visinakoraka*Kv_gn_D,Kd_gn_D,i,'D');         
                
                if(statusnaginjanjenapred<5 && ikrajnn_L==0)
                    ikrajnn_L=i;
                end
                
                if(statusgrci==3 && ikrajgd_D==0)
                    ikrajgd_D=i;
                end
                
            if(statusgrci==3)
                qvzelj(:,i) = zeros(N,1); 
                sacuvaj=1;
                faza=2;
                ikraj=i;
                resetprimitivcontrolvar;
                statusopruzi=0;
                statusgrci=0;    
                statusnanog=0;
                statusispnap=0;                
                if(ikrajnn_L==0)
                    ikrajnn_L=i;
                end                
                if(ikrajgd_D==0)
                    ikrajgd_D=i;
                end               
            end

        end
        if(faza==2)
%             if(statusnaginjanjenapred==5||statusnaginjanjenapred==0)
%                 if(skreni==1)
%                     if(broj_koraka<3)
%                         zeljeni_pravac=zeljeni_pravac+deg2rad(-10);
%                     elseif(broj_koraka<6)
%                         zeljeni_pravac=zeljeni_pravac+deg2rad(10);
%                     end
%                     skreni=0;                    
%                 end
                CMzelj(:,i)=[rprstl(1:2); 0];
                CMkr = rotx(rpetal(4))*roty(rpetal(5))*rotz(rpetal(6))*[0 -0.03*satlin(brzinakoraka-1) 0]';
                CMzelj(:,i) = CMzelj(:,i)+CMkr;                
                promenacm=3;
                [qvzelj nekoristise]= primnaginjanjenapred(qvzelj,flier,brzinakoraka*Kb_nn2_L,Kv_nn2_L,CMzelj(:,i),i,'L');
                [qvzelj statusopruzi]= primopruzinogu(qvzelj,flier,brzinakoraka*Kb_on_D,visinakoraka*Kv_on_D,duzinakoraka*Kd_on_D,Ku_on_D,0,i,'D');
                
                if(statusnaginjanjenapred<5 && ikrajnn2_L==0)
                    ikrajnn2_L=i;
                end
                
                if(statusopruzi==3 && ikrajon_D==0)
                    ikrajon_D=i;
                end
                
                %qvzelj(seglnoga(8),i)
%             end
            if(statusnaginjanjenapred==3 || statusnaginjanjenapred==4)
                qvzelj(:,i) = zeros(N,1); 
                skreni=1;
                sacuvaj=1;
                faza=2.5;
                ikraj=i;
                resetprimitivcontrolvar;
                statusopruzi=0;
                statusnaginjanjenapred=0;    
                statusispnap=0;
                
                if(ikrajnn2_L==0)
                    ikrajnn2_L=i;
                end                
                if(ikrajon_D==0)
                    ikrajon_D=i;
                end                
                
            end  
        end
        if(faza==2.5)
            CMzelj(:,i)=[rprstl(1:2)+(rpetad(1:2)-rprstl(1:2))*0.6; 0]; 
            CMkr = rotx(rpetal(4))*roty(rpetal(5))*rotz(rpetal(6))*[0 -0*satlin(brzinakoraka-1) 0]';
            CMzelj(:,i) = CMzelj(:,i)+CMkr;            
            promenacm=4;
            [qvzelj nekoristise]= primnaginjanjenapred(qvzelj,flier,brzinakoraka*Kb_nn25_L,Kv_nn25_L,CMzelj(:,i),i,'L');
            [qvzelj statusopruzi]= primopruzinogu(qvzelj,flier,brzinakoraka*Kb_on25_D,duzinakoraka^3*Kv_on25_D,duzinakoraka*Kd_on25_D,Ku_on25_D,0,i,'D');                

            if(statusnaginjanjenapred==3 && ikrajnn25_L==0)
                ikrajnn25_L=i;
            end

            if(statusopruzi==3 && ikrajon25_D==0)
                ikrajon25_D=i;
            end            
            
            if(statusnaginjanjenapred==3 && statusopruzi==3)
                qvzelj(:,i) = zeros(N,1); 
                sacuvaj=1;
                faza=3;
                ikraj=i;
                resetprimitivcontrolvar;
                statusopruzi=0;
                statusnaginjanjenapred=0;    
                statusispnap=0;
                
                if(ikrajnn25_L==0)
                    ikrajnn25_L=i;
                end
                if(ikrajon25_D==0)
                    ikrajon25_D=i;
                end 
                
            end  
        end
        if(faza==3)
            CMzelj(:,i)=[rprstl(1:2)+(rpetad(1:2)-rprstl(1:2))*0.6; 0];            
%             CMzelj(:,i)=[rpetad(1:2)-(rpetad(1:2)-rprstl(1:2))*0; 0];
            CMkr = rotx(rpetal(4))*roty(rpetal(5))*rotz(rpetal(6))*[0 -0*satlin(brzinakoraka-1) 0]';
            CMzelj(:,i) = CMzelj(:,i)+CMkr;            

            [qvzelj statusopruzinogu]= primospusti2(qvzelj,flier,brzinakoraka*Kb_sn_D,Kv_sn_D,0,i,CMzelj(:,i),'D');
            promenacm=5;
            if(statusopruzinogu==3)
                qvzelj(:,i) = zeros(N,1); 
                sacuvaj=1;
                faza=4;
%                 faza=7;                
                ikraj=i;
                resetprimitivcontrolvar;
                statusopruzinogu=0;
%                 brzinakoraka=brzinakoraka+0.25;                   
%                 duzinakoraka=duzinakoraka+0.2;
%                 if(brzinakoraka_nom<1.3)
%                     brzinakoraka_nom = brzinakoraka_nom+0.1;
%                 end                

            end
        end
        if(faza==4)
                       
%             if(i-ikraj>200)
                CMzelj(:,i)=[rpetad(1:2)+(rprstd(1:2)-rpetad(1:2))*0.25; 0];
                CMkr = rotx(rpetad(4))*roty(rpetad(5))*rotz(rpetad(6))*[0 0.03*satlin(brzinakoraka-1) 0]';
                CMzelj(:,i) = CMzelj(:,i)+CMkr;                    
%                 CMzelj(:,i)=[rpetad(1:2)-(rpetad(1:2)-rprstl(1:2))*0; 0];
%                 [qvzelj statusnanog]= primosloninanogu(qvzelj,flier,brzinakoraka*0.2,0.94,CMzelj(:,i),i,'D');
                [qvzelj statusnanog]= primosloninanogu2(qvzelj,flier,brzinakoraka*Kb_osl_D,Kv_osl_D,CMzelj(:,i),i,'D');

                promenacm=7;    

%             end
            if(statusnanog==3)
                qvzelj(:,i) = zeros(N,1); 
                sacuvaj=1;
                faza=4.5;                        
%                 faza=7;                        
                ikraj=i;
                resetprimitivcontrolvar;
                statusopruzi=0;
                statusgrci=0;    
                statusnanog=0;
                statusispnap=0;
            end

        end
        if(faza==4.5)            
            CMzelj(:,i)=[rpetad(1:2)+(rprstd(1:2)-rpetad(1:2))*0.5; 0];
                CMkr = rotx(rpetad(4))*roty(rpetad(5))*rotz(rpetad(6))*[0 0.03*satlin(brzinakoraka-1) 0]';
                CMzelj(:,i) = CMzelj(:,i)+CMkr;                    
            
            promenacm=7.5;
            [qvzelj statusispnap]= primnaginjanjenapred(qvzelj,flier,brzinakoraka*Kb_nn_D,Kv_nn_D,CMzelj(:,i),i,'D'); 
            [qvzelj statusgrci]= primgrcinogu(qvzelj,flier,brzinakoraka*Kb_gn_L,visinakoraka*Kv_gn_L,Kd_gn_L,i,'L');                    

            if(statusnaginjanjenapred<5 && ikrajnn_D==0)
                ikrajnn_D=i;
            end

            if(statusgrci==3 && ikrajgn_L==0)
                ikrajgn_L=i;
            end    
            
            if(statusgrci==3)
                qvzelj(:,i) = zeros(N,1); 
                sacuvaj=1;
                faza=5;
                ikraj=i;
                resetprimitivcontrolvar;
                statusopruzi=0;
                statusgrci=0;    
                statusnanog=0;
                statusispnap=0;
                if(ikrajnn_D==0)
                    ikrajnn_D=i;
                end

                if(ikrajgn_L==0)
                    ikrajgn_L=i;
                end                  
            end
        end        
        if(faza==5)
%             if(statusnaginjanjenapred==5||statusnaginjanjenapred==0)
%                 if(skreni==1)
%                     if(broj_koraka<3)
%                         zeljeni_pravac=zeljeni_pravac+deg2rad(-10);
%                     elseif(broj_koraka<6)
%                         zeljeni_pravac=zeljeni_pravac+deg2rad(10);
%                     end
%                     skreni=0;
%                 end
                CMzelj(:,i)=[rprstd(1:2); 0];
                CMkr = rotx(rpetad(4))*roty(rpetad(5))*rotz(rpetad(6))*[0 0.03*satlin(brzinakoraka-1) 0]';
                CMzelj(:,i) = CMzelj(:,i)+CMkr;                    
                
                promenacm=3;
                [qvzelj nekoristise]= primnaginjanjenapred(qvzelj,flier,brzinakoraka*Kb_nn5_D,Kv_nn5_D,CMzelj(:,i),i,'D');
                [qvzelj statusopruzi]= primopruzinogu(qvzelj,flier,brzinakoraka*Kb_on_L,visinakoraka*Kv_on_L,duzinakoraka*Kd_on_L,Ku_on_L,0,i,'L');

                if(statusnaginjanjenapred<5 && ikrajnn5_D==0)
                    ikrajnn5_D=i;
                end

                if(statusopruzi==3 && ikrajon_L==0)
                    ikrajon_L=i;
                end   
                %qvzelj(seglnoga(8),i)
%             end
            if(statusnaginjanjenapred==3 || statusnaginjanjenapred==4)
                qvzelj(:,i) = zeros(N,1); 
                skreni=1;
                sacuvaj=1;
                faza=5.5;
                ikraj=i;
                resetprimitivcontrolvar;
                statusopruzi=0;
                statusnaginjanjenapred=0;    
                statusispnap=0;
                if(ikrajnn5_D==0)
                    ikrajnn5_D=i;
                end

                if(ikrajon_L==0)
                    ikrajon_L=i;
                end                   
            end  
        end
        if(faza==5.5)
            CMzelj(:,i)=[rprstd(1:2)+(rpetal(1:2)-rprstd(1:2))*0.6; 0];  
            CMkr = rotx(rpetad(4))*roty(rpetad(5))*rotz(rpetad(6))*[0 0.0*satlin(brzinakoraka-1) 0]';
            CMzelj(:,i) = CMzelj(:,i)+CMkr;                    

            promenacm=9;
            [qvzelj nekoristise]= primnaginjanjenapred(qvzelj,flier,brzinakoraka*Kb_nn55_D,Kv_nn55_D,CMzelj(:,i),i,'D');
            [qvzelj statusopruzi]= primopruzinogu(qvzelj,flier,brzinakoraka*Kb_on55_L,duzinakoraka^3*Kv_on55_L,duzinakoraka*Kd_on55_L,Ku_on55_L,0,i,'L');

            if(statusnaginjanjenapred==3 && ikrajnn55_D==0)
                ikrajnn55_D=i;
            end

            if(statusopruzi==3 && ikrajon55_L==0)
                ikrajon55_L=i;
            end               
            
            if(statusnaginjanjenapred==3 && statusopruzi==3)
                qvzelj(:,i) = zeros(N,1); 
                sacuvaj=1;
                faza=6;
                ikraj=i;
                resetprimitivcontrolvar;
                statusopruzi=0;
                statusnaginjanjenapred=0;    
                statusispnap=0;

                if(ikrajnn55_D==0)
                    ikrajnn55_D=i;
                end

                if(ikrajon55_L==0)
                    ikrajon55_L=i;
                end               
                            
            end  
        end
        
        if(faza==6)           
            CMzelj(:,i)=[rprstd(1:2)+(rpetal(1:2)-rprstd(1:2))*0.6; 0];
            CMkr = rotx(rpetad(4))*roty(rpetad(5))*rotz(rpetad(6))*[0 0.0*satlin(brzinakoraka-1) 0]';
            CMzelj(:,i) = CMzelj(:,i)+CMkr;                    
            
            [qvzelj statusopruzinogu]= primospusti2(qvzelj,flier,brzinakoraka*Kb_sn_L,Kv_sn_L,0,i,CMzelj(:,i),'L');                        
            brojkoraka = brojkoraka+1;
             promenacm=10;            
            if(statusopruzinogu==3)
                sacuvaj=1;
                faza=1;
                ikraj=i;
                resetprimitivcontrolvar;
                statusopruzinogu=0;
%                 duzinakoraka=duzinakoraka+0.2;
                broj_koraka= broj_koraka+1;                                
                
                podesi_koeficijente_za_sinhronizaciju=1;
                
%                 if(brzinakoraka_nom<1.3)
%                     brzinakoraka_nom = brzinakoraka_nom+0.1;
%                 end                
            end            
        end
        
        [qvzelj statruke] = primruke(qvzelj,flier,0.4,i);

        if(podesi_koeficijente_za_sinhronizaciju==1)
            podesi_koeficijente_za_sinhronizaciju=0;
%             if(brzinakoraka<2.5)
%                 brzinakoraka = brzinakoraka+0.1;
%             end
            if(broj_koraka>2)
                d_faza15=ikrajgd_D-ikrajnn_L;
                zd_faza15=0;
                
                %Prioritetni primitiv grcenje noge i njega menjamo 10% sporije
                if((d_faza15-zd_faza15)>0)%grcenje se zavrsilo pre naginjanja                
                    %Usporavamo grcenje ubrzavamo naginjanje
    %                 if(d_faza15<10)%ukoliko je greska u sinhronizaciji mala
    %                     d_faza15=10;%postavljamo je na konacnu vrednost
    %                 end                
                    Kb_gn_D = Kb_gn_D + 0.04*dt*abs(d_faza15-zd_faza15);
                    Kb_nn_L = Kb_nn_L - 0.04*dt*abs(d_faza15-zd_faza15);

                else%%naginjanjae se zavrsilo pre grcenja
                    %Ubrzavamo grcenje usporavamo naginjanje
    %                 if(d_faza15>-10)%ukoliko je greska u sinhronizaciji mala
    %                     d_faza15=-10;%postavljamo je na konacnu vrednost
    %                 end                
                    Kb_gn_D = Kb_gn_D - 0.04*dt*abs(d_faza15-zd_faza15);
                    Kb_nn_L = Kb_nn_L + 0.04*dt*abs(d_faza15-zd_faza15);                                
                end

                d_faza2=ikrajon_D-ikrajnn2_L;
                zd_faza2=0;
                %Prioritetni primitiv opruzanje noge i njega menjamo 10% sporije
                if((d_faza2-zd_faza2)>0)%opruzanje se zavrsilo pre naginjanja                
                    %Usporavamo opruzanje ubrzavamo naginjanje
    %                 if(d_faza2<10)%ukoliko je greska u sinhronizaciji mala
    %                     d_faza2=10;%postavljamo je na konacnu vrednost
    %                 end                
                    Kb_on_D = Kb_on_D + 0.04*dt*abs(d_faza2-zd_faza2);
                    Kb_nn2_L = Kb_nn2_L - 0.04*dt*abs(d_faza2-zd_faza2);

                else%%naginjanjae se zavrsilo pre opruzanje
                    %Ubrzavamo opruzanje usporavamo naginjanje
    %                 if(d_faza2>-10)%ukoliko je greska u sinhronizaciji mala
    %                     d_faza2=-10;%postavljamo je na konacnu vrednost
    %                 end                
                    Kb_on_D = Kb_on_D - 0.04*dt*abs(d_faza2-zd_faza2);
                    Kb_nn2_L = Kb_nn2_L + 0.04*dt*abs(d_faza2-zd_faza2);                                
                end

                d_faza25=ikrajon25_D-ikrajnn25_L;
                zd_faza25=0;
                %Prioritetni primitiv opruzanje noge i njega menjamo 10% sporije
                if((d_faza25-zd_faza25)>0)%grcenje se zavrsilo pre naginjanja                
                    %Usporavamo opruzanje ubrzavamo naginjanje
    %                 if(d_faza25<10)%ukoliko je greska u sinhronizaciji mala
    %                     d_faza25=10;%postavljamo je na konacnu vrednost
    %                 end                
                    Kb_on25_D = Kb_on25_D + 0.04*dt*abs(d_faza25-zd_faza25);
                    Kb_nn25_L = Kb_nn25_L - 0.04*dt*abs(d_faza25-zd_faza25);

                else%%naginjanjae se zavrsilo pre grcenja
                    %Ubrzavamo opruzanje usporavamo naginjanje
    %                 if(d_faza25>-10)%ukoliko je greska u sinhronizaciji mala
    %                     d_faza25=-10;%postavljamo je na konacnu vrednost
    %                 end                
                    Kb_on25_D = Kb_on25_D - 0.04*dt*abs(d_faza25-zd_faza25);
                    Kb_nn25_L = Kb_nn25_L + 0.04*dt*abs(d_faza25-zd_faza25);                                
                end

                d_faza45 = ikrajgn_L-ikrajnn_D;
                zd_faza45=0;
                %Prioritetni primitiv grcenje noge i njega menjamo 10% sporije
                if((d_faza45-zd_faza45)>0)%grcenje se zavrsilo pre naginjanja                
                    %Usporavamo grcenje ubrzavamo naginjanje
    %                 if(d_faza45<10)%ukoliko je greska u sinhronizaciji mala
    %                     d_faza45=10;%postavljamo je na konacnu vrednost
    %                 end                
                    Kb_gn_L = Kb_gn_L + 0.04*dt*abs(d_faza45-zd_faza45);
                    Kb_nn_D = Kb_nn_D - 0.04*dt*abs(d_faza45-zd_faza45);

                else%%naginjanjae se zavrsilo pre grcenja
                    %Ubrzavamo grcenje usporavamo naginjanje
    %                 if(d_faza45>-10)%ukoliko je greska u sinhronizaciji mala
    %                     d_faza45=-10;%postavljamo je na konacnu vrednost
    %                 end                
                    Kb_gn_L = Kb_gn_L - 0.04*dt*abs(d_faza45-zd_faza45);
                    Kb_nn_D = Kb_nn_D + 0.04*dt*abs(d_faza45-zd_faza45);                                
                end            

                d_faza5 = ikrajnn5_D-ikrajon_L;
                zd_faza5=0;
                %Prioritetni primitiv opruzanje noge i njega menjamo 10% sporije
                if((d_faza5-zd_faza5)>0)%opruzanje se zavrsilo pre naginjanja                
                    %Usporavamo opruzanje ubrzavamo naginjanje
    %                 if(d_faza5<10)%ukoliko je greska u sinhronizaciji mala
    %                     d_faza5=10;%postavljamo je na konacnu vrednost
    %                 end                
                    Kb_on_L = Kb_on_L - 0.04*dt*abs(d_faza5-zd_faza5);
                    Kb_nn5_D = Kb_nn5_D + 0.04*dt*abs(d_faza5-zd_faza5);

                else%%naginjanjae se zavrsilo pre opruzanje
                    %Ubrzavamo opruzanje usporavamo naginjanje
    %                 if(d_faza5>-10)%ukoliko je greska u sinhronizaciji mala
    %                     d_faza5=-10;%postavljamo je na konacnu vrednost
    %                 end                
                    Kb_on_L = Kb_on_L + 0.04*dt*abs(d_faza5-zd_faza5);
                    Kb_nn5_D = Kb_nn5_D - 0.04*dt*abs(d_faza5-zd_faza5);                                
                end

                d_faza55 = ikrajnn55_D-ikrajon55_L;
                zd_faza55=0;
                %Prioritetni primitiv opruzanje noge i njega menjamo 10% sporije
                if((d_faza55-zd_faza55)>0)%grcenje se zavrsilo pre naginjanja                
                    %Usporavamo opruzanje ubrzavamo naginjanje
    %                 if(d_faza55<10)%ukoliko je greska u sinhronizaciji mala
    %                     d_faza55=10;%postavljamo je na konacnu vrednost
    %                 end                
                    Kb_on55_L = Kb_on55_L - 0.04*dt*abs(d_faza55-zd_faza55);
                    Kb_nn55_D = Kb_nn55_D + 0.04*dt*abs(d_faza55-zd_faza55);

                else%%naginjanjae se zavrsilo pre grcenja
                    %Ubrzavamo opruzanje usporavamo naginjanje
    %                 if(d_faza55>-10)%ukoliko je greska u sinhronizaciji mala
    %                     d_faza55=-10;%postavljamo je na konacnu vrednost
    %                 end                
                    Kb_on55_L = Kb_on55_L + 0.04*dt*abs(d_faza55-zd_faza55);
                    Kb_nn55_D = Kb_nn55_D - 0.04*dt*abs(d_faza55-zd_faza55);                                
                end            
            
                save(['sim1_faza' num2str(broj_koraka) '.mat'],'d_faza15','d_faza2','d_faza25','d_faza45','d_faza5','d_faza55');
                save(['podeseni_parametri' num2str(broj_koraka) '.mat'],...
                        'Kb_osl_L','Kv_osl_L',...
                        'Kb_nn_L','Kv_nn_L','ikrajnn_L','Kb_gn_D','Kv_gn_D','Kd_gn_D','ikrajgd_D',...
                        'Kb_nn2_L','Kv_nn2_L','ikrajnn2_L','Kb_on_D','Kv_on_D','Kd_on_D','Ku_on_D','ikrajon_D',...
                        'Kb_nn25_L','Kv_nn25_L','ikrajnn25_L','Kb_on25_D','Kv_on25_D','Kd_on25_D','Ku_on25_D','ikrajon25_D',...
                        'Kb_sn_D','Kv_sn_D',...
                        'Kb_osl_D','Kv_osl_D',...
                        'Kb_nn_D','Kv_nn_D','ikrajnn_D','Kb_gn_L','Kv_gn_L','Kd_gn_L','ikrajgn_L',...
                        'Kb_nn5_D','Kv_nn5_D','ikrajnn5_D','Kb_on_L','Kv_on_L','Kd_on_L','Ku_on_L','ikrajon_L',...
                        'Kb_nn55_D','Kv_nn55_D','ikrajnn55_D','Kb_on55_L','Kv_on55_L','Kd_on55_L','Ku_on55_L','ikrajon55_L',...
                        'Kb_sn_L','Kv_sn_L');
            end

            
            ikrajnn_L = 0;
            ikrajgd_D = 0;
            
            ikrajnn2_L = 0;
            ikrajon_D = 0;
            
            ikrajnn25_L = 0;
            ikrajon25_D = 0;
            
            ikrajnn_D = 0;
            ikrajgn_L = 0;
            
            ikrajnn5_D = 0;
            ikrajon_L = 0;
            
            ikrajnn55_D = 0;
            ikrajon55_L = 0;            
            
        end

        
        
    %     qvzelj(seglnoga,i) = -pinv(Jlpeta(:,seglnoga))*qvb;
    %     qvzelj(segdnoga,i) = -pinv(Jdpeta(:,segdnoga))*qvb;
    
        if(promenacm~=promenacms)
            icmzeljtr=i;
%             dCMzelj = (CMzelj(:,i)-CMzeljs)/500;
            CMzeljs=CMzelj(:,i);
            promenacms=promenacm;
        end
        
        trbrit = 500/brzinakoraka;
        
        if(i-icmzeljtr<trbrit)
            dCMzelj = (CMzelj(:,i)-CMzeljt(:,i-1))/(trbrit-(i-icmzeljtr));
            CMzeljt(:,i)=CMzeljt(:,i-1)+dCMzelj;
        else
            CMzeljt(:,i)=CMzelj(:,i);
        end
    
        if(any(isnan(ZMPpos(:,i))))
            ZMPpos(:,i)=ZMPpos(:,i-1);
        end
        
        CMerr(:,i) = CMzeljt(:,i)-[rcm(1:2,i);0];
%         CMerr(:,i) = saturacija(CMerr(:,i),0.05);
        ZMPCMerr(:,i)=[rcm(1:2,i);0]-ZMPpos(:,i);
%         ZMPCMerr(:,i)=0.5*([rcm(1:2,i);0]-ZMPpos(:,i))+0.5*([CMzeljt(1:2,i);0]-ZMPpos(:,i));
%         ZMPCMerr(:,i)=CMzeljt(:,i)-ZMPpos(:,i);
%         ZMPCMerr(:,i) = saturacija(ZMPCMerr(:,i),0.02);
        
%         %Koeficijent odnosa izmedju ZMP i CM regulatora u nogama
%         if(VecMod(ZMPCMerr(:,i))>0.02)
%             xa=xa+0.05;
%             if(xa>3)
%                 xa=3;
%             end
%         else
%             xa=xa-0.05;
%             if(xa<-3)
%                 xa=-3;
%             end
%         end
% 
%             %Koeficijent odnosa izmedju ZMP i CM regulatora u nogama
%         if(VecMod(CMerr(:,i))>0.02)
%             xb=xb+0.05;
%             if(xb>3)
%                 xb=3;
%             end
%         else
%             xb=xb-0.05;
%             if(xb<-3)
%                 xb=-3;
%             end
%         end
%         
%         a=logsig(xa);
%         b=logsig(xb);

        %%%%IZVRSITI PROJEKCIJE GRESKA NA OSE ZGLOBOVA

        %%%%ADAPTIVNI KOEFICIJENTI  !!!!!!

%         if(abs(CMerr(1,i))>0.02)
%             Kpcmxn=Kpcmxn+0.05;
%         else
%             Kpcmxn=Kpcmxn-0.05;
%         end
%         
%         if(Kpcmxn>6)
%             Kpcmxn=3;
%         end
%         if(Kpcmxn<2)
%             Kpcmxn=2;
%         end
%         
%         
%         if(abs(CMerr(2,i))>0.02)
%             Kpcmyn=Kpcmyn+0.05;
%         else
%             Kpcmyn=Kpcmyn-0.05;
%         end
%         
%         if(Kpcmyn>6)
%             Kpcmyn=6;
%         end
%         if(Kpcmyn<2)
%             Kpcmyn=2;
%         end
%         Kpzmpcmyn = 0.17;
%         Kizmpcmyn = 0.0;
%         Kdzmpcmyn = 0.12;    
%         Kpzmpcmxn = 0.17;
%         Kizmpcmxn = 0.0;
%         Kdzmpcmxn = 0.12;
% 
%         
%         Kpcmyn = 1.36;
%         Kicmyn = 0.;
%         Kdcmyn = 2.050;    
%         Kpcmxn = 1.36;
%         Kicmxn = 0.;
%         Kdcmxn = 2.050;
 
%         Kpzmpcmyn = 0;
%         Kizmpcmyn = 0;
%         Kdzmpcmyn = 0;    
%         Kpzmpcmxn = 0;
%         Kizmpcmxn = 0;
%         Kdzmpcmxn = 0;       
%         
%         Kpcmyn = 0;
%         Kicmyn = 0;
%         Kdcmyn = 0;    
%         Kpcmxn = 0;
%         Kicmxn = 0;
%         Kdcmxn = 0;

%     %     a=1;
        if size(XCL,2)>2
%     %         if(rad2deg(abs(qnew(segtrupx(1),i)))<3)
%                 %ZMP
%                 qvzelj(21,i) = qvzelj(21,i) +(Kpzmpcmxn*ZMPCMerr(1,i)+Kizmpcmxn*sum(ZMPCMerr(1,i),2)+Kdzmpcmyn*(ZMPCMerr(1,i)-ZMPCMerr(1,i-1)));
%                 qvzelj(17,i) = qvzelj(17,i) -(Kpzmpcmxn*ZMPCMerr(1,i)+Kizmpcmxn*sum(ZMPCMerr(1,i),2)+Kdzmpcmyn*(ZMPCMerr(1,i)-ZMPCMerr(1,i-1)));
% 
%                 qvzelj(19,i) = qvzelj(19,i) -(Kpzmpcmyn*ZMPCMerr(2,i)+Kizmpcmyn*sum(ZMPCMerr(2,i),2)+Kdzmpcmyn*(ZMPCMerr(2,i)-ZMPCMerr(2,i-1)));
%                 qvzelj(16,i) = qvzelj(16,i) +(Kpzmpcmyn*ZMPCMerr(2,i)+Kizmpcmyn*sum(ZMPCMerr(2,i),2)+Kdzmpcmyn*(ZMPCMerr(2,i)-ZMPCMerr(2,i-1)));
% 
%                 %CM
%                 qvzelj(21,i) = qvzelj(21,i) -(Kpcmxn*CMerr(1,i)+Kicmxn*sum(CMerr(1,i),2)+Kdcmyn*(CMerr(1,i)-CMerr(1,i-1)));
%                 qvzelj(17,i) = qvzelj(17,i) +(Kpcmxn*CMerr(1,i)+Kicmxn*sum(CMerr(1,i),2)+Kdcmyn*(CMerr(1,i)-CMerr(1,i-1)));
% 
%                 qvzelj(19,i) = qvzelj(19,i) +(Kpcmyn*CMerr(2,i)+Kicmyn*sum(CMerr(2,i),2)+Kdcmyn*(CMerr(2,i)-CMerr(2,i-1)));
%                 qvzelj(16,i) = qvzelj(16,i) -(Kpcmyn*CMerr(2,i)+Kicmyn*sum(CMerr(2,i),2)+Kdcmyn*(CMerr(2,i)-CMerr(2,i-1)));
% 
%     %         end
            LevaNogaZMPReg=1;
        else
            LevaNogaZMPReg=0;

        end
        if size(XCR,2)>2
%     %         if(rad2deg(abs(qnew(segtrupx(1),i)))<3)
% 
%                 %ZMP
%                 qvzelj(13,i) = qvzelj(13,i) +(Kpzmpcmxn*ZMPCMerr(1,i)+Kizmpcmxn*sum(ZMPCMerr(1,i),2)+Kdzmpcmxn*(ZMPCMerr(1,i)-ZMPCMerr(1,i-1)));
%                 qvzelj( 9,i) = qvzelj( 9,i) -(Kpzmpcmxn*ZMPCMerr(1,i)+Kizmpcmxn*sum(ZMPCMerr(1,i),2)+Kdzmpcmxn*(ZMPCMerr(1,i)-ZMPCMerr(1,i-1)));
% 
%                 qvzelj(11,i) = qvzelj(11,i) -(Kpzmpcmyn*ZMPCMerr(2,i)+Kizmpcmyn*sum(ZMPCMerr(2,i),2)+Kdzmpcmyn*(ZMPCMerr(2,i)-ZMPCMerr(2,i-1)));
%                 qvzelj( 8,i) = qvzelj( 8,i) +(Kpzmpcmyn*ZMPCMerr(2,i)+Kizmpcmyn*sum(ZMPCMerr(2,i),2)+Kdzmpcmyn*(ZMPCMerr(2,i)-ZMPCMerr(2,i-1)));
% 
%                 %CM
%                 qvzelj(13,i) = qvzelj(13,i) -(Kpcmxn*CMerr(1,i)+Kicmxn*sum(CMerr(1,i),2)+Kdcmxn*(CMerr(1,i)-CMerr(1,i-1)));
%                 qvzelj( 9,i) = qvzelj( 9,i) +(Kpcmxn*CMerr(1,i)+Kicmxn*sum(CMerr(1,i),2)+Kdcmxn*(CMerr(1,i)-CMerr(1,i-1)));
% 
%                 qvzelj(11,i) = qvzelj(11,i) +(Kpcmyn*CMerr(2,i)+Kicmyn*sum(CMerr(2,i),2)+Kdcmyn*(CMerr(2,i)-CMerr(2,i-1)));
%                 qvzelj( 8,i) = qvzelj( 8,i) -(Kpcmyn*CMerr(2,i)+Kicmyn*sum(CMerr(2,i),2)+Kdcmyn*(CMerr(2,i)-CMerr(2,i-1)));
%     %         end
            DesnaNogaZMPReg=1;
        else
            DesnaNogaZMPReg=0;
        end

        delKZMPCM=100;
        if(LevaNogaZMPReg==1)
            if(KpzmpcmynL<KpzmpcmynLmax)
                KpzmpcmynL = KpzmpcmynL+KpzmpcmynLmax/delKZMPCM;%4.17;
                KizmpcmynL = KizmpcmynL+KizmpcmynLmax/delKZMPCM;
                KdzmpcmynL = KdzmpcmynL+KdzmpcmynLmax/delKZMPCM;%20.02;    
                KpzmpcmxnL = KpzmpcmxnL+KpzmpcmxnLmax/delKZMPCM;%4.17;
                KizmpcmxnL = KizmpcmxnL+KizmpcmxnLmax/delKZMPCM;
                KdzmpcmxnL = KdzmpcmxnL+KdzmpcmxnLmax/delKZMPCM;%20.02;


                KpcmynL = KpcmynL+KpcmynLmax/delKZMPCM;%4;
                KicmynL = KicmynL+KicmynLmax/delKZMPCM;
                KdcmynL = KdcmynL+KdcmynLmax/delKZMPCM;%50.50;    
                KpcmxnL = KpcmxnL+KpcmxnLmax/delKZMPCM;%4;
                KicmxnL = KicmxnL+KicmxnLmax/delKZMPCM;
                KdcmxnL = KdcmxnL+KdcmxnLmax/delKZMPCM;%50.50;
            else
                KpzmpcmynL = KpzmpcmynLmax;%4.17;
                KizmpcmynL = KizmpcmynLmax;
                KdzmpcmynL = KdzmpcmynLmax;%20.02;    
                KpzmpcmxnL = KpzmpcmxnLmax;%4.17;
                KizmpcmxnL = KizmpcmxnLmax;
                KdzmpcmxnL = KdzmpcmxnLmax;%20.02;


                KpcmynL = KpcmynLmax;%4;
                KicmynL = KicmynLmax;
                KdcmynL = KdcmynLmax;%50.50;    
                KpcmxnL = KpcmxnLmax;%4;
                KicmxnL = KicmxnLmax;
                KdcmxnL = KdcmxnLmax;%50.50;                
            end       
        else
             if(KpzmpcmynL>0)
                KpzmpcmynL = KpzmpcmynL-KpzmpcmynLmax/delKZMPCM;%4.17;
                KizmpcmynL = KizmpcmynL-KizmpcmynLmax/delKZMPCM;
                KdzmpcmynL = KdzmpcmynL-KdzmpcmynLmax/delKZMPCM;%20.02;    
                KpzmpcmxnL = KpzmpcmxnL-KpzmpcmxnLmax/delKZMPCM;%4.17;
                KizmpcmxnL = KizmpcmxnL-KizmpcmxnLmax/delKZMPCM;
                KdzmpcmxnL = KdzmpcmxnL-KdzmpcmxnLmax/delKZMPCM;%20.02;


                KpcmynL = KpcmynL-KpcmynLmax/delKZMPCM;%4;
                KicmynL = KicmynL-KicmynLmax/delKZMPCM;
                KdcmynL = KdcmynL-KdcmynLmax/delKZMPCM;%50.50;    
                KpcmxnL = KpcmxnL-KpcmxnLmax/delKZMPCM;%4;
                KicmxnL = KicmxnL-KicmxnLmax/delKZMPCM;
                KdcmxnL = KdcmxnL-KdcmxnLmax/delKZMPCM;%50.50;
             else
                KpzmpcmynL = 0;
                KizmpcmynL = 0;
                KdzmpcmynL = 0;    
                KpzmpcmxnL = 0;
                KizmpcmxnL = 0;
                KdzmpcmxnL = 0;


                KpcmynL = 0;
                KicmynL = 0;
                KdcmynL = 0;    
                KpcmxnL = 0;
                KicmxnL = 0;
                KdcmxnL = 0;                 
            end   
        end
        
        if(DesnaNogaZMPReg==1)
            if(KpzmpcmynD<KpzmpcmynDmax)
                KpzmpcmynD = KpzmpcmynD+KpzmpcmynDmax/delKZMPCM;%4.17;
                KizmpcmynD = KizmpcmynD+KizmpcmynDmax/delKZMPCM;
                KdzmpcmynD = KdzmpcmynD+KdzmpcmynDmax/delKZMPCM;%20.02;    
                KpzmpcmxnD = KpzmpcmxnD+KpzmpcmxnDmax/delKZMPCM;%4.17;
                KizmpcmxnD = KizmpcmxnD+KizmpcmxnDmax/delKZMPCM;
                KdzmpcmxnD = KdzmpcmxnD+KdzmpcmxnDmax/delKZMPCM;%20.02;


                KpcmynD = KpcmynD+KpcmynDmax/delKZMPCM;%4;
                KicmynD = KicmynD+KicmynDmax/delKZMPCM;
                KdcmynD = KdcmynD+KdcmynDmax/delKZMPCM;%50.50;    
                KpcmxnD = KpcmxnD+KpcmxnDmax/delKZMPCM;%4;
                KicmxnD = KicmxnD+KicmxnDmax/delKZMPCM;
                KdcmxnD = KdcmxnD+KdcmxnDmax/delKZMPCM;%50.50;
            else
                KpzmpcmynD = KpzmpcmynDmax;%4.17;
                KizmpcmynD = KizmpcmynDmax;
                KdzmpcmynD = KdzmpcmynDmax;%20.02;    
                KpzmpcmxnD = KpzmpcmxnDmax;%4.17;
                KizmpcmxnD = KizmpcmxnDmax;
                KdzmpcmxnD = KdzmpcmxnDmax;%20.02;


                KpcmynD = KpcmynDmax;%4;
                KicmynD = KicmynDmax;
                KdcmynD = KdcmynDmax;%50.50;    
                KpcmxnD = KpcmxnDmax;%4;
                KicmxnD = KicmxnDmax;
                KdcmxnD = KdcmxnDmax;%50.50;                
            end       
        else
             if(KpzmpcmynD>0)
                KpzmpcmynD = KpzmpcmynD-KpzmpcmynDmax/delKZMPCM;%4.17;
                KizmpcmynD = KizmpcmynD-KizmpcmynDmax/delKZMPCM;
                KdzmpcmynD = KdzmpcmynD-KdzmpcmynDmax/delKZMPCM;%20.02;    
                KpzmpcmxnD = KpzmpcmxnD-KpzmpcmxnDmax/delKZMPCM;%4.17;
                KizmpcmxnD = KizmpcmxnD-KizmpcmxnDmax/delKZMPCM;
                KdzmpcmxnD = KdzmpcmxnD-KdzmpcmxnDmax/delKZMPCM;%20.02;


                KpcmynD = KpcmynD-KpcmynDmax/delKZMPCM;%4;
                KicmynD = KicmynD-KicmynDmax/delKZMPCM;
                KdcmynD = KdcmynD-KdcmynDmax/delKZMPCM;%50.50;    
                KpcmxnD = KpcmxnD-KpcmxnDmax/delKZMPCM;%4;
                KicmxnD = KicmxnD-KicmxnDmax/delKZMPCM;
                KdcmxnD = KdcmxnD-KdcmxnDmax/delKZMPCM;%50.50;
             else
                KpzmpcmynD = 0;
                KizmpcmynD = 0;
                KdzmpcmynD = 0;    
                KpzmpcmxnD = 0;
                KizmpcmxnD = 0;
                KdzmpcmxnD = 0;


                KpcmynD = 0;
                KicmynD = 0;
                KdcmynD = 0;    
                KpcmxnD = 0;
                KicmxnD = 0;
                KdcmxnD = 0;                 
            end   
        end               
        
        
%         %%%LEVA NOGA
% 
%                 ZMPCMerr(:,i) = rotZ(qsave(6))*ZMPCMerr(:,i);
%                 CMerr(:,i) = rotZ(qsave(6))*CMerr(:,i);
%                 
%                 qvzelj(21,i) = qvzelj(21,i) +(KpzmpcmxnL*ZMPCMerr(1,i)+KizmpcmxnL*sum(ZMPCMerr(1,:),2)+KdzmpcmxnL*(ZMPCMerr(1,i)-ZMPCMerr(1,i-1)));
%                 qvzelj(17,i) = qvzelj(17,i) -(KpzmpcmxnL*ZMPCMerr(1,i)+KizmpcmxnL*sum(ZMPCMerr(1,:),2)+KdzmpcmxnL*(ZMPCMerr(1,i)-ZMPCMerr(1,i-1)));
% 
%                 qvzelj(19,i) = qvzelj(19,i) -(KpzmpcmynL*ZMPCMerr(2,i)+KizmpcmynL*sum(ZMPCMerr(2,:),2)+KdzmpcmynL*(ZMPCMerr(2,i)-ZMPCMerr(2,i-1)));
%                 qvzelj(16,i) = qvzelj(16,i) +(KpzmpcmynL*ZMPCMerr(2,i)+KizmpcmynL*sum(ZMPCMerr(2,:),2)+KdzmpcmynL*(ZMPCMerr(2,i)-ZMPCMerr(2,i-1)));
% 
%                 %CM
%                 qvzelj(21,i) = qvzelj(21,i) -(KpcmxnL*CMerr(1,i)+KicmxnL*sum(CMerr(1,:),2)+KdcmxnL*(CMerr(1,i)-CMerr(1,i-1)));
%                 qvzelj(17,i) = qvzelj(17,i) +(KpcmxnL*CMerr(1,i)+KicmxnL*sum(CMerr(1,:),2)+KdcmxnL*(CMerr(1,i)-CMerr(1,i-1)));
% 
%                 qvzelj(19,i) = qvzelj(19,i) +(KpcmynL*CMerr(2,i)+KicmynL*sum(CMerr(2,:),2)+KdcmynL*(CMerr(2,i)-CMerr(2,i-1)));
%                 qvzelj(16,i) = qvzelj(16,i) -(KpcmynL*CMerr(2,i)+KicmynL*sum(CMerr(2,:),2)+KdcmynL*(CMerr(2,i)-CMerr(2,i-1)));
% 
% 
%         %%%%DESNA NOGA                
%         
%                 qvzelj(13,i) = qvzelj(13,i) +(KpzmpcmxnD*ZMPCMerr(1,i)+KizmpcmxnD*sum(ZMPCMerr(1,:),2)+KdzmpcmxnD*(ZMPCMerr(1,i)-ZMPCMerr(1,i-1)));
%                 qvzelj( 9,i) = qvzelj( 9,i) -(KpzmpcmxnD*ZMPCMerr(1,i)+KizmpcmxnD*sum(ZMPCMerr(1,:),2)+KdzmpcmxnD*(ZMPCMerr(1,i)-ZMPCMerr(1,i-1)));
% 
%                 qvzelj(11,i) = qvzelj(11,i) -(KpzmpcmynD*ZMPCMerr(2,i)+KizmpcmynD*sum(ZMPCMerr(2,:),2)+KdzmpcmynD*(ZMPCMerr(2,i)-ZMPCMerr(2,i-1)));
%                 qvzelj( 8,i) = qvzelj( 8,i) +(KpzmpcmynD*ZMPCMerr(2,i)+KizmpcmynD*sum(ZMPCMerr(2,:),2)+KdzmpcmynD*(ZMPCMerr(2,i)-ZMPCMerr(2,i-1)));
% 
%                 %CM
%                 qvzelj(13,i) = qvzelj(13,i) -(KpcmxnD*CMerr(1,i)+KicmxnD*sum(CMerr(1,:),2)+KdcmxnD*(CMerr(1,i)-CMerr(1,i-1)));
%                 qvzelj( 9,i) = qvzelj( 9,i) +(KpcmxnD*CMerr(1,i)+KicmxnD*sum(CMerr(1,:),2)+KdcmxnD*(CMerr(1,i)-CMerr(1,i-1)));
% 
%                 qvzelj(11,i) = qvzelj(11,i) +(KpcmynD*CMerr(2,i)+KicmynD*sum(CMerr(2,:),2)+KdcmynD*(CMerr(2,i)-CMerr(2,i-1)));
%                 qvzelj( 8,i) = qvzelj( 8,i) -(KpcmynD*CMerr(2,i)+KicmynD*sum(CMerr(2,:),2)+KdcmynD*(CMerr(2,i)-CMerr(2,i-1)));

                        %%%LEVA NOGA
                            
        e_ort=flier.ge;
        A_mat=flier.A;
                for ik = 1 : N
                    tteemmpp =  A_mat(:,:,ik)'*((dot(e_ort(:,ik),(rotz(-pi/2))*ZMPCMerr(:,i)).*e_ort(:,ik)));
                    ZMPCMerr_r(ik,i) = sum(tteemmpp);
                    
                    tteemmpp =  A_mat(:,:,ik)'*((dot(e_ort(:,ik),(rotz(-pi/2))*CMerr(:,i)).*e_ort(:,ik)));
                    CMerr_r(ik,i) = sum(tteemmpp);

                end
                jn=21;
                qvzelj(21,i) = qvzelj(21,i) -(KpzmpcmxnL*ZMPCMerr_r(jn,i)+KizmpcmxnL*sum(ZMPCMerr_r(jn,:),2)+KdzmpcmxnL*(ZMPCMerr_r(jn,i)-ZMPCMerr_r(jn,i-1)));
                jn=17;
                qvzelj(17,i) = qvzelj(17,i) +(KpzmpcmxnL*ZMPCMerr_r(jn,i)+KizmpcmxnL*sum(ZMPCMerr_r(jn,:),2)+KdzmpcmxnL*(ZMPCMerr_r(jn,i)-ZMPCMerr_r(jn,i-1)));

                jn=19;
                qvzelj(19,i) = qvzelj(19,i) -(KpzmpcmynL*ZMPCMerr_r(jn,i)+KizmpcmynL*sum(ZMPCMerr_r(jn,:),2)+KdzmpcmynL*(ZMPCMerr_r(jn,i)-ZMPCMerr_r(jn,i-1)));
                jn=16;
                qvzelj(16,i) = qvzelj(16,i) +(KpzmpcmynL*ZMPCMerr_r(jn,i)+KizmpcmynL*sum(ZMPCMerr_r(jn,:),2)+KdzmpcmynL*(ZMPCMerr_r(jn,i)-ZMPCMerr_r(jn,i-1)));

                %CM
                jn=21;
                qvzelj(21,i) = qvzelj(21,i) +(KpcmxnL*CMerr_r(jn,i)+KicmxnL*sum(CMerr_r(jn,:),2)+KdcmxnL*(CMerr_r(jn,i)-CMerr_r(jn,i-1)));
                jn=17;
                qvzelj(17,i) = qvzelj(17,i) -(KpcmxnL*CMerr_r(jn,i)+KicmxnL*sum(CMerr_r(jn,:),2)+KdcmxnL*(CMerr_r(jn,i)-CMerr_r(jn,i-1)));

                jn=19;
                qvzelj(19,i) = qvzelj(19,i) +(KpcmynL*CMerr_r(jn,i)+KicmynL*sum(CMerr_r(jn,:),2)+KdcmynL*(CMerr_r(jn,i)-CMerr_r(jn,i-1)));
                jn=16;
                qvzelj(16,i) = qvzelj(16,i) -(KpcmynL*CMerr_r(jn,i)+KicmynL*sum(CMerr_r(jn,:),2)+KdcmynL*(CMerr_r(jn,i)-CMerr_r(jn,i-1)));


        %%%%DESNA NOGA                
        
                jn=13;
                qvzelj(13,i) = qvzelj(13,i) -(KpzmpcmxnD*ZMPCMerr_r(jn,i)+KizmpcmxnD*sum(ZMPCMerr_r(jn,:),2)+KdzmpcmxnD*(ZMPCMerr_r(jn,i)-ZMPCMerr_r(jn,i-1)));
                jn=9;
                qvzelj( 9,i) = qvzelj( 9,i) +(KpzmpcmxnD*ZMPCMerr_r(jn,i)+KizmpcmxnD*sum(ZMPCMerr_r(jn,:),2)+KdzmpcmxnD*(ZMPCMerr_r(jn,i)-ZMPCMerr_r(jn,i-1)));

                jn=11;
                qvzelj(11,i) = qvzelj(11,i) -(KpzmpcmynD*ZMPCMerr_r(jn,i)+KizmpcmynD*sum(ZMPCMerr_r(jn,:),2)+KdzmpcmynD*(ZMPCMerr_r(jn,i)-ZMPCMerr_r(jn,i-1)));
                jn=8;
                qvzelj( 8,i) = qvzelj( 8,i) +(KpzmpcmynD*ZMPCMerr_r(jn,i)+KizmpcmynD*sum(ZMPCMerr_r(jn,:),2)+KdzmpcmynD*(ZMPCMerr_r(jn,i)-ZMPCMerr_r(jn,i-1)));

                %CM
                jn=13;
                qvzelj(13,i) = qvzelj(13,i) +(KpcmxnD*CMerr_r(jn,i)+KicmxnD*sum(CMerr_r(jn,:),2)+KdcmxnD*(CMerr_r(jn,i)-CMerr_r(jn,i-1)));
                jn=9;
                qvzelj( 9,i) = qvzelj( 9,i) -(KpcmxnD*CMerr_r(jn,i)+KicmxnD*sum(CMerr_r(jn,:),2)+KdcmxnD*(CMerr_r(jn,i)-CMerr_r(jn,i-1)));

                jn=11;
                qvzelj(11,i) = qvzelj(11,i) +(KpcmynD*CMerr_r(jn,i)+KicmynD*sum(CMerr_r(jn,:),2)+KdcmynD*(CMerr_r(jn,i)-CMerr_r(jn,i-1)));
                jn=8;
                qvzelj( 8,i) = qvzelj( 8,i) -(KpcmynD*CMerr_r(jn,i)+KicmynD*sum(CMerr_r(jn,:),2)+KdcmynD*(CMerr_r(jn,i)-CMerr_r(jn,i-1)));
                
             
%%Ogranicenje za trup
        if(rad2deg(abs(qnew(segtrupx(1),i)))<5)
            qvzelj(segtrupx,i)=qvzelj(segtrupx,i)+(Kpzmpcmyt*ZMPCMerr_r(segtrupx,i)+Kizmpcmyt*sum(ZMPCMerr_r(segtrupx,:),2)+Kdzmpcmyt*(ZMPCMerr_r(segtrupx,i)-ZMPCMerr_r(segtrupx,i-1)));
            qvzelj(segtrupx,i)=qvzelj(segtrupx,i)-(Kpcmyt*CMerr_r(segtrupx,i)+Kicmyt*sum(CMerr_r(segtrupx,:),2)+Kdcmyt*(CMerr_r(segtrupx,i)-CMerr_r(segtrupx,i-1)));        
        end
        if(rad2deg(abs(qnew(segtrupy(1),i)))<5)
            qvzelj(segtrupy,i)=qvzelj(segtrupy,i)+(Kpzmpcmxt*ZMPCMerr_r(segtrupy,i)+Kizmpcmxt*sum(ZMPCMerr_r(segtrupy,:),2)+Kdzmpcmxt*(ZMPCMerr_r(segtrupy,i)-ZMPCMerr_r(segtrupy,i-1)));   
            qvzelj(segtrupy,i)=qvzelj(segtrupy,i)-(Kpcmxt*CMerr_r(segtrupy,i)+Kicmxt*sum(CMerr_r(segtrupy,:),2)+Kdcmxt*(CMerr_r(segtrupy,i)-CMerr_r(segtrupy,i-1)));        
        end
             
% % % %         %%%LEVA NOGA
% % % %         or = flier.A;
% % % %         u21 = TrotXYZ(or(:,:,21));
% % % %         u17 = TrotXYZ(or(:,:,17));
% % % %         u19 = TrotXYZ(or(:,:,19));
% % % %         u16 = TrotXYZ(or(:,:,16));
% % % % 
% % % % %                qvzelj(21,i) = qvzelj(21,i) +(KpzmpcmxnL*ZMPCMerr(1,i)+KizmpcmxnL*sum(ZMPCMerr(1,i),2)+KdzmpcmxnL*(ZMPCMerr(1,i)-ZMPCMerr(1,i-1)));
% % % %                 
% % % %                 temp = or(:,:,21)*(KpzmpcmxnL*[ZMPCMerr(1:2,i);0]+KizmpcmxnL*[sum(ZMPCMerr(1:2,:),2);0]+KdzmpcmxnL*([ZMPCMerr(1:2,i);0]-[ZMPCMerr(1:2,i-1);0]));
% % % %                 qvzelj(21,i) = qvzelj(21,i)+ temp(1);%(KpzmpcmxnL*ZMPCMerr(1,i)+KizmpcmxnL*sum(ZMPCMerr(1,i),2)+KdzmpcmxnL*(ZMPCMerr(1,i)-ZMPCMerr(1,i-1)));
% % % %                 temp = or(:,:,17)*(KpzmpcmxnL*[ZMPCMerr(1:2,i);0]+KizmpcmxnL*[sum(ZMPCMerr(1:2,:),2);0]+KdzmpcmxnL*([ZMPCMerr(1:2,i);0]-[ZMPCMerr(1:2,i-1);0]));
% % % %                 qvzelj(17,i) = qvzelj(17,i) -temp(1);%(KpzmpcmxnL*ZMPCMerr(1,i)+KizmpcmxnL*sum(ZMPCMerr(1,i),2)+KdzmpcmxnL*(ZMPCMerr(1,i)-ZMPCMerr(1,i-1)));
% % % % 
% % % %                 temp = or(:,:,19)*(KpzmpcmynL*[ZMPCMerr(1:2,i);0]+KizmpcmynL*[sum(ZMPCMerr(1:2,:),2);0]+KdzmpcmynL*([ZMPCMerr(1:2,i);0]-[ZMPCMerr(1:2,i-1);0]));
% % % %                 qvzelj(19,i) = qvzelj(19,i) -temp(2);%(KpzmpcmynL*ZMPCMerr(2,i)+KizmpcmynL*sum(ZMPCMerr(2,i),2)+KdzmpcmynL*(ZMPCMerr(2,i)-ZMPCMerr(2,i-1)));
% % % %                 temp = or(:,:,16)*(KpzmpcmynL*[ZMPCMerr(1:2,i);0]+KizmpcmynL*[sum(ZMPCMerr(1:2,:),2);0]+KdzmpcmynL*([ZMPCMerr(1:2,i);0]-[ZMPCMerr(1:2,i-1);0]));
% % % %                 qvzelj(16,i) = qvzelj(16,i) +temp(2);%(KpzmpcmynL*ZMPCMerr(2,i)+KizmpcmynL*sum(ZMPCMerr(2,i),2)+KdzmpcmynL*(ZMPCMerr(2,i)-ZMPCMerr(2,i-1)));
% % % % 
% % % %                 %CM
% % % %                 temp = or(:,:,19)*(KpcmxnL*[CMerr(1:2,i);0]+KicmxnL*[sum(CMerr(1:2,:),2);0]+KdcmxnL*([CMerr(1:2,i);0]-[CMerr(1:2,i-1);0]));
% % % %                 qvzelj(21,i) = qvzelj(21,i) -(KpcmxnL*CMerr(1,i)+KicmxnL*sum(CMerr(1,:),2)+KdcmxnL*(CMerr(1,i)-CMerr(1,i-1)));
% % % %                 qvzelj(17,i) = qvzelj(17,i) +(KpcmxnL*CMerr(1,i)+KicmxnL*sum(CMerr(1,:),2)+KdcmxnL*(CMerr(1,i)-CMerr(1,i-1)));
% % % % 
% % % %                 qvzelj(19,i) = qvzelj(19,i) +(KpcmynL*CMerr(2,i)+KicmynL*sum(CMerr(2,:),2)+KdcmynL*(CMerr(2,i)-CMerr(2,i-1)));
% % % %                 qvzelj(16,i) = qvzelj(16,i) -(KpcmynL*CMerr(2,i)+KicmynL*sum(CMerr(2,i),2)+KdcmynL*(CMerr(2,i)-CMerr(2,i-1)));
% % % % 
% % % % 
% % % %         %%%%DESNA NOGA                
% % % %         
% % % %                 qvzelj(13,i) = qvzelj(13,i) +(KpzmpcmxnD*ZMPCMerr(1,i)+KizmpcmxnD*sum(ZMPCMerr(1,i),2)+KdzmpcmxnD*(ZMPCMerr(1,i)-ZMPCMerr(1,i-1)));
% % % %                 qvzelj( 9,i) = qvzelj( 9,i) -(KpzmpcmxnD*ZMPCMerr(1,i)+KizmpcmxnD*sum(ZMPCMerr(1,i),2)+KdzmpcmxnD*(ZMPCMerr(1,i)-ZMPCMerr(1,i-1)));
% % % % 
% % % %                 qvzelj(11,i) = qvzelj(11,i) -(KpzmpcmynD*ZMPCMerr(2,i)+KizmpcmynD*sum(ZMPCMerr(2,i),2)+KdzmpcmynD*(ZMPCMerr(2,i)-ZMPCMerr(2,i-1)));
% % % %                 qvzelj( 8,i) = qvzelj( 8,i) +(KpzmpcmynD*ZMPCMerr(2,i)+KizmpcmynD*sum(ZMPCMerr(2,i),2)+KdzmpcmynD*(ZMPCMerr(2,i)-ZMPCMerr(2,i-1)));
% % % % 
% % % %                 %CM
% % % %                 qvzelj(13,i) = qvzelj(13,i) -(KpcmxnD*CMerr(1,i)+KicmxnD*sum(CMerr(1,i),2)+KdcmxnD*(CMerr(1,i)-CMerr(1,i-1)));
% % % %                 qvzelj( 9,i) = qvzelj( 9,i) +(KpcmxnD*CMerr(1,i)+KicmxnD*sum(CMerr(1,i),2)+KdcmxnD*(CMerr(1,i)-CMerr(1,i-1)));
% % % % 
% % % %                 qvzelj(11,i) = qvzelj(11,i) +(KpcmynD*CMerr(2,i)+KicmynD*sum(CMerr(2,i),2)+KdcmynD*(CMerr(2,i)-CMerr(2,i-1)));
% % % %                 qvzelj( 8,i) = qvzelj( 8,i) -(KpcmynD*CMerr(2,i)+KicmynD*sum(CMerr(2,i),2)+KdcmynD*(CMerr(2,i)-CMerr(2,i-1)));
% % % %                 
% % % %                 
% % % %         
%         Kpzmpcmyt = 0;
%         Kizmpcmyt = 0;
%         Kdzmpcmyt = 0;    
%         Kpzmpcmxt = 0;
%         Kizmpcmxt = 0;
%         Kdzmpcmxt = 0;
% 
%         Kpcmyt = 0;
%         Kicmyt = 0;
%         Kdcmyt = 0;    
%         Kpcmxt = 0;
%         Kicmxt = 0;
%         Kdcmxt = 0;
        
%         if(abs(CMerr(1,i))>0.02)
%             Kpcmxt=Kpcmxt+0.01;
%         else
%             Kpcmxt=Kpcmxt-0.01;
%         end
%         
%         if(Kpcmxt>0.3)
%             Kpcmxt=0.3;
%         end
%         if(Kpcmxt<0.01)
%             Kpcmxt=0.01;
%         end
%         
%         
%         if(abs(CMerr(2,i))>0.02)
%             Kpcmyt=Kpcmyt+0.01;
%         else
%             Kpcmyt=Kpcmyt-0.01;
%         end
%         
%         if(Kpcmyt>0.3)
%             Kpcmyt=0.3;
%         end
%         if(Kpcmyt<0.01)
%             Kpcmyt=0.01;
%         end
        
% % %         Kpzmpcmyt = 3;
% % %         Kizmpcmyt = 0;%0.002;
% % %         Kdzmpcmyt = 0;%.22;    
% % %         Kpzmpcmxt = 3;
% % %         Kizmpcmxt = 0;%0.002;
% % %         Kdzmpcmxt = 0;%0.22;
% % % 
% % % 
% % %         Kpcmyt = 0.8;
% % %         Kicmyt = 0;%0.002;
% % %         Kdcmyt = 0;%1.12;    
% % %         Kpcmxt = 0.8;
% % %         Kicmxt = 0;%0.002;
% % %         Kdcmxt = 0;%1.12;


%%Ogranicenje za trup
%         if(rad2deg(abs(qnew(segtrupx(1),i)))<5)
%             qvzelj(segtrupx,i)=qvzelj(segtrupx,i)-(+1)*(Kpcmyt*CMerr(2,i)+Kicmyt*sum(CMerr(2,:),2)+Kdcmyt*(CMerr(2,i)-CMerr(2,i-1)));        
%             qvzelj(segtrupx,i)=qvzelj(segtrupx,i)+(Kpzmpcmyt*ZMPCMerr(2,i)+Kizmpcmyt*sum(ZMPCMerr(2,:),2)+Kdzmpcmyt*(ZMPCMerr(2,i)-ZMPCMerr(2,i-1)));        
%         end
%         if(rad2deg(abs(qnew(segtrupy(1),i)))<5)
%             qvzelj(segtrupy,i)=qvzelj(segtrupy,i)+(+1)*(Kpcmxt*CMerr(1,i)+Kicmxt*sum(CMerr(1,:),2)+Kdcmxt*(CMerr(1,i)-CMerr(1,i-1)));        
%             qvzelj(segtrupy,i)=qvzelj(segtrupy,i)-(Kpzmpcmxt*ZMPCMerr(1,i)+Kizmpcmxt*sum(ZMPCMerr(1,:),2)+Kdzmpcmxt*(ZMPCMerr(1,i)-ZMPCMerr(1,i-1)));        
%         end

%         if(rad2deg(abs(qnew(segtrupx(1),i)))<5)
%             qvzelj(segtrupx,i)=qvzelj(segtrupx,i)-(-1)*(Kpcmyt*CMerr(2,i)+Kicmyt*sum(CMerr(2,:),2)+Kdcmyt*(CMerr(2,i)-CMerr(2,i-1)));        
%             qvzelj(segtrupx,i)=qvzelj(segtrupx,i)+(Kpzmpcmyt*ZMPCMerr(2,i)+Kizmpcmyt*sum(ZMPCMerr(2,:),2)+Kdzmpcmyt*(ZMPCMerr(2,i)-ZMPCMerr(2,i-1)));        
%         end
%         if(rad2deg(abs(qnew(segtrupy(1),i)))<5)
%             qvzelj(segtrupy,i)=qvzelj(segtrupy,i)+(-1)*(Kpcmxt*CMerr(1,i)+Kicmxt*sum(CMerr(1,:),2)+Kdcmxt*(CMerr(1,i)-CMerr(1,i-1)));        
%             qvzelj(segtrupy,i)=qvzelj(segtrupy,i)-(Kpzmpcmxt*ZMPCMerr(1,i)+Kizmpcmxt*sum(ZMPCMerr(1,:),2)+Kdzmpcmxt*(ZMPCMerr(1,i)-ZMPCMerr(1,i-1)));        
%         end        
        
%         qvzelj(seglruka,i)=qvzelj(seglruka,i)-1.8*qnew(seglruka,i);
%         qvzelj(segdruka,i)=qvzelj(segdruka,i)-1.8*qnew(segdruka,i);
        qvzelj(segtrupx,i)=qvzelj(segtrupx,i)-1.8*qnew(segtrupx,i);
        qvzelj(segtrupy,i)=qvzelj(segtrupy,i)-1.8*qnew(segtrupy,i);
        qvzelj(segtrupx2,i)=qvzelj(segtrupx2,i)-1.8*qnew(segtrupx2,i);
        qvzelj(segtrupy2,i)=qvzelj(segtrupy2,i)-1.8*qnew(segtrupy2,i);
        
        qvzelj(seglnoga(1),i)=qvzelj(seglnoga(1),i)-1.8*qnew(seglnoga(1),i);
        qvzelj(segdnoga(1),i)=qvzelj(segdnoga(1),i)-1.8*qnew(segdnoga(1),i);
        
%Ogranicenje za koleno        
        if(rad2deg(qnew(seglnoga(4),i))<5)
            qvzelj(seglnoga(4),i) = 4*(deg2rad(5)-qnew(seglnoga(4),i));
        end

        if(rad2deg(qnew(segdnoga(4),i))<5)
            qvzelj(segdnoga(4),i)=  4*(deg2rad(5)-qnew(segdnoga(4),i)); 
        end
        

        %%%Filter
        ett = exp(-dt/0.0005);

        qvzelj(:,i)=qvzelj(:,i)*(1-ett)+ett*qvzelj(:,i-1);
        
%         qverr(:,i) = qvzelj(:,i)-qvnew(:,i);
        qverr(:,i) = qvnew(:,i)-qvzelj(:,i);
%         qvcerr(:,i) = qverr(:,i)-qverr(:,i-1);

        qzelj(:,i) = qnew(:,i-1)+dt*qvzelj(:,i);
        qazelj(:,i) = (qvzelj(:,i)-qvzelj(:,i-1))/dt;

%         qerr(:,i) = qzelj(:,i)-qnew(:,i);
        qerr(:,i) = qnew(:,i)-qzelj(:,i);
        
       %%% NELINEARNO UPRAVLJANJE
        if(NL_Control)

            %Disturbance estimator
        Xest(1:N,i)=qnewest(:,i-1)+(qvnewest(:,i-1)+l1*(qnew(:,i-1)-qnewest(:,i-1)))*dt;%integralim da dobijem novo estimirano stanje
        Xest(N+1:2*N,i)=qvnewest(:,i-1)+(dest(:,i-1)+v(:,i-1)+l2*(qnew(:,i-1)-qnewest(:,i-1)))*dt;
        Xest(2*N+1:3*N,i)=dest(:,i-1)+destd(:,i-1)*dt;
        Xpest(:,i)=fliergeneral_est(t(i),Xest(:,i));
        qnewest(:,i)=Xest(1:N,i);
        qvnewest(:,i)=Xest(N+1:2*N,i);
        qanewest(:,i)=Xpest(N+1:2*N,i);
        dest(:,i)=Xest(2*N+1:3*N,i);
        destd(:,i)=Xpest(2*N+1:3*N,i);

            
        eta=5;  
        KI=150;
        fi=5;%*[0.5 0.5 0.5 0.5 0.5 0.5 1.8 1.8 2 1.8 5 1.8 1.8 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5]';
        sigmaI=0;            
            
            flier_m=kk_geo(flier_m,qnew(:,i));
            flier_m=kk_kine(flier_m,qvnew(:,i));
            flier_m= kk_dyn(flier_m);
            [H_m h0_m]=kk_inemat(flier_m);
            
            taud_m=zeros(N,1);            
            for ii=1:12                
                [Jk Ak]=km_jakP(flier_m,ii);     %#ok<NASGU>
                taud_m=taud_m+Jk(1:3,:)'*FF(:,ii);
            end;


        km_geo(flier_m,qzelj(:,i-1));
        km_kine(flier_m,qvzelj(:,i-1));
        km_dyn(flier_m);
        [H h0] = km_inemat(flier_m);
        
        km_geo(flier,qzelj(:,i-1));
        km_kine(flier,qvzelj(:,i-1));
        km_dyn(flier);
        [Hmax h0max] = km_inemat(flier);
        
        
        Lambda=20;
        
        Kp_FL=4800;
        Kd_FL=3;
            
            
        s(:,i)=qverr(:,i)+Lambda*qerr(:,i);%jednacina S povrsi
        v1(:,i)=v1(:,i-1)+sign(s(:,i))*dt;
            
%         k(actseg)=abs(dest(actseg,i))+eta(actseg)';
        
        beta=abs(inv(H)*Hmax);%ovo radi
        deltaf=abs(inv(H)*h0-inv(Hmax)*h0max);
        Ks=(abs(beta-1)*abs(qazelj(:,i)-Kd_FL*qverr(:,i)-Kp_FL*qerr(:,i))+beta*(eta*ones(N,1)+deltaf));%+dest(actseg,i)


%         if flag==true%iskljucujem sliding mode u saturaciji
%             k=zeros(N,1);
%         else
%             k=100*[0 0 0 0 0 0 0.01 0.005 0.004 0.01 0.006 0.01 0.0075 0.01 0.004 0.004 0.01 0.006 0.01 0.0075 0 0]';
%         end
        
%         K=eye(N).*[k k k k k k k k k k k k k k k k k k k k k k k k k k k k k k];


    %         if flag==true%iskljucujem sliding mode u saturaciji
    %             k=zeros(N,1);
    %         else
    %             k=100*[0 0 0 0 0 0 0.01 0.005 0.004 0.01 0.006 0.01 0.0075 0.01 0.004 0.004 0.01 0.006 0.01 0.0075 0 0]';
    %         end

%             K=eye(N).*[k k k k k k k k k k k k k k k k k k k k k k k k k k k k k k];

    %         fi(:,1)=s(:,i)+0.1;
            sigma(:,i)=s(:,i)./fi(:,1);

            if sigma>1
                ro(:,i)=1;
                sigmaI=0;
            elseif sigma<-1
                ro(:,i)=-1;
                sigmaI=0;
            else
                sigmaI=sigmaI+sigma(:,i)*dt;
                ro(:,i)=sigma(:,i)+KI*sigmaI;
            end            
            
%                 tau_z(:,i) = H_m*qazelj(:,i)+h0_m-taud_m;
%                 Irzelj(:,i) = SKR*((1/Cm*eye(N))*(tau_z(:,i)-Bc*eye(N)*qvzelj(:,i)-Jr*eye(N)*qazelj(:,i)));
%                 Irdzelj(:,i) = (Irzelj(:,i)-Irnew(:,i-1))/dt;
%             
%             %%%%% SRDJANOVE JEDNACINE
%              
%             Fx = inv((Jr*Rr/Cm)*eye(N)+Rr*H)*(-(Bc*Rr/Cm+Ce)*eye(N)*qvnew(:,i)-(h0-taud_m)*Rr-Lr*PRO*Irpnew(:,i));
% %             Gx = inv((Jr*Rr/Cm)*eye(N)+Rr*H);
%             invGx = ((Jr*Rr/Cm)*eye(N)+Rr*H);
%             
%             
%             Kdfc=0;
% %             Kpfc=250;
%             Kpfc=0;
%             
%             v=qazelj(:,i)-Kdfc*(qvnew(:,i)-qvzelj(:,i))-Kpfc*(qnew(:,i)-qzelj(:,i));
%             
%             Upr(:,i)=invGx*(-Fx+v);
%             Up(:,i+1)=SKR*Upr(:,i);%saturacija(Upr(:,i),200);
%             
%             %%%%%   DOVDE
            

            %%%%%%%    MIRKOVE JEDNACIE
            AA=inv(Jr*eye(N)-H_m);
            Fx = [AA*(-h0_m+taud_m-(Cm*Ce/Rr)*qvzelj(:,i))];
            Gx = (Cm/Rr)*AA;
            v(:,i)=qazelj(:,i)-Kd_FL*(qverr(:,i))-Kp_FL*(qerr(:,i-1))-Ks.*ro(:,i)-dest(:,i);
            
            Upr(:,i) = inv(Gx)*(v(:,i)-Fx);
            Up(:,i+1)=SKR*saturacija(Upr(:,i),200);
            
%             Fx = [qvnew(:,i); inv(Jr*eye(N)-H)*(h0-taudod-Cm*PRO*Irnew(:,i)-Bc*qvnew(:,i));(1/Lr)*(Rr*PRO*Irnew(:,i)+Ce*qvnew(:,i))];
%             Gx = [zeros(N,N);zeros(N,N);1/Lr*eye(N);];
            
%             x = [qnew(:,i);qvnew(:,i);PRO*Irnew(:,i)];
%             xd = [qvnew(:,i);qanew(:,i);PRO*Irpnew(:,i)];
%             xz = [qzelj(:,i);qvzelj(:,i);PRO*Irzelj(:,i)];
%             xzd = [qvzelj(:,i);qazelj(:,i);PRO*Irdzelj(:,i)];
%             
%             xe = x-xz;
%             xed = xd-xzd;
%             xe_all(:,i)=xe;
%             xei = sum(xe_all,2);
%             Kpfc=250;
%             Kpfc=200;
%             Kpfc=[0*eye(N) zeros(N,2*N); zeros(N) 1000*eye(N) zeros(N); zeros(N,2*N) 0*eye(N)];
%             Kpfc=[zeros(N) zeros(N) zeros(N);zeros(N) zeros(N) zeros(N); zeros(N) zeros(N) 60*eye(N)];
% 
%             Kifc=[zeros(N) zeros(N) zeros(N); zeros(N) zeros(N) zeros(N); zeros(N) zeros(N) 0.2*eye(N);];

%             Kdfc=-0.05;
%             Kdfc=[0*eye(N) zeros(N,2*N); zeros(N) 40*eye(N) zeros(N); zeros(N,2*N) 0*eye(N)];
%             Kdfc=[zeros(N) zeros(N) zeros(N); zeros(N) zeros(N) zeros(N); zeros(N) zeros(N) 0.02*eye(N);];
            
            
%             v=[1*eye(N) zeros(N,2*N); zeros(N) 1*eye(N) zeros(N); zeros(N,2*N) 0*eye(N)]*xzd-Kpfc*xe-Kifc*xei-Kdfc*xed;%-Ksat(s)-d_est;
%             Upr(:,i)=pinv(Gx)*(-Fx-v);
% %             Upr(:,i)=Upr(:,i)+dUpr(:,i);
%             Up(:,i+1)=SKR*saturacija(Upr(:,i),200);
            
        else        
        
             %%%%ADAPTIVNI KOEFICIJENTI  !!!!!!


            %%%RUKA osa
            Kpruka=5.2;
            Kiruka=0.0001;
            Kdruka=64.0;

            %%%TRUP X osa
            Kptx=5.2;
            Kitx=0.0001;
            Kdtx=64.0;

            %%%TRUP Y osa
            Kpty=5.2;
            Kity=0.0001;
            Kdty=64.0;

            %%%KUK LEVA
            Kplnku=10.12;
            Kilnku=0.002;
            Kdlnku=230.0;               

            %%%KUK DESNA
            Kpdnku=10.12;
            Kidnku=0.002;
            Kddnku=230.0;

            %%%KOLENO LEVA
            Kplnko=12.12;
            Kilnko=0.003;
            Kdlnko=280.0;

            %%%KOLENO DESNA
            Kpdnko=12.12;
            Kidnko=0.003;
            Kddnko=280.0;

            %%%SKZ LEVA
            Kplnsz=15.12;
            Kilnsz=0.004;
            Kdlnsz=230.0;

            %%%SKZ DESNA
            Kpdnsz=15.12;
            Kidnsz=0.004;
            Kddnsz=230.0;

            %%%PRSTI LEVA
            Kpprstl=2.0512;
            Kiprstl=0.0003;
            Kdprstl=80.0400;

            %%%PRSTI DESNA
            Kpprstd=2.0512;
            Kiprstd=0.0003;
            Kdprstd=80.0400;        


            dUpr(seglruka,i) = Kpruka*qverr(seglruka,i)+Kiruka*sum(qverr(seglruka,:),2)+Kdruka*(qvcerr(seglruka,i));
            dUpr(segdruka,i) = Kpruka*qverr(segdruka,i)+Kiruka*sum(qverr(segdruka,:),2)+Kdruka*(qvcerr(segdruka,i));

            dUpr(segtrupx,i) = Kptx*qverr(segtrupx,i)+Kitx*sum(qverr(segtrupx,:),2)+Kdtx*(qvcerr(segtrupx,i));
            dUpr(segtrupy,i) = Kpty*qverr(segtrupy,i)+Kity*sum(qverr(segtrupy,:),2)+Kdty*(qvcerr(segtrupy,i));
            dUpr(segtrupx2,i) = Kptx*qverr(segtrupx2,i)+Kitx*sum(qverr(segtrupx2,:),2)+Kdtx*(qvcerr(segtrupx2,i));
            dUpr(segtrupy2,i) = Kpty*qverr(segtrupy2,i)+Kity*sum(qverr(segtrupy2,:),2)+Kdty*(qvcerr(segtrupy2,i));
            dUpr(seglnoga(1:3),i) = Kplnku*qverr(seglnoga(1:3),i)+Kilnku*sum(qverr(seglnoga(1:3),:),2)+Kdlnku*(qvcerr(seglnoga(1:3),i));
            dUpr(segdnoga(1:3),i) = Kpdnku*qverr(segdnoga(1:3),i)+Kidnku*sum(qverr(segdnoga(1:3),:),2)+Kddnku*(qvcerr(segdnoga(1:3),i));
            dUpr(seglnoga(4),i) = Kplnko*qverr(seglnoga(4),i)+Kilnko*sum(qverr(seglnoga(4),:),2)+Kdlnko*(qvcerr(seglnoga(4),i));
            dUpr(segdnoga(4),i) = Kpdnko*qverr(segdnoga(4),i)+Kidnko*sum(qverr(segdnoga(4),:),2)+Kddnko*(qvcerr(segdnoga(4),i));
            dUpr(seglnoga(5:7),i) = Kplnsz*qverr(seglnoga(5:7),i)+Kilnsz*sum(qverr(seglnoga(5:7),:),2)+Kdlnsz*(qvcerr(seglnoga(5:7),i));
            dUpr(segdnoga(5:7),i) = Kpdnsz*qverr(segdnoga(5:7),i)+Kidnsz*sum(qverr(segdnoga(5:7),:),2)+Kddnsz*(qvcerr(segdnoga(5:7),i));
            dUpr(seglprst,i) = Kpprstl*qverr(seglprst,i)+Kiprstl*sum(qverr(seglprst,:),2)+Kdprstl*(qvcerr(seglprst,i));
            dUpr(segdprst,i) = Kpprstd*qverr(segdprst,i)+Kiprstd*sum(qverr(segdprst,:),2)+Kdprstd*(qvcerr(segdprst,i));


    %         Upr(seglprst,i)=-taunew(seglprst-6,i);
    %         Upr(segdprst,i)=-taunew(segdprst-6,i);


            Upr(:,i)=Upr(:,i)+dUpr(:,i);
            Up(:,i+1)=SKR*Upr(:,i);
        end
        %IScrtavanje
        if(mod(i,200)==0)

            figure(1);
            subplot(4,2,[1 3],'replace');
            hold on;
            draw_fliermirko(flier,seg, 42,[13 21],[14 22],[rcm(1,i)-1 rcm(1,i)+1 -0.5 0.5 -0.1 1.9],[sin(zeljeni_pravac) -cos(zeljeni_pravac) 0],[1:12],2,'rrr');
            crtaj_zid();
%             plot3(rcm(1,i),rcm(2,i),rcm(3,i),'ko');    

            plot3(r_saka_d_p(1),r_saka_d_p(2),r_saka_d_p(3),'go');    
            plot3(r_saka_l_p(1),r_saka_l_p(2),r_saka_l_p(3),'bo');    

            plot3(rrukad(1),rrukad(2),rrukad(3),'gx');    
            plot3(rrukal(1),rrukal(2),rrukal(3),'bx');    

            title(['zp: ' num2str(rad2deg(zeljeni_pravac))])
            
            subplot(4,2,[5 7],'replace');
            hold on;
            draw_fliermirko(flier,seg, 42,[13 21],[14 22],[rcm(1,i)-1 rcm(1,i)+1 -0.5 0.5 -0.1 1.9],[0 0 1],[1:12]);
            axis([-0.1 8 -2 2 -0.1 1.9]);
            crtaj_zid();
%             plot3(rcm(1,i),rcm(2,i),rcm(3,i),'ko');    

            plot3(r_saka_d_p(1),r_saka_d_p(2),r_saka_d_p(3),'go');    
            plot3(r_saka_l_p(1),r_saka_l_p(2),r_saka_l_p(3),'bo');    

            plot3(rrukad(1),rrukad(2),rrukad(3),'gx');    
            plot3(rrukal(1),rrukal(2),rrukal(3),'bx');    

            plot(putanja(1,:),putanja(2,:))
            
            title(['Kpcmxt:' num2str(Kpcmxt) ' - Kpcmyt:' num2str(Kpcmyt)]);
            
        
            subplot(4,2,[2 4],'replace');
            hold on;          
            tac=1:nknew;        
            crtaj=1;
            if size(XC,2)>2
                try
                    tac=convhull(XC,YC);
                catch
                    crtaj=0;
                end
                if crtaj   
                [XCm YCm] = polyoffset(XC(tac),YC(tac),0.01);
                plot(XC(tac),YC(tac),'k-');  
                plot(XCm,YCm,'r-');                 
                end
            end;
            if size(XCL,2)>2 && crtaj
                tacl=convhull(XCL,YCL);      
                [XCLm YCLm] = polyoffset(XCL(tacl),YCL(tacl),0.01);
                plot(XCL(tacl),YCL(tacl),'k-');  
                plot(XCLm,YCLm,'g-'); 

            end;
            if size(XCR,2)>2 && crtaj
                tacr=convhull(XCR,YCR);
                [XCRm YCRm] = polyoffset(XCR(tacr),YCR(tacr),0.01);
                plot(XCR(tacr),YCR(tacr),'k-');  
                plot(XCRm,YCRm,'g-'); 
            end;
            if size(XCPL,2)>2 && crtaj
                tacpl=convhull(XCPL,YCPL);      
                [XCPLm YCPLm] = polyoffset(XCPL(tacpl),YCPL(tacpl),0.01);
                plot(XCPL(tacpl),YCPL(tacpl),'k-');  
%                 plot(XCPLm,YCPLm,'g-');                  
            end
            if size(XCPR,2)>2 && crtaj
                tacpr=convhull(XCPR,YCPR);
                [XCPRm YCPRm] = polyoffset(XCPR(tacpr),YCPR(tacpr),0.01);
                plot(XCPR(tacpr),YCPR(tacpr),'k-');  
%                 plot(XCPRm,YCPRm,'g-');     
            end
            plot(rpetal(1),rpetal(2),'g*');      
            plot(rpetad(1),rpetad(2),'b*');      
            plot(rprstl(1),rprstl(2),'g*');      
            plot(rprstd(1),rprstd(2),'b*');      
            plot(ZMPpos(1,i),ZMPpos(2,i),'rx');  
            plot(CMzeljt(1,i),CMzeljt(2,i),'ro');  
            
            
            plot(rcm(1,i),rcm(2,i),'k*');    
            axis equal;
            title(['faza:' num2str(faza)]);         

            subplot(4,2,[6],'replace');
            plot(ZMPCMerr(1,1:i)');       
            title([num2str(statusnanog) ':' num2str(statusnaginjanjenapred) ':' num2str(statusgrci) ':' num2str(statusopruzi)]);         
            subplot(4,2,[8],'replace');        
            plot(ZMPCMerr(2,1:i)');
%             title([num2str(a) ' : ' num2str(b) ': R:' num2str(size(XCR,2)) ': L:' num2str(size(XCL,2))]);
            
                     figure(2);
              plot(qverr(:,max(1,i-300):i)');
              title(['brzina:' num2str(brzinakoraka)]);
%             
%             subplot(4,4,[1],'replace');
%             plot(qverr(7:10,max(1,i-300):i)');       
%             title(['D kuk&kol e']);         
%                        
%             subplot(4,4,[2],'replace');
%             plot(qverr(11:14,max(1,i-300):i)');       
%             title(['D skz&prst e']);         
%             
%             subplot(4,4,[3],'replace');
%             plot(qverr(15:18,max(1,i-300):i)');       
%             title(['L kuk&kol e']);         
%                        
%             subplot(4,4,[4],'replace');
%             plot(qverr(19:22,max(1,i-300):i)');       
%             title(['L skz&prst e']);         
%             
%             subplot(4,4,[5],'replace');
%             plot(qverr(segtrupx,max(1,i-300):i)');       
%             title(['T x osa e']);         
% 
%             subplot(4,4,[6],'replace');
%             plot(qverr(segtrupy,max(1,i-300):i)');       
%             title(['T y osa e']);         
%             
%             
%             subplot(4,4,[7],'replace');
%             plot(qverr([seglruka segdruka],max(1,i-300):i)');       
%             title(['Reka e']);         
%             
%   
%             
%             subplot(4,4,[9],'replace');
%             plot(qvzelj(7:10,max(1,i-300):i)');       
%             title(['D kuk&kol ž']);         
%                        
%             subplot(4,4,[10],'replace');
%             plot(qvzelj(11:14,max(1,i-300):i)');       
%             title(['D skz&prst ž']);         
%             
%             subplot(4,4,[11],'replace');
%             plot(qvzelj(15:18,max(1,i-300):i)');       
%             title(['L kuk&kol ž']);         
%                        
%             subplot(4,4,[12],'replace');
%             plot(qvzelj(19:22,max(1,i-300):i)');       
%             title(['L skz&prst ž']);         
%             
%             subplot(4,4,[13],'replace');
%             plot(qvzelj(segtrupx,max(1,i-300):i)');       
%             title(['T x osa ž']);         
% 
%             subplot(4,4,[14],'replace');
%             plot(qvzelj(segtrupy,max(1,i-300):i)');       
%             title(['T y osa ž']);         
%             
%             
%             subplot(4,4,[15],'replace');
%             plot(qvzelj([seglruka segdruka],max(1,i-300):i)');       
%             title(['Ruka ž']);  
%             
                        
%             subplot(4,4,[16],'replace');
%             plot(rbazerr(:,ikraj:i)');       
%             title([num2str(VecMod(rbazerr(:,i)))]);  
        
             pause(0.00001);
        end


       
        if(sacuvaj)
            sacuvaj=0;            
            save(['sim1_faza' num2str(faza) 'b' num2str(brzinakoraka_nom) 'v' num2str(visinakoraka) 'd' num2str(duzinakoraka_nom)  'n.mat']);
            save(['sim1_faza' num2str(faza) 'b' num2str(brzinakoraka_nom) 'v' num2str(visinakoraka) 'd' num2str(duzinakoraka_nom)   'n_iscrtaj.mat'],...
                        'qnew','rcm','ZMPpos','dt','zp','putanja_s');
        end
    end;


    save sim1nast -v6
end
function resetprimitivcontrolvar
        global primnaginjanjnapredakt brznaginjanjelnapred brzgrcidnogur brzgrcilnogur 
        global primosloninanoguakt primgrcinoguakt 
        global brzgrcidnogul brzgrcilnogul  primopruzinoguakt  brzopruzidnogul brzopruzidnogur brzopruzidnogu primspustinoguakt 
        global brzspustidnogul brzspustidnogur
        global primlnl primlnd primlrl primlrd  primlt nlbsy ndbsy tbsy rbsy rdbsy
        global statusnaginjanjenapred
        global primodgurninoguakt 
        global brzodgurnidnogu brzodgurnilnogu brzspustidpetu brzpocetnaspustidpetu brzspustilprst brzpocetnaspustilprst
        global brzspustilnogu brzspustidnogu brzpocetnapodigninogu
        global brzspustilpetu brzpocetnaspustilpetu brzspustidprst brzpocetnaspustidprst
        global  brzosloninalnogul brzosloninadnogul brzosloninalnogur brzosloninadnogur;
        global  brzspustilnogul brzspustilnogur;
        brzgrcidnogur=0; brzgrcilnogur=0;
        brzpocetnapodigninogu=0;
        brzspustidprst=0; brzpocetnaspustidprst=0;
        brzspustilpetu=0; brzpocetnaspustilpetu=0;
        brzspustilprst=0; brzpocetnaspustilprst=0;
        brzspustidpetu=0; brzpocetnaspustidpetu=0;
        primodgurninoguakt=0; brzodgurnidnogu=0; brzodgurnilnogu=0;
        statusnaginjanjenapred=0;
        %%%Svi primitivi postaju neaktiovni
        primosloninanoguakt=0; brzosloninalnogul=0;brzosloninalnogur=0; brzosloninadnogul=0;brzosloninadnogur=0;primgrcinoguakt=0; 
        brzgrcidnogul=0; brzgrcilnogul=0; primopruzinoguakt=0; brzopruzidnogu=0;  brzopruzidnogur=0;  brzopruzidnogul=0;  primspustinoguakt=0; 
        brzspustidnogul=0; brzspustidnogur=0; brzspustilnogul=0; brzspustilnogur=0;  brzspustilnogu=0; brzspustidnogu=0; 
        primnaginjanjnapredakt=0; brznaginjanjelnapred=0;
        %%postavljamo nivo svim primititvima na 0
        primlnl=0; primlnd=0; primlrl=0; primlrd=0;  primlt=0;
        %%oslobadjamo sve elstremitete        
        nlbsy=0; ndbsy=0; tbsy=0; rbsy=0; rdbsy=0;        
end

