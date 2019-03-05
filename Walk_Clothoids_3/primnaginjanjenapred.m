function [qvzelj status] = primnaginjanjnapred(qvzelj,flier,brz,opruz,CMz,i,noga)
    global seglnoga segdnoga segtrup seglruka segdruka 
    global primlnl primlnd primlrl primlrd  primlt nlbsy ndbsy tbsy rbsy rdbsy
    global rbaz rskzl rskzd rkukl rkukd rpetal rpetad rramel rramed rrukal rrukad rtrup rprstl rprstd 
    global XC YC XCL YCL XCR YCR XCPL YCPL XCPR YCPR 
    global primnaginjanjnapredakt  ZMPtr CMtr
    global brznaginjanjelnapredl brznaginjanjelnapredr brzpocetnanagninapred itpocnagninapred
    global qnew qvnew qanew statusnaginjanjenapred
    global zeljeni_pravac dt
    priml=5;
    status=0;
        brztranz=0;
    if(noga=='L')
        if((priml<=primlnl) && (nlbsy~=0))
            if(primnaginjanjnapredakt==0 || primnaginjanjnapredakt==3)
                %%%Ako je zauzeta leva noga sa prioritetom visim od prioriteta
                %%%ovog primitiva status je jedan i vracamo se iz funkcija
                statusnaginjanjenapred=1;
                return;
            end
        else
            if(primnaginjanjnapredakt==0)
                %%%Imamo dozvolu da zauzmemo nogu i postavljamo nas prioritet
                nlbsy=1;
                primnaginjanjnapredakt=1;
                primlnl=priml;
                
                [Jkbaz Akbaz] =km_jakP(flier,13);            
                brzpocetnanagninapred = Jkbaz(:,:)*qvnew(:,i);
                brznaginjanjelnapredl = VecMod(brzpocetnanagninapred(1:3));
                brznaginjanjelnapredr = VecMod(brzpocetnanagninapred(4:6));
                itpocnagninapred = i;
                
            end
        end               
        if size(XCL,2)>0
            %%% Ovaj primitiv moze da se izvrsava samo ako su bar tri
            %%% tacke u kontaktu sa podlogom
            
            %%%Kako izvrsiti primitiv?
            %%%Zadati da se bazni segment nadje na sredini stopala
            
            %mesto gde hocu da dodje centar mase je rpetal
            %zbog toga posmatram gresku
            if(primnaginjanjnapredakt==1||primnaginjanjnapredakt==2)
                if(statusnaginjanjenapred~=4)
                    statusnaginjanjenapred=5;
                end

                
                if(statusnaginjanjenapred==5)
                    rbazzelj = [CMz(1:2);opruz;0; 0; zeljeni_pravac];
                elseif(statusnaginjanjenapred==4)
                    rbazzelj = [CMz(1:2);rbaz(3);0; 0; zeljeni_pravac];
                end                

                rbazerr = rbazzelj-rbaz;
                Cme = [CMtr(1:2);0] - [CMz(1:2);0];                               

                if((VecMod([rbazerr(1:2);0])>0.01) || rbazerr(3)>0.02)% || (VecMod([Cme(1:2);0])>0.03))% && VecMod(Cme)>0.015)     
                    % Nije se zavrsio primitiv
                    if size(XCR,2)>0
                        ork=0;
                    else
                        ork=1;
                    end
                    rbazort = [(rbazerr(1:3))/VecMod(rbazerr(1:3)) ; ork*(rbazerr(4:6))/VecMod(rbazerr(4:6))];
                                        
                    korak=0.01;
                    
                    if(brznaginjanjelnapredl>200*VecMod(rbazerr(1:3)))
                        brznaginjanjelnapredl=200*VecMod(rbazerr(1:3));
                    elseif(brznaginjanjelnapredl<(brz-korak))                        
                        brznaginjanjelnapredl=brznaginjanjelnapredl+korak*brz;
                    elseif(brznaginjanjelnapredl>(brz+korak))                                                
                        brznaginjanjelnapredl=brznaginjanjelnapredl-korak*brz;
                    else
                        brznaginjanjelnapredl=brz;
                    end

                    if(brznaginjanjelnapredr>200*VecMod(rbazerr(4:6)))
                        brznaginjanjelnapredr=200*VecMod(rbazerr(4:6));
                    elseif(brznaginjanjelnapredr<(brz-korak))                        
                        brznaginjanjelnapredr=brznaginjanjelnapredr+korak*brz;
                    elseif(brznaginjanjelnapredr>(brz+korak))                                                
                        brznaginjanjelnapredr=brznaginjanjelnapredr-korak*brz;
                    else
                        brznaginjanjelnapredr=brz;
                    end


                    qvbb(1:3) = [1 1 2]'.*(brznaginjanjelnapredl*rbazort(1:3));
                    qvbb(4:6) = brznaginjanjelnapredr*rbazort(4:6);
                    
                    brztranz = satlin(((i-itpocnagninapred)/300));
                    qvb = (1-brztranz).*brzpocetnanagninapred+brztranz.*qvbb';
                    
                    [Jkpetal Akpetal] =km_jakP(flier,18);
                    if(statusnaginjanjenapred==4)
                        qvzelj(seglnoga([2:3 5:7]),i)=pinv(Jkpetal(:,seglnoga([2:3 5:7])))*(0-qvb);
                    else                            
                        qvzelj(seglnoga(2:7),i)=pinv(Jkpetal(:,seglnoga(2:7)))*(0-qvb);
                    end
                    
                    if(size(XCPL,2)>2)                        
                        k=convhull(XCPL,YCPL);
                        [XCPLm YCPLm] = polyoffset(XCPL(k),YCPL(k),0.01);
                        if sum(inpolygon([CMtr(1)],[CMtr(2)],XCPLm,YCPLm))==1
                            statusnaginjanjenapred=4;
                            %%%Ako su ZMP i CM unutar poligona pnda moze da
                            %%%se dize peta
%                             qvzelj(seglnoga(8),i)=-abs(pinv(Jkprstl(:,seglnoga(8)))*qvb);
                        end
                    end
                    
                    if((size(XCR,2)>0) && (statusnaginjanjenapred==4))
                        % Zavrsio se primitiv
                        %%%Proveri da li je zauzeta leva noga i oslobodi je
                        if(primlnl==5 && nlbsy==1)
                            primlnl=0;
                            nlbsy=0;
                        end
                        %%%Proveri da li je zauzeta desna noga i oslobodi je
                        if(ndbsy==1 && primlnd==1)                    
                            primlnd=0;
                            ndbsy=0;
                        end
                        statusnaginjanjenapred=3;
                        primnaginjanjnapredakt=3;
                    end
%                     size(XCR,2)
                    %%%Proveravamo da li je druga noga noga u kontaktu
                    if size(XCR,2)>0
                        %%%Proveravamo da li je zauzeta nekim visim
                        %%%primitivom za ovu nogu
                        if(ndbsy==0 && primlnd<1)
                            primlnd=1;
                            ndbsy=1;
                        end
                        if(ndbsy==1 && primlnd==1)
                            [Jkpetad Akpetad] =km_jakP(flier,19);
                            qvzelj(segdnoga(2:8),i)=-pinv(Jkpetad(:,segdnoga(2:8)))*qvb;
                        end
                        
                    end
                    

                else
                    % Zavrsio se primitiv
                    %%%Proveri da li je zauzeta leva noga i oslobodi je
                    if(primlnl==5 && nlbsy==1)
                        primlnl=0;
                        nlbsy=0;
                    end
                    %%%Proveri da li je zauzeta desna noga i oslobodi je
                    if(ndbsy==1 && primlnd==1)                    
                        primlnd=0;
                        ndbsy=0;
                    end
                    statusnaginjanjenapred=3;
                    primnaginjanjnapredakt=3;
                    %obaranje lagano brzina na 0
                    ett = exp(-dt/0.005);
                    qvzelj(seglnoga(2:7),i)=qvzelj(seglnoga(2:7),i)*(1-ett)+ett*qvzelj(seglnoga(2:7),i-1);
                    
                end
                    
            elseif(primnaginjanjnapredakt==3)
                    % Zavrsio se primitiv
                ett = exp(-dt/0.005);
                qvzelj(seglnoga(2:7),i)=qvzelj(seglnoga(2:7),i)*(1-ett)+ett*qvzelj(seglnoga(2:7),i-1);

                statusnaginjanjenapred=3;
                nlbsy=0;
                primnaginjanjnapredakt=3;                
            end
            prstl3 = km_linkX(flier,3,'xyz');
            prstl5 =km_linkX(flier,5,'xyz');

            a = VecMod([prstl5(1:2)-prstl3(1:2);0]);
            b = prstl5(3)-prstl3(3);

            ugao = atan2(b,a);
            qvzelj(seglnoga(8),i) = 10*ugao;

        else
            ett = exp(-dt/0.005);
            qvzelj(seglnoga(2:7),i)=qvzelj(seglnoga(2:7),i)*(1-ett)+ett*qvzelj(seglnoga(2:7),i-1);

            statusnaginjanjenapred=3;
            nlbsy=0;
            primnaginjanjnapredakt=3; 
            
            prstl3 = km_linkX(flier,3,'xyz');
            prstl5 =km_linkX(flier,5,'xyz');

            a = VecMod([prstl5(1:2)-prstl3(1:2);0]);
            b = prstl5(3)-prstl3(3);

            ugao = atan2(b,a);
            qvzelj(seglnoga(8),i) = 10*ugao;            
            % Ako nema kontakta leve noge sa podlogom javljamo da ne moze
            % da se izvrsi ovaj primitiv i oslobadjamo nogu            
%             status=2;
%             nlbsy=0;
        end
    elseif(noga=='D')
        if((priml<=primlnd) && (ndbsy~=0))
            if(primnaginjanjnapredakt==0 || primnaginjanjnapredakt==3)
                %%%Ako je zauzeta leva noga sa prioritetom visim od prioriteta
                %%%ovog primitiva status je jedan i vracamo se iz funkcija
                statusnaginjanjenapred=1;
                return;
            end
        else
            if(primnaginjanjnapredakt==0)
                %%%Imamo dozvolu da zauzmemo nogu i postavljamo nas prioritet
                ndbsy=1;
                primnaginjanjnapredakt=1;
                primlnd=priml;
                
                [Jkbaz Akbaz] =km_jakP(flier,13);            
                brzpocetnanagninapred = Jkbaz(:,:)*qvnew(:,i);
                brznaginjanjelnapredl = VecMod(brzpocetnanagninapred(1:3));
                brznaginjanjelnapredr = VecMod(brzpocetnanagninapred(4:6));
                itpocnagninapred = i;
                
            end
        end               
        if size(XCR,2)>0
            %%% Ovaj primitiv moze da se izvrsava samo ako su bar tri
            %%% tacke u kontaktu sa podlogom
            
            %%%Kako izvrsiti primitiv?
            %%%Zadati da se bazni segment nadje na sredini stopala
            
            %mesto gde hocu da dodje centar mase je rpetal
            %zbog toga posmatram gresku
            if(primnaginjanjnapredakt==1||primnaginjanjnapredakt==2)
                if(statusnaginjanjenapred~=4)
                    statusnaginjanjenapred=5;
                end

                
                if(statusnaginjanjenapred==5)
                    rbazzelj = [CMz(1:2);opruz;0; 0; zeljeni_pravac];
                elseif(statusnaginjanjenapred==4)
                    rbazzelj = [CMz(1:2);rbaz(3);0; 0; zeljeni_pravac];
                end                
                
                rbazerr = rbazzelj-rbaz;                
                Cme = [CMtr(1:2);0] - [CMz(1:2);0];                               
                
                if((VecMod([rbazerr(1:2);0])>0.01) || rbazerr(3)>0.02)% || (VecMod([Cme(1:2);0])>0.03))% && VecMod(Cme)>0.015)     
                    % Nije se zavrsio primitiv
                    if size(XCL,2)>0
                        ork=0;
                    else
                        ork=1;
                    end
                    rbazort = [(rbazerr(1:3))/VecMod(rbazerr(1:3)) ; ork*(rbazerr(4:6))/VecMod(rbazerr(4:6))];
                                          
                    korak=0.01;
                                      
                    if(brznaginjanjelnapredl>200*VecMod(rbazerr(1:3)))
                        brznaginjanjelnapredl=200*VecMod(rbazerr(1:3));
                    elseif(brznaginjanjelnapredl<(brz-korak))                        
                        brznaginjanjelnapredl=brznaginjanjelnapredl+korak*brz;
                    elseif(brznaginjanjelnapredl>(brz+korak))                                                
                        brznaginjanjelnapredl=brznaginjanjelnapredl-korak*brz;
                    else
                        brznaginjanjelnapredl=brz;
                    end

                    if(brznaginjanjelnapredr>200*VecMod(rbazerr(4:6)))
                        brznaginjanjelnapredr=200*VecMod(rbazerr(4:6));
                    elseif(brznaginjanjelnapredr<(brz-korak))                        
                        brznaginjanjelnapredr=brznaginjanjelnapredr+korak*brz;
                    elseif(brznaginjanjelnapredr>(brz+korak))                                                
                        brznaginjanjelnapredr=brznaginjanjelnapredr-korak*brz;
                    else
                        brznaginjanjelnapredr=brz;
                    end


                    qvbb(1:3) = [1 1 2]'.*(brznaginjanjelnapredl*rbazort(1:3));
                    qvbb(4:6) = brznaginjanjelnapredr*rbazort(4:6);
                    
                    brztranz = satlin(((i-itpocnagninapred)/300));
                    qvb = (1-brztranz).*brzpocetnanagninapred+brztranz.*qvbb';
                    
                    [Jkpetad Akpetad] =km_jakP(flier,19);
                    if(statusnaginjanjenapred==4)
                        qvzelj(segdnoga([2:3 5:7]),i)=pinv(Jkpetad(:,segdnoga([2:3 5:7])))*(0-qvb);
                    else                            
                        qvzelj(segdnoga(2:7),i)=pinv(Jkpetad(:,segdnoga(2:7)))*(0-qvb);
                    end

                    if((size(XCPR,2)>2)&& (statusnaginjanjenapred~=4))
                        k=convhull(XCPR,YCPR);
                        [XCPRm YCPRm] = polyoffset(XCPR(k),YCPR(k),0.01);
                        if sum(inpolygon([CMtr(1)],[CMtr(2)],XCPRm,YCPRm))==1
                            statusnaginjanjenapred=4;
                            %%%Ako su ZMP i CM unutar poligona pnda moze da
                            %%%se dize peta
%                             qvzelj(seglnoga(8),i)=-abs(pinv(Jkprstl(:,seglnoga(8)))*qvb);
                        end
                    end
                    
                    if((size(XCL,2)>0) && (statusnaginjanjenapred==4))
                        % Zavrsio se primitiv
                        %%%Proveri da li je zauzeta leva noga i oslobodi je
                        if(primlnd==5 && ndbsy==1)
                            primlnd=0;
                            ndbsy=0;
                        end
                        %%%Proveri da li je zauzeta desna noga i oslobodi je
                        if(nlbsy==1 && primlnl==1)                    
                            primlnl=0;
                            nlbsy=0;
                        end
                        statusnaginjanjenapred=3;
                        primnaginjanjnapredakt=3;
                    end
%                     size(XCR,2)
                    %%%Proveravamo da li je druga noga noga u kontaktu
                    if size(XCL,2)>0
                        %%%Proveravamo da li je zauzeta nekim visim
                        %%%primitivom za ovu nogu
                        if(nlbsy==0 && primlnl<1)
                            primlnl=1;
                            nlbsy=1;
                        end
                        if(nlbsy==1 && primlnl==1)
                            [Jkpetal Akpetal] =km_jakP(flier,18);
                            qvzelj(seglnoga(2:8),i)=-pinv(Jkpetal(:,seglnoga(2:8)))*qvb;
                        end
                        
                    end
                    

                else
                    % Zavrsio se primitiv
                    %%%Proveri da li je zauzeta leva noga i oslobodi je
                    if(primlnd==5 && ndbsy==1)
                        primlnd=0;
                        ndbsy=0;
                    end
                    %%%Proveri da li je zauzeta desna noga i oslobodi je
                    if(nlbsy==1 && primlnl==1)                    
                        primlnl=0;
                        nlbsy=0;
                    end
                    statusnaginjanjenapred=3;
                    primnaginjanjnapredakt=3;
                    ett = exp(-dt/0.005);
                    qvzelj(segdnoga(2:7),i)=qvzelj(segdnoga(2:7),i)*(1-ett)+ett*qvzelj(segdnoga(2:7),i-1);

                            
                end

                
                
                    
            elseif(primnaginjanjnapredakt==3)
                    % Zavrsio se primitiv
                ett = exp(-dt/0.005);
                qvzelj(segdnoga(2:7),i)=qvzelj(segdnoga(2:7),i)*(1-ett)+ett*qvzelj(segdnoga(2:7),i-1);
                    
                statusnaginjanjenapred=3;
                ndbsy=0;
                primnaginjanjnapredakt=3;                
            end
            
            prstd9 = km_linkX(flier,9,'xyz');
            prstd11 =km_linkX(flier,11,'xyz');

            a = VecMod([prstd11(1:2)-prstd9(1:2);0]);
            b = prstd11(3)-prstd9(3);

            ugao = atan2(b,a);
            qvzelj(segdnoga(8),i) = 10*ugao;
        else
            ett = exp(-dt/0.005);
            qvzelj(segdnoga(2:7),i)=qvzelj(segdnoga(2:7),i)*(1-ett)+ett*qvzelj(segdnoga(2:7),i-1);

            statusnaginjanjenapred=3;
            ndbsy=0;
            primnaginjanjnapredakt=3;    
            
            prstd9 = km_linkX(flier,9,'xyz');
            prstd11 =km_linkX(flier,11,'xyz');

            a = VecMod([prstd11(1:2)-prstd9(1:2);0]);
            b = prstd11(3)-prstd9(3);

            ugao = atan2(b,a);
            qvzelj(segdnoga(8),i) = 10*ugao;
            
            % Ako nema kontakta leve noge sa podlogom javljamo da ne moze
            % da se izvrsi ovaj primitiv i oslobadjamo nogu            
%             status=2;
%             nlbsy=0;
        end
    end
%         qvzelj(:,i) = (1-brztranz)*qvzelj(:,itpocnagninapred-1)+brztranz*qvzelj(:,i); 

end