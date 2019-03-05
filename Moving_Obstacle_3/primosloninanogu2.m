function [qvzelj status] = primosloninanogu(qvzelj,flier,brz,opruz,CMz,i,noga)
    global seglnoga segdnoga segtrup seglruka segdruka 
    global primlnl primlnd primlrl primlrd  primlt nlbsy ndbsy tbsy rbsy rdbsy
    global rbaz rskzl rskzd rkukl rkukd rpetal rpetad rramel rramed rrukal rrukad rtrup rprstl rprstd 
    global XC YC XCL YCL XCR YCR
    global primosloninanoguakt  ZMPtr CMtr
    global brzpocetnapodigninogu brzosloninalnogul brzosloninalnogur brzpocetnaosloninalnogu itpocosloninalnogu
    global brzosloninadnogul brzosloninadnogur brzpocetnaosloninadnogu itpocosloninadnogu 
    global brzspustilprst brzspustidnogu brzpocetnaspustilprst brzspustilnogu brzspustidprst
    global brzpocetnaosloninadnogu brzpocetnaspustidnogu brzpocetnaspustilprst brzpocetnaosloninalnogu brzpocetnaspustilnogu brzpocetnaspustidprst
    global zeljeni_pravac
    global qnew qvnew qanew
    
    priml=5;
    status=0;
    brztranz=0;
    if(noga=='L')
        if((priml<=primlnl) && (nlbsy~=0))
            if(primosloninanoguakt==0 || primosloninanoguakt==3)
                %%%Ako je zauzeta leva noga sa prioritetom visim od prioriteta
                %%%ovog primitiva status je jedan i vracamo se iz funkcija
                status=1;
                return;
            end
        else
            if(primosloninanoguakt==0)
                %%%Imamo dozvolu da zauzmemo nogu i postavljamo nas prioritet
                nlbsy=1;
                primosloninanoguakt=1;
                primlnl=priml;
                
                [Jkbaz Akbaz] =km_jakP(flier,13);       
                brzpocetnaosloninalnogu = Jkbaz(:,:)*qvnew(:,i);
                brzosloninalnogul = VecMod(brzpocetnaosloninalnogu(1:3));
                brzosloninalnogur = VecMod(brzpocetnaosloninalnogu(4:6));

                
                [Jknogal Aknogal] =km_jakP(flier,18);    
                brzpocetnaspustilnogu = Jknogal(:,:)*qvnew(:,i);
                brzspustilnogu = VecMod(brzpocetnaspustilnogu(1:3));
                
                                
                [Jkprstd Akprstd] =km_jakP(flier,26);    
                brzpocetnaspustidprst = Jkprstd(:,:)*qvnew(:,i);
                brzspustidprst = VecMod(brzpocetnaspustidprst(1:3));
                
                itpocosloninalnogu = i;
                                
            end
        end               
%         if size(XCL,2)>0
            %%% Ovaj primitiv moze da se izvrsava samo ako su bar tri
            %%% tacke u kontaktu sa podlogom
            
            %%%Kako izvrsiti primitiv?
            %%%Zadati da se bazni segment nadje na sredini stopala
            
            %mesto gde hocu da dodje centar mase je rpetal
            %zbog toga posmatram gresku
            if(primosloninanoguakt==1)
                rz = 1*opruz;

                status=5;
                
                rbazzelj = [CMz(1:2);rz;0; 0; zeljeni_pravac];
                rbazerr = rbazzelj-rbaz;                

                rlnoga = km_linkX(flier,18,'xyz');
                rlnogazelj = [rlnoga(1:2);0;0;0;rlnoga(6)];%[rpetad2(1:2);0;0;0;rpetad2(6)];                
                rlnogaerr = rlnogazelj-rlnoga;

                rdprst = km_linkX(flier,26,'xyz');
                rdprstzelj = [rdprst(1:2);0;0;0;rdprst(6)];%[rpetad2(1:2);0;0;0;rpetad2(6)];                
                rdprsterr = rdprstzelj-rdprst;

                Cme = [CMtr(1:2);0] - [CMz(1:2);0];
%                 VecMod([rbazerr(1:2) ;0]);
                
                if((VecMod([rbazerr(1:2);0])>0.03))%&& VecMod(Cme)>0.015)     
                    % Nije se zavrsio primitiv
%                     if size(XCR,2)>0
%                         ork=0;
%                     else
                        ork=1;
%                     end
                    rbazort = [(rbazerr(1:3))/VecMod(rbazerr(1:3)) ; ork*(rbazerr(4:6))/VecMod(rbazerr(4:6))];
                    rlnogaort = [(rlnogaerr(1:3))/VecMod(rlnogaerr(1:3)) ; (rlnogaerr(4:6))/VecMod(rlnogaerr(4:6))];                                       
                    rdprstort = [(rdprsterr(1:3))/VecMod(rdprsterr(1:3)) ; (rdprsterr(4:6))/VecMod(rdprsterr(4:6))];                                       
                    
                    rlnogaort(isnan(rlnogaort))=0;
                    rdprstort(isnan(rdprstort))=0;
                    
                    korak=0.002;
                    
                    if(brzosloninalnogul>100*VecMod(rbazerr(1:3)))
                        brzosloninalnogul=100*VecMod(rbazerr(1:3));
                    elseif(brzosloninalnogul<(brz-korak))                        
                        brzosloninalnogul=brzosloninalnogul+korak*brz;
                    elseif(brzosloninalnogul>(brz+korak))                                                
                        brzosloninalnogul=brzosloninalnogul-korak*brz;
                    else
                        brzosloninalnogul=brz;
                    end

                    if(brzosloninalnogur>100*VecMod(rbazerr(4:6)))
                        brzosloninalnogur=100*VecMod(rbazerr(4:6));
                    elseif(brzosloninalnogur<(brz-korak))                        
                        brzosloninalnogur=brzosloninalnogur+korak*brz;
                    elseif(brzosloninalnogur>(brz+korak))                                                
                        brzosloninalnogur=brzosloninalnogur-korak*brz;
                    else
                        brzosloninalnogur=brz;
                    end
                    
                    brzzalepilnogul=brzosloninalnogul;
                    brzzalepilnogur=brzosloninalnogur;
                    
                    brzzalepidnogul=brzosloninalnogul;
                    brzzalepidnogur=brzosloninalnogur;
                    
                    if(brzzalepidnogul>200*VecMod(rdprsterr(1:3)))
                        brzzalepidnogul=200*VecMod(rdprsterr(1:3));
                    end
                    if(brzzalepidnogur>200*VecMod(rdprsterr(4:6)))
                        brzzalepidnogur=200*VecMod(rdprsterr(4:6));
                    end
                                        
                    if(brzzalepilnogul>200*VecMod(rlnogaerr(1:3)))
                        brzzalepilnogul=200*VecMod(rlnogaerr(1:3));
                    end
                    if(brzzalepilnogur>200*VecMod(rlnogaerr(4:6)))
                        brzzalepilnogur=200*VecMod(rlnogaerr(4:6));
                    end                    
                    
                    brztranz = satlin(((i-itpocosloninalnogu)/300));

                    qvbb(1:3) = brzosloninalnogul*rbazort(1:3);
                    qvbb(4:6) = brzosloninalnogur*rbazort(4:6);
                    qvb = 0*(1-brztranz)*brzpocetnaosloninalnogu+brztranz.*qvbb';

                    qvpp(1:3) = brzzalepilnogul*rlnogaort(1:3);
                    qvpp(4:6) = brzzalepilnogur*rlnogaort(4:6);
                    qvp = 0*(1-brztranz)*brzpocetnaspustilnogu+brztranz.*qvpp';
                    
                    qvpr(1:3) = brzzalepidnogul*rdprstort(1:3);
                    qvpr(4:6) = brzzalepidnogur*rdprstort(4:6);
                    qvr = 0*(1-brztranz)*brzpocetnaspustidprst+brztranz.*qvpr';
                    
                    [Jkpetal Akpetal] =km_jakP(flier,18);
                    qvzelj(seglnoga(2:7),i)=pinv(Jkpetal(:,seglnoga(2:7)))*(qvp-qvb);
                    
                    [Jkprstd Akprstd] =km_jakP(flier,26);
                    qvzelj(segdnoga([2:3 5:6 8]),i)=pinv(Jkprstd(:,segdnoga([2:3 5:6 8])))*(qvr-qvb);
%                     qvzelj(segdnoga([7]),i)=0.01; 
                    
                    prstd9 = km_linkX(flier,9,'xyz');
                    prstd11 =km_linkX(flier,11,'xyz');

                    a = VecMod([prstd11(1:2)-prstd9(1:2);0]);
                    b = prstd11(3)-prstd9(3);

%                     ugao = atan2(b,a);
%                     qvzelj(segdnoga(8),i) = qvzelj(segdnoga(8),i) + 10*(ugao-qnew(segdnoga(8),i-1));
%                     qvzelj(segdnoga(8),i) = qvzelj(segdnoga(8),i-1) -0.1*qnew(segdnoga(8),i);
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
                    status=3;
                    primosloninanoguakt=3;
                            
                end
            elseif(primosloninanoguakt==3)
                    % Zavrsio se primitiv
                status=3;
                nlbsy=0;
                primosloninanoguakt=3;                
            end
%         else
            % Ako nema kontakta leve noge sa podlogom javljamo da ne moze
            % da se izvrsi ovaj primitiv i oslobadjamo nogu            
%             status=2;
%             nlbsy=0;
%         end
    end
    
    if(noga=='D')
        if((priml<=primlnd) && (ndbsy~=0))
            if(primosloninanoguakt==0 || primosloninanoguakt==3)
                %%%Ako je zauzeta leva noga sa prioritetom visim od prioriteta
                %%%ovog primitiva status je jedan i vracamo se iz funkcija
                status=1;
                return;
            end
        else
            if(primosloninanoguakt==0)
                %%%Imamo dozvolu da zauzmemo nogu i postavljamo nas prioritet
                ndbsy=1;
                primosloninanoguakt=1;
                primlnd=priml;
                
                [Jkbaz Akbaz] =km_jakP(flier,13);                
                brzpocetnaosloninadnogu = Jkbaz(:,:)*qvnew(:,i);
                brzosloninadnogul = VecMod(brzpocetnaosloninadnogu(1:3));
                brzosloninadnogur = VecMod(brzpocetnaosloninadnogu(4:6));
                itpocosloninadnogu = i;
                
                [Jknogad Aknogad] =km_jakP(flier,19);    
                brzpocetnaspustidnogu = Jknogad(:,:)*qvnew(:,i);
                brzspustidnogu = VecMod(brzpocetnaspustidnogu(1:3));
                
                                
                [Jkprstl Akprstl] =km_jakP(flier,25);    
                brzpocetnaspustilprst = Jkprstl(:,:)*qvnew(:,i);
                brzspustilprst = VecMod(brzpocetnaspustilprst(1:3));
               
                itpocosloninalnogu = i;

            end
        end               
%         if size(XCR,2)>0
            %%% Ovaj primitiv moze da se izvrsava samo ako su bar tri
            %%% tacke u kontaktu sa podlogom
            
            %%%Kako izvrsiti primitiv?
            %%%Zadati da se bazni segment nadje na sredini stopala
            
            %mesto gde hocu da dodje centar mase je rpetal
            %zbog toga posmatram gresku
            if(primosloninanoguakt==1)
                rz = 1*opruz;

                status=5;
                
                rbazzelj = [CMz(1:2);rz;0; 0; zeljeni_pravac];
                rbazerr = rbazzelj-rbaz;                

                rdnoga = km_linkX(flier,19,'xyz');
                rdnogazelj = [rdnoga(1:2);0;0;0;rdnoga(6)];
                rdnogaerr = rdnogazelj-rdnoga;                
                
                rlprst = km_linkX(flier,25,'xyz');
                rlprstzelj = [rlprst(1:2);0;0;0;rlprst(6)];
                rlprsterr = rlprstzelj-rlprst;                
                
                Cme = [CMtr(1:2);0] - [CMz(1:2);0];
%                 VecMod([rbazerr(1:2) ;0])

                if((VecMod([rbazerr(1:2);0])>0.03))% && VecMod(Cme)>0.015)     
                    % Nije se zavrsio primitiv
%                     if size(XCL,2)>0
%                         ork=0;
%                     else
                        ork=1;
%                     end
                    rbazort = [(rbazerr(1:3))/VecMod(rbazerr(1:3)) ; ork*(rbazerr(4:6))/VecMod(rbazerr(4:6))];
                    rdnogaort = [(rdnogaerr(1:3))/VecMod(rdnogaerr(1:3)) ; ork*(rdnogaerr(4:6))/VecMod(rdnogaerr(4:6))];
                    rlprstort = [(rlprsterr(1:3))/VecMod(rlprsterr(1:3)) ; ork*(rlprsterr(4:6))/VecMod(rlprsterr(4:6))];
                    
                    rdnogaort(isnan(rdnogaort))=0;
                    rlprstort(isnan(rlprstort))=0;
                    
                    korak=0.002;
                    
                    if(brzosloninadnogul>100*VecMod(rbazerr(1:3)))
                        brzosloninadnogul=100*VecMod(rbazerr(1:3));
                    elseif(brzosloninadnogul<(brz-korak))                        
                        brzosloninadnogul=brzosloninadnogul+korak*brz;
                    elseif(brzosloninadnogul>(brz+korak))                                                
                        brzosloninadnogul=brzosloninadnogul-korak*brz;
                    else
                        brzosloninadnogul=brz;
                    end

                    if(brzosloninadnogur>100*VecMod(rbazerr(4:6)))
                        brzosloninadnogur=100*VecMod(rbazerr(4:6));
                    elseif(brzosloninadnogur<(brz-korak))                        
                        brzosloninadnogur=brzosloninadnogur+korak*brz;
                    elseif(brzosloninadnogur>(brz+korak))                                                
                        brzosloninadnogur=brzosloninadnogur-korak*brz;
                    else
                        brzosloninadnogur=brz;
                    end
                    
                    brzzalepilnogul=brzosloninadnogul;
                    brzzalepilnogur=brzosloninadnogur;
                    
                    brzzalepidnogul=brzosloninadnogul;
                    brzzalepidnogur=brzosloninadnogur;
                    
                    if(brzzalepilnogul>200*VecMod(rlprsterr(1:3)))
                        brzzalepilnogul=200*VecMod(rlprsterr(1:3));
                    end
                    if(brzzalepilnogur>200*VecMod(rlprsterr(4:6)))
                        brzzalepilnogur=200*VecMod(rlprsterr(4:6));
                    end
                                        
                    if(brzzalepidnogul>200*VecMod(rdnogaerr(1:3)))
                        brzzalepidnogul=200*VecMod(rdnogaerr(1:3));
                    end
                    if(brzzalepidnogur>200*VecMod(rdnogaerr(4:6)))
                        brzzalepidnogur=200*VecMod(rdnogaerr(4:6));
                    end
                    
                    
                    brztranz = satlin(((i-itpocosloninalnogu)/300));
                     
                    qvbb(1:3) = brzosloninadnogul*rbazort(1:3);
                    qvbb(4:6) = brzosloninadnogur*rbazort(4:6);
                    qvb = 0*(1-brztranz)*brzpocetnaosloninadnogu+brztranz.*qvbb';
                    
                    qvpp(1:3) = brzzalepidnogul*rdnogaort(1:3);
                    qvpp(4:6) = brzzalepidnogur*rdnogaort(4:6);
                    qvp = 0*(1-brztranz)*brzpocetnaspustidnogu+brztranz.*qvpp';
                    
                    qvpr(1:3) = brzzalepilnogul*rlprstort(1:3);
                    qvpr(4:6) = brzzalepilnogur*rlprstort(4:6);
                    qvr = 0*(1-brztranz)*brzpocetnaspustilprst+brztranz.*qvpr';
                    
                    [Jkpetad Akpetad] =km_jakP(flier,19);
                    qvzelj(segdnoga([2:7]),i)=pinv(Jkpetad(:,segdnoga([2:7])))*(qvp-qvb);

                    [Jkprstl Akprstl] =km_jakP(flier,25);
                    qvzelj(seglnoga([2:3 5:6 8]),i)=pinv(Jkprstl(:,seglnoga([2:3 5:6 8])))*(qvr-qvb);
%                     qvzelj(seglnoga([7]),i)=0.01; 

                    prstl3 = km_linkX(flier,3,'xyz');
                    prstl5 =km_linkX(flier,5,'xyz');

                    a = VecMod([prstl5(1:2)-prstl3(1:2);0]);
                    b = prstl5(3)-prstl3(3);

%                     ugao = atan2(b,a);
%                     qvzelj(seglnoga(8),i) = qvzelj(seglnoga(8),i) + 10*(ugao-qnew(seglnoga(8),i));     
%                     qvzelj(seglnoga(8),i) = qvzelj(seglnoga(8),i-1)-0.1*qnew(seglnoga(8),i);
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
                    status=3;
                    primosloninanoguakt=3;
                            
                end
            elseif(primosloninanoguakt==3)
                    % Zavrsio se primitiv
                status=3;
                ndbsy=0;
                primosloninanoguakt=3;                
            end
%         else
            % Ako nema kontakta leve noge sa podlogom javljamo da ne moze
            % da se izvrsi ovaj primitiv i oslobadjamo nogu            
%             status=2;
%             nlbsy=0;
        end
%     end
%      qvzelj(:,i) = (1-brztranz)*qvzelj(:,itpocosloninalnogu-1)+brztranz*qvzelj(:,i); 
end