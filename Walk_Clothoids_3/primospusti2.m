function [qvzelj status] = primospusti2(qvzelj,flier,brz,spust1,spust2,i,CMz,noga)
    global seglnoga segdnoga segtrup seglruka segdruka 
    global primlnl primlnd primlrl primlrd  primlt nlbsy ndbsy tbsy rbsy rdbsy
    global rbaz rskzl rskzd rkukl rkukd rpetal rpetad rramel rramed rrukal rrukad rtrup rprstl rprstd 
    global XC YC XCL YCL XCR YCR
    global primspustinoguakt  ZMPtr CMtr
    global brzspustidnogul brzspustidnogur brzspustilnogul brzspustilnogur brzpocetnaspusti itpocspusti
    global brzspustidpetu brzpocetnaspustidpetu brzspustilprst brzpocetnaspustilprst
    global brzspustilnogu
    global brzspustilpetu brzpocetnaspustilpetu brzspustidprst brzpocetnaspustidprst
    global qnew qvnew qanew rpetadzeljpoc rpetalzeljpoc rprstdzeljpoc rprstlzeljpoc
    global zeljeni_pravac
    priml=3;
    status=0;
        brztranz=0;
    if(noga=='L')
        if((priml<=primlnd) && (ndbsy~=0))
            if(primspustinoguakt==0 || primspustinoguakt==3)
                %%%Ako je zauzeta leva noga sa prioritetom visim od prioriteta
                %%%ovog primitiva status je jedan i vracamo se iz funkcija
                status=1;
                return;
            end
        else
            if(primspustinoguakt==0)
                %%%Imamo dozvolu da zauzmemo nogu i postavljamo nas prioritet
                ndbsy=1;
                primspustinoguakt=1;
                primlnd=priml;    

                [Jkbaz Akbaz] =km_jakP(flier,13);    
                brzpocetnaspusti = Jkbaz(:,:)*qvnew(:,i);
                brzspustilnogul = VecMod(brzpocetnaspusti(1:3));
                brzspustilnogur = VecMod(brzpocetnaspusti(1:3));
                itpocspusti = i;
            end
        end               
            if(primspustinoguakt==1)
                rz = spust1;
                status=5;
                rbazzelj = [CMz(1:2);rz;0; 0; zeljeni_pravac];
%                 rbazzelj = [rpetal(1:2)-(rpetal(1:2)-rprstd(1:2))*0.5; rz;0;0;0];                
                rbazerr = rbazzelj-rbaz;
                
                rpetal2 = km_linkX(flier,2,'xyz');
                rpetalzelj = [rpetal2(1:2);0;0;0;rpetal2(6)];%[rpetad2(1:2);0;0;0;rpetad2(6)];                
                rpetalerr = rpetalzelj-(rpetal2);
                
                rprstdzelj = [rprstd(1:2);0;0;0;rprstd(6)];                
                rprstderr = rprstdzelj-rprstd;
                
                if(size(XCL,2)>3)
                    status=3;
                end
                
                if(VecMod(rbazerr(1:3))>0.01 || VecMod(rbazerr(4:6))>0.01) 
                    if size(XCL,2)>0
                        ork=1;
                    else
                        ork=1;
                    end
                    
                    rpetalort = [(rpetalerr(1:3))/VecMod(rpetalerr(1:3)) ; (rpetalerr(4:6))/VecMod(rpetalerr(4:6))];
                    rprstdort = [(rprstderr(1:3))/VecMod(rprstderr(1:3)) ; (rprstderr(4:6))/VecMod(rprstderr(4:6))];
                    rbazort = [(rbazerr(1:3))/VecMod(rbazerr(1:3)) ; ork*(rbazerr(4:6))/VecMod(rbazerr(4:6))];
                    
                    rpetalort(isnan(rpetalort))=0;
                    rprstdort(isnan(rprstdort))=0;
                    rbazort(isnan(rbazort))=0;
                    
                    korak=0.01;
                    
                    if(brzspustilnogul>300*VecMod(rbazerr(1:3)))
                        brzspustilnogul=300*VecMod(rbazerr(1:3));
                    elseif(brzspustilnogul<(brz-korak))                        
                        brzspustilnogul=brzspustilnogul+korak*brz;
                    elseif(brzspustilnogul>(brz+korak))                                                
                        brzspustilnogul=brzspustilnogul-korak*brz;
                    else
                        brzspustilnogul=brz;
                    end

                    if(brzspustilnogur>300*VecMod(rbazerr(4:6)))
                        brzspustilnogur=300*VecMod(rbazerr(4:6));
                    elseif(brzspustilnogur<(brz-korak))                        
                        brzspustilnogur=brzspustilnogur+korak*brz;
                    elseif(brzspustilnogur>(brz+korak))                                                
                        brzspustilnogur=brzspustilnogur-korak*brz;
                    else
                        brzspustilnogur=brz;
                    end                   
                    
                    brzzalepilnogul=brzspustilnogul;
                    brzzalepilnogur=brzspustilnogur;
                    
                    brzzalepidnogul=brzspustilnogul;
                    brzzalepidnogur=brzspustilnogur;
                    
                    if(brzzalepidnogul>300*VecMod(rprstderr(1:3)))
                        brzzalepidnogul=300*VecMod(rprstderr(1:3));
                    end
                    if(brzzalepidnogur>300*VecMod(rprstderr(4:6)))
                        brzzalepidnogur=300*VecMod(rprstderr(4:6));
                    end
                                        
                    if(brzzalepilnogul>300*VecMod(rpetalerr(1:3)))
                        brzzalepilnogul=300*VecMod(rpetalerr(1:3));
                    end
                    if(brzzalepilnogur>300*VecMod(rpetalerr(4:6)))
                        brzzalepilnogur=300*VecMod(rpetalerr(4:6));
                    end      
                    
                    
                    brztranz = satlin(((i-itpocspusti)/300));
                                        
                    qvbb(1:3) = brzspustilnogul*rbazort(1:3);
                    qvbb(4:6) = brzspustilnogur*rbazort(4:6);
                    qvb = (1-brztranz).*brzpocetnaspusti+brztranz.*qvbb';     
%                     qvb = [10 .01 .01 .1 .1 .1]'.*qvb;
                    
                    qvpp(1:3) = brzzalepilnogul*rpetalort(1:3);
                    qvpp(4:6) = brzzalepilnogur*rpetalort(4:6);
                    qvp = brztranz.*qvpp';

                    
                    qvprr(1:3) = brzzalepidnogul*rprstdort(1:3);
                    qvprr(4:6) = brzzalepidnogur*rprstdort(4:6);
                    qvpr = brztranz.*qvprr';


    %                 qvb(3)=-0.005;
                    [Jkprstd Akprstd] =km_jakP(flier,26);
                    qvzelj(segdnoga([2:3 5:6 8]),i)=pinv(Jkprstd(:,segdnoga([2:3 5:6 8])))*(qvpr-qvb);
                    
                    [Jketal2 Aketal2] =km_jakP(flier,2);
                    qvzelj(seglnoga(2:7),i)=pinv(Jketal2(:,seglnoga(2:7)))*(qvp-qvb);
                    
                    prstd9 = km_linkX(flier,9,'xyz');
                    prstd11 =km_linkX(flier,11,'xyz');

                    a = VecMod([prstd11(1:2)-prstd9(1:2);0]);
                    b = prstd11(3)-prstd9(3);

%                     ugao = atan2(b,a);
%                     qvzelj(segdnoga(8),i) = 10*ugao;

                else
                    if(primlnl==3 && nlbsy==1)
                        primlnl=0;
                        ndbsy=0;
                    end
                    status=3;
                    primspustinoguakt=3;                            
                end
            elseif(primspustinoguakt==3)
                % Zavrsio se primitiv                    
                status=3;
                nlbsy=0;
                primspustinoguakt=3;                
            end

    elseif(noga=='D')
        if((priml<=primlnl) && (nlbsy~=0))
            if(primspustinoguakt==0 || primspustinoguakt==3)
                %%%Ako je zauzeta leva noga sa prioritetom visim od prioriteta
                %%%ovog primitiva status je jedan i vracamo se iz funkcija
                status=1;
                return;
            end
        else
            if(primspustinoguakt==0)
                %%%Imamo dozvolu da zauzmemo nogu i postavljamo nas prioritet
                nlbsy=1;
                primspustinoguakt=1;
                primlnl=priml;    

                [Jkbaz Akbaz] =km_jakP(flier,13);                
                brzpocetnaspusti = Jkbaz(:,:)*qvnew(:,i);
                brzspustidnogul = VecMod(brzpocetnaspusti(1:3));
                brzspustidnogur = VecMod(brzpocetnaspusti(4:6));
                itpocspusti = i;
            end
        end               
            if(primspustinoguakt==1)
                rz = spust1;
                status=5;
                rbazzelj = [CMz(1:2);rz;0; 0; zeljeni_pravac];
%                 rbazzelj = [rpetad(1:2)-(rpetad(1:2)-rprstl(1:2))*0.5; rbaz(3)-0.01;0;0;0];                
                rbazerr = rbazzelj-rbaz;
                
                rpetad2 = km_linkX(flier,7,'xyz');
                rpetadzelj = [rpetad2(1:2);0;0;0;rpetad2(6)];%[rpetad2(1:2);0;0;0;rpetad2(6)];                
                rpetaderr = rpetadzelj-(rpetad2);
                
                rprstlzelj = [rprstl(1:2);0;0;0;rprstl(6)];                
                rprstlerr = rprstlzelj-(rprstl);
                
                if(size(XCR,2)>3)
                    status=3;
                end
                
                if(VecMod(rbazerr(1:3))>0.01 || VecMod(rbazerr(4:6))>0.01) 
                    if size(XCR,2)>0
                        ork=1;
                    else
                        ork=1;
                    end
                    
                    rpetadort = [(rpetaderr(1:3))/VecMod(rpetaderr(1:3)) ; (rpetaderr(4:6))/VecMod(rpetaderr(4:6))];
                    rprstlort = [(rprstlerr(1:3))/VecMod(rprstlerr(1:3)) ; (rprstlerr(4:6))/VecMod(rprstlerr(4:6))];
                    rbazort = [(rbazerr(1:3))/VecMod(rbazerr(1:3)) ; ork*(rbazerr(4:6))/VecMod(rbazerr(4:6))];
                    
                    rpetadort(isnan(rpetadort))=0;
                    rprstlort(isnan(rprstlort))=0;
                    rbazort(isnan(rbazort))=0;
                    
                    korak=0.01;
                    
                    if(brzspustidnogul>300*VecMod(rbazerr(1:3)))
                        brzspustidnogul=300*VecMod(rbazerr(1:3));
                    elseif(brzspustidnogul<(brz-korak))                        
                        brzspustidnogul=brzspustidnogul+korak*brz;
                    elseif(brzspustidnogul>(brz+korak))                                                
                        brzspustidnogul=brzspustidnogul-korak*brz;
                    else
                        brzspustidnogul=brz;
                    end

                    if(brzspustidnogur>300*VecMod(rbazerr(4:6)))
                        brzspustidnogur=300*VecMod(rbazerr(4:6));
                    elseif(brzspustidnogur<(brz-korak))                        
                        brzspustidnogur=brzspustidnogur+korak*brz;
                    elseif(brzspustidnogur>(brz+korak))                                                
                        brzspustidnogur=brzspustidnogur-korak*brz;
                    else
                        brzspustidnogur=brz;
                    end
                    
                    
                    brzzalepilnogul=brzspustidnogul;
                    brzzalepilnogur=brzspustidnogur;
                    
                    brzzalepidnogul=brzspustidnogul;
                    brzzalepidnogur=brzspustidnogur;
                    
                    if(brzzalepidnogul>300*VecMod(rpetaderr(1:3)))
                        brzzalepidnogul=300*VecMod(rpetaderr(1:3));
                    end
                    if(brzzalepidnogur>300*VecMod(rpetaderr(4:6)))
                        brzzalepidnogur=300*VecMod(rpetaderr(4:6));
                    end
                                        
                    if(brzzalepilnogul>300*VecMod(rprstlerr(1:3)))
                        brzzalepilnogul=300*VecMod(rprstlerr(1:3));
                    end
                    if(brzzalepilnogur>300*VecMod(rprstlerr(4:6)))
                        brzzalepilnogur=300*VecMod(rprstlerr(4:6));
                    end      
                                        
                    
                    brztranz = satlin(((i-itpocspusti)/300));
                                        
                    qvbb(1:3) = brzspustidnogul*rbazort(1:3);
                    qvbb(4:6) = brzspustidnogur*rbazort(4:6);
                    qvb = (1-brztranz).*brzpocetnaspusti+brztranz.*qvbb';     
%                     qvb = [10 .01 .01 .1 .1 .1]'.*qvb;
                    
                    qvpp(1:3) = brzzalepidnogul*rpetadort(1:3);
                    qvpp(4:6) = brzzalepidnogur*rpetadort(4:6);
                    qvp = brztranz.*qvpp';

                    
                    qvprr(1:3) = brzzalepilnogul*rprstlort(1:3);
                    qvprr(4:6) = brzzalepilnogur*rprstlort(4:6);
                    qvpr = brztranz.*qvprr';

    %                 qvb(3)=-0.005;
                    [Jkprstl Akprstl] =km_jakP(flier,25);
                    qvzelj(seglnoga([2:3 5:6 8]),i)=pinv(Jkprstl(:,seglnoga([2:3 5:6 8])))*(qvpr-qvb);
                    
                    [Jketad2 Aketad2] =km_jakP(flier,7);
                    qvzelj(segdnoga(2:7),i)=pinv(Jketad2(:,segdnoga(2:7)))*(qvp-qvb);
%                     
                    prstl3 = km_linkX(flier,3,'xyz');
                    prstl5 =km_linkX(flier,5,'xyz');

                    a = VecMod([prstl5(1:2)-prstl3(1:2);0]);
                    b = prstl5(3)-prstl3(3);

%                     ugao = atan2(b,a);
%                     qvzelj(seglnoga(8),i) = 10*ugao; 
                else
                    % Zavrsio se primitiv
                    %%%Proveri da li je zauzeta leva noga i oslobodi je
                    if(primlnd==3 && ndbsy==1)
                        primlnd=0;
                        ndbsy=0;
                    end
                    status=3;
                    primspustinoguakt=3;                            
                end
            elseif(primspustinoguakt==3)
                % Zavrsio se primitiv                    
                status=3;
                ndbsy=0;
                primspustinoguakt=3;                
            end
% % %         else
% % %             % Ako nema kontakta leve noge sa podlogom javljamo da ne moze
% % %             % da se izvrsi ovaj primitiv i oslobadjamo nogu            
% % %             status=2;
% % %             nlbsy=0;
% % %         end
    end
end