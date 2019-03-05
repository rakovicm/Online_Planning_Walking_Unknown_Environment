function [qvzelj status] = primospusti(qvzelj,flier,brz,brzspust,spust1,spust2,i,noga)
    global seglnoga segdnoga segtrup seglruka segdruka 
    global primlnl primlnd primlrl primlrd  primlt nlbsy ndbsy tbsy rbsy rdbsy
    global rbaz rskzl rskzd rkukl rkukd rpetal rpetad rramel rramed rrukal rrukad rtrup rprstl rprstd 
    global XC YC XCL YCL XCR YCR
    global primspustinoguakt  ZMPtr CMtr
    global brzspustidnogu brzpocetnaspusti itpocspusti
    global brzspustidpetu brzpocetnaspustidpetu brzspustilprst brzpocetnaspustilprst
    global brzspustilnogu
    global brzspustilpetu brzpocetnaspustilpetu brzspustidprst brzpocetnaspustidprst
    global qnew qvnew qanew rpetadzeljpoc rpetalzeljpoc rprstdzeljpoc rprstlzeljpoc

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
                brzspustidnogu = VecMod(brzpocetnaspusti(1:3));


                [Jkpetal Akpetal] =km_jakP(flier,2);    
                brzpocetnaspustilpetu = Jkpetal(:,:)*qvnew(:,i);
                brzspustilpetu = VecMod(brzpocetnaspustilpetu(1:3));
                rpetal2 = km_linkX(flier,2,'xyz');
                rpetalzeljpoc = [rpetal2(1:2);0;0;0;0];

                [Jkprstd Akprstd] =km_jakP(flier,26);    
                brzpocetnaspustidprst = Jkprstd(:,:)*qvnew(:,i);
                brzspustidprst = VecMod(brzpocetnaspustidprst(1:3));
                itpocspusti = i;

                rprstdzeljpoc = [rprstd(1:2);0;0;0;0]; 
            end
        end               
            if(primspustinoguakt==1)
                rz = 0.94*spust1;
                status=5;

                rbazzelj = [rpetal(1:2)-(rpetal(1:2)-rprstd(1:2))*0.5; rz;0;0;0];                
                rbazerr = rbazzelj-rbaz;
                
                rpetal2 = km_linkX(flier,2,'xyz');
                rpetalzelj = [rpetal2(1:2);0;0;0;rpetal2(6)];%[rpetad2(1:2);0;0;0;rpetad2(6)];                
                rpetalerr = rpetalzelj-(rpetal2);
                
                rprstdzelj = [rprstd(1:2);0;0;0;rprstd(6)];                
                rprstderr = rprstdzelj-(rprstd);
                
                if(size(XCL,2)>3)
                    status=3;
                end
                
                if(VecMod(rbazerr(1:3))>0.01 || VecMod(rbazerr(4:6))>0.01) 
                    if size(XCL,2)>0
                        ork=1;
                    else
                        ork=1;
                    end
                    
                    rbazort = [(rbazerr(1:3))/VecMod(rbazerr(1:3)) ; ork*(rbazerr(4:6))/VecMod(rbazerr(4:6))];

                    if(brzspustilnogu>max(10*VecMod(rbazerr(1:3)),10*VecMod(rbazerr(4:6))))
                        brzspustilnogu=max(10*VecMod(rbazerr(1:3)),10*VecMod(rbazerr(4:6)));
                    elseif(brzspustilnogu<brz)                        
                        brzspustilnogu=brzspustilnogu+0.002*brz;
                    else                        
                        brzspustilnogu=brz;
                    end
                    
                    rpetalort = [(rpetalerr(1:3))/VecMod(rpetalerr(1:3)) ; (rpetalerr(4:6))/VecMod(rpetalerr(4:6))];

                    if(brzspustilpetu>max(10*VecMod(rpetalerr(1:3)),10*VecMod(rpetalerr(4:6))))
                        brzspustilpetu=max(10*VecMod(rpetalerr(1:3)),10*VecMod(rpetalerr(4:6)));
                    elseif(brzspustilpetu<brzspust)                        
                        brzspustilpetu=brzspustilpetu+0.002*brzspust;
                    else                        
                        brzspustilpetu=brzspust;
                    end
                    
                    rprstdort = [(rprstderr(1:3))/VecMod(rprstderr(1:3)) ; (rprstderr(4:6))/VecMod(rprstderr(4:6))];

                    if(brzspustidprst>max(10*VecMod(rprstderr(1:3)),10*VecMod(rprstderr(4:6))))
                        brzspustidprst=max(10*VecMod(rprstderr(1:3)),10*VecMod(rprstderr(4:6)));
                    elseif(brzspustidprst<brz)                        
                        brzspustidprst=brzspustidprst+0.002*brz;
                    else                        
                        brzspustidprst=brz;
                    end
                    
                    qvbb = brzspustilnogu*rbazort;
                    brztranz = satlin(((i-itpocspusti)/200));
                    qvb = (1-brztranz).*brzpocetnaspusti+brztranz.*qvbb;     
%                     qvb = [10 .01 .01 .1 .1 .1]'.*qvb;
                    
                    qvpp = 1*brzspustilpetu*rpetalort;
                    brztranz = satlin(((i-itpocspusti)/200));
                    qvp = (1-brztranz).*brzpocetnaspustilpetu+brztranz.*qvpp;

                    
                    qvprr = brzspustidprst*rprstdort;
                    brztranz = satlin(((i-itpocspusti)/200));
                    qvpr = (1-brztranz).*brzpocetnaspustidprst+brztranz.*qvprr;

    %                 qvb(3)=-0.005;
                    [Jkprstd Akprstd] =km_jakP(flier,26);
                    qvzelj(segdnoga(2:8),i)=-pinv(Jkprstd(:,segdnoga(2:8)))*(qvb-qvpr);
                    
                    [Jketal2 Aketal2] =km_jakP(flier,2);
                    qvzelj(seglnoga(2:7),i)=-pinv(Jketal2(:,seglnoga(2:7)))*(qvb-qvp);
                    
                    prstd9 = km_linkX(flier,9,'xyz');
                    prstd11 =km_linkX(flier,11,'xyz');

                    a = VecMod([prstd11(1:2)-prstd9(1:2);0]);
                    b = prstd11(3)-prstd9(3);

                    ugao = atan2(b,a);
                    qvzelj(segdnoga(8),i) = 2*ugao;

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
                brzspustidnogu = VecMod(brzpocetnaspusti(1:3));


                [Jkpetad Akpetad] =km_jakP(flier,7);    
                brzpocetnaspustidpetu = Jkpetad(:,:)*qvnew(:,i);
                brzspustidpetu = VecMod(brzpocetnaspustidpetu(1:3));
                rpetad2 = km_linkX(flier,7,'xyz');
                rpetadzeljpoc = [rpetad2(1:2);0;0;0;0];

                [Jkprstl Akprstl] =km_jakP(flier,25);    
                brzpocetnaspustilprst = Jkprstl(:,:)*qvnew(:,i);
                brzspustilprst = VecMod(brzpocetnaspustilprst(1:3));
                itpocspusti = i;

                rprstlzeljpoc = [rprstl(1:2);0;0;0;0]; 
            end
        end               
% % %         if sum(inpolygon([CMtr(1) ZMPtr(1)],[CMtr(2) ZMPtr(2)],XCL,YCL))==2
            %%% Ovaj primitiv moze da se izvrsava samo ako se ZMP i CM
            %%% nalaze ispod druge noge
            
            %%%Kako izvrsiti primitiv?
            %%%Zadati da se stopalo dize na odgovarajucu koordinatu
            
            if(primspustinoguakt==1)
                rz = 0.94*spust1;
                status=5;

                rbazzelj = [rpetad(1:2)-(rpetad(1:2)-rprstl(1:2))*0.5; rz;0;0;0];                
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
                    
                    rbazort = [(rbazerr(1:3))/VecMod(rbazerr(1:3)) ; ork*(rbazerr(4:6))/VecMod(rbazerr(4:6))];

                    if(brzspustidnogu>max(10*VecMod(rbazerr(1:3)),10*VecMod(rbazerr(4:6))))
                        brzspustidnogu=max(10*VecMod(rbazerr(1:3)),10*VecMod(rbazerr(4:6)));
                    elseif(brzspustidnogu<brz)                        
                        brzspustidnogu=brzspustidnogu+0.002*brz;
                    else                        
                        brzspustidnogu=brz;
                    end
                    
                    rpetadort = [(rpetaderr(1:3))/VecMod(rpetaderr(1:3)) ; (rpetaderr(4:6))/VecMod(rpetaderr(4:6))];

                    if(brzspustidpetu>max(10*VecMod(rpetaderr(1:3)),10*VecMod(rpetaderr(4:6))))
                        brzspustidpetu=max(10*VecMod(rpetaderr(1:3)),10*VecMod(rpetaderr(4:6)));
                    elseif(brzspustidpetu<brzspust)                        
                        brzspustidpetu=brzspustidpetu+0.002*brzspust;
                    else                        
                        brzspustidpetu=brzspust;
                    end
                    
                    rprstlort = [(rprstlerr(1:3))/VecMod(rprstlerr(1:3)) ; (rprstlerr(4:6))/VecMod(rprstlerr(4:6))];

                    if(brzspustilprst>max(10*VecMod(rprstlerr(1:3)),10*VecMod(rprstlerr(4:6))))
                        brzspustilprst=max(10*VecMod(rprstlerr(1:3)),10*VecMod(rprstlerr(4:6)));
                    elseif(brzspustilprst<brz)                        
                        brzspustilprst=brzspustilprst+0.002*brz;
                    else                        
                        brzspustilprst=brz;
                    end
                    
                    qvbb = brzspustidnogu*rbazort;
                    brztranz = satlin(((i-itpocspusti)/300));
                    qvb = (1-brztranz).*brzpocetnaspusti+brztranz.*qvbb;     
%                     qvb = [10 .01 .01 .1 .1 .1]'.*qvb;
                    
                    qvpp = 1*brzspustidpetu*rpetadort;
                    brztranz = satlin(((i-itpocspusti)/300));
                    qvp = (1-brztranz).*brzpocetnaspustidpetu+brztranz.*qvpp;

                    
                    qvprr = brzspustilprst*rprstlort;
                    brztranz = satlin(((i-itpocspusti)/300));
                    qvpr = (1-brztranz).*brzpocetnaspustilprst+brztranz.*qvprr;

    %                 qvb(3)=-0.005;
                    [Jkprstl Akprstl] =km_jakP(flier,18);
                    qvzelj(seglnoga(2:7),i)=-pinv(Jkprstl(:,seglnoga(2:7)))*(qvb-qvpr);
                    
                    [Jketad2 Aketad2] =km_jakP(flier,7);
                    qvzelj(segdnoga(2:7),i)=-pinv(Jketad2(:,segdnoga(2:7)))*(qvb+qvp);
                    
                    prstl3 = km_linkX(flier,3,'xyz');
                    prstl5 =km_linkX(flier,5,'xyz');

                    a = VecMod([prstl5(1:2)-prstl3(1:2);0]);
                    b = prstl5(3)-prstl3(3);

                    ugao = atan2(b,a);
                    qvzelj(seglnoga(8),i) = 20*ugao;
                
                    
                    %%%Proveravamo da li je druga noga noga u kontaktu
%                     if size(XCR,2)>0
%                         %%%Proveravamo da li je zauzeta nekim visim
%                         %%%primitivom za ovu nogu
%                         if(ndbsy==0 && primlnd<1)
%                             primlnd=1;
%                             ndbsy=1;
%                         end
%                         if(ndbsy==1 && primlnd==1)
%                             [Jkpetad Akpetad] =km_jakP(flier,19);
%                             qvzelj(segdnoga,i)=-pinv(Jkpetad(:,segdnoga))*qvb;
%                         end
%                         
%                     end
%                     rpetaerrort = [0; 0; 0; (rpetaerr(4:6))/VecMod(rpetaerr(4:6))];
%                     qvp = 0.5*brztranz*rpetaerrort;
% 
%                     [Jkpetad Akpetad] =km_jakP(flier,19);
%                     
%                     qvzelj(segdnoga(5:7),i)=pinv(Jkpetad(4:6,segdnoga(5:7)))*qvp(4:6);
                    
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