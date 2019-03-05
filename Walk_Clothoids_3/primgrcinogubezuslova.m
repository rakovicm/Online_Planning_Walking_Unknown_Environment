function [qvzelj status] = primgrcinogu(qvzelj,flier,brz,grci1,grci2,i,noga)
    global seglnoga segdnoga segtrup seglruka segdruka 
    global primlnl primlnd primlrl primlrd  primlt nlbsy ndbsy tbsy rbsy rdbsy
    global rbaz rskzl rskzd rkukl rkukd rpetal rpetad rramel rramed rrukal rrukad rtrup rprstl rprstd 
    global XC YC XCL YCL XCR YCR
    global primgrcinoguakt  ZMPtr CMtr
    global brzgrcidnogu brzgrcilnogu brzpocetnagrci itpocgrci
    global qnew qvnew qanew

    priml=3;
    status=0;
    if(noga=='L')
        if((priml<=primlnl) && (nlbsy~=0))
            if(primgrcinoguakt==0 || primgrcinoguakt==3)
                %%%Ako je zauzeta leva noga sa prioritetom visim od prioriteta
                %%%ovog primitiva status je jedan i vracamo se iz funkcija
                status=1;
                return;
            end
        else
            if(primgrcinoguakt==0)
%                 if(size(XCR,2)>2)
%                     k=convhull(XCR,YCR);
%                     [XCRm YCRm] = polyoffset(XCR(k),YCR(k),0.01);
%                     if sum(inpolygon([CMtr(1) ZMPtr(1)],[CMtr(2) ZMPtr(2)],XCRm,YCRm))==2
                        %%%Imamo dozvolu da zauzmemo nogu i postavljamo nas prioritet
                        nlbsy=1;
                        primgrcinoguakt=1;
                        primlnl=priml;
                        [Jkpetal Akpetal] =km_jakP(flier,18);
                        brzpocetnagrci = Jkpetal(:,:)*qvnew(:,i);
                        brzgrcilnogu = VecMod(brzpocetnagrci(1:3));
                        itpocgrci = i;
%                     end
%                 end
            end
        end               
% % %         if sum(inpolygon([CMtr(1) ZMPtr(1)],[CMtr(2) ZMPtr(2)],XCL,YCL))==2
            %%% Ovaj primitiv moze da se izvrsava samo ako se ZMP i CM
            %%% nalaze ispod druge noge
            
            %%%Kako izvrsiti primitiv?
            %%%Zadati da se stopalo dize na odgovarajucu koordinatu
            
            if(primgrcinoguakt==1)
                rz = 0.9*(1-grci1);
                rx = 0.9*grci2;
                status=5;
                rpetazelj = [rx;0;-rz;0;0;0];
                rpetatr = rpetal-rkukl;
                rpetaerr = rpetazelj-rpetatr;
                if(VecMod(rpetaerr(1:3))>0.1 || VecMod(rpetaerr(4:6))>0.01)                
                    % Nije se zavrsio primitiv
                    rpetaerrort = [(rpetaerr(1:3))/VecMod(rpetaerr(1:3)) ; (rpetaerr(4:6))/VecMod(rpetaerr(4:6))];
                    if(brzgrcilnogu>max(3*VecMod(rpetaerr(1:3)),3*VecMod(rpetaerr(4:6))))
                        brzgrcilnogu=max(3*VecMod(rpetaerr(1:3)),3*VecMod(rpetaerr(4:6)));
                    elseif(brzgrcilnogu<brz)                        
                        brzgrcilnogu=brzgrcilnogu+0.0005*brz;
                    else                        
                        brzgrcilnogu=brz;
                    end
%                     qvp = -brzgrcidnogu*rpetaerrort;

                    if(mod(i-itpocgrci,20)==0)
                        mirko=1;
                    end

                    qvpp = brzgrcilnogu*rpetaerrort;
                    brztranz = satlin(((i-itpocgrci)/200));
                    qvp = (1-brztranz).*brzpocetnagrci+brztranz.*qvpp;
                    
                    
                    
    %                 qvb(3)=-0.005;
                    [Jkpetal Akpetal] =km_jakP(flier,18);
                    qvzelj(seglnoga(2:7),i)=pinv(Jkpetal(:,seglnoga(2:7)))*qvp;
                else
                    % Zavrsio se primitiv
                    %%%Proveri da li je zauzeta leva noga i oslobodi je
                    if(primlnl==3 && nlbsy==1)
                        primlnl=0;
                        nlbsy=0;
                    end
                    status=3;
                    primgrcinoguakt=3;                            
                end
            elseif(primgrcinoguakt==3)
                % Zavrsio se primitiv                    
                status=3;
                nlbsy=0;
                primgrcinoguakt=3;                
            end

    elseif(noga=='D')
        if((priml<=primlnd) && (ndbsy~=0))
            if(primgrcinoguakt==0 || primgrcinoguakt==3)
                %%%Ako je zauzeta leva noga sa prioritetom visim od prioriteta
                %%%ovog primitiva status je jedan i vracamo se iz funkcija
                status=1;
                return;
            end
        else
            if(primgrcinoguakt==0)
%                 if(size(XCL,2)>2)
%                     k=convhull(XCL,YCL);
%                     [XCLm YCLm] = polyoffset(XCL(k),YCL(k),0.01);
%                     if sum(inpolygon([CMtr(1) ZMPtr(1)],[CMtr(2) ZMPtr(2)],XCLm,YCLm))==2
                        %%%Imamo dozvolu da zauzmemo nogu i postavljamo nas prioritet
                        ndbsy=1;
                        primgrcinoguakt=1;
                        primlnd=priml;
                        [Jkpetad Akpetad] =km_jakP(flier,19);
                        brzpocetnagrci = Jkpetad(:,:)*qvnew(:,i);
                        brzgrcidnogu = VecMod(brzpocetnagrci(1:3));
                        itpocgrci = i;
%                     end
%                 end
            end
        end               
% % %         if sum(inpolygon([CMtr(1) ZMPtr(1)],[CMtr(2) ZMPtr(2)],XCL,YCL))==2
            %%% Ovaj primitiv moze da se izvrsava samo ako se ZMP i CM
            %%% nalaze ispod druge noge
            
            %%%Kako izvrsiti primitiv?
            %%%Zadati da se stopalo dize na odgovarajucu koordinatu
            
            if(primgrcinoguakt==1)
                rz = 0.9*(1-grci1);
                rx = 0.9*grci2;
                status=5;
                rpetazelj = [rx;0;-rz;0;0;0];
                rpetatr = rpetad-rkukd;
                rpetaerr = rpetazelj-rpetatr;
                if(VecMod(rpetaerr(1:3))>0.1 || VecMod(rpetaerr(4:6))>0.01)                
                    % Nije se zavrsio primitiv
                    rpetaerrort = [(rpetaerr(1:3))/VecMod(rpetaerr(1:3)) ; (rpetaerr(4:6))/VecMod(rpetaerr(4:6))];
                    if(brzgrcidnogu>max(3*VecMod(rpetaerr(1:3)),3*VecMod(rpetaerr(4:6))))
                        brzgrcidnogu=max(3*VecMod(rpetaerr(1:3)),3*VecMod(rpetaerr(4:6)));
                    elseif(brzgrcidnogu<brz)                        
                        brzgrcidnogu=brzgrcidnogu+0.0005*brz;
                    else                        
                        brzgrcidnogu=brz;
                    end
%                     qvp = -brzgrcidnogu*rpetaerrort;

                    if(mod(i-itpocgrci,20)==0)
                        mirko=1;
                    end

                    qvpp = brzgrcidnogu*rpetaerrort;
                    brztranz = satlin(((i-itpocgrci)/200));
                    qvp = (1-brztranz).*brzpocetnagrci+brztranz.*qvpp;
                    
                    
                    
    %                 qvb(3)=-0.005;
                    [Jkpetad Akpetad] =km_jakP(flier,19);
                    qvzelj(segdnoga(2:7),i)=pinv(Jkpetad(:,segdnoga(2:7)))*qvp;
                else
                    % Zavrsio se primitiv
                    %%%Proveri da li je zauzeta leva noga i oslobodi je
                    if(primlnd==3 && ndbsy==1)
                        primlnd=0;
                        ndbsy=0;
                    end
                    status=3;
                    primgrcinoguakt=3;                            
                end
            elseif(primgrcinoguakt==3)
                % Zavrsio se primitiv                    
                status=3;
                ndbsy=0;
                primgrcinoguakt=3;                
            end
% % %         else
% % %             % Ako nema kontakta leve noge sa podlogom javljamo da ne moze
% % %             % da se izvrsi ovaj primitiv i oslobadjamo nogu            
% % %             status=2;
% % %             nlbsy=0;
% % %         end
    end
    
end