function [qvzelj status] = primgrcinogu(qvzelj,flier,brz,odgurni1,odgurni2,i,noga)
    global seglnoga segdnoga segtrup seglruka segdruka 
    global primlnl primlnd primlrl primlrd  primlt nlbsy ndbsy tbsy rbsy rdbsy
    global rbaz rskzl rskzd rkukl rkukd rpetal rpetad rramel rramed rrukal rrukad rtrup rprstl rprstd 
    global XC YC XCL YCL XCR YCR
    global primodgurninoguakt  ZMPtr CMtr
    global brzodgurnidnogu brzodgurnilnogu brzpocetnaodgurni itpocodgurni
    global qnew qvnew qanew

    priml=3;
    status=0;
    if(noga=='L')
        if((priml<=primlnl) && (nlbsy~=0))
            if(primodgurninoguakt==0 || primodgurninoguakt==3)
                %%%Ako je zauzeta leva noga sa prioritetom visim od prioriteta
                %%%ovog primitiva status je jedan i vracamo se iz funkcija
                status=1;
                return;
            end
        else
            if(primodgurninoguakt==0)
%                 if(size(XCR,2)>2)
%                     k=convhull(XCR,YCR);
%                     [XCRm YCRm] = polyoffset(XCR(k),YCR(k),0.01);
%                     if sum(inpolygon([CMtr(1) ZMPtr(1)],[CMtr(2) ZMPtr(2)],XCRm,YCRm))==2
                        %%%Imamo dozvolu da zauzmemo nogu i postavljamo nas prioritet
                        nlbsy=1;
                        primodgurninoguakt=1;
                        primlnl=priml;
                        [Jkpetal Akpetal] =km_jakP(flier,18);
                        brzpocetnaodgurni = Jkpetal(:,:)*qvnew(:,i);
                        brzodgurnilnogu = VecMod(brzpocetnaodgurni(1:3));
                        itpocodgurni = i;
%                     end
%                 end
            end
        end               
% % %         if sum(inpolygon([CMtr(1) ZMPtr(1)],[CMtr(2) ZMPtr(2)],XCL,YCL))==2
            %%% Ovaj primitiv moze da se izvrsava samo ako se ZMP i CM
            %%% nalaze ispod druge noge
            
            %%%Kako izvrsiti primitiv?
            %%%Zadati da se stopalo dize na odgovarajucu koordinatu
            
            if(primodgurninoguakt==1)
                rz = 0.9*(1-odgurni1);
                rx = 0.9*odgurni2;
                status=5;
                rpetazelj = [rpetal(1)-rkukl(1)-0.1;rpetal(2)-rkukl(2);-rkukl(3)-0.1;0;pi/4;0];
                rpetatr = rpetal-rkukl;
                rpetaerr =rpetazelj-rpetatr;
                if(VecMod(rpetaerr(1:3))>0.1 || VecMod(rpetaerr(4:6))>0.01)                
                    % Nije se zavrsio primitiv
                    rpetaerrort = [(rpetaerr(1:3))/VecMod(rpetaerr(1:3)) ; (rpetaerr(4:6))/VecMod(rpetaerr(4:6))];
                    if(brzodgurnilnogu>max(3*VecMod(rpetaerr(1:3)),3*VecMod(rpetaerr(4:6))))
                        brzodgurnilnogu=max(3*VecMod(rpetaerr(1:3)),3*VecMod(rpetaerr(4:6)));
                    elseif(brzodgurnilnogu<brz)                        
                        brzodgurnilnogu=brzodgurnilnogu+0.0005*brz;
                    else                        
                        brzodgurnilnogu=brz;
                    end
%                     qvp = -brzodgurnidnogu*rpetaerrort;

                    if(mod(i-itpocodgurni,20)==0)
                        mirko=1;
                    end

                    qvpp = brzodgurnilnogu*rpetaerrort;
                    brztranz = satlin(((i-itpocodgurni)/200));
                    qvp = (1-brztranz).*brzpocetnaodgurni+brztranz.*qvpp;
                    
                    
                    
    %                 qvb(3)=-0.005;
                    [Jkpetal Akpetal] =km_jakP(flier,18);
                    qvzelj(seglnoga(5:7),i)=pinv(Jkpetal(:,seglnoga(5:7)))*qvp(:);
                    
                    qvzelj(seglnoga(7),i)=0.5;
                    qvzelj(seglnoga(8),i)=-0.5;
                    
                    
                else
                    % Zavrsio se primitiv
                    %%%Proveri da li je zauzeta leva noga i oslobodi je
                    if(primlnl==3 && nlbsy==1)
                        primlnl=0;
                        nlbsy=0;
                    end
                    status=3;
                    primodgurninoguakt=3;                            
                end
            elseif(primodgurninoguakt==3)
                % Zavrsio se primitiv                    
                status=3;
                nlbsy=0;
                primodgurninoguakt=3;                
            end

    elseif(noga=='D')
        if((priml<=primlnd) && (ndbsy~=0))
            if(primodgurninoguakt==0 || primodgurninoguakt==3)
                %%%Ako je zauzeta leva noga sa prioritetom visim od prioriteta
                %%%ovog primitiva status je jedan i vracamo se iz funkcija
                status=1;
                return;
            end
        else
            if(primodgurninoguakt==0)
%                 if(size(XCL,2)>2)
%                     k=convhull(XCL,YCL);
%                     [XCLm YCLm] = polyoffset(XCL(k),YCL(k),0.01);
%                     if sum(inpolygon([CMtr(1) ZMPtr(1)],[CMtr(2) ZMPtr(2)],XCLm,YCLm))==2
                        %%%Imamo dozvolu da zauzmemo nogu i postavljamo nas prioritet
                        ndbsy=1;
                        primodgurninoguakt=1;
                        primlnd=priml;
                        [Jkpetad Akpetad] =km_jakP(flier,19);
                        brzpocetnaodgurni = Jkpetad(:,:)*qvnew(:,i);
                        brzodgurnidnogu = VecMod(brzpocetnaodgurni(1:3));
                        itpocodgurni = i;
%                     end
%                 end
            end
        end               
% % %         if sum(inpolygon([CMtr(1) ZMPtr(1)],[CMtr(2) ZMPtr(2)],XCL,YCL))==2
            %%% Ovaj primitiv moze da se izvrsava samo ako se ZMP i CM
            %%% nalaze ispod druge noge
            
            %%%Kako izvrsiti primitiv?
            %%%Zadati da se stopalo dize na odgovarajucu koordinatu
            
            if(primodgurninoguakt==1)
                rz = 0.9*(1-odgurni1);
                rx = 0.9*odgurni2;
                status=5;
                rpetazelj = [rx;0;-rz;0;0;0];
                rpetatr = rpetad-rkukd;
                rpetaerr = rpetazelj-rpetatr;
                if(VecMod(rpetaerr(1:3))>0.1 || VecMod(rpetaerr(4:6))>0.01)                
                    % Nije se zavrsio primitiv
                    rpetaerrort = [(rpetaerr(1:3))/VecMod(rpetaerr(1:3)) ; (rpetaerr(4:6))/VecMod(rpetaerr(4:6))];
                    if(brzodgurnidnogu>max(3*VecMod(rpetaerr(1:3)),3*VecMod(rpetaerr(4:6))))
                        brzodgurnidnogu=max(3*VecMod(rpetaerr(1:3)),3*VecMod(rpetaerr(4:6)));
                    elseif(brzodgurnidnogu<brz)                        
                        brzodgurnidnogu=brzodgurnidnogu+0.0005*brz;
                    else                        
                        brzodgurnidnogu=brz;
                    end
%                     qvp = -brzodgurnidnogu*rpetaerrort;

                    if(mod(i-itpocodgurni,20)==0)
                        mirko=1;
                    end

                    qvpp = brzodgurnidnogu*rpetaerrort;
                    brztranz = satlin(((i-itpocodgurni)/200));
                    qvp = (1-brztranz).*brzpocetnaodgurni+brztranz.*qvpp;
                    
                    
                    
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
                    primodgurninoguakt=3;                            
                end
            elseif(primodgurninoguakt==3)
                % Zavrsio se primitiv                    
                status=3;
                ndbsy=0;
                primodgurninoguakt=3;                
            end
% % %         else
% % %             % Ako nema kontakta leve noge sa podlogom javljamo da ne moze
% % %             % da se izvrsi ovaj primitiv i oslobadjamo nogu            
% % %             status=2;
% % %             nlbsy=0;
% % %         end
    end
    
end