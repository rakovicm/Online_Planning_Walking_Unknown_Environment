function [qvzelj status] = primgrcinogu(qvzelj,flier,brz,grci1,grci2,i,noga)
    global seglnoga segdnoga segtrup seglruka segdruka 
    global primlnl primlnd primlrl primlrd  primlt nlbsy ndbsy tbsy rbsy rdbsy
    global rbaz rskzl rskzd rkukl rkukd rpetal rpetad rramel rramed rrukal rrukad rtrup rprstl rprstd 
    global XC YC XCL YCL XCR YCR
    global primgrcinoguakt  ZMPtr CMtr CMzelj
    global brzgrcidnogul brzgrcilnogul brzgrcidnogur brzgrcilnogur brzpocetnagrci itpocgrci
    global qnew qvnew qanew
    global zeljeni_pravac
    global zeljeni_pravac_L
    global zeljeni_pravac_D
    priml=3;
    status=0;
    brztranz=0;
    if(noga=='L')
    priml=5;
        if((priml<=primlnl) && (nlbsy~=0))
            if(primgrcinoguakt==0 || primgrcinoguakt==3)
                %%%Ako je zauzeta leva noga sa prioritetom visim od prioriteta
                %%%ovog primitiva status je jedan i vracamo se iz funkcija
                status=1;
                return;
            end
        else
            if(primgrcinoguakt==0)
                if(size(XCR,2)>2)
%                     k=convhull(XCR,YCR);
%                     [XCRm YCRm] = polyoffset(XCR(k),YCR(k),0.02);
%                     if (sum(inpolygon([CMtr(1) ZMPtr(1)],[CMtr(2) ZMPtr(2)],XCRm,YCRm))==2) || (VecMod([CMzelj(1:2,i);0]-[CMtr(1:2);0])<0.05)
                        %%%Imamo dozvolu da zauzmemo nogu i postavljamo nas
                        %%%prioritet
                        nlbsy=1;
                        primgrcinoguakt=1;
                        primlnl=priml;
                        brzgrcilnogul=0;
                        [Jkpetal Akpetal] =km_jakP(flier,18);
                        brzpocetnagrci = 0*Jkpetal(:,:)*qvnew(:,i);
                        brzgrcilnogul = VecMod(brzpocetnagrci(1:3));
                        brzgrcilnogur = VecMod(brzpocetnagrci(4:6));
                        itpocgrci = i;
%                     end
                end
            end
        end               
% % %         if sum(inpolygon([CMtr(1) ZMPtr(1)],[CMtr(2) ZMPtr(2)],XCL,YCL))==2
            %%% Ovaj primitiv moze da se izvrsava samo ako se ZMP i CM
            %%% nalaze ispod druge noge
            
            %%%Kako izvrsiti primitiv?
            %%%Zadati da se stopalo dize na odgovarajucu koordinatu
            
            if(primgrcinoguakt==1)
%           STARO               
%                 rz = 0.9*(1-grci1);
%                 rx = 0.9*grci2;
                status=5;
%                 rpetatr = rpetal-rkukl;
%                        
% %                 if(rx<0)
% %                     rpetazelj = [rpetatr(1:2);-rz;rpetatr(4:6)];
% %                 else
%                     rpetazelj = [rx;0;-(1-grci1)*rbaz(3);0;0;zeljeni_pravac];
% %                 end
%                 
%                 
%                 rpetaerr = rpetazelj-rpetatr;
                
                %NOVO
                A=flier.A;
                rpetazelj(1:3,1) = [rbaz(1:2);0] + A(:,:,6)*[grci2;0.1019;grci1];
                rpetazelj(3) = grci1;
                rpetazelj(4:6) = TrotXYZ(A(:,:,6)*(RotX(0)*RotY(0)*RotZ(zeljeni_pravac_L))*A(:,:,6)')';
                rpetaerr = rpetazelj-rpetal;                
                
                if(VecMod(rpetaerr(1:3))>0.04)                
                    % Nije se zavrsio primitiv
                    rpetaerrort = [(rpetaerr(1:3))/VecMod(rpetaerr(1:3)) ; (rpetaerr(4:6))/VecMod(rpetaerr(4:6))];
                    rpetaerrort(isnan(rpetaerrort))=0;
                    
                    korak=0.01;

                    if(brzgrcilnogul>100*VecMod(rpetaerr(1:3)))
                        brzgrcilnogul=100*VecMod(rpetaerr(1:3));
                    elseif(brzgrcilnogul<(brz-korak))                        
                        brzgrcilnogul=brzgrcilnogul+korak*brz;
                    elseif(brzgrcilnogul>(brz+korak))                                                
                        brzgrcilnogul=brzgrcilnogul-korak*brz;
                    else
                        brzgrcilnogul=brz;
                    end


                    if(brzgrcilnogur>100*VecMod(rpetaerr(4:6)))
                        brzgrcilnogur=100*VecMod(rpetaerr(4:6));
                    elseif(brzgrcilnogur<(brz-korak))                        
                        brzgrcilnogur=brzgrcilnogur+korak*brz;
                    elseif(brzgrcilnogur>(brz+korak))                                                
                        brzgrcilnogur=brzgrcilnogur-korak*brz;
                    else
                        brzgrcilnogur=brz;
                    end

                    brztranz = satlin(((i-itpocgrci)/300));

                    qvpp(1:3) = brzgrcilnogul*rpetaerrort(1:3);
                    qvpp(4:6) = brzgrcilnogur*rpetaerrort(4:6);
                    qvp = (1-brztranz).*brzpocetnagrci+brztranz.*qvpp';
                    
                    %ako je bilo koja tacka u kontaktu samo idemo na gore i nista drugo
                    for kt = 1:6 
                        ktln(:,kt) = km_linkX(flier,kt,'xyz');
                    end
                    if(any(ktln(3,:)<0.03))
%                         qvp = qvp.*[0 0 1 0 0 0]';
                        qvp = [0 0 VecMod(qvp(1:3)) 0 0 0]';
                    end

    %                 qvb(3)=-0.005;
                    [Jkpetal Akpetal] =km_jakP(flier,18);
                    
                    qvzelj(seglnoga(1:7),i)=pinv(Jkpetal(:,seglnoga(1:7)))*qvp;
                    
                    if(any(ktln(3,:)<0.03))
                        
                    else
                        qvzelj(seglnoga([7 8]),i) = -0.01*qnew(seglnoga([7 8]),i);
                    end
                    
                    
                else
                    % Zavrsio se primitiv
                    %%%Proveri da li je zauzeta leva noga i oslobodi je
                    nlbsy=0;
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
    priml=5;    
        if((priml<primlnd) && (ndbsy~=0))
            if(primgrcinoguakt==0 || primgrcinoguakt==3)
                %%%Ako je zauzeta leva noga sa prioritetom visim od prioriteta
                %%%ovog primitiva status je jedan i vracamo se iz funkcija
                status=1;
                return;
            end
        else
            if(primgrcinoguakt==0)
                if(size(XCL,2)>2)
%                     k=convhull(XCL,YCL);
%                     [XCLm YCLm] = polyoffset(XCL(k),YCL(k),0.01);
%                     if sum(inpolygon([CMtr(1) ZMPtr(1)],[CMtr(2) ZMPtr(2)],XCLm,YCLm))==2
                        %%%Imamo dozvolu da zauzmemo nogu i postavljamo nas prioritet
                        ndbsy=1;
                        primgrcinoguakt=1;
                        primlnd=priml;
                        brzgrcidnogul=0;
                        [Jkpetad Akpetad] =km_jakP(flier,19);
                        brzpocetnagrci = Jkpetad(:,:)*qvnew(:,i);
                        brzgrcidnogul = VecMod(brzpocetnagrci(1:3));
                        brzgrcidnogur = VecMod(brzpocetnagrci(4:6));
                        itpocgrci = i;
%                     end
                end
            end
        end               
% % %         if sum(inpolygon([CMtr(1) ZMPtr(1)],[CMtr(2) ZMPtr(2)],XCL,YCL))==2
            %%% Ovaj primitiv moze da se izvrsava samo ako se ZMP i CM
            %%% nalaze ispod druge noge
            
            %%%Kako izvrsiti primitiv?
            %%%Zadati da se stopalo dize na odgovarajucu koordinatu
            
            if(primgrcinoguakt==1)
                %STARO
% %                 rz = 0.9*(1-grci1);
%                 rx = 0.9*grci2;                
%                 status=5;
% 
%                 rpetatr = rpetad-rkukd;
%                 
% %                 if(rx<0)
% %                     rpetazelj = [rpetatr(1:2);-rz;rpetatr(4:6)];
% %                 else
%                     rpetazelj = [rx;0;-(1-grci1)*rbaz(3);0;0;zeljeni_pravac];
% %                 end
% 
%                 
%                 rpetaerr = rpetazelj-rpetatr;

%               NOVO
                A=flier.A;
                rpetazelj(1:3,1) = [rbaz(1:2);0] + A(:,:,6)*[grci2;-0.1019;grci1];
                rpetazelj(3) = grci1;
                rpetazelj(4:6) = TrotXYZ(A(:,:,6)*(RotX(0)*RotY(0)*RotZ(zeljeni_pravac_D))*A(:,:,6)')';
                rpetaerr = rpetazelj-rpetad;   
                
                if(VecMod(rpetaerr(1:3))>0.04)                
                    % Nije se zavrsio primitiv
                    rpetaerrort = [(rpetaerr(1:3))/VecMod(rpetaerr(1:3)) ; (rpetaerr(4:6))/VecMod(rpetaerr(4:6))];
                    rpetaerrort(isnan(rpetaerrort))=0;
                    
                    korak=0.01;
                    
                    if(brzgrcidnogul>100*VecMod(rpetaerr(1:3)))
                        brzgrcidnogul=100*VecMod(rpetaerr(1:3));
                    elseif(brzgrcidnogul<(brz-korak))                        
                        brzgrcidnogul=brzgrcidnogul+korak*brz;
                    elseif(brzgrcidnogul>(brz+korak))                                                
                        brzgrcidnogul=brzgrcidnogul-korak*brz;
                    else
                        brzgrcidnogul=brz;
                    end
                    
                    if(brzgrcidnogur>100*VecMod(rpetaerr(4:6)))
                        brzgrcidnogur=100*VecMod(rpetaerr(4:6));
                    elseif(brzgrcidnogur<(brz-korak))                        
                        brzgrcidnogur=brzgrcidnogur+korak*brz;
                    elseif(brzgrcidnogur>(brz+korak))                                                
                        brzgrcidnogur=brzgrcidnogur-korak*brz;
                    else
                        brzgrcidnogur=brz;
                    end

                    brztranz = satlin(((i-itpocgrci)/300));

                    qvpp(1:3) = brzgrcidnogul*rpetaerrort(1:3);
                    qvpp(4:6) = brzgrcidnogur*rpetaerrort(4:6);
                    qvp = (1-brztranz).*brzpocetnagrci+brztranz.*qvpp';

                    for kt = 1:6 
                        ktln(:,kt) = km_linkX(flier,kt+6,'xyz');
                    end
                    if(any(ktln(3,:)<0.03))
%                         qvp = qvp.*[0 0 1 0 0 0]';
                        qvp = [0 0 VecMod(qvp(1:3)) 0 0 0]';
                    end
                    
                    
                    [Jkpetad Akpetad] =km_jakP(flier,19);
                    qvzelj(segdnoga(1:7),i)=pinv(Jkpetad(:,segdnoga(1:7)))*qvp;

                    if(any(ktln(3,:)<0.03))
                        
                    else
                         qvzelj(segdnoga([7 8]),i) = -0.01*qnew(segdnoga([7 8]),i); 
                    end
                    

                    
                else
                    % Zavrsio se primitiv
                    %%%Proveri da li je zauzeta leva noga i oslobodi je
                    status=3;
                    ndbsy=0;
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
%     qvzelj(:,i) = (1-brztranz)*qvzelj(:,itpocgrci-1)+brztranz*qvzelj(:,i); 
end