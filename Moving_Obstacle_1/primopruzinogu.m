function [qvzelj status] = primopruzi(qvzelj,flier,brz,opruz1,opruz2,sugao,skrugao,i,noga)
    global seglnoga segdnoga segtrup seglruka segdruka 
    global primlnl primlnd primlrl primlrd  primlt nlbsy ndbsy tbsy rbsy rdbsy
    global rbaz rskzl rskzd rkukl rkukd rpetal rpetad rramel rramed rrukal rrukad rtrup rprstl rprstd 
    global XC YC XCL YCL XCR YCR
    global primopruzinoguakt  ZMPtr CMtr
    global brzopruzidnogul brzopruzidnogur brzpocetnaopruzi itpocopruzi
    global qnew qvnew qanew
    global zeljeni_pravac
    global zeljeni_pravac_L
    global zeljeni_pravac_D
        
    priml=3;
    status=0;
        brztranz=0;
    if(noga=='L')
        if((priml<=primlnl) && (nlbsy~=0))
            if(primopruzinoguakt==0 || primopruzinoguakt==3)
                %%%Ako je zauzeta leva noga sa prioritetom visim od prioriteta
                %%%ovog primitiva status je jedan i vracamo se iz funkcija
                status=1;
                return;
            end
        else
            if(primopruzinoguakt==0)
%                 if(size(XCR,2)>2)
%                     k=convhull(XCR,YCR);
%                     [XCRm YCRm] = polyoffset(XCR(k),YCR(k),0.01);
%                     if sum(inpolygon([CMtr(1) ZMPtr(1)],[CMtr(2) ZMPtr(2)],XCRm,YCRm))==2
                        %%%Imamo dozvolu da zauzmemo nogu i postavljamo nas prioritet
                        nlbsy=1;
                        primopruzinoguakt=1;
                        primlnl=priml;    
                        
                        [Jkpetal Akpetal] =km_jakP(flier,18);
                        brzpocetnaopruzi = Jkpetal(:,:)*qvnew(:,i);
                        brzopruzidnogul = VecMod(brzpocetnaopruzi(1:3));
                        brzopruzidnogur = VecMod(brzpocetnaopruzi(4:6));
                        itpocopruzi = i;
%                     end
%                 end
            end
        end               
% % %         if sum(inpolygon([CMtr(1) ZMPtr(1)],[CMtr(2) ZMPtr(2)],XCL,YCL))==2
            %%% Ovaj primitiv moze da se izvrsava samo ako se ZMP i CM
            %%% nalaze ispod druge noge
            
            %%%Kako izvrsiti primitiv?
            %%%Zadati da se stopalo dize na odgovarajucu koordinatu
            
            if(primopruzinoguakt==1)
%                 rz = 1*(1-opruz1);
                rx = 1*opruz2;
                status=5;
                %STARO
%                 rpetazelj = [rx;-0.06;-(1-opruz1)*rbaz(3);0;sugao;skrugao];
%                 rpetatr = rpetal-rkukl;
%                 rpetaerr = rpetazelj-rpetatr;
                %NOVO
                A=flier.A;
                rpetazelj(1:3,1) = [rbaz(1:2);0] + A(:,:,6)*[rx;0.08;opruz1];
                rpetazelj(3)=opruz1;
                rpetazelj(4:6) = TRotXYZ(A(:,:,6)*(rotx(0)*roty(sugao)*rotz(zeljeni_pravac_L))*A(:,:,6)')';
                rpetaerr = rpetazelj-rpetal;
                
                if(VecMod(rpetaerr(1:3))>0.01 || VecMod(rpetaerr(4:6))>0.01)                
                    % Nije se zavrsio primitiv
                    rpetaerrort = [1 1 1 3 3 3]'.*[(rpetaerr(1:3))/VecMod(rpetaerr(1:3)) ; (rpetaerr(4:6))/VecMod(rpetaerr(4:6))];
                    
                    korak=0.01;
                    
                    if(brzopruzidnogul>100*VecMod(rpetaerr(1:3)))
                        brzopruzidnogul=100*VecMod(rpetaerr(1:3));
                    elseif(brzopruzidnogul<(brz-korak))                        
                        brzopruzidnogul=brzopruzidnogul+korak*brz;
                    elseif(brzopruzidnogul>(brz+korak))                                                
                        brzopruzidnogul=brzopruzidnogul-korak*brz;
                    else
                        brzopruzidnogul=brz;
                    end

                    if(brzopruzidnogur>100*VecMod(rpetaerr(4:6)))
                        brzopruzidnogur=100*VecMod(rpetaerr(4:6));
                    elseif(brzopruzidnogur<(brz-korak))                        
                        brzopruzidnogur=brzopruzidnogur+korak*brz;
                    elseif(brzopruzidnogur>(brz+korak))                                                
                        brzopruzidnogur=brzopruzidnogur-korak*brz;
                    else
                        brzopruzidnogur=brz;
                    end

                    qvpp(1:3) = brzopruzidnogul*rpetaerrort(1:3);
                    qvpp(4:6) = brzopruzidnogur*rpetaerrort(4:6);
                    
                    brztranz = satlin(((i-itpocopruzi)/300));
                    
                    qvp = (1-brztranz).*brzpocetnaopruzi+brztranz.*qvpp';
                    
    %                 qvb(3)=-0.005;
                    [Jkpetal Akpetal] =km_jakP(flier,18);
                    qvzelj(seglnoga(1:7),i)=[1 1 1 1 1 1 1]'.*(pinv(Jkpetal(:,seglnoga(1:7)))*qvp);

%                     qvzelj(seglnoga(5),i) = -10*qnew(seglnoga(5),i);
%                     qvzelj(seglnoga(6),i) = -10*qnew(seglnoga(6),i);
                    qvzelj(seglnoga(8),i) = -20*qnew(seglnoga(8),i);
  
                else
                    % Zavrsio se primitiv
                    %%%Proveri da li je zauzeta leva noga i oslobodi je
                    if(primlnl==3 && nlbsy==1)
                        primlnl=0;
                        nlbsy=0;
                    end
                    status=3;
                    primopruzinoguakt=3;                            
                end
                if(size(XCL,2)>0)
                    % Zavrsio se primitiv
                    % Nastala je dvooslonacka faza
                    if(primlnl==3 && nlbsy==1)
                        primlnl=0;
                        nlbsy=0;
                    end
                    status=3;
                    primopruzinoguakt=3;                            

                end
            elseif(primopruzinoguakt==3)
                % Zavrsio se primitiv                    
                status=3;
                nlbsy=0;
                primopruzinoguakt=3;                
            end
% % %         else
% % %             % Ako nema kontakta leve noge sa podlogom javljamo da ne moze
% % %             % da se izvrsi ovaj primitiv i oslobadjamo nogu            
% % %             status=2;
% % %             nlbsy=0;
% % %         end
    elseif(noga=='D')
        if((priml<=primlnd) && (ndbsy~=0))
            if(primopruzinoguakt==0 || primopruzinoguakt==3)
                %%%Ako je zauzeta leva noga sa prioritetom visim od prioriteta
                %%%ovog primitiva status je jedan i vracamo se iz funkcija
                status=1;
                return;
            end
        else
            if(primopruzinoguakt==0)
%                 if(size(XCL,2)>2)
%                     k=convhull(XCL,YCL);
%                     [XCLm YCLm] = polyoffset(XCL(k),YCL(k),0.01);
%                     if sum(inpolygon([CMtr(1) ZMPtr(1)],[CMtr(2) ZMPtr(2)],XCLm,YCLm))==2
                        %%%Imamo dozvolu da zauzmemo nogu i postavljamo nas prioritet
                        ndbsy=1;
                        primopruzinoguakt=1;
                        primlnd=priml;    
                        
                        [Jkpetad Akpetad] =km_jakP(flier,19);
                        brzpocetnaopruzi = Jkpetad(:,:)*qvnew(:,i);
                        brzopruzidnogul = VecMod(brzpocetnaopruzi(1:3));
                        brzopruzidnogur = VecMod(brzpocetnaopruzi(4:6));
                        itpocopruzi = i;
%                     end
%                 end
            end
        end               
% % %         if sum(inpolygon([CMtr(1) ZMPtr(1)],[CMtr(2) ZMPtr(2)],XCL,YCL))==2
            %%% Ovaj primitiv moze da se izvrsava samo ako se ZMP i CM
            %%% nalaze ispod druge noge
            
            %%%Kako izvrsiti primitiv?
            %%%Zadati da se stopalo dize na odgovarajucu koordinatu
            
            if(primopruzinoguakt==1)
%                 rz = 1*(1-opruz1);
                rx = 1*opruz2;
                status=5;
%                 rpetazelj = [rx;0.06;-(1-opruz1)*rbaz(3);0;sugao;skrugao];
%                 rpetatr = rpetad-rkukd;
%                 rpetaerr = rpetazelj-rpetatr;

                %NOVO
                A=flier.A;
                rpetazelj(1:3,1) = [rbaz(1:2);0] + A(:,:,6)*[rx;-0.08;opruz1];
                rpetazelj(3) = opruz1;
                rpetazelj(4:6) = TRotXYZ(A(:,:,6)*(rotx(0)*roty(sugao)*rotz(zeljeni_pravac_D))*A(:,:,6)')';
                rpetaerr = rpetazelj-rpetad;
                

                if(VecMod(rpetaerr(1:3))>0.01 || VecMod(rpetaerr(4:6))>0.01)                
                    % Nije se zavrsio primitiv
                    rpetaerrort = [1 1 1 3 3 3]'.*[(rpetaerr(1:3))/VecMod(rpetaerr(1:3)) ; (rpetaerr(4:6))/VecMod(rpetaerr(4:6))];
                    
                    korak=0.01;
                    
                    if(brzopruzidnogul>100*VecMod(rpetaerr(1:3)))
                        brzopruzidnogul=100*VecMod(rpetaerr(1:3));
                    elseif(brzopruzidnogul<(brz-korak))                        
                        brzopruzidnogul=brzopruzidnogul+korak*brz;
                    elseif(brzopruzidnogul>(brz+korak))                                                
                        brzopruzidnogul=brzopruzidnogul-korak*brz;
                    else
                        brzopruzidnogul=brz;
                    end

                    if(brzopruzidnogur>100*VecMod(rpetaerr(4:6)))
                        brzopruzidnogur=100*VecMod(rpetaerr(4:6));
                    elseif(brzopruzidnogur<(brz-korak))                        
                        brzopruzidnogur=brzopruzidnogur+korak*brz;
                    elseif(brzopruzidnogur>(brz+korak))                                                
                        brzopruzidnogur=brzopruzidnogur-korak*brz;
                    else
                        brzopruzidnogur=brz;
                    end

                    qvpp(1:3) = brzopruzidnogul*rpetaerrort(1:3);
                    qvpp(4:6) = brzopruzidnogur*rpetaerrort(4:6);
                    
                    brztranz = satlin(((i-itpocopruzi)/300));
                    
                    qvp = (1-brztranz).*brzpocetnaopruzi+brztranz.*qvpp';
                    
    %                 qvb(3)=-0.005;
                    [Jkpetad Akpetad] =km_jakP(flier,19);
                    qvzelj(segdnoga(1:7),i)= [1 1 1 1 1 1 1]'.*(pinv(Jkpetad(:,segdnoga(1:7)))*qvp);


%                     qvzelj(segdnoga(5),i) = -10*qnew(segdnoga(5),i);
%                     qvzelj(segdnoga(6),i) = -10*qnew(segdnoga(6),i);
                    qvzelj(segdnoga(8),i) = -20*qnew(segdnoga(8),i);
                else
                    % Zavrsio se primitiv
                    %%%Proveri da li je zauzeta leva noga i oslobodi je
                    if(primlnd==3 && ndbsy==1)
                        primlnd=0;
                        ndbsy=0;
                    end
                    status=3;
                    primopruzinoguakt=3;                            
                end
                if(size(XCR,2)>0)
                    % Zavrsio se primitiv
                    % Nastala je dvooslonacka faza
                    if(primlnd==3 && ndbsy==1)
                        primlnd=0;
                        ndbsy=0;
                    end
                    status=3;
                    primopruzinoguakt=3;                            

                end
            elseif(primopruzinoguakt==3)
                % Zavrsio se primitiv                    
                status=3;
                ndbsy=0;
                primopruzinoguakt=3;                
            end
% % %         else
% % %             % Ako nema kontakta leve noge sa podlogom javljamo da ne moze
% % %             % da se izvrsi ovaj primitiv i oslobadjamo nogu            
% % %             status=2;
% % %             nlbsy=0;
% % %         end
    end
    
%     qvzelj(:,i) = (1-brztranz)*qvzelj(:,itpocopruzi-1)+brztranz*qvzelj(:,i); 
    
end