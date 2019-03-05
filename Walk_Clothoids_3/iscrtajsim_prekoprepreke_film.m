function iscrtajsim(fps,FileName,Output2)


global Up t Fuks ZMPtr CMtr
global XC YC XCL YCL XCR YCR XCPL YCPL XCPR YCPR TAC
global segdkuk segdkol segdskz segdprst seglkuk seglkol seglskz seglprst 
global segtrupx segtrupy segdrame segdlakat seglrame segllakat 
global seglnoga segdnoga segtrup seglruka segdruka 

load(FileName)

if(nargin==1)
    Output2=0;
end

%Inicijalizacija promenljivih
letac_motor_kontakti;

if(Output2~=0)
     fig=gcf;
     set(fig,'DoubleBuffer','on');
     set(gcf,'NextPlot','replace','Visible','off');
     mov = avifile(Output2);
end


zelj_p = zp(:,:);

ikraj = sum(qnew(3,:)>0);
incr = round(1/fps/dt);
incr = 100;

i=1;
c=1;
cc='rrr';
% while(i<ikraj)
for i = 1 : incr : ikraj
        i
        km_geo(flier,qnew(:,i));
        
        rbaz = km_linkX(flier,13,'xyz');
        rskzl = km_linkX(flier,14,'xyz');
        rskzd = km_linkX(flier,15,'xyz');
        rkukl = km_linkX(flier,16,'xyz');
        rkukd = km_linkX(flier,17,'xyz');
        rpetal = km_linkX(flier,18,'xyz');
        rpetad = km_linkX(flier,19,'xyz');
        rramel = km_linkX(flier,20,'xyz');
        rramed = km_linkX(flier,21,'xyz');
        rrukal = km_linkX(flier,22,'xyz');
        rrukad = km_linkX(flier,23,'xyz');
        rtrup = km_linkX(flier,24,'xyz');
        rprstl = km_linkX(flier,25,'xyz');
        rprstd = km_linkX(flier,26,'xyz');


        XC=[];
        YC=[];
        XCL=[];
        YCL=[];
        XCR=[];
        YCR=[];
        XCPL=[];
        YCPL=[];
        XCPR=[];
        YCPR=[];

        for ik=1:12
            Rk=km_linkX(flier,ik,'xyz');
            if Rk(3)<0
                XC=[XC Rk(1)];
                YC=[YC Rk(2)];
                if ik>0 && ik<7 
                    XCL=[XCL Rk(1)];
                    YCL=[YCL Rk(2)];
                    if ik>2
                        XCPL=[XCPL Rk(1)];
                        YCPL=[YCPL Rk(2)];
                    end            
                end 
                if ik>6 && ik<13
                    XCR=[XCR Rk(1)];
                    YCR=[YCR Rk(2)];
                    if ik>8
                        XCPR=[XCPR Rk(1)];
                        YCPR=[YCPR Rk(2)];
                    end
                end   
            end
        end
        
        xprep=0.58;
        yprep=0.09;
        zprep=0.2;

        if(mod(i,1000)==1)
            if(c==1)
                c=2;
                cc='rrr';
            else
                c=1;
                cc='bbb';
            end

            %         flier=kk_kinedyn(flier,qvnew(:,i));
    %         subplot(4,2,[1 3],'replace');
%             subplot(4,2,[1 2 3 4]);
% 
%             hold on;
%             draw_fliermirko(flier,seg, 42,[13 21],[14 22],[rcm(1,i)-1 rcm(1,i)+1 -0.5 0.5 -0.1 1.9],[sin(zelj_p(1,i)) -cos(zelj_p(1,i)) 0],[],2,cc);
%     %             for pri = 0:0.01:0.1
%     %                 line('XData',[xprep+0.5*pri xprep+0.5*pri],'YData',[yprep+pri yprep+pri],'ZData',[0 zprep],'Color','m','LineWidth',4);  
%     %             end
%     %         plot3(rcm(1,i),rcm(2,i),rcm(3,i),'ko');    
%             plot(putanja_s(1,:,1),putanja_s(2,:,1),'b');
%             plot(putanja_s(1,:,2),putanja_s(2,:,2),'r');
%             plot(putanja_s(1,:,3),putanja_s(2,:,3),'m');
%             crtaj_zid();
% 
%     %         plot(putanja_s(1,:,4),putanja_s(2,:,4),'c');
%     %         plot(putanja_s(1,:,5),putanja_s(2,:,5),'k');
%             axis([rcm(1,i)-1 rcm(1,i)+1 -0.5 0.5 -0.1 1.9]);
% 
    %         subplot(4,2,[2 4] ,'replace');
            subplot(4,2,[1 2 3 4] );
            hold on;
            draw_fliermirko(flier,seg, 42,[13 21],[14 22],[-0.1 6 -0.3 3.5 -0.1 1.9],[0 0 1],[],2,cc);
    %             draw_fliermirko(flier,seg, 42,[13 21],[14 22],[rcm(1,i)-1 rcm(1,i)+1 rcm(2,i)-1 rcm(2,i)+1 -0.1 1.9],[0 0 1],[],2);
    %             for pri = 0:0.01:0.1
    %                 line('XData',[xprep+0.5*pri xprep+0.5*pri],'YData',[yprep+pri yprep+pri],'ZData',[0 zprep],'Color','m','LineWidth',4);  
    %             end
    %         plot3(rcm(1,i),rcm(2,i),rcm(3,i),'ko'); 
            axis([-0.1 8 -2 2 -0.1 1.9]);

            plot(putanja_s(1,:,1),putanja_s(2,:,1),'b');
            plot(putanja_s(1,:,2),putanja_s(2,:,2),'r');
            plot(putanja_s(1,:,3),putanja_s(2,:,3),'m');
            crtaj_zid();
    %         plot(putanja_s(1,:,4),putanja_s(2,:,4),'c');
    %         plot(putanja_s(1,:,5),putanja_s(2,:,5),'k');

    %         subplot(4,2,[5 7],'replace');
%             subplot(4,2,[5 7]);
%             hold on;
%             draw_fliermirko(flier,seg, 42,[13 21],[14 22],[rcm(1,i)-1 rcm(1,i)+1 -0.5 0.5 -0.1 1.9],[cos(zelj_p(1,i)) sin(zelj_p(1,i)) 0],[],2,cc);
%     %             draw_fliermirko(flier,seg, 42,[13 21],[14 22],[rcm(1,i)-1 rcm(1,i)+1 rcm(2,i)-1 rcm(2,i)+1 -0.1 1.9],[cos(zp(i)) sin(zp(i)) 0],[],2);
%     %             for pri = 0:0.01:0.1
%     %                 line('XData',[xprep+0.5*pri xprep+0.5*pri],'YData',[yprep+pri yprep+pri],'ZData',[0 zprep],'Color','m','LineWidth',4);  
%     %             end
%     %         plot3(rcm(1,i),rcm(2,i),rcm(3,i),'ko');    
% 
%              axis([rcm(1,i)-1 rcm(1,i)+1 rcm(2,i)-1 rcm(2,i)+1 -0.1 1.9]);
%             plot(putanja_s(1,:,1),putanja_s(2,:,1),'b');
%             plot(putanja_s(1,:,2),putanja_s(2,:,2),'r');
%             plot(putanja_s(1,:,3),putanja_s(2,:,3),'m');
%             crtaj_zid();
        end
%         plot(putanja_s(1,:,4),putanja_s(2,:,4),'c');
%         plot(putanja_s(1,:,5),putanja_s(2,:,5),'k');

        
        subplot(4,2,[5 6 7 8]);
%         subplot(4,2,[2 4],'replace');
        hold on;          
        tac=1:12;        
        crtaj=1;
        if(mod(i,1000)==1)
            if size(XC,2)>2
                try
                    tac=convhull(XC,YC);
                catch
                    crtaj=0;
                end
                if crtaj   
                [XCm YCm] = polyoffset(XC(tac),YC(tac),0.01);
    %             plot(XC(tac),YC(tac),'k-');  
    %             plot(XCm,YCm,'r-');                 
                end
            end;
            if size(XCL,2)>2 && crtaj
                tacl=convhull(XCL,YCL);      
                [XCLm YCLm] = polyoffset(XCL(tacl),YCL(tacl),0.01);
                plot(XCL(tacl),YCL(tacl),'k-');  
    %             plot(XCLm,YCLm,'g-'); 

            end;
            if size(XCR,2)>2 && crtaj
                tacr=convhull(XCR,YCR);
                [XCRm YCRm] = polyoffset(XCR(tacr),YCR(tacr),0.01);
                plot(XCR(tacr),YCR(tacr),'k-');  
    %             plot(XCRm,YCRm,'g-'); 
            end;
            if size(XCPL,2)>2 && crtaj
                tacpl=convhull(XCPL,YCPL);      
                [XCPLm YCPLm] = polyoffset(XCPL(tacpl),YCPL(tacpl),0.01);
                plot(XCPL(tacpl),YCPL(tacpl),'k-');  
    %             plot(XCPLm,YCPLm,'g-');                  
            end
            if size(XCPR,2)>2 && crtaj
                tacpr=convhull(XCPR,YCPR);
                [XCPRm YCPRm] = polyoffset(XCPR(tacpr),YCPR(tacpr),0.01);
                plot(XCPR(tacpr),YCPR(tacpr),'k-');  
    %             plot(XCPRm,YCPRm,'g-');     
            end
        end
%         plot(rpetal(1),rpetal(2),'g*');      
%         plot(rpetad(1),rpetad(2),'b*');      
        plot(ZMPpos(1,i),ZMPpos(2,i),'rx');
        plot(rcm(1,i),rcm(2,i),'b*');  
        if(i==1)
            plot(putanja_s(1,:,1),putanja_s(2,:,1),'b');
            plot(putanja_s(1,:,2),putanja_s(2,:,2),'r');
            plot(putanja_s(1,:,3),putanja_s(2,:,3),'m');
            crtaj_zid();

%             plot(putanja_s(1,:,4),putanja_s(2,:,4),'c');
%             plot(putanja_s(1,:,5),putanja_s(2,:,5),'k');
        end
%         axis([rcm(1,i)-0.5 rcm(1,i)+0.5 rcm(2,i)-0.4 rcm(2,i)+0.4 -0.1 1.9]);
       
        axis equal;

%         subplot(4,2,[6],'replace');
%         plot(rcm(1:2,1:i)');
%         subplot(4,2,[8],'replace');        
%         plot(ZMPCMerr(1:2,1:i)');
        
        if(any(Output2~=0))
            Fa = getframe(gcf);
            mov = addframe(mov,Fa);
        end         
        
        pause(0.00001);
        
%         if(i==ikraj-1)
%             i=ikraj;
%         else            
%             i=i+incr;
%             if( i>ikraj)
%                 i=ikraj-1;
%             end
%         end
end
% 
% % save trainingsetall6_2.mat 'taus' 'dtaus' 'qs' 'qvs' 'qas' 'spds' 'skpoz' 'skpozt' 'skpozz' 'skbrzt' 'ZMPs' 'CMs' 'segnogad' 'skzzp' 'skzzv' 'skpozerr';
% i=ikraj;
%         km_geo(flier,qnew(:,i));
%         
%         rbaz = km_linkX(flier,13,'xyz');
%         rskzl = km_linkX(flier,14,'xyz');
%         rskzd = km_linkX(flier,15,'xyz');
%         rkukl = km_linkX(flier,16,'xyz');
%         rkukd = km_linkX(flier,17,'xyz');
%         rpetal = km_linkX(flier,18,'xyz');
%         rpetad = km_linkX(flier,19,'xyz');
%         rramel = km_linkX(flier,20,'xyz');
%         rramed = km_linkX(flier,21,'xyz');
%         rrukal = km_linkX(flier,22,'xyz');
%         rrukad = km_linkX(flier,23,'xyz');
%         rtrup = km_linkX(flier,24,'xyz');
%         rprstl = km_linkX(flier,25,'xyz');
%         rprstd = km_linkX(flier,26,'xyz');
% 
% 
%         XC=[];
%         YC=[];
%         XCL=[];
%         YCL=[];
%         XCR=[];
%         YCR=[];
%         XCPL=[];
%         YCPL=[];
%         XCPR=[];
%         YCPR=[];
% 
%         for ik=1:12
%             Rk=km_linkX(flier,ik,'xyz');
%             if Rk(3)<0
%                 XC=[XC Rk(1)];
%                 YC=[YC Rk(2)];
%                 if ik>0 && ik<7 
%                     XCL=[XCL Rk(1)];
%                     YCL=[YCL Rk(2)];
%                     if ik>2
%                         XCPL=[XCPL Rk(1)];
%                         YCPL=[YCPL Rk(2)];
%                     end            
%                 end 
%                 if ik>6 && ik<13
%                     XCR=[XCR Rk(1)];
%                     YCR=[YCR Rk(2)];
%                     if ik>8
%                         XCPR=[XCPR Rk(1)];
%                         YCPR=[YCPR Rk(2)];
%                     end
%                 end   
%             end
%         end
%         
%         
%         
%         %         flier=kk_kinedyn(flier,qvnew(:,i));
%         subplot(4,2,[1 3],'replace');
%         hold on;
%             draw_fliermirko(flier,seg, 42,[13 21],[14 22],[rcm(1,i)-1 rcm(1,i)+1 rcm(2,i)-1 rcm(2,i)+1 -0.1 1.9],[1 0 0],[]);
% %             for pri = 0:0.01:0.1
% %                 line('XData',[xprep+0.5*pri xprep+0.5*pri],'YData',[yprep+pri yprep+pri],'ZData',[0 zprep],'Color','m','LineWidth',4);  
% %             end
% %         plot3(rcm(1,i),rcm(2,i),rcm(3,i),'ko');    
%         subplot(4,2,[5 7],'replace');
% 
%         hold on;
%             draw_fliermirko(flier,seg, 42,[13 21],[14 22],[rcm(1,i)-1 rcm(1,i)+1 rcm(2,i)-1 rcm(2,i)+1 -0.1 1.9],[0 0 1],[]);
% %             for pri = 0:0.01:0.1
% %                 line('XData',[xprep+0.5*pri xprep+0.5*pri],'YData',[yprep+pri yprep+pri],'ZData',[0 zprep],'Color','m','LineWidth',4);  
% %             end
% %         plot3(rcm(1,i),rcm(2,i),rcm(3,i),'ko');    
%         subplot(4,2,[6 8],'replace');
%         hold on;
%             draw_fliermirko(flier,seg, 42,[13 21],[14 22],[rcm(1,i)-1 rcm(1,i)+1 rcm(2,i)-1 rcm(2,i)+1 -0.1 1.9],[0 -1 0],[]);
% %             for pri = 0:0.01:0.1
% %                 line('XData',[xprep+0.5*pri xprep+0.5*pri],'YData',[yprep+pri yprep+pri],'ZData',[0 zprep],'Color','m','LineWidth',4);  
% %             end
% %         plot3(rcm(1,i),rcm(2,i),rcm(3,i),'ko');    
% 
%         
%         subplot(4,2,[2 4]);
% %         subplot(4,2,[2 4],'replace');
%         hold on;          
%         tac=1:nknew;        
%         crtaj=1;
%         if size(XC,2)>2
%             try
%                 tac=convhull(XC,YC);
%             catch
%                 crtaj=0;
%             end
%             if crtaj   
%             [XCm YCm] = polyoffset(XC(tac),YC(tac),0.01);
% %             plot(XC(tac),YC(tac),'k-');  
% %             plot(XCm,YCm,'r-');                 
%             end
%         end;
%         if size(XCL,2)>2 && crtaj
%             tacl=convhull(XCL,YCL);      
%             [XCLm YCLm] = polyoffset(XCL(tacl),YCL(tacl),0.01);
%             plot(XCL(tacl),YCL(tacl),'k-');  
% %             plot(XCLm,YCLm,'g-'); 
%             
%         end;
%         if size(XCR,2)>2 && crtaj
%             tacr=convhull(XCR,YCR);
%             [XCRm YCRm] = polyoffset(XCR(tacr),YCR(tacr),0.01);
%             plot(XCR(tacr),YCR(tacr),'k-');  
% %             plot(XCRm,YCRm,'g-'); 
%         end;
%         if size(XCPL,2)>2 && crtaj
%             tacpl=convhull(XCPL,YCPL);      
%             [XCPLm YCPLm] = polyoffset(XCPL(tacpl),YCPL(tacpl),0.01);
%             plot(XCPL(tacpl),YCPL(tacpl),'k-');  
% %             plot(XCPLm,YCPLm,'g-');                  
%         end
%         if size(XCPR,2)>2 && crtaj
%             tacpr=convhull(XCPR,YCPR);
%             [XCPRm YCPRm] = polyoffset(XCPR(tacpr),YCPR(tacpr),0.01);
%             plot(XCPR(tacpr),YCPR(tacpr),'k-');  
% %             plot(XCPRm,YCPRm,'g-');     
%         end
% %         plot(rpetal(1),rpetal(2),'g*');      
% %         plot(rpetad(1),rpetad(2),'b*');      
%         plot(ZMPpos(1,i),ZMPpos(2,i),'rx');
%         plot(rcm(1,i),rcm(2,i),'b*');  
%         
%         axis([rcm(1,i)-0.5 rcm(1,i)+0.5 -0.2 +0.2 -0.1 1.9]);
%         axis equal;
%         
% %         subplot(4,2,[6],'replace');
% %         plot(rcm(1:2,1:i)');
% %         subplot(4,2,[8],'replace');        
% %         plot(ZMPCMerr(1:2,1:i)');
%         
        if(any(Output2~=0))
            Fa = getframe(gcf);
            mov = addframe(mov,Fa);
        end         
        
        pause(0.00001);
        
if(Output2~=0)
    mov = close(mov);   
end
