function iscrtajsim(fps,FileName,Output2)


global Up t Fuks ZMPtr CMtr
global XC YC XCL YCL XCR YCR XCPL YCPL XCPR YCPR TAC
global segdkuk segdkol segdskz segdprst seglkuk seglkol seglskz seglprst 
global segtrupx segtrupy segdrame segdlakat seglrame segllakat 
global seglnoga segdnoga segtrup seglruka segdruka 


if(nargin==0)
    fps=15;
    FileName= 'sim1nast1.mat';
    Output2=0;
end

if(nargin==1)
    Output2=0;
end

load(FileName)

%Inicijalizacija promenljivih
letac_motor_kontakti;

%   mrle	   |  Microsoft - Run Length Encoding(msrle32.dll)
%   msvc	   |  Microsoft - Video 1(msvidc32.dll)
%   i420	   |  Intel - Indeo 4 codec(iyuv_32.dll)
%   uyvy	   |  Microsoft - UYVY(msyuv.dll)
%   yuy2	   |  Microsoft - YUY2(msyuv.dll)
%   yvyu	   |  Microsoft - YVYU(msyuv.dll)
%   iyuv	   |  Intel - Indeo(iyuv_32.dll)
%   yvu9	   |  Toshiba Video Codec(tsbyuv.dll)
%   cvid	   |  Supermac - Cinepak(iccvid.dll)
%   XVID	   |  Xvid MPEG-4 Video Codec(xvidvfw.dll)


figure('Color',[1 1 1],'Position',+[200 200  1200 800]);
H=subplot(1,1,1);
set(H,'LineWidth',1,'FontSize',10,'FontName','Times New Roman')   


if(Output2~=0)
     fig=gcf;
     set(fig,'DoubleBuffer','on');
     set(gcf,'NextPlot','replace','Visible','off');
     mov = avifile(Output2,'COMPRESSION','FFDS');
end



col='bbb';
ikraj = 45000;
ipoc=1;
ikorak=99;%round(1/fps/dt)
for i = ipoc : ikorak : ikraj
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
        if(mod(i,ikorak*10)==1)
        %         flier=kk_kinedyn(flier,qvnew(:,i));
        subplot(3,1,[1 2]);
        if(col=='bbb')
            col='ggg';
            Lw=2.7;
        else
            col='bbb';
            Lw=2;
        end
        hold on;
            draw_fliermirko(flier,seg, 42,[13 21],[14 22],[-0.4 3+0.7 -0.5 +0.5 -0.01 1.9],[0 -1 0],[],Lw,col);
            %             for pri = 0:0.01:0.1
%                 line('XData',[xprep+0.5*pri xprep+0.5*pri],'YData',[yprep+pri yprep+pri],'ZData',[0 zprep],'Color','m','LineWidth',4);  
%             end
%         plot3(rcm(1,i),rcm(2,i),rcm(3,i),'ko');    
%         subplot(1,3,[2]);
% 
%         hold on;
%             draw_fliermirko(flier,seg, 42,[13 21],[14 22],[-0.4 1+0.6 -0.5 +0.5 -0.01 1.9],[0 -1 0],[],2,'rrr');
%         grid on;
%             for pri = 0:0.01:0.1
%                 line('XData',[xprep+0.5*pri xprep+0.5*pri],'YData',[yprep+pri yprep+pri],'ZData',[0 zprep],'Color','m','LineWidth',4);  
%             end
%         plot3(rcm(1,i),rcm(2,i),rcm(3,i),'ko');    
%         subplot(4,2,[6 8],'replace');
%         hold on;
%             draw_fliermirko(flier,seg, 42,[13 21],[14 22],[rcm(1,i)-1 rcm(1,i)+1 rcm(2,i)-1 rcm(2,i)+1 -0.1 1.9],[0 -1 0],[],2,'rrr');
% %             for pri = 0:0.01:0.1
% %                 line('XData',[xprep+0.5*pri xprep+0.5*pri],'YData',[yprep+pri yprep+pri],'ZData',[0 zprep],'Color','m','LineWidth',4);  
% %             end
% %         plot3(rcm(1,i),rcm(2,i),rcm(3,i),'ko');    


        subplot(3,1,[3]);
%         subplot(4,2,[2 4],'replace');
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
        plot(rcm(1,i),rcm(2,i),'go');  
        
        axis([-0.4 5+0.6 -0.3 +0.3 -0.01 1.9]);
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
end

% save trainingsetall6_2.mat 'taus' 'dtaus' 'qs' 'qvs' 'qas' 'spds' 'skpoz' 'skpozt' 'skpozz' 'skbrzt' 'ZMPs' 'CMs' 'segnogad' 'skzzp' 'skzzv' 'skpozerr';
i=ikraj;
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
        
        
        subplot(3,1,[1 2]);
        if(col=='bbb')
            col='ggg';
            Lw=2.7;
        else
            col='bbb';
            Lw=2;
        end
        hold on;
            draw_fliermirko(flier,seg, 42,[13 21],[14 22],[-0.4 3+0.7 -0.5 +0.5 -0.01 1.9],[0 -1 0],[],Lw,col);
%             for pri = 0:0.01:0.1
%                 line('XData',[xprep+0.5*pri xprep+0.5*pri],'YData',[yprep+pri yprep+pri],'ZData',[0 zprep],'Color','m','LineWidth',4);  
%             end
%         plot3(rcm(1,i),rcm(2,i),rcm(3,i),'ko');    
%         subplot(1,3,[2]);

%         hold on;
%             draw_fliermirko(flier,seg, 42,[13 21],[14 22],[-0.4 1+0.6 -0.5 +0.5 -0.01 1.9],[0 -1 0],[],2,'bbb');
% %             for pri = 0:0.01:0.1
% %                 line('XData',[xprep+0.5*pri xprep+0.5*pri],'YData',[yprep+pri yprep+pri],'ZData',[0 zprep],'Color','m','LineWidth',4);  
% %             end
% %         plot3(rcm(1,i),rcm(2,i),rcm(3,i),'ko');    
% %         subplot(4,2,[6 8]);
% %         hold on;
% %             draw_fliermirko(flier,seg, 42,[13 21],[14 22],[rcm(1,i)-1 rcm(1,i)+1 rcm(2,i)-1 rcm(2,i)+1 -0.1 1.9],[0 -1 0],[],2,'bbb');
% % %             for pri = 0:0.01:0.1
% % %                 line('XData',[xprep+0.5*pri xprep+0.5*pri],'YData',[yprep+pri yprep+pri],'ZData',[0 zprep],'Color','m','LineWidth',4);  
% % %             end
% % %         plot3(rcm(1,i),rcm(2,i),rcm(3,i),'ko');    
% 
%         
        subplot(3,1,[3]);
%         subplot(4,2,[2 4],'replace');
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
%         plot(rpetal(1),rpetal(2),'g*');      
%         plot(rpetad(1),rpetad(2),'b*');      
        plot(ZMPpos(1,i),ZMPpos(2,i),'rx');
        plot(rcm(1,i),rcm(2,i),'go');  
        
        axis([-0.4 5+0.6 -0.3 +0.3 -0.01 1.9]);
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
        
if(Output2~=0)
    mov = close(mov);   
end
