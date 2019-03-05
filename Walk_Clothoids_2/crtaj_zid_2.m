function crtaj_zid(sto_poz)


if(nargin==0)
    sto_poz=[0 0 0];
end
% 
    moguce_boje={['r']; ['b']; ['g']};
    moguci_oblici = {['prizma_4s'];['cilindar'];['prizma_3s']};
% 
%     zid1.boja = moguce_boje{1};
%     zid1.oblik = moguci_oblici{1};
%     zid1.dim(1) = 0.2;
%     zid1.dim(2) = 3.5;
%     zid1.dim(3) = 2.7;
%     zid1.poz(1)=2.1;%0.01*randi([2 5]);
%     zid1.poz(2)=1.25;%0.01*randi([-5 5]);
%     zid1.poz(3)=0;%0.01*randi([0 10]);
%     zid1.or(1)=0;%randi([-180 180])*pi/180;
%     zid1.or(2)=0;%randi([-180 180])*pi/180;
%     zid1.or(3)=0;%randi([-180 180])*pi/180;
%     
% 
%     zid2.boja = moguce_boje{2};
%     zid2.oblik = moguci_oblici{1};
%     zid2.dim(1) = 0.2;
%     zid2.dim(2) = 4.1;
%     zid2.dim(3) = 2.7;
%     zid2.poz(1)=4.1;%0.01*randi([2 5]);
%     zid2.poz(2)=-2;%0.01*randi([-5 5]);
%     zid2.poz(3)=0;%0.01*randi([0 10]);
%     zid2.or(1)=0;%randi([-180 180])*pi/180;
%     zid2.or(2)=0;%randi([-180 180])*pi/180;
%     zid2.or(3)=0;%randi([-180 180])*pi/180;
% 
% 
%     zid3.boja = moguce_boje{3};
%     zid3.oblik = moguci_oblici{1};
%     zid3.dim(1) = 0.2;
%     zid3.dim(2) = 4.5;
%     zid3.dim(3) = 2.7;
%     zid3.poz(1)=6.1;%0.01*randi([2 5]);
%     zid3.poz(2)=1.75;%0.01*randi([-5 5]);
%     zid3.poz(3)=0;%0.01*randi([0 10]);
%     zid3.or(1)=0;%randi([-180 180])*pi/180;
%     zid3.or(2)=0;%randi([-180 180])*pi/180;
%     zid3.or(3)=0;%randi([-180 180])*pi/180;
% 
%     zid1.poz=zid1.poz+sto_poz;
%     zid2.poz=zid2.poz+sto_poz;
%     zid3.poz=zid3.poz+sto_poz;
%     
%     hold on;


    walls = [   2 -0.5 2.2 -0.5 2.2 3 2 3 2 -0.5 ;
                5 -3 5.2 -3 5.2 -1.5 5 -1.5 5 -3 ];

% %     walls = [   2 -0.5 2.2 -0.5 2.2 3 2 3 2 -0.5 ;
% %                 4 -3 4.2 -3 4.2 0.3 4 0.3 4 -3 ;
% %                 6 -1.5 6.2 -1.5 6.2 3 6 3 6 -1.5 ]; 


    zid1.boja = moguce_boje{2};
    zid1.oblik = moguci_oblici{1};
    zid1.dim(1) = 0.2;
    zid1.dim(2) = 3.5;
    zid1.dim(3) = 1.5;
    zid1.poz(1)=2.1;%0.01*randi([2 5]);
    zid1.poz(2)=3.5/2-0.5;%0.01*randi([-5 5]);
    zid1.poz(3)=0;%0.01*randi([0 10]);
    zid1.or(1)=0;%randi([-180 180])*pi/180;
    zid1.or(2)=0;%randi([-180 180])*pi/180;
    zid1.or(3)=0;%randi([-180 180])*pi/180;
            
    zid2.boja = moguce_boje{3};
    zid2.oblik = moguci_oblici{1};
    zid2.dim(1) = 0.2;
    zid2.dim(2) = 1.5;
    zid2.dim(3) = 1.5;
    zid2.poz(1)=5.1;%0.01*randi([2 5]);
    zid2.poz(2)=1.5/2-1.5-1.5;%0.01*randi([-5 5]);
    zid2.poz(3)=0;%0.01*randi([0 10]);
    zid2.or(1)=0;%randi([-180 180])*pi/180;
    zid2.or(2)=0;%randi([-180 180])*pi/180;
    zid2.or(3)=0;%randi([-180 180])*pi/180;
            
    wi=0;
    hold on;
    while(wi<size(walls,1))
        wi=wi+1;
        fill(walls(wi,1:2:end),walls(wi,2:2:end),moguce_boje{mod(wi,3)+1});
    end
     crtaj_predmet(zid1);
     crtaj_predmet(zid2);

end