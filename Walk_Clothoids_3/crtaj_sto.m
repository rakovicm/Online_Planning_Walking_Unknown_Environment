function crtaj_sto(sto_poz)

if(nargin==0)
    sto_poz=[0.3 -0.3 0];
end

    moguce_boje={['r']; ['b']; ['g']};
    moguci_oblici = {['prizma_4s'];['cilindar'];['prizma_3s']};

    noga1.boja = moguce_boje{3};
    noga1.oblik = moguci_oblici{1};
    noga1.dim(1) = 0.15;
    noga1.dim(2) = 0.15;
    noga1.dim(3) = 0.7;
    noga1.poz(1)=0.05;%0.01*randi([2 5]);
    noga1.poz(2)=0.05;%0.01*randi([-5 5]);
    noga1.poz(3)=0;%0.01*randi([0 10]);
    noga1.or(1)=0;%randi([-180 180])*pi/180;
    noga1.or(2)=0;%randi([-180 180])*pi/180;
    noga1.or(3)=0;%randi([-180 180])*pi/180;
    

    noga2.boja = moguce_boje{3};
    noga2.oblik = moguci_oblici{1};
    noga2.dim(1) = 0.15;
    noga2.dim(2) = 0.15;
    noga2.dim(3) = 0.7;
    noga2.poz(1)=1.15;%0.01*randi([2 5]);
    noga2.poz(2)=0.05;%0.01*randi([-5 5]);
    noga2.poz(3)=0;%0.01*randi([0 10]);
    noga2.or(1)=0;%randi([-180 180])*pi/180;
    noga2.or(2)=0;%randi([-180 180])*pi/180;
    noga2.or(3)=0;%randi([-180 180])*pi/180;


    noga3.boja = moguce_boje{3};
    noga3.oblik = moguci_oblici{1};
    noga3.dim(1) = 0.15;
    noga3.dim(2) = 0.15;
    noga3.dim(3) = 0.7;
    noga3.poz(1)=0.05;%0.01*randi([2 5]);
    noga3.poz(2)=0.55;%0.01*randi([-5 5]);
    noga3.poz(3)=0;%0.01*randi([0 10]);
    noga3.or(1)=0;%randi([-180 180])*pi/180;
    noga3.or(2)=0;%randi([-180 180])*pi/180;
    noga3.or(3)=0;%randi([-180 180])*pi/180;


    noga4.boja = moguce_boje{3};
    noga4.oblik = moguci_oblici{1};
    noga4.dim(1) = 0.15;
    noga4.dim(2) = 0.15;
    noga4.dim(3) = 0.7;
    noga4.poz(1)=1.15;%0.01*randi([2 5]);
    noga4.poz(2)=0.55;%0.01*randi([-5 5]);
    noga4.poz(3)=0;%0.01*randi([0 10]);
    noga4.or(1)=0;%randi([-180 180])*pi/180;
    noga4.or(2)=0;%randi([-180 180])*pi/180;
    noga4.or(3)=0;%randi([-180 180])*pi/180;

    radnapovrs.boja = moguce_boje{3};
    radnapovrs.oblik = moguci_oblici{1};
    radnapovrs.dim(1) = 1.3;
    radnapovrs.dim(2) = 0.7;
    radnapovrs.dim(3) = 0.02;
    radnapovrs.poz(1)=0.6;%0.01*randi([2 5]);
    radnapovrs.poz(2)=0.3;%0.01*randi([-5 5]);
    radnapovrs.poz(3)=0.7;%0.01*randi([0 10]);
    radnapovrs.or(1)=0;%randi([-180 180])*pi/180;
    radnapovrs.or(2)=0;%randi([-180 180])*pi/180;
    radnapovrs.or(3)=0;%randi([-180 180])*pi/180;

    noga1.poz=noga1.poz+sto_poz;
    noga2.poz=noga2.poz+sto_poz;
    noga3.poz=noga3.poz+sto_poz;
    noga4.poz=noga4.poz+sto_poz;
    radnapovrs.poz=radnapovrs.poz+sto_poz;
    
    hold on;
    crtaj_predmet(noga1);
    crtaj_predmet(noga2);
    crtaj_predmet(noga3);
    crtaj_predmet(noga4);
    crtaj_predmet(radnapovrs);

end