function crtaj_predmet(predmet)

    switch(predmet.oblik)
        case 'cilindar'
            crtaj_cilindar(predmet);
        case 'prizma_3s'
            crtaj_prizma_3s(predmet);
        case 'prizma_4s'
            crtaj_prizma_4s(predmet);
    end
end