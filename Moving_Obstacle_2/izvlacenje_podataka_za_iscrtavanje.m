function izvlacenje_podataka_za_iscrtavanje(FileName,Output2)


    load(FileName);

    save([Output2 '.mat'],...
                        'qnew','rcm','dt','zp');
    
end