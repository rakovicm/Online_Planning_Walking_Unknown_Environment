function izlaz = saturacija(ulaz,mm)

vece = ulaz>mm;
vm1 = ulaz<mm;
vm2 = ulaz>-mm;
vm=vm1.*vm2;
manje = ulaz<-mm;

if(sum(vece)>1 || sum(manje)>1)
    disp('sat');
end

vece = mm.*vece;
manje = -mm.*manje;
vm = vm.*ulaz;


izlaz = vece+manje+vm;
end

