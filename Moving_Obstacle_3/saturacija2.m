function izlaz = saturacija(ulaz,min,max)

izlaz = ulaz;
if(izlaz>max)
    izlaz=max;
end
if(izlaz<min)
    izlaz=min;
end
    
end

