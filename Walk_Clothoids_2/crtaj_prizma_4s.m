function []=crtaj_prizma_4s(predmet)

% funkcija za crtanje cevorostane prizme
% prosledjuje se struktura predmet sa poljima
% 'boja'
% 'oblik'
% 'dim'
% 'poz'
% 'or'
try
    if(strcmp(predmet.oblik,'prizma_4s'))
        alpha=0.8;

        d=[predmet.dim(1)/2; predmet.dim(2)/2; predmet.dim(3)/2];

        X = predmet.dim(1)*[0 0 0 0 0 1; 
               1 0 1 1 1 1; 
               1 0 1 1 1 1; 
               0 0 0 0 0 1];

        Y = predmet.dim(2)*[0 0 0 0 1 0; 
               0 1 0 0 1 1; 
               0 1 1 1 1 1; 
               0 0 1 1 1 0];

        Z = predmet.dim(3)*[0 0 1 0 0 0; 
               0 0 1 0 0 0; 
               1 1 1 0 1 1; 
               1 1 1 0 1 1];

        matrica_rotacije=rotx(predmet.or(1))*roty(predmet.or(2))*rotz(predmet.or(3));

        for i=1:24
                temena(:,i)=[X(i); Y(i); Z(i)];
                temena_rotirano(:,i)=matrica_rotacije*temena(:,i);
        end

        D=matrica_rotacije*d;

        for i=1:24
            X(i)=temena_rotirano(1,i)-matrica_rotacije(1,:)*d+predmet.poz(1);
            Y(i)=temena_rotirano(2,i)-matrica_rotacije(2,:)*d+predmet.poz(2);
            Z(i)=temena_rotirano(3,i)+predmet.poz(3);%-matrica_rotacije(3,:)*d;
        end

        fill3(X,Y,Z,predmet.boja,'FaceAlpha',alpha);  
        
        axis equal;

    end
catch
    error('Prosledjeni oblik ne moze da se iscta');
    return;
end

end