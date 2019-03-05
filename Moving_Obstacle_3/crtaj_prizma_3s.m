function []=crtaj_prizma_3s(predmet)

% funkcija za crtanje cevorostane prizme
% prosledjuje se struktura predmet sa poljima
% 'boja'
% 'oblik'
% 'dim'
% 'poz'
% 'or

try
    alpha=0.8;
    h=predmet.dim(1)*0.5*sqrt(3);
    T=[predmet.dim(1)/2; h/3; predmet.dim(2)/2];

    X = predmet.dim(1)*[0 0           0.5         0           0; 
           1 0           0.5         1           1; 
           1 0.5         1           0.5         0.5; 
           0 0.5         1           0           0];

    Y =   [0 0           h           0           0; 
           0 0           h           0           0; 
           0 h           0           h           h; 
           0 h           0           0           0];

    Z = predmet.dim(2)*[0 0           1           0           1; 
           0 1           0           0           1; 
           1 1           0           0           1; 
           1 0           1           0           1];

    matrica_rotacije=rotx(predmet.or(1))*roty(predmet.or(2))*rotz(predmet.or(3));

    for i=1:20
            temena(:,i)=[X(i); Y(i); Z(i)];
            temena_rotirano(:,i)=matrica_rotacije*temena(:,i);
    end

    for i=1:20
        X(i)=temena_rotirano(1,i)-matrica_rotacije(1,:)*T+predmet.poz(1);
        Y(i)=temena_rotirano(2,i)-matrica_rotacije(2,:)*T+predmet.poz(2);
        Z(i)=temena_rotirano(3,i)-matrica_rotacije(3,:)*T+predmet.poz(3);
    end

    fill3(X,Y,Z,predmet.boja,'FaceAlpha',alpha);    

    axis equal;
catch
    error('Prosledjeni oblik ne moze da se iscta');
    return;
end


end