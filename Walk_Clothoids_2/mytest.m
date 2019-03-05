%=============================================================================%
%                                                                             %
%  Autors: Enrico Bertolazzi and Marco Frego                                  %
%          Department of Industrial Engineering                               %
%          University of Trento                                               %
%          enrico.bertolazzi@unitn.it                                         %
%          m.fregox@gmail.com                                                 %
%                                                                             %
%=============================================================================%
% Driver test program to check clothoid computation                           %
%=============================================================================%

function [zp,duz_k,putanja] = mytest(flier, rbaz,cilj)
    global dt t;

    %t=linspace(0,100,1000);
    %dt=1e-3*1;%perioda odabiranja
    %t = 0:dt:15;
    brit=size(t,2);%broj iteracija

    tol  = 1e-12 ;
    npts = 10000;%brit;

    % Robot base position

    %x_robot=zeros(size(t));
    %y_robot=zeros(size(t));

    x_robot=rbaz(1);
    y_robot=rbaz(2);
    A=flier.A;
    rbazor = TrotXYZ(A(:,:,6));
                
    theta_robot=rbazor(3);

    % initial point with angle direction
    x0     = x_robot ;
    y0     = y_robot   ;
    theta0 = theta_robot  ;

    % final point with angle direction
    x1     = cilj(1) ;
    y1     = cilj(2) ;
    theta1 = cilj(3)   ;

    % fprintf('Testing G1 Clothoid interpolation\n') ;
    % fprintf('initial point (%g,%g) initial angle = %g\n', x0, y0, theta0) ;
    % fprintf('final point (%g,%g) final angle = %g\n', x1, y1, theta1) ;
    % fprintf('tolerance = %g\n', tol ) ;

    % compute clothoid parameters
%     [k,dk,L] = buildClothoid( x0, y0, theta0, x1, y1, theta1 ) ;%, tol

%     fprintf('Computed parameters: k = %g, k'' = %g, L = %g\n', k, dk, L ) ;

    % compute points on clothoid
%     XY = pointsOnClothoid( x0, y0, theta0, k, dk, L, npts ) ;

    % compute minimal distance from robot base to sampled points on
    % curve (npts)

%     Kp=0.2;
%     brz_k=1;
    duz_k=1;
    zp = rbazor(3);


%     while (j<(brit-1)&&(min_distance_index<400))%size(t,2)

%         j=j+1

%         if (j==500)
                % initial point with angle direction
%             x0     = x_robot ;
%             y0     = y_robot   ;
%             theta0 = ugao  ;

            % final point with angle direction
%             x1     = 30 ;
%             y1     = 110 ;
%             theta1 = pi   ;

            % fprintf('Testing G1 Clothoid interpolation\n') ;
            % fprintf('initial point (%g,%g) initial angle = %g\n', x0, y0, theta0) ;
            % fprintf('final point (%g,%g) final angle = %g\n', x1, y1, theta1) ;
            % fprintf('tolerance = %g\n', tol ) ;

            % compute clothoid parameters
            [k,dk,L] = buildClothoid( x0, y0, theta0, x1, y1, theta1 ) ;%, tol

%             fprintf('Computed parameters: k = %g, k'' = %g, L = %g\n', k, dk, L ) ;

            % compute points on clothoid
            XY = pointsOnClothoid( x0, y0, theta0, k, dk, L, npts ) ;
%         end

%         if (j==1000)
            % initial point with angle direction
%             x0     = x_robot ;
%             y0     = y_robot   ;
%             theta0 = ugao  ;

            % final point with angle direction
%             x1     = 30 ;
%             y1     = 80 ;
%             theta1 = -pi/2   ;

            % fprintf('Testing G1 Clothoid interpolation\n') ;
            % fprintf('initial point (%g,%g) initial angle = %g\n', x0, y0, theta0) ;
            % fprintf('final point (%g,%g) final angle = %g\n', x1, y1, theta1) ;
            % fprintf('tolerance = %g\n', tol ) ;

            % compute clothoid parameters
%             [k,dk,L] = buildClothoid( x0, y0, theta0, x1, y1, theta1 ) ;%, tol

%             fprintf('Computed parameters: k = %g, k'' = %g, L = %g\n', k, dk, L ) ;

            % compute points on clothoid
%             XY = pointsOnClothoid( x0, y0, theta0, k, dk, L, npts ) ;
%         end

%         min_distance=sqrt((XY(1,1)-x_robot)^2+(XY(2,1)-y_robot)^2);
%         min_distance_index=1;
% 
%         for i=2:npts+1
%            distance=sqrt((XY(1,i)-x_robot)^2+( XY(2,i)-y_robot)^2);
%            l=[x_robot-XY(1,i),y_robot-XY(2,i)]';
%            vec_distance(:,i)=[x_robot-XY(1,i),y_robot-XY(2,i)]';
%            poz_ort(:,i) = vec_distance(:,i)/VecMod2(vec_distance(:,i));
%            if (distance<min_distance)
%                min_distance=distance;
%                min_distance_index=i;
%            end
%         end
% 
%         ugao=XY(3,min_distance_index);
% 
%         vektor=[x_robot-XY(1,min_distance_index) y_robot-XY(2,min_distance_index) 0];
%         tangenta=[cos(ugao) sin(ugao) 1]';
%         vek_proizvod=cross(vektor,tangenta);
%         smer=sign(vek_proizvod(1,3));
% 
%         if (smer>0)
%             Kp=0.1;
%         else
%             Kp=-0.1;
%         end
% 
%         zp=ugao-Kp*(0-atan2(min_distance,duz_k));
        putanja = XY;
        %rad2deg(-Kp*(0-atan2(min_distance,duzina_koraka)))
%         y_robot=y_robot+brzina*sin(skretanje)*dt;
%         x_robot=x_robot+brzina*cos(skretanje)*dt;
% 
%         if(mod(j,67)==0)
% 
%             %min_distance
%             %smer
% 
%             figure(1)
%             plot( XY(1,:), XY(2,:), '-r' );
%             hold on
%             plot(x_robot, y_robot, 'r.', 'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',20)
%             line([x_robot XY(1,min_distance_index)],[y_robot XY(2,min_distance_index)],[0 0],'Marker','.','LineStyle','-')
%             line([XY(1,min_distance_index) XY(1,min_distance_index)+15*cos(ugao) ],[XY(2,min_distance_index) XY(2,min_distance_index)+15*sin(ugao)],[0 0],'Marker','.','LineStyle','-')
% 
%             grid on ;
%             axis([-20,100,0,140])
%             %axis equal
%             hold off
% 
% %             pause(0.5);
%         end
%     end
end
