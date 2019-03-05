function [x y XY sece] = search_point_for_clatoid(c_pose, goal, walls)
%     
%     c_pose = [1.9 -0.75 -pi/8];
%     c_pose = [4.8 0.5 pi/8];
%     c_pose = [0 0 0];
%     goal = [8 0 0];
    x = 0;
    y = 0;
    npts = 2000;
    x_robot=c_pose(1);
    y_robot=c_pose(2);
    theta_robot=c_pose(3);
    
    x_e = goal(1);
    y_e = goal(2);
    theta_e = goal(3);
    
    sece = 0;
    
    if(VecMod([(goal(1:2)-c_pose(1:2)) 0])>5)
        dir = [goal(1:2)-c_pose(1:2) 0 ]./VecMod([goal(1:2)-c_pose(1:2) 0 ]);
        x_e = x_robot+5*dir(1);
        y_e = y_robot+5*dir(2);
        theta_e = atan2((y_e-y_robot),(x_e-x_robot));
    end
    
    
% % %     walls = [   2 -0.5 2.2 -0.5 2.2 3 2 3 2 -0.5 ;
% % %                 4 -3 4.2 -3 4.2 0.3 4 0.3 4 -3 ;
% % %                 6 -1.5 6.2 -1.5 6.2 3 6 3 6 -1.5 ]; 
    
    
    [k,dk,L] = buildClothoid( x_robot, y_robot, theta_robot, x_e, y_e, theta_e ) ;%, tol
    XY = pointsOnClothoid( x_robot, y_robot, theta_robot, k, dk, L, npts ) ;
    i=0;
    
    
    
    while(i<size(walls,1))
        i=i+1;
        d_min = p_poly_dist(walls(i,1:2:end),walls(i,2:2:end),XY(1,:),XY(2,:));
        [xi,yi] = polyxpoly(walls(i,1:2:end),walls(i,2:2:end),XY(1,:),XY(2,:));
        [x_e,i_min] = min(xi);
        y_e_1 = yi(i_min);
        y_e_2 = yi(i_min);
        wall_i = i;
        if(~isempty(x_e) || d_min(1)<0.02)%Proveriti ovaj uslov
            sece = 1;
            point_found = 0;
            XC = [c_pose(1) walls(wall_i,1:2:end)];
            YC = [c_pose(2) walls(wall_i,2:2:end)];
            tac1=convhull(XC,YC);
            XCw = XC(tac1(2:end-1));
            YCw = YC(tac1(2:end-1));
            tac2=convhull(XCw,YCw);
            min_dist = Inf;
            min_p1 = tac2(1);
            min_p2 = tac2(end-1);

            x_e_1 = XCw(min_p1);
            x_e_2 = XCw(min_p2);
            y_e_1 = YCw(min_p1);
            y_e_2 = YCw(min_p2);
            
            P1=[x_e_1;y_e_1;0;];
            P2=[x_e_2;y_e_2;0;];
            
            dir_e1 = [x_e_1-c_pose(1);y_e_1-c_pose(2);0];
            dir_e2 = [x_e_2-c_pose(1);y_e_2-c_pose(2);0];
            dir_e1 = dir_e1./VecMod(dir_e1);
            dir_e2 = dir_e2./VecMod(dir_e2);
            c=cross(dir_e1,dir_e2);
            if(c(3)>0)
                dir_e1 = rotz(-pi/2)*dir_e1;
                dir_e2 = rotz(pi/2)*dir_e2;                
            else
                dir_e1 = rotz(pi/2)*dir_e1;
                dir_e2 = rotz(-pi/2)*dir_e2;
            end
            
            P1_d = p_poly_dist(P1(1),P1(2),XY(1,:),XY(2,:));
            P2_d = p_poly_dist(P2(1),P2(2),XY(1,:),XY(2,:));
            
            while(~point_found)

                npts = 200;
                delta = 0.01;
                
                if(P1_d<P2_d)
                P1= P1 + delta*dir_e1;
                theta_e_1 = atan2((P1(2)-c_pose(2)),(P1(1)-c_pose(1)));
                else
                P2= P2 + delta*dir_e2;
                theta_e_2 = atan2((P2(2)-c_pose(2)),(P2(1)-c_pose(1)));
                end
                        
                
                
                if(P1_d<P2_d)
                    [k,dk,L] = buildClothoid( x_robot, y_robot, theta_robot, P1(1), P1(2), theta_e_1 ) ;%, tol
                    XY_1 = pointsOnClothoid( x_robot, y_robot, theta_robot, k, dk, L, npts );
                    [xi_1,yi_1] = polyxpoly(walls(wall_i,1:2:end),walls(wall_i,2:2:end),XY_1(1,:),XY_1(2,:));
%                     plot(XY_1(1,:),XY_1(2,:));
                    if(isempty(xi_1))                    
                        d_min_1 = p_poly_dist(walls(wall_i,1:2:end),walls(wall_i,2:2:end),XY_1(1,:),XY_1(2,:));
                        if(min(d_min_1)>0.03)
                            x = x_e;
                            y = y_e_1;
                            point_found = 1;
                        end
                    end
                else
                    [k,dk,L] = buildClothoid( x_robot, y_robot, theta_robot, P2(1), P2(2), theta_e_2 ) ;%, tol
                    XY_2 = pointsOnClothoid( x_robot, y_robot, theta_robot, k, dk, L, npts );
                    [xi_2,yi_2] = polyxpoly(walls(wall_i,1:2:end),walls(wall_i,2:2:end),XY_2(1,:),XY_2(2,:));
%                     plot(XY_2(1,:),XY_2(2,:));
                    if(isempty(xi_2))                    
                        d_min_2 = p_poly_dist(walls(wall_i,1:2:end),walls(wall_i,2:2:end),XY_2(1,:),XY_2(2,:));                    
                        if(min(d_min_2)>0.03)
                            x = x_e;
                            y = y_e_2;
                            point_found = 1;
                        end                    
                    end
                end
                %search on the wall 
                i = size(walls,1)+1;
            end
        end
    end
    npts = 10000;
    XY = pointsOnClothoid( x_robot, y_robot, theta_robot, k, dk, L, npts ) ;
    crtaj_zid([0 0 0],walls);
    hold on; 
    plot(XY(1,:),XY(2,:));
end
