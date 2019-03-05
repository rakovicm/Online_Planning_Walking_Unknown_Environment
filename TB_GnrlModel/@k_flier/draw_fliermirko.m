function draw_flier(a, ds, varargin)
% Draws the flier object i 3D (stick diagram)
% Usage : draw_flier(flier_object, seg_display_scheme, head segment, [lfoot rfoot] )
% seg_display_scheme : [ [seg_nr chain]; [seg_nr chain]; ... [seg_nr chain] ]
%
sc = 0.035;
[X, Y, Z] = GenSphr(sc,5);
%newplot;
%clf;
hold on;

LineW=1;
Col = 'bgr';

[Xm,Ym] = meshgrid(-5:.2:5);
Zm=0*Xm;
mesh(Xm,Ym,Zm,0*Zm+0);

axis equal;

switch nargin
	case 2
		view(51,41);
		axis([-2 2 -2 2 -1 3]);
	case 3
		dh = varargin{1};
	case 4
		dh = varargin{1};
		df = varargin{2};
		view(51,41);
		axis([-2 2 -2 2 -1 3]);
    case 5
		dh = varargin{1};
		df = varargin{2};
		dtoe = varargin{3};
		view(51,41);
		axis([-2 2 -2 2 -1 3]);
    case 6
		dh = varargin{1};
		df = varargin{2};
		dtoe = varargin{3};
		view(51,41);
		axis(varargin{4});
    case 7
		dh = varargin{1};
		df = varargin{2};
		dtoe = varargin{3};
		view(varargin{5});
		axis(varargin{4});
    case 8
        dh = varargin{1};
		df = varargin{2};
		dtoe = varargin{3};
		view(varargin{5});
		axis(varargin{4});
        konttackedisp = varargin{6};
    case 9
        dh = varargin{1};
		df = varargin{2};
		dtoe = varargin{3};
		view(varargin{5});
		axis(varargin{4});
        konttackedisp = varargin{6};
        LineW=varargin{7};
    case 10
        dh = varargin{1};
		df = varargin{2};
		dtoe = varargin{3};
		view(varargin{5});
		axis(varargin{4});
        konttackedisp = varargin{6};        
        LineW=varargin{7};
        Col = varargin{8};
	otherwise
		error('Invalid number of arguments!');
end



for i = konttackedisp  
    ln = a.Con.lnr(i);    
    rkt = a.rc(:,ln,1) + a.A(:,:,ln)*a.Con.cvec(:,i);
    orseg =a.A(:,:,ln); 
%     or=eye(3);
    xosa = rkt(1:3)+orseg*[.1 0 0]';
    yosa = rkt(1:3)+orseg*[0 .1 0]';
    zosa = rkt(1:3)+orseg*[0 0 .1]';    

    line('XData',[rkt(1) xosa(1)],'YData',[rkt(2) xosa(2)],...
		'ZData',[rkt(3) xosa(3)],'Color',Col(3),'LineWidth',LineW); 

    line('XData',[rkt(1) yosa(1)],'YData',[rkt(2) yosa(2)],...
		'ZData',[rkt(3) yosa(3)],'Color','g','LineWidth',LineW); 

    line('XData',[rkt(1) zosa(1)],'YData',[rkt(2) zosa(2)],...
		'ZData',[rkt(3) zosa(3)],'Color','b','LineWidth',LineW); 

    plot3(rkt(1),rkt(2),rkt(3),'rx');
end

for i = 1:size(ds,1)
	TV = a.rc(:,ds(i,1),1);
	DV = TV - a.gdr(:,ds(i,1));
	RV = TV - a.gpr(:,ds(i,1),ds(i,2));
% 	surf(X+TV(1), Y+TV(2), Z+TV(3));	% sphere in the COM
% 	surf(X+DV(1), Y+DV(2), Z+DV(3));	% sphere in the entrance joint
% 	surf(X+RV(1), Y+RV(2), Z+RV(3));	% sphere in the output joint
	line('XData',[TV(1) DV(1)],'YData',[TV(2) DV(2)],...
		'ZData',[TV(3) DV(3)],'Color','b','LineWidth',LineW);    
	line('XData',[TV(1) RV(1)],'YData',[TV(2) RV(2)],...
		'ZData',[TV(3) RV(3)],'Color','r','LineWidth',LineW);    
    if(i==3)
        RVskzr=RV;
    end
    if(i==8)
        RVskzl=RV;
    end
end

% - ADDITIONAL DECORATIONS: HEAD, FEET -
% --------------------------------------
% FIXED SIZE HEAD, FIXED RELATIVE POSITION TO THE TRUNK
if(exist('dh'))
	if(dh)
		[Xh, Yh, Zh] = GenSphr(0.1,10);
		HV = a.rc(:,dh,1) + a.A(:,:,dh)*[0 0 0.25]';
		mesh(Xh+HV(1), Yh+HV(2), Zh+HV(3),'EdgeColor','m');	% draw head        
	end
end
% FIXED SIZE FEET, DRAWN RELATIVE TO THE FOOT SEGMENTS
if(exist('df'))
	if(df(1))
        rpeta_1 = km_linkX(a,7,'xyz');
        rpeta_2 = km_linkX(a,8,'xyz');
        lpeta_1 = km_linkX(a,1,'xyz');
        lpeta_2 = km_linkX(a,2,'xyz');
        X = [0.0828 -0.1128; 0.0828 -0.1128; 0.0828 -0.1128; 0.0828 -0.1128; 0.0828 -0.1128;];
		Y = [-1 -1; 1 1; 1 1; -1 -1; -1 -1] .* 0.07;
		Z = [-1 -1; -1 -1; 1 1; 1 1; -1 -1] .* 0.000000001;
		HV = a.rc(:,df(1),1) + a.A(:,:,df(1))*[0 0 -0.0414]';
		for i=1:5
			for j=1:2
				V = a.A(:,:,df(1)) * [X(i,j); Y(i,j); Z(i,j)];
				Xf(i,j) = V(1);	Yf(i,j) = V(2); Zf(i,j) = V(3);
			end
        end
    	line('XData',[RVskzr(1) (rpeta_1(1)+rpeta_2(1))/2],'YData',[RVskzr(2) (rpeta_1(2)+rpeta_2(2))/2],...
        	'ZData',[RVskzr(3) (rpeta_1(3)+rpeta_2(3))/2],'Color','b','LineWidth',LineW);        
		mesh(Xf+HV(1), Yf+HV(2), Zf+HV(3),'EdgeColor','k');	% draw one foot
		HV = a.rc(:,df(2),1) + a.A(:,:,df(2))*[0 0 -0.0414]';
		for i=1:5
			for j=1:2
				V = a.A(:,:,df(2)) * [X(i,j); Y(i,j); Z(i,j)];
				Xf(i,j) = V(1);	Yf(i,j) = V(2); Zf(i,j) = V(3);
			end
		end
		mesh(Xf+HV(1), Yf+HV(2), Zf+HV(3),'EdgeColor','k');	% draw the other foot
    	line('XData',[RVskzl(1) (lpeta_1(1)+lpeta_2(1))/2],'YData',[RVskzl(2) (lpeta_1(2)+lpeta_2(2))/2],...
        	'ZData',[RVskzl(3) (lpeta_1(3)+lpeta_2(3))/2],'Color','b','LineWidth',LineW);              
	end
end

if(exist('dtoe'))
	if(dtoe(1))
		X = [1 -1; 1 -1; 1 -1; 1 -1; 1 -1] .* (0.03);
		Y = [-1 -1; 1 1; 1 1; -1 -1; -1 -1] .* 0.07;
		Z = [-1 -1; -1 -1; 1 1; 1 1; -1 -1] .* (0.000000001);
		HV = a.rc(:,dtoe(1),1) + a.A(:,:,dtoe(1))*[0.03  0 0]';
		for i=1:5
			for j=1:2
				V = a.A(:,:,dtoe(1)) * [X(i,j); Y(i,j); Z(i,j)];
				Xf(i,j) = V(1);	Yf(i,j) = V(2); Zf(i,j) = V(3);
			end
		end
		mesh(Xf+HV(1), Yf+HV(2), Zf+HV(3),'EdgeColor','k');	% draw one foot
		HV = a.rc(:,dtoe(2),1) + a.A(:,:,dtoe(2))*[0.03  0 0]';
		for i=1:5
			for j=1:2
				V = a.A(:,:,dtoe(2)) * [X(i,j); Y(i,j); Z(i,j)];
				Xf(i,j) = V(1);	Yf(i,j) = V(2); Zf(i,j) = V(3);
			end
		end
		mesh(Xf+HV(1), Yf+HV(2), Zf+HV(3),'EdgeColor','k');	% draw the other foot
	end
end

