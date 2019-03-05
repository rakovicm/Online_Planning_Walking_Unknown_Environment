function draw_fliermirkospec(a, ds,Clr, varargin)
% Draws the flier object i 3D (stick diagram)
% Usage : draw_flier(flier_object, seg_display_scheme, head segment, [lfoot rfoot] )
% seg_display_scheme : [ [seg_nr chain]; [seg_nr chain]; ... [seg_nr chain] ]
%
sc = 0.035;
[X, Y, Z] = GenSphr(sc,5);
%newplot;
%clf;
hold on;
switch nargin
	case 3
		view(51,41);
		axis([-2 2 -2 2 -1 3]);
	case 4
		dh = varargin{1};
	case 5
		dh = varargin{1};
		df = varargin{2};
		view(51,41);
		axis([-2 2 -2 2 -1 3]);
    case 6
		dh = varargin{1};
		df = varargin{2};
		dtoe = varargin{3};
		view(51,41);
		axis([-2 2 -2 2 -1 3]);
    case 7
		dh = varargin{1};
		df = varargin{2};
		dtoe = varargin{3};
		view(51,41);
		axis(varargin{4});
    case 8
		dh = varargin{1};
		df = varargin{2};
		dtoe = varargin{3};
		view(varargin{5});
		axis(varargin{4});
% % % %     case 8
% % % % 		dh = varargin{1};
% % % % 		df = varargin{2};
% % % % 		dtoe = varargin{3};
% % % % 		view(varargin{5});
% % % % 		axis(varargin{4});
% % % %         Clr=varargin(6);
	otherwise
		error('Invalid number of arguments!');
end

axis equal;

for i = 1:size(ds,1)
	TV = a.rc(:,ds(i,1),1);
	DV = TV - a.gdr(:,ds(i,1));
	RV = TV - a.gpr(:,ds(i,1),ds(i,2));
% 	surf(X+TV(1), Y+TV(2), Z+TV(3));	% sphere in the COM
% 	surf(X+DV(1), Y+DV(2), Z+DV(3));	% sphere in the entrance joint
% 	surf(X+RV(1), Y+RV(2), Z+RV(3));	% sphere in the output joint
	line('XData',[TV(1) DV(1)],'YData',[TV(2) DV(2)],...
		'ZData',[TV(3) DV(3)],'Color',Clr,'LineWidth',2);    
	line('XData',[TV(1) RV(1)],'YData',[TV(2) RV(2)],...
		'ZData',[TV(3) RV(3)],'Color',Clr,'LineWidth',2);    
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
		HV = a.rc(:,dh,1) + a.A(:,:,dh)*[0 0 0.20]';
		mesh(Xh+HV(1), Yh+HV(2), Zh+HV(3),'EdgeColor','m');	% draw head        
	end
end
% FIXED SIZE FEET, DRAWN RELATIVE TO THE FOOT SEGMENTS
if(exist('df'))
	if(df(1))
		X = [1 -1; 1 -1; 1 -1; 1 -1; 1 -1] .* 0.072;
		Y = [-1 -1; 1 1; 1 1; -1 -1; -1 -1] .* 0.05;
		Z = [-1 -1; -1 -1; 1 1; 1 1; -1 -1] .* 0.005;
		HV = a.rc(:,df(1),1) + a.A(:,:,df(1))*[0 0 -0.041]';
		for i=1:5
			for j=1:2
				V = a.A(:,:,df(1)) * [X(i,j); Y(i,j); Z(i,j)];
				Xf(i,j) = V(1);	Yf(i,j) = V(2); Zf(i,j) = V(3);
			end
        end
    	line('XData',[RVskzr(1) HV(1)+V(1)],'YData',[RVskzr(2) HV(2)],...
        	'ZData',[RVskzr(3) HV(3)+V(3)],'Color','b','LineWidth',2);        
		mesh(Xf+HV(1), Yf+HV(2), Zf+HV(3),'EdgeColor','b');	% draw one foot
		HV = a.rc(:,df(2),1) + a.A(:,:,df(2))*[0 0 -0.041]';
		for i=1:5
			for j=1:2
				V = a.A(:,:,df(2)) * [X(i,j); Y(i,j); Z(i,j)];
				Xf(i,j) = V(1);	Yf(i,j) = V(2); Zf(i,j) = V(3);
			end
		end
		mesh(Xf+HV(1), Yf+HV(2), Zf+HV(3),'EdgeColor','k');	% draw the other foot
    	line('XData',[RVskzl(1) HV(1)+V(1)],'YData',[RVskzl(2) HV(2)],...
        	'ZData',[RVskzl(3) HV(3)+V(3)],'Color','b','LineWidth',2);              
	end
end

if(exist('dtoe'))
	if(dtoe(1))
		X = [1 -1; 1 -1; 1 -1; 1 -1; 1 -1] .* 0.025;
		Y = [-1 -1; 1 1; 1 1; -1 -1; -1 -1] .* 0.05;
		Z = [-1 -1; -1 -1; 1 1; 1 1; -1 -1] .* 0.005;
		HV = a.rc(:,dtoe(1),1) + a.A(:,:,dtoe(1))*[0.025 0 0]';
		for i=1:5
			for j=1:2
				V = a.A(:,:,dtoe(1)) * [X(i,j); Y(i,j); Z(i,j)];
				Xf(i,j) = V(1);	Yf(i,j) = V(2); Zf(i,j) = V(3);
			end
		end
		mesh(Xf+HV(1), Yf+HV(2), Zf+HV(3),'EdgeColor','k');	% draw one foot
		HV = a.rc(:,dtoe(2),1) + a.A(:,:,dtoe(2))*[0.025 0 0]';
		for i=1:5
			for j=1:2
				V = a.A(:,:,dtoe(2)) * [X(i,j); Y(i,j); Z(i,j)];
				Xf(i,j) = V(1);	Yf(i,j) = V(2); Zf(i,j) = V(3);
			end
		end
		mesh(Xf+HV(1), Yf+HV(2), Zf+HV(3),'EdgeColor','k');	% draw the other foot
	end
end
