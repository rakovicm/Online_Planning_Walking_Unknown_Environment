function draw_flierkrst(a, ds, varargin)
% Draws the flier object i 3D (stick diagram)
% Usage : draw_flier(flier_object, seg_display_scheme, head segment, [lfoot rfoot] )
% seg_display_scheme : [ [seg_nr chain]; [seg_nr chain]; ... [seg_nr chain] ]
%
sc = 0.035;
hold on;
[X, Y, Z] = GenSphr(sc,5);
% % % newplot;
% % % clf;
hold on;
lw=2;
rcol=[1 0 0];
bcol=[0 0 1];
switch nargin
	case 2
		view(51,41);
		axis([-2 2 -2 2 -1 3]);
	case 3
		dh = varargin{1};
	case 4
		dh = varargin{1};
		dtoe = varargin{2};
		view(51,41);
		axis([-2 2 -2 2 -1 3]);
	case 5
		dh = varargin{1};
		dtoe = varargin{2};
		view(51,41);
		axis(varargin{3});
	case 6
		dh = varargin{1};
		dtoe = varargin{2};
		view(varargin{4});
		axis(varargin{3});
    case 9
        dh = varargin{1};
		dtoe = varargin{2};
		view(varargin{4});
		axis(varargin{3});
        rcol=varargin{5};
        bcol=varargin{6};
        lw=varargin{7};
	otherwise
		error('Invalid number of arguments!');
end

axis equal;
for i = 1:size(ds,1)
    if ds(i,1)~=6
        TV = a.rc(:,ds(i,1),1);
        DV = TV - a.gdr(:,ds(i,1));
    	RV = TV - a.gpr(:,ds(i,1),ds(i,2));
%	surf(X+TV(1), Y+TV(2), Z+TV(3));	% sphere in the COM
%	surf(X+DV(1), Y+DV(2), Z+DV(3));	% sphere in the entrance joint
%	surf(X+RV(1), Y+RV(2), Z+RV(3));	% sphere in the output joint
        line('XData',[TV(1) DV(1)],'YData',[TV(2) DV(2)],...
        	'ZData',[TV(3) DV(3)],'Color',bcol,'LineWidth',lw);
    	line('XData',[TV(1) RV(1)],'YData',[TV(2) RV(2)],...
        	'ZData',[TV(3) RV(3)],'Color',rcol,'LineWidth',lw);
    else
        TV = a.rc(:,ds(i,1),1);
       	RV1 = TV - a.A(:,:,6)*[0;0;a.lpr(3,ds(i,1),ds(i,2))];
        RV2 = TV - a.A(:,:,6)*[0;a.lpr(2:3,ds(i,1),ds(i,2))];
        RV3 = TV - a.A(:,:,6)*a.lpr(:,ds(i,1),ds(i,2));
        line('XData',[TV(1) RV1(1)],'YData',[TV(2) RV1(2)],...
        	'ZData',[TV(3) RV1(3)],'Color','k','LineWidth',lw);
        line('XData',[RV1(1) RV2(1)],'YData',[RV1(2) RV2(2)],...
        	'ZData',[RV1(3) RV2(3)],'Color','k','LineWidth',lw);
                line('XData',[RV2(1) RV3(1)],'YData',[RV2(2) RV3(2)],...
        	'ZData',[RV2(3) RV3(3)],'Color','k','LineWidth',lw);
    end;
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
if(exist('dtoe'))
	dtoe.lfc=[dtoe.lfc dtoe.lfc(1)];
    dtoe.rfc=[dtoe.rfc dtoe.rfc(1)];
	ll=a.Con.lnr(dtoe.lfc);
    for i=1:length(ll)-1
        lc=[a.rc(:,ll(i))-a.gdr(:,ll(i))...
            a.rc(:,ll(i))+a.A(:,:,ll(i))*a.Con.cvec(:,dtoe.lfc(i))...
            a.rc(:,ll(i+1))+a.A(:,:,ll(i+1))*a.Con.cvec(:,dtoe.lfc(i+1))];
        fill3(lc(1,:),lc(2,:),lc(3,:), [0.25 0 0]');
    end;
    
	ll=a.Con.lnr(dtoe.lfc);
    lr=a.Con.lnr(dtoe.rfc);
    for i=1:length(lr)-1
        lc=[a.rc(:,lr(i))-a.gdr(:,lr(i))...
            a.rc(:,lr(i))+a.A(:,:,lr(i))*a.Con.cvec(:,dtoe.rfc(i))...
            a.rc(:,lr(i+1))+a.A(:,:,lr(i+1))*a.Con.cvec(:,dtoe.rfc(i+1))];
        fill3(lc(1,:),lc(2,:),lc(3,:), [0.25 0.25 0]');
    end;
end
