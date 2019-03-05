function draw_flierW(a, ds, varargin)
% Draws the flier object i 3D (stick diagram)
% Usage : draw_flier(flier_object, seg_display_scheme [,extra] [ axes, view] )
% seg_display_scheme : [ [seg_nr chain]; [seg_nr chain]; ... [seg_nr chain] ]
%
newplot;
clf;
hold on;
switch nargin
	case 2
		view(51,41);
		axis([-2 2 -2 2 -1 3]);
		sr = 0.01;
	case 3
		view(51,41);
		axis([-2 2 -2 2 -1 3]);
		sr = varargin{1};
	case 4
		view(51,41);
		axis([-2 2 -2 2 -1 3]);
		sr = varargin{1};
		sc = varargin{2};
	case 6
		sr = varargin{1};
		sc = varargin{2};
		axis(varargin{3});
		view(varargin{4});
	otherwise
		error('Invalid number of arguments!');
end

[X, Y, Z] = GenSphr(sr,5);
axis equal;

for i = 1:size(ds,1)
	TV = a.rc(:,ds(i,1),1);
	DV = TV - a.gdr(:,ds(i,1));
	RV = TV - a.gpr(:,ds(i,1),ds(i,2));
%	surf(X+TV(1), Y+TV(2), Z+TV(3));	% sphere in the COM
	surf(X+DV(1), Y+DV(2), Z+DV(3));	% sphere in the entrance joint
%	surf(X+RV(1), Y+RV(2), Z+RV(3));	% sphere in the output joint
	line('XData',[TV(1) DV(1)],'YData',[TV(2) DV(2)],...
		'ZData',[TV(3) DV(3)],'Color','b','LineWidth',2);
	line('XData',[TV(1) RV(1)],'YData',[TV(2) RV(2)],...
		'ZData',[TV(3) RV(3)],'Color','r','LineWidth',2);
end

if exist('sc')
	for i=1:length(sc)
		TV = a.rc(:,sc(i).seg,1);
		switch sc(i).type
			case 'l'
				P(:,1) = a.A(:,:,sc(i).seg) * sc(i).matrix(:,1) + TV;
				P(:,2) = a.A(:,:,sc(i).seg) * sc(i).matrix(:,2) + TV;
				te = line(P(1,:),P(2,:),P(3,:));
				set(te, 'Color', sc(i).color);
			otherwise
				error('Illegal element type!');
		end
	end
end
