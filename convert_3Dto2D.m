%% convert 3D electrode positions into 2D layout

% configs
cfg.elec.pnt(:,1) = [ ALLEEG(1).chanlocs.X ];
cfg.elec.pnt(:,2) = [ ALLEEG(1).chanlocs.Y ];
cfg.elec.pnt(:,3) = [ ALLEEG(1).chanlocs.Z ];
fprintf('creating layout for %s system\n', limo_ft_senstype(cfg.elec));

cfg.elec.label = { ALLEEG(1).chanlocs.labels };
cfg.rotate = [];
cfg.style = '2d';
cfg.projection = 'polar';
cfg.layout = [];
cfg.grad = [];
cfg.gradfile = [];
cfg.elecfile = [];
cfg.output = [];
cfg.feedback = 'no';
cfg.montage = 'no';
cfg.image = [];
cfg.bw = 0;
cfg.channel = 'all';
cfg.skipscale = 'no';
cfg.skipcomnt = 'no';
skipscale = strcmp(cfg.skipscale, 'yes'); % in general a scale is desired
skipcomnt = strcmp(cfg.skipcomnt, 'yes'); % in general a comment desired

% apply rotation
if isempty(cfg.rotate)
    switch limo_ft_senstype(cfg.elec)
        case {'ctf151', 'ctf275', 'bti148', 'bti248', 'ctf151_planar', 'ctf275_planar', ...
                'bti148_planar', 'bti248_planar', 'electrode', 'ext1020'}
            rz = 90;
        otherwise
            rz = 0;
    end
end
cfg.elec.pnt = limo_ft_warp_apply(limo_ft_rotate([0 0 rz]), cfg.elec.pnt, 'homogenous');

% use helper function for 3D layout
[pnt, label] = limo_ft_channelposition(cfg.elec);

if ~strcmpi(cfg.style, '3d')
    prj = limo_ft_elproj(pnt, cfg.projection);
    d = limo_ft_dist(prj');
    d(find(eye(size(d)))) = inf;
    mindist = min(d(:));
    X = prj(:,1);
    Y = prj(:,2);
    Width  = ones(size(X)) * mindist * 0.8;
    Height = ones(size(X)) * mindist * 0.6;
    lay.pos    = [X Y];
    lay.width  = Width;
    lay.height = Height;
    lay.label  = label;
else
    lay.pos   = pnt;
    lay.label = label;
end

%% check whether outline and mask are available. if not, add default "circle with triangle"
% to resemble the head in case of "circle with triangle", the electrode positions should also be
% scaled
if ~strcmp(cfg.style, '3d')
    if ~isfield(lay, 'outline') || ~isfield(lay, 'mask')
        rmax  = 0.5;
        l     = 0:2*pi/100:2*pi;
        HeadX = cos(l).*rmax;
        HeadY = sin(l).*rmax;
        NoseX = [0.18*rmax 0 -0.18*rmax];
        NoseY = [rmax-.004 rmax*1.15 rmax-.004];
        EarX  = [.497 .510 .518 .5299 .5419 .54 .547 .532 .510 .489];
        EarY  = [.0555 .0775 .0783 .0746 .0555 -.0055 -.0932 -.1313 -.1384 -.1199];

        % Scale the electrode positions to fit within a unit circle, i.e. electrode radius = 0.45
        ind_scale = strmatch('SCALE', lay.label);
        ind_comnt = strmatch('COMNT', lay.label);
        sel = setdiff(1:length(lay.label), [ind_scale ind_comnt]); % these are excluded for scaling
        x = lay.pos(sel,1);
        y = lay.pos(sel,2);
        xrange = limo_ft_range(x);
        yrange = limo_ft_range(y);

        % First scale the width and height of the box for multiplotting
        lay.width  = lay.width./xrange;
        lay.height = lay.height./yrange;

        % Then shift and scale the electrode positions
        lay.pos(:,1) = 0.9*((lay.pos(:,1)-min(x))/xrange-0.5);
        lay.pos(:,2) = 0.9*((lay.pos(:,2)-min(y))/yrange-0.5);

        % Define the outline of the head, ears and nose
        lay.outline{1} = [HeadX(:) HeadY(:)];
        lay.outline{2} = [NoseX(:) NoseY(:)];
        lay.outline{3} = [ EarX(:)  EarY(:)];
        lay.outline{4} = [-EarX(:)  EarY(:)];

        % Define the anatomical mask based on a circular head
        lay.mask{1} = [HeadX(:) HeadY(:)];
    end
end

%% add axes positions for comments and scale information if required
% add a placeholder for the comment in the upper left corner
lay.label{end+1}  = 'COMNT';
lay.width(end+1)  = mean(lay.width);
lay.height(end+1) = mean(lay.height);
X                 = min(lay.pos(:,1));
Y                 = max(lay.pos(:,2));
Y                 = min(lay.pos(:,2));
lay.pos(end+1,:)  = [X Y];

% add a placeholder for the scale in the upper right corner
lay.label{end+1}  = 'SCALE';
lay.width(end+1)  = mean(lay.width);
lay.height(end+1) = mean(lay.height);
X                 = max(lay.pos(:,1));
Y                 = max(lay.pos(:,2));
Y                 = min(lay.pos(:,2));
lay.pos(end+1,:)  = [X Y];

%% plot the layout for debugging (ONLY FOR 2D)
% tmpcfg = [];
% tmpcfg.layout = lay;
% ft_layoutplot(tmpcfg);

%% to write the layout to a text file, you can use this code snippet
% fprintf('writing layout to ''%s''\n', cfg.output);
% fid = fopen(cfg.output, 'wt');
% for i=1:numel(lay.label)
%     fprintf(fid, '%d %f %f %f %f %s\n', i, lay.pos(i,1), lay.pos(i,2), lay.width(i), lay.height(i), lay.label{i});
% end
% fclose(fid);
