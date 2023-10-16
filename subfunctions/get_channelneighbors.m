% Get channel neighors matrix from EEG channel locations
% using 2D polar (angular) projection and then 2D Delaunay triangulation.
% 
% inputs:
%       chanlocs - structure with EEG channel XYZ coordinates
%       vis - plot visualization (1) or not (0)
% 
% Cedric Cannard, Sep 2022

function [neighbors, channeighbstructmat] = get_channelneighbors(chanlocs,vis)

tmppath = fileparts(which('compute_mcc.m'));
addpath(fullfile(tmppath, 'subfunctions'))

cfg.elec.elecpos(:,1) = [ chanlocs.X ];
cfg.elec.elecpos(:,2) = [ chanlocs.Y ];
cfg.elec.elecpos(:,3) = [ chanlocs.Z ];
cfg.elec.label = { chanlocs.labels };
cfg.label = { chanlocs.labels };
cfg.method = 'distance';

% find channel neighbors
data = cfg;
cfg  = rmfield(cfg, 'label');     % first input must not be data
data = rmfield(data, 'method');   % second input must not be method
[neighbours, cfg] = ft_prepare_neighbours(cfg);

% neighbors = ft_prepare_neighbours(cfg, data); %%% TO EXTRACT
% cfg.feedback = 'no';
cfg.channel  = 'all';
cfg.compress = 'yes';
% cfg.parcellation = 'parcellation';

% check if the input data is valid for this function
% data = ft_checkdata(data);

% set the default for senstype depending on the data
cfg.senstype = 'eeg';

% get 3D positions from the sensor description
% addpath('C:\Users\Tracy\Documents\MATLAB\eeglab\plugins\Fieldtrip-lite20231015\external\artinis\private')
sens = ft_fetch_sens(cfg, data);
chanpos = sens.chanpos;
label   = sens.label;

% remove channels that are not in data
[dum, sensidx] = match_str(data.label, label);
chanpos = chanpos(sensidx, :);
label   = label(sensidx);

% select the desired channels
desired = ft_channelselection(cfg.channel, label);
[sensidx] = match_str(label, desired);
chanpos = chanpos(sensidx, :);
label   = label(sensidx);

% Project sensor positions on 2D plane if not already the case
if size(chanpos, 2) == 2 || all( chanpos(:,3) == 0 )
    proj = chanpos(:,1:2); % already on a 2D plane
else
    % project sensor on a 2D plane (from function elproj)
    x = chanpos(:,1);
    y = chanpos(:,2);
    if size(chanpos, 2) == 3
        z = chanpos(:,3);
    end

    % use default polar (angular) projection
    [az, el, r] = cart2sph(x, y, z);
    [x, y] = pol2cart(az, pi/2 - el);
    proj = [x, y];
end

% 2D Delaunay triangulation of the projected points
tri = delaunay(proj(:,1), proj(:,2));
if strcmp(cfg.compress, 'yes')
    tri_x = delaunay(proj(:,1)./2, proj(:,2)); % compress in the x-direction
    tri_y = delaunay(proj(:,1), proj(:,2)./2); % compress in the y-direction
    tri = [tri; tri_x; tri_y];
end

% Compute the neighbourhood geometry from the gradiometer/electrode positions
% neighbors = compneighbstructfromtri(chanpos, label, tri);

% mark neighbors according to triangulation
nchan = length(label);
channeighbstructmat = zeros(nchan,nchan);
for i = 1:size(tri, 1)
    channeighbstructmat(tri(i, 1), tri(i, 2)) = 1;
    channeighbstructmat(tri(i, 1), tri(i, 3)) = 1;
    channeighbstructmat(tri(i, 2), tri(i, 1)) = 1;
    channeighbstructmat(tri(i, 3), tri(i, 1)) = 1;
    channeighbstructmat(tri(i, 2), tri(i, 3)) = 1;
    channeighbstructmat(tri(i, 3), tri(i, 2)) = 1;
end

% construct a structured cell-array with all neighbors
neighbors = struct;
alldist = [];
for i = 1:nchan
    neighbors(i).label          = label{i};
    neighbidx                   = find(channeighbstructmat(i,:));
    neighbors(i).dist           = sqrt(sum((repmat(chanpos(i, :), numel(neighbidx), 1) - chanpos(neighbidx, :)).^2, 2));
    alldist                     = [alldist; neighbors(i).dist];
    neighbors(i).neighblabel    = label(neighbidx);
end

% remove neighbouring channels that are too far away (IMPORTANT if missing sensors)
neighbdist = mean(alldist)+3*std(alldist);
for i=1:nchan
    idx = neighbors(i).dist > neighbdist;
    neighbors(i).dist(idx)         = [];
    neighbors(i).neighblabel(idx)  = [];
end
neighbors = rmfield(neighbors, 'dist');

% Only select channels that are in the data
% if isfield(cfg, 'channel') && ~isempty(cfg.channel)
% %     desired = ft_channelselection(cfg.channel, data.label);
%     desired = data.label;
% end
% complete = struct;
% for i = 1:numel(desired)
% complete(i).label = desired{i};
% sel = find(strcmp({neighbors(:).label}, desired{i}));
% if numel(sel)==1
%   % take the set of neighbors from the definition
%   complete(i).neighblabel = neighbors(sel).neighblabel;
% else
%   % there are no neighbors defined for this channel
%   complete(i).neighblabel = {};
% end
% end
% neighbors = complete;

% Convert neighbors into row-arrays for a nicer code representation
% for i = 1:length(neighbors)
%   neighbors(i).neighblabel = neighbors(i).neighblabel(:)';
% end

% check that all chans have neighbors
k = 0;
for i = 1:length(neighbors)
    if isempty(neighbors(i).neighblabel)
        warning('no neighbors found for %s', neighbors(i).label);
    end
    k = k + length(neighbors(i).neighblabel);
end
if k==0
    warning('No neighbouring channels were specified or found');
else
    fprintf('there are on average %.1f neighbors per channel\n', k/length(neighbors));
end

% Visual feedback (or try convert_3dto2d for eeglab topo with nose and labels)
if vis
    tmpcfg = keepfields(cfg, {'layout', 'rows', 'columns', 'commentpos', 'skipcomnt', ...
        'scalepos', 'skipscale', 'projection', 'viewpoint', 'rotate', 'width', 'height', ...
        'elec', 'grad', 'opto', 'showcallinfo', 'trackcallinfo', 'trackconfig', ...
        'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo'});
    tmpcfg.neighbours = neighbors;
    tmpcfg.senstype = cfg.senstype;
    ft_neighbourplot(tmpcfg, data);
end


%% 
function [sens] = ft_fetch_sens(cfg, data)

% FT_FETCH_SENS mimics the behavior of FT_READ_SENS, but for a FieldTrip
% data structure or a FieldTrip configuration instead of a file on disk.
%
% Use as
%   [sens] = ft_fetch_sens(cfg)
% or as
%   [sens] = ft_fetch_sens(cfg, data)
%
% The sensor configuration can be passed into this function in four ways:
%  (1) in a configuration field
%  (2) in a file whose name is passed in a configuration field, see FT_READ_SENS
%  (3) in a layout file, see FT_PREPARE_LAYOUT
%  (4) in a data field
%
% The following fields are used from the configuration:
%   cfg.elec     = structure with electrode positions or filename, see FT_READ_SENS
%   cfg.grad     = structure with gradiometer definition or filename, see FT_READ_SENS
%   cfg.opto     = structure with optode definition or filename, see FT_READ_SENS
%   cfg.layout   = structure with layout definition or filename, see FT_PREPARE_LAYOUT
%   cfg.senstype = string, can be 'meg', 'eeg', or 'nirs', this is used to choose in combined data (default = 'eeg')
%
% When the sensors are not specified in the configuration, this function will
% fetch the grad, elec or opto field from the data.
%
% See also FT_READ_SENS, FT_DATATYPE_SENS, FT_FETCH_DATA, FT_PREPARE_LAYOUT

% Copyright (C) 2011-2016, Jorn M. Horschig
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without evft_neighbourploten the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

if nargin > 1 && ~isempty(data)
  data = ft_checkdata(data);
else
  data = struct; % initialize as empty struct
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed', {'elecfile', 'elec'});
cfg = ft_checkconfig(cfg, 'renamed', {'gradfile', 'grad'});
cfg = ft_checkconfig(cfg, 'renamed', {'optofile', 'opto'});

% set the defaults
cfg.senstype = ft_getopt(cfg, 'senstype');

% meg booleans
hasgradfile = isfield(cfg, 'grad') && ischar(cfg.grad);
hascfggrad  = isfield(cfg, 'grad') && isstruct(cfg.grad);
hasdatagrad = isfield(data, 'grad');

% eeg booleans
haselecfile = isfield(cfg, 'elec') && ischar(cfg.elec);
hascfgelec  = isfield(cfg, 'elec') && isstruct(cfg.elec);
hasdataelec = isfield(data, 'elec');

% nirs booleans
hasoptofile = isfield(cfg, 'opto') && ischar(cfg.opto);
hascfgopto  = isfield(cfg, 'opto') && isstruct(cfg.opto);
hasdataopto = isfield(data, 'opto');

% other
haslayout   = isfield(cfg, 'layout');
iscfgsens   = isfield(cfg, 'pnt')  || isfield(cfg, 'chanpos');
isdatasens  = isfield(data, 'pnt') || isfield(data, 'chanpos');

if isempty(cfg.senstype) && ((hasgradfile || hascfggrad || hasdatagrad) + (haselecfile || hascfgelec || hasdataelec) + (hasoptofile || hascfgopto || hasdataopto))>1
  ft_error('Cannot determine which sensors you want to work on. Specify cfg.senstype as ''meg'', ''eeg'' or ''nirs''');

elseif ~isempty(cfg.senstype)
  if iscell(cfg.senstype)
    % this represents combined EEG and MEG sensors, where each modality has its own sensor definition
    % use recursion to fetch all sensor descriptions
    sens = cell(size(cfg.senstype));
    for i=1:numel(cfg.senstype)
      tmpcfg = cfg;
      tmpcfg.senstype = cfg.senstype{i};
      sens{i} = ft_fetch_sens(tmpcfg, data);
    end
    return

  else
    switch lower(cfg.senstype)
      case 'meg'
        haselecfile = false;
        hascfgelec  = false;
        hasdataelec = false;
        hasoptofile = false;
        hascfgopto  = false;
        hasdataopto = false;
      case 'eeg'
        hasgradfile = false;
        hascfggrad  = false;
        hasdatagrad = false;
        hasoptofile = false;
        hascfgopto  = false;
        hasdataopto = false;
      case 'nirs'
        haselecfile = false;
        hascfgelec  = false;
        hasdataelec = false;
        hasgradfile = false;
        hascfggrad  = false;
        hasdatagrad = false;
      otherwise
        ft_error('unsupported specification of cfg.senstype as "%s"', cfg.senstype);
    end
  end
end

if (hasgradfile + hascfggrad + hasdatagrad + ...
    haselecfile + hascfgelec + hasdataelec + ...
    hasoptofile + hascfgopto + hasdataopto + ...
    haslayout + iscfgsens + isdatasens) > 1
  fprintf('Your data and configuration allow for multiple sensor definitions.\n');
  display = @warning;
else
  display = @fprintf;
end

if hasgradfile
  display('reading gradiometers from file ''%s''\n', cfg.grad);
  sens = ft_read_sens(cfg.grad, 'senstype', 'meg');
elseif hascfggrad
  display('using gradiometers specified in the configuration\n');
  sens = cfg.grad;
elseif hasdatagrad
  display('using gradiometers specified in the data\n');
  sens = data.grad;

elseif haselecfile
  display('reading electrodes from file ''%s''\n', cfg.elec);
  sens = ft_read_sens(cfg.elec, 'senstype', 'eeg');
  % only keep positions and labels in case of EEG electrodes
  sens = keepfields(sens, {'elecpos', 'chanpos', 'unit', 'coordsys', 'label','tra'});
elseif hascfgelec
  display('using electrodes specified in the configuration\n');
  sens = cfg.elec;
  % only keep positions and labels in case of EEG electrodes
  sens = keepfields(sens, {'elecpos', 'chanpos', 'unit', 'coordsys', 'label','tra'});
elseif hasdataelec
  display('using electrodes specified in the data\n');
  sens = data.elec;
  % only keep positions and labels in case of EEG electrodes
  sens = keepfields(sens, {'elecpos', 'chanpos', 'unit', 'coordsys', 'label','tra'});

elseif hasoptofile
  display('reading optodes from file ''%s''\n', cfg.opto);
  sens = ft_read_sens(cfg.opto, 'senstype', 'nirs');
  % only keep known fields in case of NIRS optodes
  sens = keepfields(sens, {'chanpos', 'label', 'optopos', 'optotype', 'optolabel', 'wavelength', 'tra', 'unit', 'coordsys'});
elseif hascfgopto
  display('using optodes specified in the configuration\n');
  sens = cfg.opto;
  % only keep known fields in case of NIRS optodes
  sens = keepfields(sens, {'chanpos', 'label', 'optopos', 'optotype', 'optolabel', 'wavelength', 'tra', 'unit', 'coordsys'});
elseif hasdataopto
  display('using optodes specified in the data\n');
  sens = data.opto;
  % only keep known fields in case of NIRS optodes
  sens = keepfields(sens, {'chanpos', 'label', 'optopos', 'optotype', 'optolabel', 'wavelength', 'tra', 'unit', 'coordsys'});

elseif haslayout
  display('Using the 2-D layout to determine the sensor position\n');
  lay = ft_prepare_layout(cfg);

  % remove the COMNT and SCALE labels
  sel = ~ismember(lay.label, {'COMNT' 'SCALE'});

  sens = [];
  sens.label = lay.label(sel);
  sens.chanpos = lay.pos(sel,:);
  sens.chanpos(:,3) = 0;

elseif iscfgsens
  % could be a sensor description
  display('The configuration input might already be a sensor description.\n');
  sens = cfg;

elseif isdatasens
  % could be a sensor description
  display('The data input might already be a sensor description.\n');
  sens = data;

else
  ft_error('no electrodes, gradiometers or optodes specified.');
end

% ensure that the sensor description is up-to-date
if (hasgradfile + hascfggrad + hasdatagrad + ...
    haselecfile + hascfgelec + hasdataelec + ...
    hasoptofile + hascfgopto + hasdataopto)
  % this should only be called if the sensor definition is a complete one, and not constructed from a layout
  % see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=3143#c9
  sens = ft_datatype_sens(sens);
end
