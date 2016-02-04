% [cfg] = ft_SASICA(cfg,comp,data)
%
% This is a wrapper to eeg_SASICA. See help eeg_SASICA for more info.


function [cfg] = ft_SASICA(cfg,comp,data)

% create a minimal EEG structure to pass to eeg_SASICA.
EEG = eeg_emptyset;
EEG.setname = 'internal';
EEG.nbchan = numel(comp.topolabel);
EEG.trials = numel(comp.trial);
EEG.pnts = size(comp.trial{1},2);
EEG.srate = comp.fsample;
EEG.xmin = comp.time{1}(1);
EEG.xmax = comp.time{end}(end);
EEG.times = comp.time{1}*1000;
EEG.icaact = cat(3,comp.trial{:});
EEG.icawinv = comp.topo;
EEG.icachansind = 1:size(EEG.icaact,1);

% attempt to create a chanlocs
if exist('data','var') && isfield(data,'elec')
    EEG.chanlocs = struct('labels',data.elec.label(:)');
    [EEG.chanlocs.X] = rep2struct(data.elec.pnt(:,1));
    [EEG.chanlocs.Y] = rep2struct(data.elec.pnt(:,2));
    [EEG.chanlocs.Z] = rep2struct(data.elec.pnt(:,3));
else
    cfg.layout = ft_prepare_layout(cfg);
    EEG.chanlocs = struct('labels',cfg.layout.label(:)');
    [EEG.chanlocs.X] = rep2struct(cfg.layout.pos(:,1));
    [EEG.chanlocs.Y] = rep2struct(cfg.layout.pos(:,2));
    [EEG.chanlocs.Z] = rep2struct(zeros(size(EEG.chanlocs)));
end
EEG.chanlocs = convertlocs(EEG.chanlocs,'cart2all');
EEG.chaninfo = [];

if not(exist('data','var'))
    EEG.data = reshape(EEG.icawinv * EEG.icaact(:,:),EEG.nbchan,EEG.pnts,EEG.trials);
else
    EEG.data = cat(3,data.trial{:});
end

if not(cfg.opts.noplot)
    try
        evalin('base','EEG;');
        try
            rep = uigetpref('SASICA','overwriteEEG','Overwrite another EEG in memory?',...
                {'I need to create an EEG variable in your base workspace.' 'Currently, there is another one can I overwrite it?'},...
                'Yes|Cancel');
            if strcmpi(rep,'cancel')
                try
                    rmpref('SASICA','overwriteEEG')
                end
                error()
            end
        catch
            disp('Rename the EEG variable in your base workspace to continue')
            return
        end
    end
end


EEG = eeg_SASICA(EEG,cfg);

if ~cfg.opts.noplot
    assignin('base','EEG',EEG);
end

cfg.reject = EEG.reject;


function [varargout] = rep2struct(varargin)

% [s.target] = rep2struct(dat)
% replicate the value dat into each element of structure s in field target.
% if dat has same nb or elements as s, each element of dat goes into one
% element of s. if dat is more dimensional and doesn't have the same number
% of elements as s, and has same size along dimension 1, then pass each
% slice into s.target.

if numel(varargin) == 1
    dat = varargin{1};
    if numel(dat) == nargout
        for i = 1:nargout
            varargout{i} = dat(i);
        end
    elseif size(dat,1) == nargout
        for i = 1:nargout
            varargout{i} = dat(i,:);
        end
    else
        for i = 1:nargout
            varargout{i} = dat;
        end
    end
elseif numel(varargin) == nargout
    for i = 1:nargout
        varargout{i} = varargin{i};
    end
else
    error('Wrong number of arguments');
end

% convertlocs() - Convert electrode locations between coordinate systems
%                 using the EEG.chanlocs structure.
%
% Usage: >> newchans = convertlocs( EEG, 'command');
%
% Input:
%   chanlocs  - An EEGLAB EEG dataset OR a EEG.chanlocs channel locations structure
%   'command' - ['cart2topo'|'sph2topo'|'sphbesa2topo'| 'sph2cart'|'topo2cart'|'sphbesa2cart'|
%               'cart2sph'|'sphbesa2sph'|'topo2sph'| 'cart2sphbesa'|'sph2sphbesa'|'topo2sphbesa'|
%               'cart2all'|'sph2all'|'sphbesa2all'|'topo2all']
%                These command modes convert between four coordinate frames: 3-D Cartesian
%                (cart), Matlab spherical (sph), Besa spherical (sphbesa), and 2-D polar (topo)
%               'auto' -- Here, the function finds the most complex coordinate frame
%                 and constrains all the others to this one. It searches first for Cartesian
%                 coordinates, then for spherical and finally for polar. Default is 'auto'.
%
% Optional input
%   'verbose' - ['on'|'off'] default is 'off'.
%
% Outputs:
%   newchans - new EEGLAB channel locations structure
%
% Ex:  CHANSTRUCT = convertlocs( CHANSTRUCT, 'cart2topo');
%      % Convert Cartesian coordinates to 2-D polar (topographic).
%
% Author: Arnaud Delorme, CNL / Salk Institute, 22 Dec 2002
%
% See also: readlocs()

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 22 Dec 2002, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function chans = convertlocs(chans, command, varargin);

if nargin < 1
    help convertlocs;
    return;
end;

if nargin < 2
    command = 'auto';
end;
if nargin == 4 && strcmpi(varargin{2}, 'on')
    verbose = 1;
else
    verbose = 0; % off
end;

% test if value exists for default
% --------------------------------
if strcmp(command, 'auto')
    if isfield(chans, 'X') && ~isempty(chans(1).X)
        command = 'cart2all';
        if verbose
            disp('Make all coordinate frames uniform using Cartesian coords');
        end;
    else
        if isfield(chans, 'sph_theta') && ~isempty(chans(1).sph_theta)
            command = 'sph2all';
            if verbose
                disp('Make all coordinate frames uniform using spherical coords');
            end;
        else
            if isfield(chans, 'sph_theta_besa') && ~isempty(chans(1).sph_theta_besa)
                command = 'sphbesa2all';
                if verbose
                    disp('Make all coordinate frames uniform using BESA spherical coords');
                end;
            else
                command = 'topo2all';
                if verbose
                    disp('Make all coordinate frames uniform using polar coords');
                end;
            end;
        end;
    end;
end;

% convert
% -------
switch command
    case 'topo2sph',
        theta  = {chans.theta};
        radius = {chans.radius};
        indices = find(~cellfun('isempty', theta));
        [sph_phi sph_theta] = topo2sph( [ [ theta{indices} ]' [ radius{indices}]' ] );
        if verbose
            disp('Warning: electrodes forced to lie on a sphere for polar to 3-D conversion');
        end;
        for index = 1:length(indices)
            chans(indices(index)).sph_theta  = sph_theta(index);
            chans(indices(index)).sph_phi    = sph_phi  (index);
        end;
        if isfield(chans, 'sph_radius'),
            meanrad = mean([ chans(indices).sph_radius ]);
            if isempty(meanrad), meanrad = 1; end;
        else
            meanrad = 1;
        end;
        sph_radius(1:length(indices)) = {meanrad};
    case 'topo2sphbesa',
        chans = convertlocs(chans, 'topo2sph', varargin{:}); % search for spherical coords
        chans = convertlocs(chans, 'sph2sphbesa', varargin{:}); % search for spherical coords
    case 'topo2cart'
        chans = convertlocs(chans, 'topo2sph', varargin{:}); % search for spherical coords
        if verbose
            disp('Warning: spherical coordinates automatically updated');
        end;
        chans = convertlocs(chans, 'sph2cart', varargin{:}); % search for spherical coords
    case 'topo2all',
        chans = convertlocs(chans, 'topo2sph', varargin{:}); % search for spherical coords
        chans = convertlocs(chans, 'sph2sphbesa', varargin{:}); % search for spherical coords
        chans = convertlocs(chans, 'sph2cart', varargin{:}); % search for spherical coords
    case 'sph2cart',
        sph_theta  = {chans.sph_theta};
        sph_phi    = {chans.sph_phi};
        indices = find(~cellfun('isempty', sph_theta));
        if ~isfield(chans, 'sph_radius'), sph_radius(1:length(indices)) = {1};
        else                              sph_radius = {chans.sph_radius};
        end;
        inde = find(cellfun('isempty', sph_radius));
        if ~isempty(inde)
            meanrad = mean( [ sph_radius{:} ]);
            sph_radius(inde) = { meanrad };
        end;
        [x y z] = sph2cart([ sph_theta{indices} ]'/180*pi, [ sph_phi{indices} ]'/180*pi, [ sph_radius{indices} ]');
        for index = 1:length(indices)
            chans(indices(index)).X = x(index);
            chans(indices(index)).Y = y(index);
            chans(indices(index)).Z = z(index);
        end;
    case 'sph2topo',
        if verbose
            % disp('Warning: all radii constrained to one for spherical to topo transformation');
        end;
        sph_theta  = {chans.sph_theta};
        sph_phi    = {chans.sph_phi};
        indices = find(~cellfun('isempty', sph_theta));
        [chan_num,angle,radius] = sph2topo([ ones(length(indices),1)  [ sph_phi{indices} ]' [ sph_theta{indices} ]' ], 1, 2); % using method 2
        for index = 1:length(indices)
            chans(indices(index)).theta  = angle(index);
            chans(indices(index)).radius = radius(index);
            if ~isfield(chans, 'sph_radius') || isempty(chans(indices(index)).sph_radius)
                chans(indices(index)).sph_radius = 1;
            end;
        end;
    case 'sph2sphbesa',
        % using polar coordinates
        sph_theta  = {chans.sph_theta};
        sph_phi    = {chans.sph_phi};
        indices = find(~cellfun('isempty', sph_theta));
        [chan_num,angle,radius] = sph2topo([ones(length(indices),1)  [ sph_phi{indices} ]' [ sph_theta{indices} ]' ], 1, 2);
        [sph_theta_besa sph_phi_besa] = topo2sph([angle radius], 1, 1);
        for index = 1:length(indices)
            chans(indices(index)).sph_theta_besa  = sph_theta_besa(index);
            chans(indices(index)).sph_phi_besa    = sph_phi_besa(index);
        end;
    case 'sph2all',
        chans = convertlocs(chans, 'sph2topo', varargin{:}); % search for spherical coords
        chans = convertlocs(chans, 'sph2sphbesa', varargin{:}); % search for spherical coords
        chans = convertlocs(chans, 'sph2cart', varargin{:}); % search for spherical coords
    case 'sphbesa2sph',
        % using polar coordinates
        sph_theta_besa  = {chans.sph_theta_besa};
        sph_phi_besa    = {chans.sph_phi_besa};
        indices = find(~cellfun('isempty', sph_theta_besa));
        [chan_num,angle,radius] = sph2topo([ones(length(indices),1)  [ sph_theta_besa{indices} ]' [ sph_phi_besa{indices} ]' ], 1, 1);
        %for index = 1:length(chans)
        %   chans(indices(index)).theta  = angle(index);
        %   chans(indices(index)).radius = radius(index);
        %   chans(indices(index)).labels = int2str(index);
        %end;
        %figure; topoplot([],chans, 'style', 'blank', 'electrodes', 'labelpoint');
        
        [sph_phi sph_theta] = topo2sph([angle radius], 2);
        for index = 1:length(indices)
            chans(indices(index)).sph_theta  = sph_theta(index);
            chans(indices(index)).sph_phi    = sph_phi  (index);
        end;
    case 'sphbesa2topo',
        chans = convertlocs(chans, 'sphbesa2sph', varargin{:}); % search for spherical coords
        chans = convertlocs(chans, 'sph2topo', varargin{:}); % search for spherical coords
    case 'sphbesa2cart',
        chans = convertlocs(chans, 'sphbesa2sph', varargin{:}); % search for spherical coords
        chans = convertlocs(chans, 'sph2cart', varargin{:}); % search for spherical coords
    case 'sphbesa2all',
        chans = convertlocs(chans, 'sphbesa2sph', varargin{:}); % search for spherical coords
        chans = convertlocs(chans, 'sph2all', varargin{:}); % search for spherical coords
    case 'cart2topo',
        chans = convertlocs(chans, 'cart2sph', varargin{:}); % search for spherical coords
        chans = convertlocs(chans, 'sph2topo', varargin{:}); % search for spherical coords
    case 'cart2sphbesa',
        chans = convertlocs(chans, 'cart2sph', varargin{:}); % search for spherical coords
        chans = convertlocs(chans, 'sph2sphbesa', varargin{:}); % search for spherical coords
    case 'cart2sph',
        if verbose
            disp('WARNING: If XYZ center has not been optimized, optimize it using Edit > Channel Locations');
        end;
        X  = {chans.X};
        Y  = {chans.Y};
        Z  = {chans.Z};
        indices = find(~cellfun('isempty', X));
        [th phi radius] = cart2sph( [ X{indices} ], [ Y{indices} ], [ Z{indices} ]);
        for index = 1:length(indices)
            chans(indices(index)).sph_theta     = th(index)/pi*180;
            chans(indices(index)).sph_phi       = phi(index)/pi*180;
            chans(indices(index)).sph_radius    = radius(index);
        end;
    case 'cart2all',
        chans = convertlocs(chans, 'cart2sph', varargin{:}); % search for spherical coords
        chans = convertlocs(chans, 'sph2all', varargin{:}); % search for spherical coords
end;
% sph2topo() - Convert from a 3-column headplot file in spherical coordinates
%              to 3-column topoplot() locs file in polar (not cylindrical) coords.
%              Used for topoplot() and other 2-D topographic plotting programs.
%              Assumes a spherical coordinate system in which horizontal angles
%              have a range [-180,180] deg,  with zero pointing to the right ear.
%              In the output polar coordinate system, zero points to the nose.
%              See  >> help readlocs
% Usage:
%          >> [chan_num,angle,radius] = sph2topo(input,shrink_factor,method);
%
% Inputs:
%   input         = [channo,az,horiz] = chan_number, azumith (deg), horiz. angle (deg)
%                   When az>0, horiz=0 -> right ear, 90 -> nose
%                   When az<0, horiz=0 -> left ear, -90 -> nose
%   shrink_factor = arc_length shrinking factor>=1 (deprecated).
%                   1 -> plot edge is 90 deg azimuth {default};
%                   1.5 -> plot edge is +/-135 deg azimuth See
%                   >> help topoplot().
%   method        = [1|2], optional. 1 is for Besa compatibility, 2 is for
%                   compatibility with Matlab function cart2sph(). Default is 2
%
% Outputs:
%   channo  = channel number (as in input)
%   angle   = horizontal angle (0 -> nose; 90 -> right ear; -90 -> left ear)
%   radius  = arc_lengrh from vertex (Note: 90 deg az -> 0.5/shrink_factor);
%             By topoplot() convention, radius=0.5 is the nasion-ear_canal plane.
%             Use topoplot() 'plotrad' to plot chans with abs(az) > 90 deg.
%
% Author: Scott Makeig & Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 6/12/98
%
% See also: cart2topo(), topo2sph()

% Copyright (C) 6/12/98 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% corrected left/right orientation mismatch, Blair Hicks 6/20/98
% changed name sph2pol() -> sph2topo() for compatibility -sm
% 01-25-02 reformated help & license -ad
% 01-25-02 changed computation so that it works with sph2topo -ad

function [channo,angle,radius] = sph2topo(input,factor, method)

chans = size(input,1);
angle = zeros(chans,1);
radius = zeros(chans,1);

if nargin < 1
    help sph2topo
    return
end

if nargin< 2
    factor = 0;
end
if factor==0
    factor = 1;
end
if factor < 1
    help sph2topo
    return
end

if size(input,2) ~= 3
    help sph2topo
    return
end

channo = input(:,1);
az = input(:,2);
horiz = input(:,3);

if exist('method')== 1 & method == 1
    radius = abs(az/180)/factor;
    i = find(az>=0);
    angle(i) = 90-horiz(i);
    i = find(az<0);
    angle(i) = -90-horiz(i);
else
    angle  = -horiz;
    radius = 0.5 - az/180;
end;
% topo2sph() - convert a topoplot() style 2-D polar-coordinate
%              channel locations file to a 3-D spherical-angle
%              file for use with headplot()
% Usage:
%   >> [c h] = topo2sph('eloc_file','eloc_outfile', method, unshrink);
%   >> [c h] = topo2sph( topoarray, method, unshrink );
%
% Inputs:
%   'eloc_file'    = filename of polar 2-D electrode locations file used by
%                    topoplot(). See >> topoplot example or cart2topo()
%   'eloc_outfile' = output file of 3-D electrode locations in spherical angle
%                    coords. for use in headplot().
%   topoarray      = polar array of 2-D electrode locations, with polar angle
%                    in the first column and radius in the second one.
%   method         = [1|2] 1 is for Besa compatibility, 2 is for
%                    compatibility with Matlab function cart2sph(). {default: 2}
%   unshrink       = [0<real<1] unshrink factor. Enter a shrink factor used
%                    to convert spherical to topo (see sph2topo()). Only
%                    implemented for 'method' 1 (above). Electrode 'shrink'
%                    is now deprecated. See >> help topoplot
% Outputs:
%   c = coronal rotation
%   h = horizontal rotation
%
% Author: Scott Makeig & Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 1999
%
% See also: sph2topo(), cart2topo()

% Copyright (C) 1999 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% 3-16-00 changed name to topo2sph() for compatibility with cart2topo() -sm
% 01-25-02 reformated help & license -ad
% 03-22-02 complete remodeling for returning arguments and taking arrays -ad

function [c, h] = topo2sph(eloc_locs,eloc_angles, method, unshrink)

MAXCHANS = 1024;

if nargin < 1
    help topo2sph;
    return;
end;
if nargin > 1 && ~isstr(eloc_angles)
    if nargin > 2
        unshrink = method;
    end;
    method = eloc_angles;
else
    method = 2;
end;

if isstr(eloc_locs)
    fid = fopen(eloc_locs);
    if fid<1,
        fprintf('topo2sph()^G: cannot open eloc_loc file (%s)\n',eloc_locs)
        return
    end
    E = fscanf(fid,'%d %f %f  %s',[7 MAXCHANS]);
    E = E';
    fclose(fid);
else
    E = eloc_locs;
    E = [ ones(size(E,1),1) E ];
end;

if nargin > 1 & isstr(eloc_angles)
    if exist(eloc_angles)==2,
        fprintf('topo2sph: eloc_angles file (%s) already exists and will be erased.\n',eloc_angles);
    end
    
    fid = fopen(eloc_angles,'a');
    if fid<1,
        fprintf('topo2sph()^G: cannot open eloc_angles file (%s)\n',eloc_angles)
        return
    end
end;

if method == 2
    t = E(:,2); % theta
    r = E(:,3); % radius
    h = -t;  % horizontal rotation
    c = (0.5-r)*180;
else
    for e=1:size(E,1)
        % (t,r) -> (c,h)
        
        t = E(e,2); % theta
        r = E(e,3); % radius
        r = r*unshrink;
        if t>=0
            h(e) = 90-t; % horizontal rotation
        else
            h(e) = -(90+t);
        end
        if t~=0
            c(e) = sign(t)*180*r; % coronal rotation
        else
            c(e) = 180*r;
        end
    end;
    t = t';
    r = r';
end;

for e=1:size(E,1)
    if nargin > 1 & isstr(eloc_angles)
        chan = E(e,4:7);
        fprintf('%d	%g	%g	%s\n',E(e,1),c(e),h(e),chan);
        fprintf(fid,'%d	%g	%g	%s\n',E(e,1),c(e),h(e),chan);
    end;
end

% eeg_emptyset() - Initialize an EEG dataset structure with default values.
%
% Usage:
%   >> EEG = eeg_emptyset();
%
% Outputs:
%   EEG    - empty dataset structure with default values.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeglab()

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% 01-25-02 reformated help & license -ad

function EEG = eeg_emptyset();

EEG.setname     = '';
EEG.filename    = '';
EEG.filepath    = '';
EEG.subject     = '';
EEG.group       = '';
EEG.condition   = '';
EEG.session     = [];
EEG.comments    = '';
EEG.nbchan      = 0;
EEG.trials      = 0;
EEG.pnts        = 0;
EEG.srate       = 1;
EEG.xmin        = 0;
EEG.xmax        = 0;
EEG.times       = [];
EEG.data        = [];
EEG.icaact      = [];
EEG.icawinv     = [];
EEG.icasphere   = [];
EEG.icaweights  = [];
EEG.icachansind = [];
EEG.chanlocs    = [];
EEG.urchanlocs  = [];
EEG.chaninfo    = [];
EEG.ref         = [];
EEG.event       = [];
EEG.urevent     = [];
EEG.eventdescription = {};
EEG.epoch       = [];
EEG.epochdescription = {};
EEG.reject      = [];
EEG.stats       = [];
EEG.specdata    = [];
EEG.specicaact  = [];
EEG.splinefile  = '';
EEG.icasplinefile = '';
EEG.dipfit      = [];
EEG.history     = '';
EEG.saved       = 'no';
EEG.etc         = [];

%EEG.reject.threshold  = [1 0.8 0.85];
%EEG.reject.icareject  = [];
%EEG.reject.compreject = [];
%EEG.reject.gcompreject= [];
%EEG.reject.comptrial  = [];
%EEG.reject.sigreject  = [];
%EEG.reject.elecreject = [];

%EEG.stats.kurta      = [];
%EEG.stats.kurtr      = [];
%EEG.stats.kurtd      = [];
%EEG.stats.eegentropy = [];
%EEG.stats.eegkurt    = [];
%EEG.stats.eegkurtg   = [];
%EEG.stats.entropy    = [];
%EEG.stats.kurtc      = [];
%EEG.stats.kurtt      = [];
%EEG.stats.entropyc   = [];

return;
