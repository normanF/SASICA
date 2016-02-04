% [EEG, cfg] = eeg_SASICA(EEG,cfg)
%
% Suggest components to reject from an EEG dataset with ICA decomposition.
%
% Inputs: EEG: EEGlab structure with ICA fields.
%         cfg: structure describing which methods are to use for suggesting
%              bad components (see structure called def, in the code below)
%              Available methods are:
%              Autocorrelation: detects noisy components with weak
%                               autocorrelation (muscle artifacts usually)
%              Focal components: detects components that are too focal and
%                               thus unlikely to correspond to neural
%                               activity (bad channel or muscle usually).
%              Focal trial activity: detects components with focal trial
%                               activity, with same algorhithm as focal
%                               components above. Results similar to trial
%                               variability.
%              Signal to noise ratio: detects components with weak signal
%                               to noise ratio between arbitrary baseline
%                               and interest time windows.
%              Dipole fit residual variance: detects components with high
%                               residual variance after subtraction of the
%                               forward dipole model. Note that the inverse
%                               dipole modeling using DIPFIT2 in EEGLAB
%                               must have been computed to use this
%                               measure.
%              EOG correlation: detects components whose time course
%                               correlates with EOG channels.
%              Bad channel correlation: detects components whose time course
%                               correlates with any channel(s).
%              ADJUST selection: use ADJUST routines to select components
%                               (see Mognon, A., Jovicich, J., Bruzzone,
%                               L., & Buiatti, M. (2011). ADJUST: An
%                               automatic EEG artifact detector based on
%                               the joint use of spatial and temporal
%                               features. Psychophysiology, 48(2), 229-240.
%                               doi:10.1111/j.1469-8986.2010.01061.x)
%              FASTER selection: use FASTER routines to select components
%                               (see Nolan, H., Whelan, R., & Reilly, R. B.
%                               (2010). FASTER: Fully Automated Statistical
%                               Thresholding for EEG artifact Rejection.
%                               Journal of Neuroscience Methods, 192(1),
%                               152-162. doi:16/j.jneumeth.2010.07.015)
%              MARA selection:  use MARA classification engine to select components
%                               (see Winkler I, Haufe S, Tangermann M.
%                               2011. Automatic Classification of
%                               Artifactual ICA-Components for Artifact
%                               Removal in EEG Signals. Behavioral and
%                               Brain Functions. 7:30.)
%
%              Options: noplot: just compute and store result in EEG. Do
%                           not make any plots.
%
% If you use this program in your research, please cite the following
% article:
%   Chaumon M, Bishop DV, Busch NA. A Practical Guide to the Selection of
%   Independent Components of the Electroencephalogram for Artifact
%   Correction. Journal of neuroscience methods. 2015
%
%   SASICA is a software that helps select independent components of
%   the electroencephalogram based on various signal measures.
%     Copyright (C) 2014  Maximilien Chaumon
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.


function [EEG, cfg] = eeg_SASICA(EEG,cfg)

if nargin < 1
    error('Need at least one input argument')
end
if ~exist('cfg','var')
    cfg = struct;
end
% deal with calling pop_prop_ADJ or pop_prop_FST here
if ischar(cfg) && strncmp(cfg,'pop_',4)
    eval(cfg);
    return
end
%
PLOTPERFIG = 35;
def = SASICA('getdefs');

cfg = setdef(cfg,def);
v = regexp(version,'^\d+\.\d+','match');
if num2str(v{1}) >= 8.4
    mkersize = 25;
else
    mkersize = 20;
end

if isempty(EEG.icaact) && isempty(EEG.icawinv)
    errordlg('No ica weights in the current EEG dataset! Compute ICA on your data first.')
    error('No ica weights! Compute ICA on your data first.')
end
struct2ws(cfg.opts);

rejfields = {'icarejautocorr' 'Autocorrelation' [         0         0    1.0000]
    'icarejfocalcomp' 'Focal components' [         0    0.5000         0]
    'icarejtrialfoc' 'Focal trial activity' [    0.7500         0    0.7500]
    'icarejSNR' 'Signal to noise ' [    0.8000         0         0]
    'icarejresvar' 'Residual variance' [    0     0.7500    0.7500]
    'icarejchancorr' 'Correlation with channels' [    0.7500    0.7500         0]
    'icarejADJUST' 'ADJUST selections' [    .3 .3 .3]
    'icarejFASTER' 'FASTER selections' [    0 .7 0]
    'icarejMARA' 'MARA selections' [    .5 .5 0]
    };

rejects = zeros(size(rejfields,1),1);

if numel(noplot) == 1
    noplot = noplot * ones(1,size(rejfields,1));
end

if any(~noplot)
    figure(321541);clf;% just a random number so we always work in the same figure
    BACKCOLOR           = [.93 .96 1];
    set(gcf,'numbertitle', 'off','name','Automatic component rejection measures','color',BACKCOLOR)
    isubplot = 1;
end

ncomp= size(EEG.icaact,1); % ncomp is number of components
if ncomp == 0
    ncomp = size(EEG.icaweights,1);
end
if ~nocompute
    EEG.icaact = eeg_getica(EEG,1:ncomp);
    EEG.reject.SASICA = [];
    for ifield = 1:size(rejfields,1)
        %     EEG.reject.SASICA.(rejfields{ifield}) = false(1,ncomp);
        EEG.reject.SASICA.([rejfields{ifield} 'col']) = rejfields{ifield,3};
    end
    fprintf('Computing selection methods...\n')
end
if cfg.autocorr.enable
    rejects(1) = 1;
    disp('Autocorrelation.')
    %% Autocorrelation
    % Identifying noisy components
    %----------------------------------------------------------------
    struct2ws(cfg.autocorr);
    
    if ~nocompute
        Ncorrint=round(autocorrint/(1000/EEG.srate)); % number of samples for lag
        rej = false(1,ncomp);
        for k=1:ncomp
            y=EEG.icaact(k,:,:);
            yy=xcorr(mean(y,3),Ncorrint,'coeff');
            autocorr(k) = yy(1);
        end
        dropautocorr = readauto(dropautocorr,autocorr,'-');
        for k = 1:ncomp
            if autocorr(k) < dropautocorr
                rej(k)=true;
            end
        end
        EEG.reject.SASICA.(strrep(rejfields{1,1},'rej','')) = autocorr;
        EEG.reject.SASICA.(strrep(rejfields{1,1},'rej','thresh')) = dropautocorr;
        EEG.reject.SASICA.(rejfields{1,1}) = logical(rej);
    else
        autocorr = EEG.reject.SASICA.(strrep(rejfields{1,1},'rej',''));
        dropautocorr = EEG.reject.SASICA.(strrep(rejfields{1,1},'rej','thresh'));
        rej = EEG.reject.SASICA.(rejfields{1,1});
    end
    %----------------------------------------------------------------
    if ~noplot(1)
        subplot(2,3,isubplot);cla;isubplot = isubplot+1;
        set(gca,'fontsize',FontSize)
        plot(autocorr,'k','linestyle','none');
        
        hold on
        xlim([0 ncomp+1]);
        s = std(autocorr);
        m = mean(autocorr);
        yl = ylim;xl = xlim;
        [x,y] = meshgrid(xl(1):.1:xl(2),yl(1):.1:yl(2));
        galpha = 1./(s*(2*pi)^.5).*exp(-(y-m).^2./(2.*s^2));
        %     h = surf(x,y,-ones(size(y)));shading flat
        %     color = [ 0 0 0]';
        %     C = repmat(color,[1,size(y)]);
        %     C = permute(C,[2 3 1]);
        %     set(h,'alphadata',1-galpha,'alphadatamapping','scaled','facealpha','interp',...
        %         'CData',C,'CDataMapping','direct')
        hline(dropautocorr,'r');
        
        plot(xl(2)-diff(xl)/20,yl(2)-diff(yl)/20,'marker','.','color',rejfields{1,3},'markersize',40)
        xlabel('Components')
        ylabel('Autocorrelation')
        title(['Autocorrelation at ' num2str(autocorrint) ' ms.'])
        toplot = autocorr;
        toplot(toplot > dropautocorr) = NaN;
        plot(toplot,'o','color',rejfields{1,3})
        for i = 1:numel(autocorr)
            h = scatter(i,autocorr(i),mkersize,'k','filled');
            cb = sprintf('eeg_SASICA(EEG, ''pop_prop( %s, 0, %d, findobj(''''tag'''',''''comp%d''''), { ''''freqrange'''', [1 50] })'');', inputname(1), i, i);
            set(h,'buttondownfcn',cb);
        end
    end
    
end
if cfg.focalcomp.enable
    rejects(2) = 1;
    disp('Focal components.')
    %% Focal activity
    %----------------------------------------------------------------
    struct2ws(cfg.focalcomp);
    if ~nocompute
        rej = false(1,ncomp);
        clear mywt
        for k=1:ncomp
            mywt(:,k) = sort(abs(zscore(EEG.icawinv(:,k))),'descend'); %sorts standardized weights in descending order
        end
        focalICAout = readauto(focalICAout,mywt(1,:),'+');
        for k = 1:ncomp
            if mywt(1,k) > focalICAout
                rej(k)=true;
            end
        end
        EEG.reject.SASICA.(strrep(rejfields{2,1},'rej','')) = mywt(1,:);
        EEG.reject.SASICA.(strrep(rejfields{2,1},'rej','thresh')) = focalICAout;
        EEG.reject.SASICA.(rejfields{2,1}) = logical(rej);
    else
        mywt(1,:) = EEG.reject.SASICA.(strrep(rejfields{2,1},'rej',''));
        focalICAout = EEG.reject.SASICA.(strrep(rejfields{2,1},'rej','thresh'));
        rej = EEG.reject.SASICA.(rejfields{2,1});
    end
    %----------------------------------------------------------------
    if ~noplot(2)
        subplot(2,3,isubplot);cla;isubplot = isubplot+1;
        set(gca,'fontsize',FontSize)
        toplot = mywt(1,:);
        plot(toplot,'k','linestyle','none');
        toplot(toplot < focalICAout) = NaN;
        hold on
        hline(focalICAout,'r');
        plot(toplot,'o','color',rejfields{2,3});
        xlim([0 ncomp+1]);
        xl = xlim;yl = ylim;
        plot(xl(2)-diff(xl)/20,yl(2)-diff(yl)/20,'marker','.','color',rejfields{2,3},'markersize',40)
        xlabel('Components')
        ylabel('Standardized weights')
        title('Components with focal activity')
        for i = 1:numel(mywt(1,:))
            h = scatter(i,mywt(1,i),mkersize,'k','filled');
            cb = sprintf('eeg_SASICA(EEG, ''pop_prop( %s, 0, %d, findobj(''''tag'''',''''comp%d''''), { ''''freqrange'''', [1 50] })'');', inputname(1), i, i);
            set(h,'buttondownfcn',cb);
        end
    end
    
end

if cfg.trialfoc.enable
    rejects(3) = 1;
    disp('Focal trial activity.');
    %% Focal trial activity
    struct2ws(cfg.trialfoc);
    if ~nocompute
        % Find components with focal trial activity (those that have activity
        % on just a few trials and are almost zero on others)
        %----------------------------------------------------------------
        if ndims(EEG.icaact) < 3
            error('This method cannot be used on continuous data (no ''trials''!)');
        end
        myact =sort(abs(zscore(range(EEG.icaact,2),[],3)),3,'descend'); % sorts standardized range of trial activity
        focaltrialout = readauto(focaltrialout,myact(:,:,1)','+');
        % in descending order
        rej = myact(:,:,1) > focaltrialout;
        EEG.reject.SASICA.(strrep(rejfields{3,1},'rej','')) = myact(:,:,1)';
        EEG.reject.SASICA.(strrep(rejfields{3,1},'rej','thresh')) = focaltrialout;
        EEG.reject.SASICA.(rejfields{3,1}) = rej';
    else
        myact = EEG.reject.SASICA.(strrep(rejfields{3,1},'rej',''))';
        focaltrialout = EEG.reject.SASICA.(strrep(rejfields{3,1},'rej','thresh'));
        rej = EEG.reject.SASICA.(rejfields{3,1})';
    end
    %----------------------------------------------------------------
    if ~noplot(3)
        subplot(2,3,isubplot);cla;isubplot = isubplot+1;
        if EEG.trials > 1
            set(gca,'fontsize',FontSize)
            toplot = myact(:,:,1);
            plot(toplot,'k','linestyle','none')
            hold on
            toplot(toplot < focaltrialout) = NaN;
            plot(1:ncomp,toplot,'o','color',rejfields{3,3});
            xlim([0 ncomp+1])
            hline(focaltrialout,'r');
            xl = xlim;yl =ylim;
            xlabel('Components')
            ylabel('Standardized peak trial activity')
            plot(xl(2)-diff(xl)/20,yl(2)-diff(yl)/20,'marker','.','color',rejfields{3,3},'markersize',40)
            for i = 1:numel(myact(:,:,1))
                h = scatter(i,myact(i),mkersize,'k','filled');
                cb = sprintf('eeg_SASICA(EEG, ''pop_prop( %s, 0, %d, findobj(''''tag'''',''''comp%d''''), { ''''freqrange'''', [1 50] })'');', inputname(1), i, i);
                set(h,'buttondownfcn',cb);
            end
            
            title(['Focal trial activity'])
        else
            xl = xlim;yl = ylim;
            text(xl(1)+diff(xl)/2,yl(1)+diff(yl)/2,{'Only one trial.' 'Focal trial' 'activity method'  'is inadequate.'},'horizontalalignment','center');
            axis off
        end
    end
    %----------------------------------------------------------------
end

if cfg.SNR.enable
    rejects(4) = 1;
    disp('Signal to noise ratio.')
    %% Low Signal to noise components
    struct2ws(cfg.SNR);
    if ~nocompute
        rejfields{4,2} = ['Signal to noise Time of interest ' num2str(snrPOI,'%g ') ' and Baseline ' num2str(snrBL,'%g ') ' ms.'];
        
        POIpts = timepts(snrPOI);
        BLpts = timepts(snrBL);
        
        zz = zscore(EEG.icaact,[],2);% zscore along time
        av1 = mean(zz(:,POIpts,:),3); % average activity in POI across trials
        av2 = mean(zz(:,BLpts,:),3); % activity in baseline acros trials
        SNR = std(av1,[],2)./std(av2,[],2); % ratio of the standard deviations of activity and baseline
        snrcut = readauto(snrcut,SNR,'-');
        rej = SNR < snrcut;
        EEG.reject.SASICA.(strrep(rejfields{4,1},'rej','')) = SNR';
        EEG.reject.SASICA.(strrep(rejfields{4,1},'rej','thresh')) = snrcut;
        EEG.reject.SASICA.(rejfields{4,1}) = rej';
    else
        SNR = EEG.reject.SASICA.(strrep(rejfields{4,1},'rej',''))';
        snrcut = EEG.reject.SASICA.(strrep(rejfields{4,1},'rej','thresh'));
        rej = EEG.reject.SASICA.(rejfields{4,1})';
    end
    %----------------------------------------------------------------
    if ~noplot(4)
        subplot(2,3,isubplot);cla;isubplot = isubplot+1;
        set(gca,'fontsize',FontSize)
        plot(SNR,'k','linestyle','none');
        hold on
        xlim([0 ncomp+1]);
        xl = xlim; yl = ylim;
        hline(snrcut,'r');
        toplot = SNR;
        toplot(toplot > snrcut) = NaN;
        plot(toplot,'o','color',rejfields{4,3})
        plot(xl(2)-diff(xl)/20,yl(2)-diff(yl)/20,'marker','.','color',rejfields{4,3},'markersize',40)
        for i = 1:numel(SNR)
            h = scatter(i,SNR(i),mkersize,'k','filled');
            cb = sprintf('eeg_SASICA(EEG, ''pop_prop( %s, 0, %d, findobj(''''tag'''',''''comp%d''''), { ''''freqrange'''', [1 50] })'');', inputname(1), i, i);
            set(h,'buttondownfcn',cb);
        end
        title({'Signal to noise ratio between' ['Time of interest ' num2str(snrPOI,'%g ') ' and Baseline ' num2str(snrBL,'%g ') ' ms.']})
        xlabel('Components')
        ylabel('SNR')
    end
    
    %----------------------------------------------------------------
end

if cfg.resvar.enable
    rejects(5) = 1;
    disp('Residual variance thresholding.')
    %% High residual variance
    struct2ws(cfg.resvar);
    if ~nocompute
        resvar = 100*[EEG.dipfit.model.rv];
        rej = resvar > thresh;
        
        EEG.reject.SASICA.(strrep(rejfields{5,1},'rej','')) = resvar;
        EEG.reject.SASICA.(strrep(rejfields{5,1},'rej','thresh')) = thresh;
        EEG.reject.SASICA.(rejfields{5,1}) = rej;
    else
        resvar = EEG.reject.SASICA.(strrep(rejfields{5,1},'rej',''));
        thresh = EEG.reject.SASICA.(strrep(rejfields{5,1},'rej','thresh'));
        rej = EEG.reject.SASICA.(rejfields{5,1});
    end
    %----------------------------------------------------------------
    if ~noplot(5)
        subplot(2,3,isubplot);cla;isubplot = isubplot+1;
        set(gca,'fontsize',FontSize)
        plot(resvar,'k','linestyle','none');
        hold on
        xlim([0 ncomp+1]);
        ylim([0 100]);
        xl = xlim; yl = ylim;
        hline(thresh,'r');
        toplot = resvar;
        toplot(toplot < thresh) = NaN;
        plot(toplot,'o','color',rejfields{5,3})
        plot(xl(2)-diff(xl)/20,yl(2)-diff(yl)/20,'marker','.','color',rejfields{5,3},'markersize',40)
        for i = 1:numel(resvar)
            h = scatter(i,resvar(i),mkersize,'k','filled');
            cb = sprintf('eeg_SASICA(EEG, ''pop_prop( %s, 0, %d, findobj(''''tag'''',''''comp%d''''), { ''''freqrange'''', [1 50] })'');', inputname(1), i, i);
            set(h,'buttondownfcn',cb);
        end
        title({'Residual variance of dipole fit'})
        xlabel('Components')
        ylabel('RV (%)')
    end
    
    %----------------------------------------------------------------
end

if cfg.EOGcorr.enable
    rejects(6) = 1;
    disp('Correlation with EOGs.');
    %% Correlation with EOG
    struct2ws(cfg.EOGcorr);
    if ~nocompute
        noV = 0;noH = 0;
        try
            Veogchan = chnb(Veogchannames);
        catch
            Veogchan = [];
        end
        try
            Heogchan = chnb(Heogchannames);
        catch
            Heogchan = [];
        end
        if numel(Veogchan) == 1
            VEOG = EEG.data(Veogchan,:,:);
        elseif numel(Veogchan) == 2
            VEOG = EEG.data(Veogchan(1),:,:) - EEG.data(Veogchan(2),:,:);
        else
            disp('no Vertical EOG channels...');
            noV = 1;
        end
        if numel(Heogchan) == 1
            HEOG = EEG.data(Heogchan,:,:);
        elseif numel(Heogchan) == 2
            HEOG = EEG.data(Heogchan(1),:,:) - EEG.data(Heogchan(2),:,:);
        else
            disp('no Horizontal EOG channels...');
            noH = 1;
        end
        ICs = EEG.icaact(:,:)';
        if ~noV
            VEOG = VEOG(:);
            cV  = abs(corr(ICs,VEOG))';
            corthreshV = readauto(corthreshV,cV,'+');
            rejV = cV > corthreshV ;
        else
            cV = NaN(1,size(ICs,2));
            corthreshV = 0;
            rejV = false(size(cV));
        end
        if ~noH
            HEOG = HEOG(:);
            cH  = abs(corr(ICs,HEOG))';
            corthreshH = readauto(corthreshH,cH,'+');
            rejH = cH > corthreshH;
        else
            cH = NaN(1,size(ICs,2));
            corthreshH = 0;
            rejH = false(size(cH));
        end
        
        EEG.reject.SASICA.([strrep(rejfields{6,1},'rej','') 'VEOG']) = cV;
        EEG.reject.SASICA.([strrep(rejfields{6,1},'rej','thresh') 'VEOG']) = corthreshV;
        EEG.reject.SASICA.([strrep(rejfields{6,1},'rej','') 'HEOG']) = cH;
        EEG.reject.SASICA.([strrep(rejfields{6,1},'rej','thresh') 'HEOG']) = corthreshH;
        EEG.reject.SASICA.(rejfields{6,1}) = [rejV|rejH];
    else
        if existnotempty(EEG.reject.SASICA,[strrep(rejfields{6,1},'rej','') 'VEOG'])
            cV = EEG.reject.SASICA.([strrep(rejfields{6,1},'rej','') 'VEOG']);
            corthreshV = EEG.reject.SASICA.([strrep(rejfields{6,1},'rej','thresh') 'VEOG']);
        end
        if existnotempty(EEG.reject.SASICA,[strrep(rejfields{6,1},'rej','') 'HEOG'])
            cH = EEG.reject.SASICA.([strrep(rejfields{6,1},'rej','') 'HEOG']);
            corthreshH = EEG.reject.SASICA.([strrep(rejfields{6,1},'rej','thresh') 'HEOG']);
        end
    end
    %----------------------------------------------------------------
    if ~noplot(6)
        subplot(2,3,isubplot);cla;isubplot = isubplot+1;
        set(gca,'fontsize',FontSize)
        cols = get(gca,'colororder');
        [hplotcorr] = plot([cV;cH]','.','linestyle','none');
        icol = 2;
        hold all
        xlim([0 ncomp+1]);
        xl = xlim;yl = ylim;
        hline(corthreshV,'color',cols(1,:));
        hline(corthreshH,'color',cols(2,:));
        
        title(['Correlation with EOG'])
        legstr = {'VEOG' 'HEOG'};
        ylabel('Correlation coef (r)');
        xlabel('Components');
        toplot = cV;
        toplot(toplot < corthreshV) = NaN;
        plot(1:ncomp,toplot,'o','color',rejfields{6,3})
        toplot = cH;
        toplot(toplot < corthreshH) = NaN;
        plot(1:ncomp,toplot,'o','color',rejfields{6,3})
        plot(xl(2)-diff(xl)/20,yl(2)-diff(yl)/20,'marker','.','color',rejfields{6,3},'markersize',40)
        legend(legstr,'fontsize',10, 'location', 'best');
        for i = 1:numel(cH)
            h(1) = scatter(i,cH(i),mkersize,cols(1,:),'filled');
            h(2) = scatter(i,cV(i),mkersize,cols(2,:),'filled');
            cb = sprintf('eeg_SASICA(EEG, ''pop_prop( %s, 0, %d, findobj(''''tag'''',''''comp%d''''), { ''''freqrange'''', [1 50] })'');', inputname(1), i, i);
            set(h,'buttondownfcn',cb);
        end
    end
    %----------------------------------------------------------------
end

if cfg.chancorr.enable
    rejects(6) = 1;
    disp('Correlation with other channels.')
    %% Correlation with other channels
    struct2ws(cfg.chancorr);
    if ~nocompute
        if ~cfg.EOGcorr.enable
            rejH = false(1,ncomp);
            rejV = false(1,ncomp);
        end
        if ~isempty(channames)
            try
                [chan cellchannames channames] = chnb(channames);
            end
            chanEEG = EEG.data(chan,:)';
            ICs = EEG.icaact(:,:)';
            c  = abs(corr(ICs,chanEEG))';
            corthresh = mean(readauto(corthresh,c,'+'));
            rej = c > corthresh ;
            if size(rej,1) > 1
                rej = sum(rej)>=1;
            end
            EEG.reject.SASICA.([strrep(rejfields{6,1},'rej','') 'chans']) = c;
            EEG.reject.SASICA.([strrep(rejfields{6,1},'rej','thresh') 'chans']) = corthresh;
            EEG.reject.SASICA.(rejfields{6,1}) = [rej|rejH|rejV];
        else
            noplot(6) = 1;
            disp('Could not find the channels to compute correlation.');
            c = NaN(1,ncomp);
            EEG.reject.SASICA.([strrep(rejfields{6,1},'rej','') 'chans']) = c;
            EEG.reject.SASICA.([strrep(rejfields{6,1},'rej','thresh') 'chans']) = corthresh;
            rej = false(1,ncomp);
            EEG.reject.SASICA.(rejfields{6,1}) = [rej|rejV|rejH];
        end
    else
        c = EEG.reject.SASICA.([strrep(rejfields{6,1},'rej','') 'chans']);
        corthresh = EEG.reject.SASICA.([strrep(rejfields{6,1},'rej','thresh') 'chans']);
    end
    %----------------------------------------------------------------
    if ~noplot(6);
        if exist('hplotcorr','var')
            isubplot = isubplot-1;
        end
        subplot(2,3,isubplot);
        if ~cfg.EOGcorr.enable
            cla;
            set(gca,'fontsize',FontSize);
            cols = get(gca,'colororder');
        end
        hold all
        if not(exist('hplotcorr','var'))
            hplotcorr = [];
        end
        icol = numel(hplotcorr);
        for ichan = 1:numel(chan)
            [hplotcorr(end+1)] = plot([c(ichan,:)]','.','linestyle','none','color',cols(rem(icol+ichan-1,size(cols,1))+1,:));
        end
        xlim([0 ncomp+1]);
        xl = xlim;yl = ylim;
        hline(corthresh,'r');
        title(['Correlation with channels'])
        if cfg.EOGcorr.enable
            legstr = {'VEOG' 'HEOG' cellchannames{:}};
        else
            legstr = {cellchannames{:}};
        end
        ylabel('Correlation coef (r)');
        xlabel('Components');
        toplot = c;
        for i = 1:size(toplot,1)
            toplot(i,toplot(i,:) < corthresh) = NaN;
        end
        plot(1:ncomp,toplot,'o','color',rejfields{6,3})
        plot(xl(2)-diff(xl)/20,yl(2)-diff(yl)/20,'marker','.','color',rejfields{6,3},'markersize',40)
        legend(hplotcorr,legstr,'fontsize',10, 'location', 'best');
        for ichan = 1:size(c,1)
            for i = 1:size(c,2)
                h = scatter(i,c(ichan,i),mkersize,cols(rem(icol+ichan-1,size(cols,1))+1,:),'filled');
                cb = sprintf('eeg_SASICA(EEG, ''pop_prop( %s, 0, %d, findobj(''''tag'''',''''comp%d''''), { ''''freqrange'''', [1 50] })'');', inputname(1), i, i);
                set(h,'buttondownfcn',cb);
            end
        end
        
    end
    %----------------------------------------------------------------
end
if cfg.ADJUST.enable
    rejects(7) = 1;
    disp('ADJUST methods selection')
    %% ADJUST
    struct2ws(cfg.ADJUST);
    if ~nocompute
        [art, horiz, vert, blink, disc,...
            soglia_DV, diff_var, soglia_K, med2_K, meanK, soglia_SED, med2_SED, SED, soglia_SAD, med2_SAD, SAD, ...
            soglia_GDSF, med2_GDSF, GDSF, soglia_V, med2_V, nuovaV, soglia_D, maxdin] = ADJUST (EEG);
        
        ADJ.art = art;ADJ.horiz = horiz;ADJ.vert = vert;ADJ.blink = blink;ADJ.disc = disc;
        
        ADJ.soglia_DV = soglia_DV; ADJ.diff_var = diff_var;
        ADJ.soglia_K = soglia_K;ADJ.med2_K = med2_K; ADJ.meanK = meanK;
        ADJ.soglia_SED = soglia_SED; ADJ.med2_SED = med2_SED;ADJ.SED = SED;
        ADJ.med2_SAD = med2_SAD;ADJ.soglia_SAD = soglia_SAD;ADJ.SAD = SAD;
        ADJ.soglia_GDSF = soglia_GDSF; ADJ.med2_GDSF = med2_GDSF;ADJ.GDSF = GDSF;
        ADJ.soglia_V = soglia_V;ADJ.med2_V = med2_V;ADJ.nuovaV = nuovaV;
        ADJ.soglia_D = soglia_D; ADJ.maxdin = maxdin;
        
        rej = false(1,size(EEG.icaact,1));
        rej([ADJ.art ADJ.horiz ADJ.vert ADJ.blink ADJ.disc]) = true;
        
        EEG.reject.SASICA.(strrep(rejfields{7,1},'rej','')) = ADJ;
        EEG.reject.SASICA.(rejfields{7,1}) = rej;
    else
        ADJ = EEG.reject.SASICA.(strrep(rejfields{7,1},'rej',''));
    end
    %----------------------------------------------------------------
end
if cfg.FASTER.enable
    rejects(8) = 1;
    disp('FASTER methods selection')
    %% FASTER
    struct2ws(cfg.FASTER);
    if ~nocompute
        blinkchans = chnb(blinkchans);
        listprops = component_properties(EEG,blinkchans);
        FST.rej = min_z(listprops)' ~= 0;
        FST.listprops = listprops;
        
        EEG.reject.SASICA.(strrep(rejfields{8,1},'rej','')) = FST;
        EEG.reject.SASICA.(rejfields{8,1}) = FST.rej;
    else
        FST = EEG.reject.SASICA.(strrep(rejfields{8,1},'rej',''));
    end
    %----------------------------------------------------------------
end
if cfg.MARA.enable
    rejects(9) = 1;
    disp('MARA methods selection')
    %% MARA
    struct2ws(cfg.MARA);
    if ~nocompute
        [rej info] = MARA(EEG);
        MR.rej = false(1,size(EEG.icaact,1));
        MR.rej(rej) = true;
        MR.info = info;
        
        EEG.reject.SASICA.(strrep(rejfields{9,1},'rej','')) = MR;
        EEG.reject.SASICA.(rejfields{9,1}) = MR.rej;
    else
        MR = EEG.reject.SASICA.(strrep(rejfields{9,1},'rej',''));
    end
    %----------------------------------------------------------------
end

EEG.reject.SASICA.var = var(EEG.icaact(:,:),[],2);% variance of each component

if (cfg.ADJUST.enable||cfg.FASTER.enable) && any(~noplot)
    h = uicontrol('style','text','string','for ADJUST or FASTER results, click on component buttons in the other window(s)','units','normalized','position',[0 0 1 .05],'backgroundcolor',get(gcf,'color'));
    uistack(h,'bottom')
end
fprintf('... Done.\n')

drawnow

%% Final computations
% combine in gcompreject field and pass to pop_selectcomps
EEG.reject.gcompreject = false(1,ncomp);
for ifield = 1:size(rejfields,1)
    if rejects(ifield)
        EEG.reject.gcompreject = [EEG.reject.gcompreject ; EEG.reject.SASICA.(rejfields{ifield})];
    end
end
EEG.reject.gcompreject = sum(EEG.reject.gcompreject) >= 1;

%% plotting
try
    delete(findobj('-regexp','name','pop_selectcomps'))
    drawnow
end
if any(~noplot)
    if ~isempty([EEG.chanlocs.radius])% assume we have sensor locations...
        clear hfig
        delete(findobj('tag','waitcomp'))
        textprogressbar('Drawing topos...');
        for ifig = 1:ceil((ncomp)/PLOTPERFIG)
            cmps = [1+(ifig-1)*PLOTPERFIG:min([ncomp,ifig*PLOTPERFIG])];
            eeg_SASICA(EEG,['pop_selectcomps(EEG, [' num2str(cmps) '],' num2str(ncomp) ');']);
            hfig(ifig) = gcf;
            set(hfig(ifig),'name',[get(hfig(ifig),'name') ' -- SASICA ' num2str(ifig)]);
            % find the ok button and change its callback fcn
            okbutt = findobj(hfig(ifig),'string','OK');
            set(okbutt,'callback',['delete(findobj(''-regexp'',''name'',''pop_selectcomps.* -- SASICA''));delete(findobj(''-regexp'',''name'',''Automatic component rejection measures''));' ...
                'if exist(''ALLEEG'',''var'') && exist(''EEG'',''var'') && exist(''CURRENTSET'',''var''); [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG,CURRENTSET); if not(isempty(findobj(''-regexp'',''name'',''^EEGLAB''))); eeglab(''redraw'');end;end;' ...
                'warndlg({''Remember you need to now subtract the marked components.'' ''Use Tools > Remove components''});']);
            % find the cancel button and change its callback fcn
            cancelbutt = findobj(hfig(ifig),'string','Cancel');
            closecallback = ['try; delete(findobj(''-regexp'',''name'',''pop_selectcomps''));delete(findobj(''-regexp'',''name'',''Automatic component rejection measures''));end;'];
            set(cancelbutt,'callback',[closecallback 'EEG.reject.gcompreject = false(size(EEG.reject.gcompreject));disp(''Operation cancelled. No component is selected.'');']);
            set(hfig(ifig),'closerequestfcn',closecallback)
            % crazy thing to find and order the axes for the topos.
            ax{ifig} = findobj(hfig(ifig),'type','Axes');
            ax{ifig} = ax{ifig}(end-1:-1:1);% erase pointer to the big axis behind all others and reorder the axes handles.
        end;
        ax = vertcat(ax{:});
        
        if not(numel(ax) == ncomp) || isempty(okbutt) || ~ishandle(okbutt)
            errordlg('Please do not click while I''m drawing these topos, it''s disturbing. Start over again...')
            error('Please do not click while I''m drawing these topos, it''s disturbing. Start over again...')
        end
        
        % create markers next to each topoplot showing which threshold has been
        % passed.
        for i_comp = 1:ncomp
            if EEG.reject.gcompreject(i_comp)
                %                 axes(ax(i_comp))
                f = get(ax(i_comp),'parent');
                set(0,'currentFigure',f);
                set(f,'CurrentAxes',ax(i_comp));
                drawnow;
                hold on
                for irej = 1:size(rejfields,1)
                    if isfield(EEG.reject.SASICA,rejfields{irej,1}) && ...
                            EEG.reject.SASICA.(rejfields{irej,1})(i_comp)
                        x = -.5 + (irej > 6);
                        y = .5 - .1*irej-.3*(rem(irej-1,6)+1>3);
                        h = plot(x,y,'markerfacecolor',EEG.reject.SASICA.([rejfields{irej} 'col']),'markeredgecolor',EEG.reject.SASICA.([rejfields{irej} 'col']),'marker','o');
                    end
                end
            end
        end
        set(hfig,'visible','on');
        try
            eeg_SASICA(EEG,['pop_selectcomps(EEG, [' num2str(ncomp+1) ']);']);
        end
        textprogressbar;
        hlastfig = gcf;
        set(hlastfig,'name',[get(hlastfig,'name') ' -- SASICA']);
        lastax = findobj(hlastfig,'type','Axes');
        set(lastax,'visible','off');
        axes(lastax);
        hold on
        for irej = 1:numel(rejects)
            set(gca,'xlimmode','manual');
            if rejects(irej)
                x = 0;
                y = .5 - .2*irej;
                
                scatter(x,y,'markerfacecolor',EEG.reject.SASICA.([rejfields{irej} 'col']),'markeredgecolor',EEG.reject.SASICA.([rejfields{irej} 'col']));
                text(x+.1,y,[rejfields{irej,2} ' (' num2str(sum(EEG.reject.SASICA.(rejfields{irej,1}))) ')']);
            end
        end
        for i = numel(hfig):-1:1
            figure(hfig(i));
            setctxt(hfig(i),EEG,cfg);
        end
        figure(hlastfig);
    else
        disp('No channel locations. I''m not plotting.');
    end
end
if nargout == 0
    assignin('caller','EEG',EEG);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = existnotempty(s,f)
res = isfield(s,f) && not(isempty(s.(f)));

function yticklabel(labels)

set(gca,'yticklabel',labels)

function alpha = msk2alpha(msk,m,M)

% replace in msk all 0 by m and all 1 by M
if not(islogical(msk)) || isequal(unique(msk),[0,1])
    error('Cannot deal with non binary msk')
end
alpha = msk;
alpha(alpha==0) = m;
alpha(alpha==1) = M;

function setctxt(hfig,EEG,cfg)
COLREJ = '[1 0.6 0.6]';
COLACC = '[0.75 1 0.75]';
buttons = findobj(hfig,'-regexp','tag','^comp\d{1,3}$');
buttonnums = regexp(get(buttons,'tag'),'comp(\d{1,3})','tokens');
if numel(buttonnums)>1
    buttonnums = cellfun(@(x)(str2num(x{1}{1})),buttonnums);
else
    buttonnums = str2num(buttonnums{1}{1});
end
for i = 1:numel(buttonnums)
    hcmenu = uicontextmenu;
    
    if ~isempty(EEG.reject.gcompreject)
        status = EEG.reject.gcompreject(buttonnums(i));
    else
        status = 0;
    end;
    
    hcb1 = ['EEG.reject.gcompreject(' num2str(buttonnums(i)) ') = ~EEG.reject.gcompreject(' num2str(buttonnums(i)) ');'...
        'set(gco,''backgroundcolor'',fastif(EEG.reject.gcompreject(' num2str(buttonnums(i)) '), ' COLREJ ',' COLACC '));'...
        'set(findobj(''tag'',''ctxt' num2str(buttonnums(i)) '''), ''Label'',fastif(EEG.reject.gcompreject(' num2str(buttonnums(i)) '),''ACCEPT'',''REJECT''));' ];
    uimenu(hcmenu, 'Label', fastif(status,'ACCEPT','REJECT'), 'Callback', hcb1,'tag',['ctxt' num2str(buttonnums(i))]);
    
    mycb = strrep(get(buttons(i),'Callback'),'''','''''');
    mycb = regexprep(mycb,'pop_prop','eeg_SASICA(EEG,''pop_prop');
    mycb = [mycb ''');'];
    set(buttons(i),'CallBack',mycb)
    set(buttons(i),'uicontextmenu',hcmenu)
end

% pop_prop() - plot the properties of a channel or of an independent
%              component.
% Usage:
%   >> pop_prop( EEG);           % pops up a query window
%   >> pop_prop( EEG, typecomp); % pops up a query window
%   >> pop_prop( EEG, typecomp, chanorcomp, winhandle,spectopo_options);
%
% Inputs:
%   EEG        - EEGLAB dataset structure (see EEGGLOBAL)
%
% Optional inputs:
%   typecomp   - [0|1] 1 -> display channel properties
%                0 -> component properties {default: 1 = channel}
%   chanorcomp - channel or component number[s] to display {default: 1}
%
%   winhandle  - if this parameter is present or non-NaN, buttons
%                allowing the rejection of the component are drawn.
%                If non-zero, this parameter is used to back-propagate
%                the color of the rejection button.
%   spectopo_options - [cell array] optional cell arry of options for
%                the spectopo() function.
%                For example { 'freqrange' [2 50] }
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: pop_runica(), eeglab()

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

% hidden parameter winhandle

% 01-25-02 reformated help & license -ad
% 02-17-02 removed event index option -ad
% 03-17-02 debugging -ad & sm
% 03-18-02 text settings -ad & sm
% 03-18-02 added title -ad & sm

function pop_prop(EEG, typecomp, chanorcomp, winhandle, spec_opt)


% assumed input is chanorcomp
% -------------------------
try, icadefs;
catch,
    BACKCOLOR = [0.8 0.8 0.8];
    GUIBUTTONCOLOR   = [0.8 0.8 0.8];
end;
basename = ['Component ' int2str(chanorcomp) ];

fh = figure('name', ['pop_prop() - ' basename ' properties'], 'color', BACKCOLOR, 'numbertitle', 'off', 'visible', 'on');
pos = get(gcf,'position');
set(gcf,'Position', [pos(1) pos(2)-700+pos(4) 500 700], 'visible', 'on');
pos = get(gca,'position'); % plot relative to current axes
hh = gca;
q = [pos(1) pos(2) 0 0];
s = [pos(3) pos(4) pos(3) pos(4)]./100;
delete(gca);
p = panel();
p.margin = [10 10 10 10];
p.pack('v',{.35 []});
p(1).margin = [0 0 0 0];
p(1).pack('h',{.4 [] .01});


% plotting topoplot
p(1,1).select()
topoplot( EEG.icawinv(:,chanorcomp), EEG.chanlocs, 'chaninfo', EEG.chaninfo, ...
    'shading', 'interp', 'numcontour', 3); axis square;
title(basename, 'fontsize', 14);

% plotting erpimage
p(1,2).margin = [15 15 5 15];
p(1,2).select();
eeglab_options;
if EEG.trials > 1
    % put title at top of erpimage
    axis off
    EEG.times = linspace(EEG.xmin, EEG.xmax, EEG.pnts);
    if EEG.trials < 6
        ei_smooth = 1;
    else
        ei_smooth = 3;
    end
    icaacttmp = eeg_getica(EEG, chanorcomp);
    offset = nan_mean(icaacttmp(:));
    era    = nan_mean(squeeze(icaacttmp)')-offset;
    era_limits=get_era_limits(era);
    erpimage( icaacttmp-offset, ones(1,EEG.trials)*10000, EEG.times*1000, ...
        '', ei_smooth, 1, 'caxis', 2/3, 'cbar','erp', 'yerplabel', '','erp_vltg_ticks',era_limits);
    title(sprintf('%s activity \\fontsize{10}(global offset %3.3f)', basename, offset));
else
    % put title at top of erpimage
    EI_TITLE = 'Continous data';
    ERPIMAGELINES = 200; % show 200-line erpimage
    while size(EEG.data,2) < ERPIMAGELINES*EEG.srate
        ERPIMAGELINES = 0.9 * ERPIMAGELINES;
    end
    ERPIMAGELINES = round(ERPIMAGELINES);
    if ERPIMAGELINES > 2   % give up if data too small
        if ERPIMAGELINES < 10
            ei_smooth = 1;
        else
            ei_smooth = 3;
        end
        erpimageframes = floor(size(EEG.data,2)/ERPIMAGELINES);
        erpimageframestot = erpimageframes*ERPIMAGELINES;
        eegtimes = linspace(0, erpimageframes-1, EEG.srate/1000);
        if typecomp == 1 % plot channel
            offset = nan_mean(EEG.data(chanorcomp,:));
            % Note: we don't need to worry about ERP limits, since ERPs
            % aren't visualized for continuous data
            erpimage( reshape(EEG.data(chanorcomp,1:erpimageframestot),erpimageframes,ERPIMAGELINES)-offset, ones(1,ERPIMAGELINES)*10000, eegtimes , ...
                EI_TITLE, ei_smooth, 1, 'caxis', 2/3, 'cbar');
        else % plot component
            icaacttmp = eeg_getdatact(EEG, 'component', chanorcomp);
            offset = nan_mean(icaacttmp(:));
            erpimage(reshape(icaacttmp(:,1:erpimageframestot),erpimageframes,ERPIMAGELINES)-offset,ones(1,ERPIMAGELINES)*10000, eegtimes , ...
                EI_TITLE, ei_smooth, 1, 'caxis', 2/3, 'cbar','yerplabel', '');
        end
    else
        axis off;
        text(0.1, 0.3, [ 'No erpimage plotted' 10 'for small continuous data']);
    end;
end;

% plotting spectrum
% -----------------
if ~exist('winhandle') || isempty(winhandle) || ~ishandle(winhandle)
    winhandle = NaN;
    p(2).pack('v',{.3 [] })
else
    p(2).pack('v',{.3 [] .1})
end;
p(2,1).pack('h',{.01,[],.01});
p(2,1).margin = [15 15 0 55];
p(2,1,1).margin = 0;
p(2,1,3).margin = 0;
p(2,1,2).pack('v',{.01 []});
p(2,1,2,1).margin = 0;
p(2,1,2,2).margintop = 5;
p(2,1,2,2).select();
try
    spectopo( EEG.icaact(chanorcomp,:), EEG.pnts, EEG.srate, 'mapnorm', EEG.icawinv(:,chanorcomp), spec_opt{:} );
    % set( get(gca, 'ylabel'), 'string', 'Power 10*log_{10}(\muV^{2}/Hz)', 'fontsize', 12);
    % set( get(gca, 'xlabel'), 'string', 'Frequency (Hz)', 'fontsize', 12);
    axis on
    xlabel('Frequency (Hz)')
    h = title('Activity power spectrum', 'fontsize', 10);
    %     set(h,'position',get(h,'position')+[-15 -7 0]);
    set(gca,'fontSize',10)
catch
    axis off;
    lasterror
    text(0.1, 0.3, [ 'Error: no spectrum plotted' 10 ' make sure you have the ' 10 'signal processing toolbox']);
end;

%%%% Add SASICA measures.
%               eye      muscle/noise    channel     ~ok
colors = { [0 .75 .75]      [0 0 1]      [0 .5 0] [.2 .2 .2]};
% C={[1 0 0],[.6 0 .2],[1 1 0],[0 1 0], [0 1 1]};% colors used in ADJ
computed = fieldnames(EEG.reject.SASICA);
computed = computed(regexpcell(computed,'rej|thresh|^var$','inv'));
computedthresh = regexprep(computed,'ica','icathresh');
computedrej = regexprep(computed,'ica','icarej');
toPlot = {};
toPlot_axprops = {};
toPlot_title = {}; SXticks = {};co = [];
for i = 1:numel(computed)
    if strcmp(computed{i},'icaADJUST')
        struct2ws(EEG.reject.SASICA.icaADJUST)
        toPlot{end+1}{1} = (SAD(chanorcomp)-med2_SAD)/(soglia_SAD-med2_SAD);
        toPlot{end}{2} = (SED(chanorcomp)-med2_SED)/(soglia_SED-med2_SED);
        toPlot{end}{3} = (GDSF(chanorcomp)-med2_GDSF)/(soglia_GDSF-med2_GDSF);
        toPlot{end}{4} = (nuovaV(chanorcomp)-med2_V)/(soglia_V-med2_V);
        toPlot{end}{5} = (meanK(chanorcomp)-med2_K)/(soglia_K-med2_K);
        ADJis = '';
        aco = repmat(colors{4},numel(toPlot{end}),1);
        if ismember(chanorcomp,horiz)
            ADJis = [ADJis 'HEM/'];
            aco([2 4],:) = repmat(colors{1},2,1);
        end
        if ismember(chanorcomp,vert)
            ADJis = [ADJis 'VEM/'];
            aco([1 4],:) = repmat(colors{1},2,1);
        end
        if ismember(chanorcomp,blink)
            ADJis = [ADJis 'Blink/'];
            aco([1 4 5],:) = repmat(colors{1},3,1);
        end
        if ismember(chanorcomp,disc)
            ADJis = [ADJis 'Disc/'];
            aco([3 4],:) = repmat(colors{3},2,1);
        end
        if isempty(ADJis)
            ADJis = 'OK';
        else
            ADJis(end) = [];
        end
        toPlot_title{end+1} = ['ADJUST: ' ADJis];
        toPlot_axprops{end+1} = {'ColorOrder' aco,...
            'ylim' [0 2],...
            'ytick' [1 2],...
            'yticklabel' {'Th' '2*Th'},...
            'xtick' 1:numel(toPlot{end}),...
            'xticklabel' {'SAD' 'SED' 'GDSF' 'MEV' 'TK'}};
    elseif strcmp(computed{i},'icaFASTER')
        listprops = EEG.reject.SASICA.icaFASTER.listprops;
        str='FASTER: ';
        FASTER_reasons = {'HighFreq ' 'FlatSpectrum ' 'SpatialKurtosis ' 'HurstExponent ' 'EOGCorrel '};
        %                     1 Median gradient value, for high frequency stuff
        %                     2 Mean slope around the LPF band (spectral)
        %                     3 Kurtosis of spatial map
        %                     4 Hurst exponent
        %                     5 Eyeblink correlations
        zlist = zscore(listprops);
        for i = 1:size(listprops,2)
            fst(:,i) = min_z(listprops(:,i));
        end
        reasons = FASTER_reasons(fst(chanorcomp,:));
        if isempty(reasons)
            str = [str 'OK'];
        else
            str = [str reasons{:}];
        end
        FSTis = str;
        toPlot{end+1} = {};
        for ip = 1:numel(zlist(chanorcomp,:))
            toPlot{end}{ip} = abs(zlist(chanorcomp,ip))/3;% normalized by threshold
        end
        toPlot_title{end+1} = FSTis;
        toPlot_axprops{end+1} = {'ColorOrder' [colors{2};colors{2};colors{3};colors{2};colors{1}],...
            'ylim' [0 2],...
            'ytick' [1 2],...
            'yticklabel' {'Th' '2*Th'},...
            'xtick',1:numel(toPlot{end}),...
            'xticklabel',{'MedGrad' 'SpecSl' 'SK' 'HE' 'EOGCorr'}};
    elseif strcmp(computed{i},'icaMARA')
        info = EEG.reject.SASICA.icaMARA.info;
        str='MARA: ';
        MARA_meas = {'CurrDensNorm ' 'SpatRange ' 'AvgLocSkew ' '\lambda ' '8-13 Pow' '1/F Fit '};
        %                     1 Current Density Norm
        %                     2 Spatial Range
        %                     3 Average Local Skewness
        %                     4 lambda
        %                     5 Band Power (8-13 Hz)
        %                     6 Fit Error
        if ~ EEG.reject.SASICA.icarejMARA(chanorcomp)
            str = [str 'OK       '];
        else
            str = [str 'Reject    '];
        end
        MARAis = [str '(' num2str(round(100*info.posterior_artefactprob(chanorcomp)),'%g') '%)'];
        toPlot{end+1} = {};
        for ip = 1:numel(info.normfeats(:,chanorcomp))
            toPlot{end}{ip} = info.normfeats(ip,chanorcomp) ;
        end
        toPlot_title{end+1} = MARAis;
        toPlot_axprops{end+1} = {'ColorOrder' repmat(colors{4},numel(MARA_meas),1),...
            'ylimmode' 'auto',...
            'xtick',1:numel(toPlot{end}),...
            'xticklabel',{'CDN' 'SpRg' 'AvLocSkw' 'lambda' '8-13 Hz' '1/F Fit'}
            };
    else
        rejfields = {
            'icaautocorr'       'LoAC'   colors{2}
            'icafocalcomp'      'FocCh'       colors{3}
            'icatrialfoc'       'FocTr'        colors{3}
            'icaSNR'            'LoSNR'        colors{2}
            'icaresvar'         'ResV'          colors{2}
            'icachancorrVEOG'   'CorrV'         colors{1}
            'icachancorrHEOG'   'CorrH'         colors{1}
            'icachancorrchans'  'CorrC'         colors{3}
            };
        if isempty(toPlot)
            toPlot{1} = {};
            toPlot_axprops{1} = {};
            toPlot_title{1} = 'SASICA';
        end
        switch computed{i}
            case 'icaautocorr'
                toPlot{1}{end+1} = 2 - (EEG.reject.SASICA.(computed{i})(chanorcomp) +1)/(EEG.reject.SASICA.(computedthresh{i}) +1);
            case 'icaSNR'
                toPlot{1}{end+1} = EEG.reject.SASICA.(computedthresh{i})/EEG.reject.SASICA.(computed{i})(chanorcomp);
            otherwise
                toPlot{1}{end+1} = EEG.reject.SASICA.(computed{i})(:,chanorcomp)/EEG.reject.SASICA.(computedthresh{i});
        end
        SXticks{end+1} = rejfields{strcmp(computed{i},rejfields(:,1)),2};
        co(end+1,:) = rejfields{strcmp(computed{i},rejfields(:,1)),3};
    end
end
if not(isempty(SXticks))
    toPlot_axprops{1} = {toPlot_axprops{1}{:} 'ylim' [0 2]...
        'ytick' [1 2] ...
        'yticklabel' {'Th' '2*Th'} ...
        'xtick' 1:numel(SXticks) ...
        'Xticklabel' SXticks...
        'xlim',[.5 numel(SXticks)+.5],...
        'colororder',co};
end

p(2,2).pack('v',numel(toPlot));
p(2,2).de.margintop = 0;
for i = 1:numel(toPlot)
    p(2,2,i).pack('h',{.2 []});
    p(2,2,i,1).select();
    text(1.1,0.5,strjust(strwrap(toPlot_title{i},15),'right'),'horizontalalignment','right');
    axis off
    p(2,2,i,2).select()
    hold on
    set(gca,toPlot_axprops{i}{:});
    cs = get(gca,'colorOrder');
    for j = 1:numel(toPlot{i})
        xj = linspace(j-(numel(toPlot{i}{j})>1)*.3,j+(numel(toPlot{i}{j})>1)*.3,numel(toPlot{i}{j}));
        bar(xj,toPlot{i}{j},'facecolor',cs(rem(j-1,size(cs,1))+1,:));
    end
    hline(1,':k')
end


% display buttons
% ---------------
if ishandle(winhandle)
    COLREJ = '[1 0.6 0.6]';
    COLACC = '[0.75 1 0.75]';
    % CANCEL button
    % -------------
    h  = uicontrol(gcf, 'Style', 'pushbutton', 'backgroundcolor', GUIBUTTONCOLOR, 'string', 'Cancel', 'Units','Normalized','Position',[-10 -10 15 6].*s+q, 'callback', 'close(gcf);');
    
    %     % VALUE button
    %     % -------------
    %     hval  = uicontrol(gcf, 'Style', 'pushbutton', 'backgroundcolor', GUIBUTTONCOLOR, 'string', 'Values', 'Units','Normalized', 'Position', [15 -10 15 6].*s+q);
    
    % REJECT button
    % -------------
    if ~isempty(EEG.reject.gcompreject)
        status = EEG.reject.gcompreject(chanorcomp);
    else
        status = 0;
    end;
    hr = uicontrol(gcf, 'Style', 'pushbutton', 'backgroundcolor', eval(fastif(status,COLREJ,COLACC)), ...
        'string', fastif(status, 'REJECT', 'ACCEPT'), 'Units','Normalized', 'Position', [40 -10 15 6].*s+q, 'userdata', status, 'tag', 'rejstatus');
    command = [ 'set(gcbo, ''userdata'', ~get(gcbo, ''userdata''));' ...
        'if get(gcbo, ''userdata''),' ...
        '     set( gcbo, ''backgroundcolor'',' COLREJ ', ''string'', ''REJECT'');' ...
        'else ' ...
        '     set( gcbo, ''backgroundcolor'',' COLACC ', ''string'', ''ACCEPT'');' ...
        'end;' ];
    set( hr, 'callback', command);
    
    %     % HELP button
    %     % -------------
    %     h  = uicontrol(gcf, 'Style', 'pushbutton', 'backgroundcolor', GUIBUTTONCOLOR, 'string', 'HELP', 'Units','Normalized', 'Position', [65 -10 15 6].*s+q, 'callback', 'pophelp(''pop_prop'');');
    
    % OK button
    % ---------
    command = [ 'global EEG;' ...
        'tmpstatus = get( findobj(''parent'', gcbf, ''tag'', ''rejstatus''), ''userdata'');' ...
        'EEG.reject.gcompreject(' num2str(chanorcomp) ') = tmpstatus;' ];
    if winhandle ~= 0
        command = [ command ...
            sprintf('if tmpstatus set(gcbo, ''backgroundcolor'', %s); else set(gcbo, ''backgroundcolor'', %s); end;', ...
            COLREJ, COLACC) ...
            ['obj = findobj(''-regexp'',''name'',''pop_selectcomps.* -- SASICA''); obj = fastif(isempty(obj),[],findobj(obj,''tag'',''comp' num2str(chanorcomp) '''));'] ...
            sprintf('if ~isempty(obj) && tmpstatus set(obj, ''backgroundcolor'', %s); else set(obj, ''backgroundcolor'', %s); end;', ...
            COLREJ, COLACC)];
    end;
    command = [ command 'close(gcf); clear tmpstatus' ];
    h  = uicontrol(gcf, 'Style', 'pushbutton', 'string', 'OK', 'backgroundcolor', GUIBUTTONCOLOR, 'Units','Normalized', 'Position',[90 -10 15 6].*s+q, 'callback', command);
    
    %     % draw the figure for statistical values
    %     % --------------------------------------
    %     index = num2str( chanorcomp );
    %     command = [ ...
    %         'figure(''MenuBar'', ''none'', ''name'', ''Statistics of the component'', ''numbertitle'', ''off'');' ...
    %         '' ...
    %         'pos = get(gcf,''Position'');' ...
    %         'set(gcf,''Position'', [pos(1) pos(2) 340 340]);' ...
    %         'pos = get(gca,''position'');' ...
    %         'q = [pos(1) pos(2) 0 0];' ...
    %         's = [pos(3) pos(4) pos(3) pos(4)]./100;' ...
    %         'axis off;' ...
    %         ''  ...
    %         'txt1 = sprintf(''(\n' ...
    %         'Entropy of component activity\t\t%2.2f\n' ...
    %         '> Rejection threshold \t\t%2.2f\n\n' ...
    %         ' AND                 \t\t\t----\n\n' ...
    %         'Kurtosis of component activity\t\t%2.2f\n' ...
    %         '> Rejection threshold \t\t%2.2f\n\n' ...
    %         ') OR                  \t\t\t----\n\n' ...
    %         'Kurtosis distibution \t\t\t%2.2f\n' ...
    %         '> Rejection threhold\t\t\t%2.2f\n\n' ...
    %         '\n' ...
    %         'Current thesholds sujest to %s the component\n\n' ...
    %         '(after manually accepting/rejecting the component, you may recalibrate thresholds for future automatic rejection on other datasets)'',' ...
    %         'EEG.stats.compenta(' index '), EEG.reject.threshentropy, EEG.stats.compkurta(' index '), ' ...
    %         'EEG.reject.threshkurtact, EEG.stats.compkurtdist(' index '), EEG.reject.threshkurtdist, fastif(EEG.reject.gcompreject(' index '), ''REJECT'', ''ACCEPT''));' ...
    %         '' ...
    %         'uicontrol(gcf, ''Units'',''Normalized'', ''Position'',[-11 4 117 100].*s+q, ''Style'', ''frame'' );' ...
    %         'uicontrol(gcf, ''Units'',''Normalized'', ''Position'',[-5 5 100 95].*s+q, ''String'', txt1, ''Style'',''text'', ''HorizontalAlignment'', ''left'' );' ...
    %         'h = uicontrol(gcf, ''Style'', ''pushbutton'', ''string'', ''Close'', ''Units'',''Normalized'', ''Position'', [35 -10 25 10].*s+q, ''callback'', ''close(gcf);'');' ...
    %         'clear txt1 q s h pos;' ];
    %     set( hval, 'callback', command);
    %     if isempty( EEG.stats.compenta )
    %         set(hval, 'enable', 'off');
    %     end;
    
    com = sprintf('pop_prop( %s, %d, %d, 0, %s);', inputname(1), typecomp, chanorcomp, vararg2str( { spec_opt } ) );
else
    com = sprintf('pop_prop( %s, %d, %d, NaN, %s);', inputname(1), typecomp, chanorcomp, vararg2str( { spec_opt } ) );
end;

return;

% pop_selectcomps() - Display components with button to vizualize their
%                  properties and label them for rejection.
% Usage:
%       >> OUTEEG = pop_selectcomps( INEEG, compnum );
%
% Inputs:
%   INEEG    - Input dataset
%   compnum  - vector of component numbers
%
% Output:
%   OUTEEG - Output dataset with updated rejected components
%
% Note:
%   if the function POP_REJCOMP is ran prior to this function, some
%   fields of the EEG datasets will be present and the current function
%   will have some more button active to tune up the automatic rejection.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: pop_prop(), eeglab()

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

function [EEG, com] = pop_selectcomps( EEG, compnum, comptot );
if not(exist('comptot','var'))
    comptot = max(compnum);
end
COLREJ = '[1 0.6 0.6]';
COLACC = '[0.75 1 0.75]';
PLOTPERFIG = 35;

com = '';
if nargin < 1
    help pop_selectcomps;
    return;
end;

if nargin < 2
    promptstr = { 'Components to plot:' };
    initstr   = { [ '1:' int2str(size(EEG.icaweights,1)) ] };
    
    result = inputdlg2(promptstr, 'Reject comp. by map -- pop_selectcomps',1, initstr);
    if isempty(result), return; end;
    compnum = eval( [ '[' result{1} ']' ]);
    
    if length(compnum) > PLOTPERFIG
        ButtonName=questdlg2(strvcat(['More than ' int2str(PLOTPERFIG) ' components so'],'this function will pop-up several windows'), ...
            'Confirmation', 'Cancel', 'OK','OK');
        if ~isempty( strmatch(lower(ButtonName), 'cancel')), return; end;
    end;
    
end;
% fprintf('Drawing figure...\n');
currentfigtag = ['selcomp' num2str(rand)]; % generate a random figure tag

if length(compnum) > PLOTPERFIG
    for index = 1:PLOTPERFIG:length(compnum)
        pop_selectcomps(EEG, compnum([index:min(length(compnum),index+PLOTPERFIG-1)]));
    end;
    
    com = [ 'pop_selectcomps(' inputname(1) ', ' vararg2str(compnum) ');' ];
    return;
end;

if isempty(EEG.reject.gcompreject)
    EEG.reject.gcompreject = zeros( size(EEG.icawinv,2));
end;
try, icadefs;
catch,
    BACKCOLOR = [0.8 0.8 0.8];
    GUIBUTTONCOLOR   = [0.8 0.8 0.8];
end;

% set up the figure
% -----------------
column =ceil(sqrt( length(compnum) ))+1;
rows = ceil(length(compnum)/column);
if ~exist('fig')
    figure('name', [ 'Reject components by map - pop_selectcomps() (dataset: ' EEG.setname ')'], 'tag', currentfigtag, ...
        'numbertitle', 'off', 'color', BACKCOLOR,'visible','off');
    set(gcf,'MenuBar', 'none');
    pos = get(gcf,'Position');
    set(gcf,'Position', [pos(1) 20 800/7*column 600/5*rows]);
    incx = 120;
    incy = 110;
    sizewx = 100/column;
    if rows > 2
        sizewy = 90/rows;
    else
        sizewy = 80/rows;
    end;
    pos = get(gca,'position'); % plot relative to current axes
    hh = gca;
    q = [pos(1) pos(2) 0 0];
    s = [pos(3) pos(4) pos(3) pos(4)]./100;
    axis off;
end;

% figure rows and columns
% -----------------------
if EEG.nbchan > 64
    %     disp('More than 64 electrodes: electrode locations not shown');
    plotelec = 0;
else
    plotelec = 1;
end;
count = 1;
for ri = compnum
    if ri > numel(EEG.icachansind)
        error('don''t panic')
    end
    textprogressbar(ri/comptot*100);
    if exist('fig')
        button = findobj('parent', fig, 'tag', ['comp' num2str(ri)]);
        if isempty(button)
            error( 'pop_selectcomps(): figure does not contain the component button');
        end;
    else
        button = [];
    end;
    
    if isempty( button )
        % compute coordinates
        % -------------------
        X = mod(count-1, column)/column * incx-10;
        Y = (rows-floor((count-1)/column))/rows * incy - sizewy*1.3;
        
        % plot the head
        % -------------
        if ~strcmp(get(gcf, 'tag'), currentfigtag);
            disp('Aborting plot');
            return;
        end;
        ha = axes('Units','Normalized', 'Position',[X Y sizewx sizewy].*s+q);
        if plotelec
            topoplot( EEG.icawinv(:,ri), EEG.chanlocs, 'verbose', ...
                'off', 'style' , 'fill', 'chaninfo', EEG.chaninfo, 'numcontour', 8);
        else
            topoplot( EEG.icawinv(:,ri), EEG.chanlocs, 'verbose', ...
                'off', 'style' , 'fill','electrodes','off', 'chaninfo', EEG.chaninfo, 'numcontour', 8);
        end;
        axis square;
        
        % plot the button
        % ---------------
        button = uicontrol(gcf, 'Style', 'pushbutton', 'Units','Normalized', 'Position',...
            [X Y+sizewy sizewx sizewy*0.25].*s+q, 'tag', ['comp' num2str(ri)]);
        command = sprintf('pop_prop( %s, 0, %d, gcbo, { ''freqrange'', [1 50] });', inputname(1), ri); %RMC command = sprintf('pop_prop( %s, 0, %d, %3.15f, { ''freqrange'', [1 50] });', inputname(1), ri, button);
        set( button, 'callback', command );
    end;
    set( button, 'backgroundcolor', eval(fastif(EEG.reject.gcompreject(ri), COLREJ,COLACC)), 'string', int2str(ri));
    drawnow;
    count = count +1;
end;

% draw the bottom button
% ----------------------
if ~exist('fig')
    hh = uicontrol(gcf, 'Style', 'pushbutton', 'string', 'Cancel', 'Units','Normalized', 'backgroundcolor', GUIBUTTONCOLOR, ...
        'Position',[-10 -10  15 sizewy*0.25].*s+q, 'callback', 'close(gcf); fprintf(''Operation cancelled\n'')' );
    hh = uicontrol(gcf, 'Style', 'pushbutton', 'string', 'See projection', 'Units','Normalized', 'backgroundcolor', GUIBUTTONCOLOR, ...
        'Position',[50 -10  15 sizewy*0.25].*s+q, 'callback', ' ', 'enable', 'off'  );
    hh = uicontrol(gcf, 'Style', 'pushbutton', 'string', 'Help', 'Units','Normalized', 'backgroundcolor', GUIBUTTONCOLOR, ...
        'Position',[70 -10  15 sizewy*0.25].*s+q, 'callback', 'pophelp(''pop_selectcomps'');' );
    command = '[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET); eegh(''[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);''); close(gcf)';
    hh = uicontrol(gcf, 'Style', 'pushbutton', 'string', 'OK', 'Units','Normalized', 'backgroundcolor', GUIBUTTONCOLOR, ...
        'Position',[90 -10  15 sizewy*0.25].*s+q, 'callback',  command);
    % sprintf(['eeg_global; if %d pop_rejepoch(%d, %d, find(EEG.reject.sigreject > 0), EEG.reject.elecreject, 0, 1);' ...
    %		' end; pop_compproj(%d,%d,1); close(gcf); eeg_retrieve(%d); eeg_updatemenu; '], rejtrials, set_in, set_out, fastif(rejtrials, set_out, set_in), set_out, set_in));
end;

com = [ 'pop_selectcomps(' inputname(1) ', ' vararg2str(compnum) ');' ];
return;

function out = nan_mean(in)

nans = find(isnan(in));
in(nans) = 0;
sums = sum(in);
nonnans = ones(size(in));
nonnans(nans) = 0;
nonnans = sum(nonnans);
nononnans = find(nonnans==0);
nonnans(nononnans) = 1;
out = sum(in)./nonnans;
out(nononnans) = NaN;


function era_limits=get_era_limits(era)
%function era_limits=get_era_limits(era)
%
% Returns the minimum and maximum value of an event-related
% activation/potential waveform (after rounding according to the order of
% magnitude of the ERA/ERP)
%
% Inputs:
% era - [vector] Event related activation or potential
%
% Output:
% era_limits - [min max] minimum and maximum value of an event-related
% activation/potential waveform (after rounding according to the order of
% magnitude of the ERA/ERP)

mn=min(era);
mx=max(era);
mn=orderofmag(mn)*round(mn/orderofmag(mn));
mx=orderofmag(mx)*round(mx/orderofmag(mx));
era_limits=[mn mx];


function ord=orderofmag(val)
%function ord=orderofmag(val)
%
% Returns the order of magnitude of the value of 'val' in multiples of 10
% (e.g., 10^-1, 10^0, 10^1, 10^2, etc ...)
% used for computing erpimage trial axis tick labels as an alternative for
% plotting sorting variable

val=abs(val);
if val>=1
    ord=1;
    val=floor(val/10);
    while val>=1,
        ord=ord*10;
        val=floor(val/10);
    end
    return;
else
    ord=1/10;
    val=val*10;
    while val<1,
        ord=ord/10;
        val=val*10;
    end
    return;
end

function thresh = readauto(thresh,dat,comp)
% if thresh starts with 'auto'
% compute auto threshold as mean(dat) +/- N std(dat)
% with N read in the string thresh = 'auto N'
% if not, use thresh as a value
if isstr(thresh) && strncmp(thresh,'auto',4)
    if numel(thresh) > 4
        threshsigma = str2num(thresh(5:end));
    else
        threshsigma = 2;
    end
    thresh = eval(['mean(dat,2)' comp 'threshsigma * std(dat,[],2)']);
end



function [nb,channame,strnames] = chnb(channame, varargin)

% chnb() - return channel number corresponding to channel names in an EEG
%           structure
%
% Usage:
%   >> [nb]                 = chnb(channameornb);
%   >> [nb,names]           = chnb(channameornb,...);
%   >> [nb,names,strnames]  = chnb(channameornb,...);
%   >> [nb]                 = chnb(channameornb, labels);
%
% Input:
%   channameornb  - If a string or cell array of strings, it is assumed to
%                   be (part of) the name of channels to search. Either a
%                   string with space separated channel names, or a cell
%                   array of strings.
%                   Note that regular expressions can be used to match
%                   several channels. See regexp.
%                   If only one channame pattern is given and the string
%                   'inv' is attached to it, the channels NOT matching the
%                   pattern are returned.
%   labels        - Channel names as found in {EEG.chanlocs.labels}.
%
% Output:
%   nb            - Channel numbers in labels, or in the EEG structure
%                   found in the caller workspace (i.e. where the function
%                   is called from) or in the base workspace, if no EEG
%                   structure exists in the caller workspace.
%   names         - Channel names, cell array of strings.
%   strnames      - Channel names, one line character array.
error(nargchk(1,2,nargin));
if nargin == 2
    labels = varargin{1};
else
    
    try
        EEG = evalin('caller','EEG');
    catch
        try
            EEG = evalin('base','EEG');
        catch
            error('Could not find EEG structure');
        end
    end
    if not(isfield(EEG,'chanlocs'))
        error('No channel list found');
    end
    EEG = EEG(1);
    labels = {EEG.chanlocs.labels};
end
if iscell(channame) || ischar(channame)
    
    if ischar(channame) || iscellstr(channame)
        if iscellstr(channame) && numel(channame) == 1 && isempty(channame{1})
            channame = '';
        end
        tmp = regexp(channame,'(\S*) ?','tokens');
        channame = {};
        for i = 1:numel(tmp)
            if iscellstr(tmp{i}{1})
                channame{i} = tmp{i}{1}{1};
            else
                channame{i} = tmp{i}{1};
            end
        end
        if isempty(channame)
            nb = [];
            return
        end
    end
    if numel(channame) == 1 && not(isempty(strmatch('inv',channame{1})))
        cmd = 'exactinv';
        channame{1} = strrep(channame{1},'inv','');
    else
        channame{1} = channame{1};
        cmd = 'exact';
    end
    nb = regexpcell(labels,channame,[cmd 'ignorecase']);
    
elseif isnumeric(channame)
    nb = channame;
    if nb > numel(labels)
        nb = [];
    end
end
channame = labels(nb);
strnames = sprintf('%s ',channame{:});
if not(isempty(strnames))
    strnames(end) = [];
end

function idx = regexpcell(c,pat, cmds)

% idx = regexpcell(c,pat, cmds)
%
% Return indices idx of cells in c that match pattern(s) pat (regular expression).
% Pattern pat can be char or cellstr. In the later case regexpcell returns
% indexes of cells that match any pattern in pat.
%
% cmds is a string that can contain one or several of these commands:
% 'inv' return indexes that do not match the pattern.
% 'ignorecase' will use regexpi instead of regexp
% 'exact' performs an exact match (regular expression should match the whole strings in c).
% 'all' (default) returns all indices, including repeats (if several pat match a single cell in c).
% 'unique' will return unique sorted indices.
% 'intersect' will return only indices in c that match ALL the patterns in pat.
%
% v1 Maximilien Chaumon 01/05/09
% v1.1 Maximilien Chaumon 24/05/09 - added ignorecase
% v2 Maximilien Chaumon 02/03/2010 changed input method.
%       inv,ignorecase,exact,combine are replaced by cmds

error(nargchk(2,3,nargin))
if not(iscellstr(c))
    error('input c must be a cell array of strings');
end
if nargin == 2
    cmds = '';
end
if not(isempty(regexpi(cmds,'inv', 'once' )))
    inv = true;
else
    inv = false;
end
if not(isempty(regexpi(cmds,'ignorecase', 'once' )))
    ignorecase = true;
else
    ignorecase = false;
end
if not(isempty(regexpi(cmds,'exact', 'once' )))
    exact = true;
else
    exact = false;
end
if not(isempty(regexpi(cmds,'unique', 'once' )))
    combine = 2;
elseif not(isempty(regexpi(cmds,'intersect', 'once' )))
    combine = 3;
else
    combine = 1;
end

if ischar(pat)
    pat = cellstr(pat);
end

if exact
    for i_pat = 1:numel(pat)
        pat{i_pat} = ['^' pat{i_pat} '$'];
    end
end

for i_pat = 1:length(pat)
    if ignorecase
        trouv = regexpi(c,pat{i_pat}); % apply regexp on each pattern
    else
        trouv = regexp(c,pat{i_pat}); % apply regexp on each pattern
    end
    idx{i_pat} = [];
    for i = 1:numel(trouv)
        if not(isempty(trouv{i}))% if there is a match, store index
            idx{i_pat}(end+1) = i;
        end
    end
end
switch combine
    case 1
        idx = [idx{:}];
    case 2
        idx = unique([idx{:}]);
    case 3
        for i_pat = 2:length(pat)
            idx{1} = intersect(idx{1},idx{i_pat});
        end
        idx = idx{1};
end
if inv % if we want to invert result, then do so.
    others = 1:numel(trouv);
    others(idx) = [];
    idx = others;
end

function s = setdef(s,d)
% s = setdef(s,d)
% Merges the two structures s and d recursively.
% Adding the default field values from d into s when not present or empty.

if isstruct(s) && not(isempty(s))
    fields = fieldnames(d);
    for i_f = 1:numel(fields)
        if isfield(s,fields{i_f})
            s.(fields{i_f}) = setdef(s.(fields{i_f}),d.(fields{i_f}));
        else
            s.(fields{i_f}) = d.(fields{i_f});
        end
    end
elseif not(isempty(s))
    s = s;
elseif isempty(s);
    s = d;
end

function struct2ws(s,varargin)

% struct2ws(s,varargin)
%
% Description : This function returns fields of scalar structure s in the
% current workspace
% __________________________________
% Inputs :
%   s (scalar structure array) :    a structure that you want to throw in
%                                   your current workspace.
%   re (string optional) :          a regular expression. Only fields
%                                   matching re will be returned
% Outputs :
%   No output : variables are thrown directly in the caller workspace.
%
%
% _____________________________________
% See also : ws2struct ; regexp
%
% Maximilien Chaumon v1.0 02/2007


if nargin == 0
    cd('d:\Bureau\work')
    s = dir('pathdef.m');
end
if length(s) > 1
    error('Structure should be scalar.');
end
if not(isempty(varargin))
    re = varargin{1};
else
    re = '.*';
end

vars = fieldnames(s);
vmatch = regexp(vars,re);
varsmatch = [];
for i = 1:length(vmatch)
    if isempty(vmatch{i})
        continue
    end
    varsmatch(end+1) = i;
end
for i = varsmatch
    assignin('caller',vars{i},s.(vars{i}));
end

function [sortie] = ws2struct(varargin)

% [s] = ws2struct(varargin)
%
% Description : This function returns a structure containing variables
% of the current workspace.
% __________________________________
% Inputs :
%   re (string optional) :  a regular expression matching the variables to
%                           be returned.
% Outputs :
%   s (structure array) :   a structure containing all variables of the
%                           calling workspace. If re input is specified,
%                           only variables matching re are returned.
% _____________________________________
% See also : struct2ws ; regexp
%
% Maximilien Chaumon v1.0 02/2007


if not(isempty(varargin))
    re = varargin{1};
else
    re = '.*';
end

vars = evalin('caller','who');
vmatch = regexp(vars,re);
varsmatch = [];
for i = 1:length(vmatch)
    if isempty(vmatch{i}) || not(vmatch{i} == 1)
        continue
    end
    varsmatch{end+1} = vars{i};
end

for i = 1:length(varsmatch)
    dat{i} = evalin('caller',varsmatch{i});
end

sortie = cell2struct(dat,varsmatch,2);

function [tpts tvals] = timepts(timein, varargin)

% timepts() - return time points corresponding to a certain latency range
%             in an EEG structure.
%
% Usage:
%   >> [tpts] = timepts(timein);
%   >> [tpts tvals] = timepts(timein, times);
%               Note: this last method also works with any type of numeric
%               data entered under times (ex. frequencies, trials...)
%
% Input:
%   timein        - latency range [start stop] (boundaries included). If
%                   second argument 'times' is not provided, EEG.times will
%                   be evaluated from the EEG structure found in the caller
%                   workspace (or base if not in caller).
%   times         - time vector as found in EEG.times
%
% Output:
%   tpts          - index numbers corresponding to the time range.
%   tvals         - values of EEG.times at points tpts
%

error(nargchk(1,2,nargin));
if nargin == 2
    times = varargin{1};
else
    
    try
        EEG = evalin('caller','EEG');
    catch
        try
            EEG = evalin('base','EEG');
        catch
            error('Could not find EEG structure');
        end
    end
    if not(isfield(EEG,'times'))
        error('No time list found');
    end
    times = EEG.times;
    if isempty(times)
        times = EEG.xmin:1/EEG.srate:EEG.xmax;
    end
end
if isempty(times)
    error('could not find times');
end
if numel(timein) == 1
    [dum tpts] = min(abs(times - timein));% find the closest one
    if tpts == numel(times)
        warning('Strange time is last index of times')
    end
elseif numel(timein) == 2
    tpts = find(times >= timein(1) & times <= timein(2));% find times within bounds
else
    error('timein should be a scalar or a 2 elements vector');
end
tvals = times(tpts);


function tw = strwrap(t,n)

% tw = strwrap(t,n)
%
% wrap text array t at n characters taking non alphanumeric characters as
% breaking characters (i.e. not cutting words strangely).

t = deblank(t(:)');
seps = '[\s-]';
tw = '';
while not(isempty(t))
    breaks = regexp(t,seps);
    breaks(end+1) = numel(t);
    idx = 1:min(n,breaks(find(breaks < n, 1,'last')));
    if isempty(idx)
        idx = 1:min(n,numel(t));
    end
    tw(end+1,:) = char(padarray(double(t(idx)),[0 n-numel(idx)],32,'post'));
    t(idx)= [];
    t = strtrim(t);
end


function [z,mu,sigma] = zscore(x,flag,dim)
%ZSCORE Standardized z score.
%   Z = ZSCORE(X) returns a centered, scaled version of X, the same size as X.
%   For vector input X, Z is the vector of z-scores (X-MEAN(X)) ./ STD(X). For
%   matrix X, z-scores are computed using the mean and standard deviation
%   along each column of X.  For higher-dimensional arrays, z-scores are
%   computed using the mean and standard deviation along the first
%   non-singleton dimension.
%
%   The columns of Z have sample mean zero and sample standard deviation one
%   (unless a column of X is constant, in which case that column of Z is
%   constant at 0).
%
%   [Z,MU,SIGMA] = ZSCORE(X) also returns MEAN(X) in MU and STD(X) in SIGMA.
%
%   [...] = ZSCORE(X,1) normalizes X using STD(X,1), i.e., by computing the
%   standard deviation(s) using N rather than N-1, where N is the length of
%   the dimension along which ZSCORE works.  ZSCORE(X,0) is the same as
%   ZSCORE(X).
%
%   [...] = ZSCORE(X,FLAG,DIM) standardizes X by working along the dimension
%   DIM of X. Pass in FLAG==0 to use the default normalization by N-1, or 1
%   to use N.
%
%   See also MEAN, STD.

%   Copyright 1993-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2010/03/16 00:18:33 $

% [] is a special case for std and mean, just handle it out here.
if isequal(x,[]), z = []; return; end

if nargin < 2
    flag = 0;
end
if nargin < 3
    % Figure out which dimension to work along.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end

% Compute X's mean and sd, and standardize it
mu = mean(x,dim);
sigma = std(x,flag,dim);
sigma0 = sigma;
sigma0(sigma0==0) = 1;
z = bsxfun(@minus,x, mu);
z = bsxfun(@rdivide, z, sigma0);

function narginchk(min,max)

n = evalin('caller','nargin');
if  n < min || n > max
    error('number of arguments')
end
function h = vline(x,varargin)

% h = vline(x,varargin)
% add vertical line(s) on the current axes at x
% all varargin arguments are passed to plot...

x = x(:);
ho = ishold;
hold on
h = plot([x x]',repmat(ylim,numel(x),1)',varargin{:});
if not(ho)
    hold off
end
if nargout == 0
    clear h
end
function h = hline(y,varargin)

% h = hline(y,varargin)
% add horizontal line(s) on the current axes at y
% all varargin arguments are passed to plot...

y = y(:);
ho = ishold;
hold on
h = plot(repmat(xlim,numel(y),1)',[y y]',varargin{:});
if not(ho)
    hold off
end
if nargout == 0
    clear h
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% BELOW IS ADJUST CODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ADJUST() - Automatic EEG artifact Detector
% with Joint Use of Spatial and Temporal features
%
% Usage:
%   >> [art, horiz, vert, blink, disc,...
%         soglia_DV, diff_var, soglia_K, med2_K, meanK, soglia_SED, med2_SED, SED, soglia_SAD, med2_SAD, SAD, ...
%         soglia_GDSF, med2_GDSF, GDSF, soglia_V, med2_V, nuovaV, soglia_D, maxdin]=ADJUST (EEG,out)
%
% Inputs:
%   EEG        - current dataset structure or structure array (has to be epoched)
%   out        - (string) report file name
%
% Outputs:
%   art        - List of artifacted ICs
%   horiz      - List of HEM ICs
%   vert       - List of VEM ICs
%   blink      - List of EB ICs
%   disc       - List of GD ICs
%   soglia_DV  - SVD threshold
%   diff_var   - SVD feature values
%   soglia_K   - TK threshold
%   meanK      - TK feature values
%   soglia_SED - SED threshold
%   SED        - SED feature values
%   soglia_SAD - SAD threshold
%   SAD        - SAD feature values
%   soglia_GDSF- GDSF threshold
%   GDSF       - GDSF feature values
%   soglia_V   - MEV threshold
%   nuovaV     - MEV feature values
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ADJUST
% Automatic EEG artifact Detector based on the Joint Use of Spatial and Temporal features
%
% Developed 2007-2014
% Andrea Mognon (1) and Marco Buiatti (2),
% (1) Center for Mind/Brain Sciences, University of Trento, Italy
% (2) INSERM U992 - Cognitive Neuroimaging Unit, Gif sur Yvette, France
%
% Last update: 02/05/2014 by Marco Buiatti
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reference paper:
% Mognon A, Bruzzone L, Jovicich J, Buiatti M,
% ADJUST: An Automatic EEG artifact Detector based on the Joint Use of Spatial and Temporal features.
% Psychophysiology 48 (2), 229-240 (2011).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2009-2014 Andrea Mognon (1) and Marco Buiatti (2),
% (1) Center for Mind/Brain Sciences, University of Trento, Italy
% (2) INSERM U992 - Cognitive Neuroimaging Unit, Gif sur Yvette, France
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% VERSIONS LOG
%
% 02/05/14: Modified text in Report.txt (MB).
%
% 30/03/14: Removed 'message to the user' (redundant). (MB)
%
% 22/03/14: kurtosis is replaced by kurt for compatibility if signal processing
%           toolbox is missing (MB).
%
% V2 (07 OCTOBER 2010) - by Andrea Mognon
% Added input 'nchannels' to compute_SAD and compute_SED_NOnorm;
% this is useful to differentiate the number of ICs (n) and the number of
% sensors (nchannels);
% bug reported by Guido Hesselman on October, 1 2010.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% function [art, horiz, vert, blink, disc,...
%         soglia_DV, diff_var, soglia_K, meanK, soglia_SED, SED, soglia_SAD, SAD, ...
%         soglia_GDSF, GDSF, soglia_V, nuovaV, soglia_D, maxdin]=ADJUST (EEG,out)
function [art, horiz, vert, blink, disc,...
    soglia_DV, diff_var, soglia_K, med2_K, meanK, soglia_SED, med2_SED, SED, soglia_SAD, med2_SAD, SAD, ...
    soglia_GDSF, med2_GDSF, GDSF, soglia_V, med2_V, nuovaV, soglia_D, maxdin]=ADJUST (EEG)


%% Settings

% ----------------------------------------------------
% |  Change experimental settings in this section    |
% ----------------------------------------------------

% ----------------------------------------------------
% |  Initial message to user:                        |
% ----------------------------------------------------
%
% disp(' ')
% disp('Detects Horizontal and Vertical eye movements,')
% disp('Blinks and Discontinuities in dataset:')
% disp([EEG.filename])
% disp(' ')

% ----------------------------------------------------
% |  Collect useful data from EEG structure          |
% ----------------------------------------------------

%number of ICs=size(EEG.icawinv,1);

%number of time points=size(EEG.data,2);

if length(size(EEG.data))==3
    
    num_epoch=size(EEG.data,3);
    
else
    
    num_epoch=1;
    
end

% Check the presence of ICA activations
EEG.icaact = eeg_getica(EEG);

topografie=EEG.icawinv'; %computes IC topographies

% Topographies and time courses normalization
%
% disp(' ');
% disp('Normalizing topographies...')
% disp('Scaling time courses...')

for i=1:size(EEG.icawinv,2) % number of ICs
    
    ScalingFactor=norm(topografie(i,:));
    
    topografie(i,:)=topografie(i,:)/ScalingFactor;
    
    if length(size(EEG.data))==3
        EEG.icaact(i,:,:)=ScalingFactor*EEG.icaact(i,:,:);
    else
        EEG.icaact(i,:)=ScalingFactor*EEG.icaact(i,:);
    end
    
end
%
% disp('Done.')
% disp(' ')

% Variables memorizing artifacted ICs indexes

blink=[];

horiz=[];

vert=[];

disc=[];

%% Check EEG channel position information
nopos_channels=[];
for el=1:length(EEG.chanlocs)
    if isempty(EEG.chanlocs(1,el).X)
        nopos_channels=[nopos_channels el];
    end;
end

if ~isempty(nopos_channels)
    disp(['Warning : Channels ' num2str(nopos_channels) ' have incomplete location information. They will NOT be used to compute ADJUST spatial features']);
    disp(' ');
end;

pos_channels=setdiff(1:length(EEG.chanlocs),nopos_channels);

%% Feature extraction

disp(' ')
disp('Features Extraction:')

%GDSF - General Discontinuity Spatial Feature

disp('GDSF - General Discontinuity Spatial Feature...')

GDSF = compute_GD_feat(topografie,EEG.chanlocs(1,pos_channels),size(EEG.icawinv,2));


%SED - Spatial Eye Difference

disp('SED - Spatial Eye Difference...')

[SED,medie_left,medie_right]=computeSED_NOnorm(topografie,EEG.chanlocs(1,pos_channels),size(EEG.icawinv,2));


%SAD - Spatial Average Difference

disp('SAD - Spatial Average Difference...')

[SAD,var_front,var_back,mean_front,mean_back]=computeSAD(topografie,EEG.chanlocs(1,pos_channels),size(EEG.icawinv,2));


%SVD - Spatial Variance Difference between front zone and back zone

diff_var=var_front-var_back;

%epoch dynamic range, variance and kurtosis

K=zeros(num_epoch,size(EEG.icawinv,2)); %kurtosis
Kloc=K;

Vmax=zeros(num_epoch,size(EEG.icawinv,2)); %variance

% disp('Computing variance and kurtosis of all epochs...')

for i=1:size(EEG.icawinv,2) % number of ICs
    
    for j=1:num_epoch
        Vmax(j,i)=var(EEG.icaact(i,:,j));
        %         Kloc(j,i)=kurtosis(EEG.icaact(i,:,j));
        K(j,i)=kurt(EEG.icaact(i,:,j));
    end
end

% check that kurt and kurtosis give the same values:
% [a,b]=max(abs(Kloc(:)-K(:)))

%TK - Temporal Kurtosis

disp('Temporal Kurtosis...')

meanK=zeros(1,size(EEG.icawinv,2));

for i=1:size(EEG.icawinv,2)
    if num_epoch>100
        meanK(1,i)=trim_and_mean(K(:,i));
    else meanK(1,i)=mean(K(:,i));
    end
    
end


%MEV - Maximum Epoch Variance

disp('Maximum epoch variance...')

maxvar=zeros(1,size(EEG.icawinv,2));
meanvar=zeros(1,size(EEG.icawinv,2));


for i=1:size(EEG.icawinv,2)
    if num_epoch>100
        maxvar(1,i)=trim_and_max(Vmax(:,i)');
        meanvar(1,i)=trim_and_mean(Vmax(:,i)');
    else
        maxvar(1,i)=max(Vmax(:,i));
        meanvar(1,i)=mean(Vmax(:,i));
    end
end

% MEV in reviewed formulation:

nuovaV=maxvar./meanvar;



%% Thresholds computation

disp('Computing EM thresholds...')

% soglia_K=EM(meanK);
%
% soglia_SED=EM(SED);
%
% soglia_SAD=EM(SAD);
%
% soglia_GDSF=EM(GDSF);
%
% soglia_V=EM(nuovaV);
[soglia_K,med1_K,med2_K]=EM(meanK);

[soglia_SED,med1_SED,med2_SED]=EM(SED);

[soglia_SAD,med1_SAD,med2_SAD]=EM(SAD);

[soglia_GDSF,med1_GDSF,med2_GDSF]=EM(GDSF);

[soglia_V,med1_V,med2_V]=EM(nuovaV);



%% Horizontal eye movements (HEM)


horiz=intersect(intersect(find(SED>=soglia_SED),find(medie_left.*medie_right<0)),...
    (find(nuovaV>=soglia_V)));




%% Vertical eye movements (VEM)



vert=intersect(intersect(find(SAD>=soglia_SAD),find(medie_left.*medie_right>0)),...
    intersect(find(diff_var>0),find(nuovaV>=soglia_V)));



%% Eye Blink (EB)


blink=intersect ( intersect( find(SAD>=soglia_SAD),find(medie_left.*medie_right>0) ) ,...
    intersect ( find(meanK>=soglia_K),find(diff_var>0) ));



%% Generic Discontinuities (GD)



disc=intersect(find(GDSF>=soglia_GDSF),find(nuovaV>=soglia_V));


%compute output variable
art = nonzeros( union (union(blink,horiz) , union(vert,disc)) )'; %artifact ICs

% these three are old outputs which are no more necessary in latest ADJUST version.
soglia_D=0;
soglia_DV=0;
maxdin=zeros(1,size(EEG.icawinv,2));

return

% compute_GD_feat() - Computes Generic Discontinuity spatial feature
%
% Usage:
%   >> res = compute_GD_feat(topografie,canali,num_componenti);
%
% Inputs:
%   topografie - topographies vector
%   canali     - EEG.chanlocs struct
%   num_componenti  - number of components
%
% Outputs:
%   res       - GDSF values

% Copyright (C) 2009-2014 Andrea Mognon (1) and Marco Buiatti (2),
% (1) Center for Mind/Brain Sciences, University of Trento, Italy
% (2) INSERM U992 - Cognitive Neuroimaging Unit, Gif sur Yvette, France
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


function res = compute_GD_feat(topografie,canali,num_componenti)

% Computes GDSF, discontinuity spatial feature
% topografie is the topography weights matrix
% canali is the structure EEG.chanlocs
% num_componenti is the number of ICs
% res is GDSF values

xpos=[canali.X];ypos=[canali.Y];zpos=[canali.Z];
pos=[xpos',ypos',zpos'];

res=zeros(1,num_componenti);

for ic=1:num_componenti
    
    % consider the vector topografie(ic,:)
    
    aux=[];
    
    for el=1:length(canali)-1
        
        P=pos(el,:); %position of current electrode
        d=pos-repmat(P,length(canali),1);
        %d=pos-repmat(P,62,1);
        dist=sqrt(sum((d.*d),2));
        
        [y,I]=sort(dist);
        repchas=I(2:11); % list of 10 nearest channels to el
        weightchas=exp(-y(2:11)); % respective weights, computed wrt distance
        
        aux=[aux abs(topografie(ic,el)-mean(weightchas.*topografie(ic,repchas)'))];
        % difference between el and the average of 10 neighbors
        % weighted according to weightchas
    end
    
    res(ic)=max(aux);
    
end


% computeSAD() - Computes Spatial Average Difference feature
%
% Usage:
%   >> [rapp,var_front,var_back,mean_front,mean_back]=computeSAD(topog,chanlocs,n);
%
% Inputs:
%   topog      - topographies vector
%   chanlocs   - EEG.chanlocs struct
%   n          - number of ICs
%   nchannels  - number of channels
%
% Outputs:
%   rapp       - SAD values
%   var_front  - Frontal Area variance values
%   var_back   - Posterior Area variance values
%   mean_front - Frontal Area average values
%   mean_back  - Posterior Area average values
%
%
% Copyright (C) 2009-2014 Andrea Mognon (1) and Marco Buiatti (2),
% (1) Center for Mind/Brain Sciences, University of Trento, Italy
% (2) INSERM U992 - Cognitive Neuroimaging Unit, Gif sur Yvette, France
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


function [rapp,var_front,var_back,mean_front,mean_back]=computeSAD(topog,chanlocs,n)

nchannels=length(chanlocs);

%% Define scalp zones

% Find electrodes in Frontal Area (FA)
dimfront=0; %number of FA electrodes
index1=zeros(1,nchannels); %indexes of FA electrodes

for k=1:nchannels
    if (abs(chanlocs(1,k).theta)<60) && (chanlocs(1,k).radius>0.40) %electrodes are in FA
        dimfront=dimfront+1; %count electrodes
        index1(1,dimfront)=k;
    end
end

% Find electrodes in Posterior Area (PA)
dimback=0;
index3=zeros(1,nchannels);
for h=1:nchannels
    if (abs(chanlocs(1,h).theta)>110)
        dimback=dimback+1;
        index3(1,dimback)=h;
    end
end

if dimfront*dimback==0
    disp('ERROR: no channels included in some scalp areas.')
    disp('Check channels distribution and/or change scalp areas definitions in computeSAD.m and computeSED_NOnorm.m')
    disp('ADJUST session aborted.')
    return
end

%% Outputs

rapp=zeros(1,n); % SAD
mean_front=zeros(1,n); % FA electrodes mean value
mean_back=zeros(1,n); % PA electrodes mean value
var_front=zeros(1,n); % FA electrodes variance value
var_back=zeros(1,n); % PA electrodes variance value

%% Output computation

for i=1:n % for each topography
    
    %create FA electrodes vector
    front=zeros(1,dimfront);
    for h=1:dimfront
        front(1,h)=topog(i,index1(1,h));
    end
    
    %create PA electrodes vector
    back=zeros(1,dimback);
    for h=1:dimback
        back(1,h)=topog(i,index3(1,h));
    end
    
    
    
    %compute features
    
    rapp(1,i)=abs(mean(front))-abs(mean(back)); % SAD
    mean_front(1,i)=mean(front);
    mean_back(1,i)=mean(back);
    var_back(1,i)=var(back);
    var_front(1,i)=var(front);
    
end


% computeSED_NOnorm() - Computes Spatial Eye Difference feature
% without normalization
%
% Usage:
%   >> [out,medie_left,medie_right]=computeSED_NOnorm(topog,chanlocs,n);
%
% Inputs:
%   topog      - topographies vector
%   chanlocs   - EEG.chanlocs struct
%   n          - number of ICs
%   nchannels  - number of channels
%
% Outputs:
%   out        - SED values
%   medie_left - Left Eye area average values
%   medie_right- Right Eye area average values
%
%
% Copyright (C) 2009-2014 Andrea Mognon (1) and Marco Buiatti (2),
% (1) Center for Mind/Brain Sciences, University of Trento, Italy
% (2) INSERM U992 - Cognitive Neuroimaging Unit, Gif sur Yvette, France
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

function [out,medie_left,medie_right]=computeSED_NOnorm(topog,chanlocs,n)

nchannels=length(chanlocs);

%% Define scalp zones

% Find electrodes in Left Eye area (LE)
dimleft=0; %number of LE electrodes
index1=zeros(1,nchannels); %indexes of LE electrodes

for k=1:nchannels
    if (-61<chanlocs(1,k).theta) && (chanlocs(1,k).theta<-35) && (chanlocs(1,k).radius>0.30) %electrodes are in LE
        dimleft=dimleft+1; %count electrodes
        index1(1,dimleft)=k;
    end
end

% Find electrodes in Right Eye area (RE)
dimright=0; %number of RE electrodes
index2=zeros(1,nchannels); %indexes of RE electrodes
for g=1:nchannels
    if (34<chanlocs(1,g).theta) && (chanlocs(1,g).theta<61) && (chanlocs(1,g).radius>0.30) %electrodes are in RE
        dimright=dimright+1; %count electrodes
        index2(1,dimright)=g;
    end
end

% Find electrodes in Posterior Area (PA)
dimback=0;
index3=zeros(1,nchannels);
for h=1:nchannels
    if (abs(chanlocs(1,h).theta)>110)
        dimback=dimback+1;
        index3(1,dimback)=h;
    end
end

if dimleft*dimright*dimback==0
    disp('ERROR: no channels included in some scalp areas.')
    disp('Check channels distribution and/or change scalp areas definitions in computeSAD.m and computeSED_NOnorm.m')
    disp('ADJUST session aborted.')
    return
end

%% Outputs

out=zeros(1,n); %memorizes SED
medie_left=zeros(1,n); %memorizes LE mean value
medie_right=zeros(1,n); %memorizes RE mean value

%% Output computation

for i=1:n  % for each topography
    %create LE electrodes vector
    left=zeros(1,dimleft);
    for h=1:dimleft
        left(1,h)=topog(i,index1(1,h));
    end
    
    %create RE electrodes vector
    right=zeros(1,dimright);
    for h=1:dimright
        right(1,h)=topog(i,index2(1,h));
    end
    
    %create PA electrodes vector
    back=zeros(1,dimback);
    for h=1:dimback
        back(1,h)=topog(i,index3(1,h));
    end
    
    
    
    %compute features
    out1=abs(mean(left)-mean(right));
    out2=var(back);
    out(1,i)=out1; % SED not notmalized
    medie_left(1,i)=mean(left);
    medie_right(1,i)=mean(right);
    
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EM - ADJUST package
%
% Performs automatic threshold on the digital numbers
% of the input vector 'vec'; based on Expectation - Maximization algorithm

% Reference paper:
% Bruzzone, L., Prieto, D.F., 2000. Automatic analysis of the difference image
% for unsupervised change detection.
% IEEE Trans. Geosci. Remote Sensing 38, 1171:1182

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Usage:
%   >> [last,med1,med2,var1,var2,prior1,prior2]=EM(vec);
%
% Input: vec (row vector, to be thresholded)
%
% Outputs: last (threshold value)
%          med1,med2 (mean values of the Gaussian-distributed classes 1,2)
%          var1,var2 (variance of the Gaussian-distributed classes 1,2)
%          prior1,prior2 (prior probabilities of the Gaussian-distributed classes 1,2)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2009-2014 Andrea Mognon (1) and Marco Buiatti (2),
% (1) Center for Mind/Brain Sciences, University of Trento, Italy
% (2) INSERM U992 - Cognitive Neuroimaging Unit, Gif sur Yvette, France
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


function [last,med1,med2,var1,var2,prior1,prior2]=EM(vec)

if size(vec,2)>1
    len=size(vec,2); %number of elements
else
    vec=vec';
    len=size(vec,2);
end

c_FA=1; % False Alarm cost
c_MA=1; % Missed Alarm cost

med=mean(vec);
standard=std(vec);
mediana=(max(vec)+min(vec))/2;

alpha1=0.01*(max(vec)-mediana); % initialization parameter/ righthand side
alpha2=0.01*(mediana-min(vec)); % initialization parameter/ lefthand side

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPECTATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

train1=[]; % Expectation of class 1
train2=[];
train=[]; % Expectation of 'unlabeled' samples

for i=1:(len)
    if (vec(i)<(mediana-alpha2))
        train2=[train2 vec(i)];
    elseif (vec(i)>(mediana+alpha1))
        train1=[train1 vec(i)];
    else
        train=[train vec(i)];
    end
end

n1=length(train1);
n2=length(train2);

med1=mean(train1);
med2=mean(train2);
prior1=n1/(n1+n2);
prior2=n2/(n1+n2);
var1=var(train1);
var2=var(train2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAXIMIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count=0;
dif_med_1=1; % difference between current and previous mean
dif_med_2=1;
dif_var_1=1; % difference between current and previous variance
dif_var_2=1;
dif_prior_1=1; % difference between current and previous prior
dif_prior_2=1;
stop=0.0001;

while((dif_med_1>stop)&&(dif_med_2>stop)&&(dif_var_1>stop)&&(dif_var_2>stop)&&(dif_prior_1>stop)&&(dif_prior_2>stop))
    
    count=count+1;
    
    med1_old=med1;
    med2_old=med2;
    var1_old=var1;
    var2_old=var2;
    prior1_old=prior1;
    prior2_old=prior2;
    prior1_i=[];
    prior2_i=[];
    
    % FOLLOWING FORMULATION IS ACCORDING TO REFERENCE PAPER:
    
    for i=1:len
        prior1_i=[prior1_i prior1_old*Bayes(med1_old,var1_old,vec(i))/...
            (prior1_old*Bayes(med1_old,var1_old,vec(i))+prior2_old*Bayes(med2_old,var2_old,vec(i)))];
        prior2_i=[prior2_i prior2_old*Bayes(med2_old,var2_old,vec(i))/...
            (prior1_old*Bayes(med1_old,var1_old,vec(i))+prior2_old*Bayes(med2_old,var2_old,vec(i)))];
    end
    
    
    prior1=sum(prior1_i)/len;
    prior2=sum(prior2_i)/len;
    med1=sum(prior1_i.*vec)/(prior1*len);
    med2=sum(prior2_i.*vec)/(prior2*len);
    var1=sum(prior1_i.*((vec-med1_old).^2))/(prior1*len);
    var2=sum(prior2_i.*((vec-med2_old).^2))/(prior2*len);
    
    dif_med_1=abs(med1-med1_old);
    dif_med_2=abs(med2-med2_old);
    dif_var_1=abs(var1-var1_old);
    dif_var_2=abs(var2-var2_old);
    dif_prior_1=abs(prior1-prior1_old);
    dif_prior_2=abs(prior2-prior2_old);
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THRESHOLDING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=c_MA/c_FA;
a=(var1-var2)/2;
b= ((var2*med1)-(var1*med2));
c=(log((k*prior1*sqrt(var2))/(prior2*sqrt(var1)))*(var2*var1))+(((((med2)^2)*var1)-(((med1)^2)*var2))/2);
rad=(b^2)-(4*a*c);
if rad<0
    disp('Negative Discriminant!');
    return;
end

soglia1=(-b+sqrt(rad))/(2*a);
soglia2=(-b-sqrt(rad))/(2*a);

if ((soglia1<med2)||(soglia1>med1))
    last=soglia2;
else
    last=soglia1;
end

if isnan(last) % TO PREVENT CRASHES
    last=mediana;
end

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function prob=Bayes(med,var,point)
if var==0
    prob=1;
else
    prob=((1/(sqrt(2*pi*var)))*exp((-1)*((point-med)^2)/(2*var)));
end





% trim_and_max() - Computes maximum value from vector 'vettore'
% after removing the top 1% of the values
% (to be outlier resistant)
%
% Usage:
%   >> valore=trim_and_max(vettore);
%
% Inputs:
%   vettore    - row vector
%
% Outputs:
%   valore     - result
%
%
% Author: Andrea Mognon, Center for Mind/Brain Sciences, University of
% Trento, 2009

% Motivation taken from the following comment to our paper:
% "On page 11 the authors motivate the use of the max5 function when computing
% Maximum Epoch Variance because the simple maximum would be too sensitive
% to spurious outliers. This is a good concern, however the max5 function would
% still be sensitive to spurious outliers for very large data sets. In other words, if
% the data set is large enough, one will be very likely to record more than five
% outliers. The authors should use a trimmed max function that computes the
% simple maximum after the top say .1% of the values have been removed from
% consideration. This rejection criteria scales appropriately with the size of the data
% set."

% Copyright (C) 2009-2014 Andrea Mognon (1) and Marco Buiatti (2),
% (1) Center for Mind/Brain Sciences, University of Trento, Italy
% (2) INSERM U992 - Cognitive Neuroimaging Unit, Gif sur Yvette, France
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


function valore=trim_and_max(vettore)


dim=floor(.01*size(vettore,2)); % = 1% of vector length

tmp=sort(vettore);
valore= tmp(length(vettore)-dim);



% trim_and_mean() - Computes average value from vector 'vettore'
% after removing the top .1% of the values
% (to be outlier resistant)
%
% Usage:
%   >> valore=trim_and_mean(vettore);
%
% Inputs:
%   vettore    - row vector
%
% Outputs:
%   valore     - result
%
% Copyright (C) 2009-2014 Andrea Mognon (1) and Marco Buiatti (2),
% (1) Center for Mind/Brain Sciences, University of Trento, Italy
% (2) INSERM U992 - Cognitive Neuroimaging Unit, Gif sur Yvette, France
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


function valore=trim_and_mean(vettore)


dim=floor(.01*size(vettore,2)); % = 1% of vector length

tmp=sort(vettore);
valore= mean (tmp(1:(length(vettore)-dim)));





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% END  ADJUST   CODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%   BELOW IS FASTER CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function list_properties = component_properties(EEG,blink_chans,lpf_band)

% Copyright (C) 2010 Hugh Nolan, Robert Whelan and Richard Reilly, Trinity College Dublin,
% Ireland
% nolanhu@tcd.ie, robert.whelan@tcd.ie
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


list_properties = [];
%

if ~exist('lpf_band','var') || length(lpf_band)~=2 || ~any(lpf_band)
    ignore_lpf=1;
else
    ignore_lpf=0;
end

delete_activations_after=0;
if ~isfield(EEG,'icaact') || isempty(EEG.icaact)
    delete_activations_after=1;
    EEG.icaact = eeg_getica(EEG);
end
try
    checkfunctionmatlab('pwelch', 'signal_toolbox')
end
for u = 1:size(EEG.icaact,1)
    [spectra(u,:) freqs] = pwelch(EEG.icaact(u,:),[],[],(EEG.srate),EEG.srate);
end

list_properties = zeros(size(EEG.icaact,1),5); %This 5 corresponds to number of measurements made.

for u=1:size(EEG.icaact,1)
    measure = 1;
    % TEMPORAL PROPERTIES
    
    % 1 Median gradient value, for high frequency stuff
    list_properties(u,measure) = median(diff(EEG.icaact(u,:)));
    measure = measure + 1;
    
    % 2 Mean slope around the LPF band (spectral)
    if ignore_lpf
        list_properties(u,measure) = 0;
    else
        list_properties(u,measure) = mean(diff(10*log10(spectra(u,find(freqs>=lpf_band(1),1):find(freqs<=lpf_band(2),1,'last')))));
    end
    measure = measure + 1;
    
    % SPATIAL PROPERTIES
    
    % 3 Kurtosis of spatial map (if v peaky, i.e. one or two points high
    % and everywhere else low, then it's probably noise on a single
    % channel)
    list_properties(u,measure) = kurt(EEG.icawinv(:,u));
    measure = measure + 1;
    
    % OTHER PROPERTIES
    
    % 4 Hurst exponent
    list_properties(u,measure) = hurst_exponent(EEG.icaact(u,:));
    measure = measure + 1;
    
    % 10 Eyeblink correlations
    if (exist('blink_chans','var') && ~isempty(blink_chans))
        for v = 1:length(blink_chans)
            if ~(max(EEG.data(blink_chans(v),:))==0 && min(EEG.data(blink_chans(v),:))==0);
                f = corrcoef(EEG.icaact(u,:),EEG.data(blink_chans(v),:));
                x(v) = abs(f(1,2));
            else
                x(v) = v;
            end
        end
        list_properties(u,measure) = max(x);
        measure = measure + 1;
    end
end

for u = 1:size(list_properties,2)
    list_properties(isnan(list_properties(:,u)),u)=nanmean(list_properties(:,u));
    list_properties(:,u) = list_properties(:,u) - median(list_properties(:,u));
end

if delete_activations_after
    EEG.icaact=[];
end


% The Hurst exponent
%--------------------------------------------------------------------------
% This function does dispersional analysis on a data series, then does a
% Matlab polyfit to a log-log plot to estimate the Hurst exponent of the
% series.
%
% This algorithm is far faster than a full-blown implementation of Hurst's
% algorithm.  I got the idea from a 2000 PhD dissertation by Hendrik J
% Blok, and I make no guarantees whatsoever about the rigor of this approach
% or the accuracy of results.  Use it at your own risk.
%
% Bill Davidson
% 21 Oct 2003

function [hurst] = hurst_exponent(data0)   % data set

data=data0;         % make a local copy

[M,npoints]=size(data0);

yvals=zeros(1,npoints);
xvals=zeros(1,npoints);
data2=zeros(1,npoints);

index=0;
binsize=1;

while npoints>4
    
    y=std(data);
    index=index+1;
    xvals(index)=binsize;
    yvals(index)=binsize*y;
    
    npoints=fix(npoints/2);
    binsize=binsize*2;
    for ipoints=1:npoints % average adjacent points in pairs
        data2(ipoints)=(data(2*ipoints)+data((2*ipoints)-1))*0.5;
    end
    data=data2(1:npoints);
    
end % while

xvals=xvals(1:index);
yvals=yvals(1:index);

logx=log(xvals);
logy=log(yvals);

p2=polyfit(logx,logy,1);
hurst=p2(1); % Hurst exponent is the slope of the linear fit of log-log plot

return;


function [lengths]  =  min_z(list_properties, rejection_options)
if (~exist('rejection_options', 'var'))
    rejection_options.measure = ones(1, size(list_properties, 2));
    rejection_options.z = 3*ones(1, size(list_properties, 2));
end

rejection_options.measure = logical(rejection_options.measure);
zs = list_properties - repmat(mean(list_properties, 1), size(list_properties, 1), 1);
zs = zs./repmat(std(zs, [], 1), size(list_properties, 1), 1);
zs(isnan(zs)) = 0;
all_l  =  abs(zs) > repmat(rejection_options.z, size(list_properties, 1), 1);
lengths  =  any(all_l(:, rejection_options.measure), 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%   END FASTER CODE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%   BEGIN MARA CODE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MARA() - Automatic classification of multiple artifact components
%          Classies artifactual ICs based on 6 features from the time domain,
%           the frequency domain, and the pattern
%
% Usage:
%   >> [artcomps, info] = MARA(EEG);
%
% Inputs:
%   EEG         - input EEG structure
%
% Outputs:
%   artcomps    - array containing the numbers of the artifactual
%                 components
%   info        - struct containing more information about MARA classification
%                   .posterior_artefactprob : posterior probability for each
%                            IC of being an artefact
%                   .normfeats : <6 x nIC > features computed by MARA for each IC,
%                            normalized by the training data
%                      The features are: (1) Current Density Norm, (2) Range
%                      in Pattern, (3) Local Skewness of the Time Series,
%                      (4) Lambda, (5) 8-13 Hz, (6) FitError.
%
%  For more information see:
%  I. Winkler, S. Haufe, and M. Tangermann, Automatic classification of artifactual ICA-components
%  for artifact removal in EEG signals, Behavioral and Brain Functions, 7, 2011.
%
% See also: processMARA()

% Copyright (C) 2013 Irene Winkler and Eric Waldburger
% Berlin Institute of Technology, Germany
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
function [artcomps, info] = MARA(EEG)
%%%%%%%%%%%%%%%%%%%%
%%  Calculate features from the pattern (component map)
%%%%%%%%%%%%%%%%%%%%
% extract channel labels
clab = {};
for i=1:length(EEG.chanlocs)
    clab{i} = EEG.chanlocs(i).labels;
end

% cut to channel labels common with training data
load('fv_training_MARA'); %load struct fv_tr
[clab_common i_te i_tr ] = intersect(upper(clab), upper(fv_tr.clab));
clab_common = fv_tr.clab(i_tr);
if length(clab_common) == 0
    error(['There were no matching channeldescriptions found.' , ...
        'MARA needs channel labels of the form Cz, Oz, F3, F4, Fz, etc. Aborting.'])
end
patterns = (EEG.icawinv(i_te,:));
[M100 idx] = get_M100_ADE(clab_common); %needed for Current Density Norm

disp('MARA is computing features. Please wait');
%standardize patterns
patterns = patterns./repmat(std(patterns,0,1),length(patterns(:,1)),1);

%compute current density norm
feats(1,:) = log(sqrt(sum((M100*patterns(idx,:)).^2)));
%compute spatial range
feats(2,:) = log(max(patterns) - min(patterns));

%%%%%%%%%%%%%%%%%%%%
%%  Calculate time and frequency features
%%%%%%%%%%%%%%%%%%%%
%compute time and frequency features (Current Density Norm, Range Within Pattern,
%Average Local Skewness, Band Power 8 - 13 Hz)
feats(3:6,:) = extract_time_freq_features(EEG);
disp('Features ready');


%%%%%%%%%%%%%%%%%%%%%%
%%  Adapt train features to clab
%%%%%%%%%%%%%%%%%%%%
fv_tr.pattern = fv_tr.pattern(i_tr, :);
fv_tr.pattern = fv_tr.pattern./repmat(std(fv_tr.pattern,0,1),length(fv_tr.pattern(:,1)),1);
fv_tr.x(2,:) = log(max(fv_tr.pattern) - min(fv_tr.pattern));
fv_tr.x(1,:) = log(sqrt(sum((M100 * fv_tr.pattern).^2)));

%%%%%%%%%%%%%%%%%%%%
%%  Classification
%%%%%%%%%%%%%%%%%%%%
[C, foo, posterior] = classify(feats',fv_tr.x',fv_tr.labels(1,:));
artcomps = find(C == 0)';
info.posterior_artefactprob = posterior(:, 1)';
info.normfeats = (feats - repmat(mean(fv_tr.x, 2), 1, size(feats, 2)))./ ...
    repmat(std(fv_tr.x,0, 2), 1, size(feats, 2));

function features = extract_time_freq_features(EEG)
%                             - 1st row: Average Local Skewness
%                             - 2nd row: lambda
%                             - 3rd row: Band Power 8 - 13 Hz
%                             - 4rd row: Fit Error
%
data = EEG.data(EEG.icachansind,:,:);
fs = EEG.srate; %sampling frequency

% transform epoched data into continous data
data = data(:,:);

%downsample (to 100-200Hz)
factor = max(floor(EEG.srate/100),1);
data = data(:, 1:factor:end);
fs = round(fs/factor);

%compute icaactivation and standardise variance to 1
icacomps = EEG.icaact(:,:)';%(EEG.icaweights * EEG.icasphere * data)';
icacomps = icacomps./repmat(std(icacomps,0,1),length(icacomps(:,1)),1);
icacomps = icacomps';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate featues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ic=1:size(icacomps,1)  %for each component
    fprintf('.');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Proc Spectrum for Channel
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [pxx, freq] = pwelch(icacomps(ic,:), ones(1, fs), [], fs, fs);
    pxx = 10*log10(pxx * fs/2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The average log band power between 8 and 13 Hz
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p = 0;
    for i = 8:13
        p = p + pxx(find(freq == i,1));
    end
    Hz8_13 = p / (13-8+1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % lambda and FitError: deviation of a component's spectrum from
    % a protoptypical 1/frequency curve
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p1.x = 2; %first point: value at 2 Hz
    p1.y = pxx(find(freq == p1.x,1));
    
    p2.x = 3; %second point: value at 3 Hz
    p2.y = pxx(find(freq == p2.x,1));
    
    %third point: local minimum in the band 5-13 Hz
    p3.y = min(pxx(find(freq == 5,1):find(freq == 13,1)));
    p3.x = freq(find(pxx == p3.y,1));
    
    %fourth point: min - 1 in band 5-13 Hz
    p4.x = p3.x - 1;
    p4.y = pxx(find(freq == p4.x,1));
    
    %fifth point: local minimum in the band 33-39 Hz
    p5.y = min(pxx(find(freq == 33,1):find(freq == 39,1)));
    p5.x = freq(find(pxx == p5.y,1));
    
    %sixth point: min + 1 in band 33-39 Hz
    p6.x = p5.x + 1;
    p6.y = pxx(find(freq == p6.x,1));
    
    pX = [p1.x; p2.x; p3.x; p4.x; p5.x; p6.x];
    pY = [p1.y; p2.y; p3.y; p4.y; p5.y; p6.y];
    
    myfun = @(x,xdata)(exp(x(1))./ xdata.^exp(x(2))) - x(3);
    xstart = [4, -2, 54];
    fittedmodel = lsqcurvefit(myfun,xstart,double(pX),double(pY), [], [], optimset('Display', 'off'));
    
    %FitError: mean squared error of the fit to the real spectrum in the band 2-40 Hz.
    ts_8to15 = freq(find(freq == 8) : find(freq == 15));
    fs_8to15 = pxx(find(freq == 8) : find(freq == 15));
    fiterror = log(norm(myfun(fittedmodel, ts_8to15)-fs_8to15)^2);
    
    %lambda: parameter of the fit
    lambda = fittedmodel(2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Averaged local skewness 15s
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    interval = 15;
    abs_local_scewness = [];
    for i=1:interval:length(icacomps(ic,:))/fs-interval
        abs_local_scewness = [abs_local_scewness, abs(skewness(icacomps(ic, i * fs:(i+interval) * fs)))];
    end
    
    if isempty(abs_local_scewness)
        error('MARA needs at least 15ms long ICs to compute its features.')
    else
        mean_abs_local_scewness_15 = log(mean(abs_local_scewness));
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Append Features
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    features(:,ic)= [mean_abs_local_scewness_15, lambda, Hz8_13, fiterror];
end
disp('.');


function [M100, idx_clab_desired] = get_M100_ADE(clab_desired)
% [M100, idx_clab_desired] = get_M100_ADEC(clab_desired)
%
% IN  clab_desired - channel setup for which M100 should be calculated
% OUT M100
%     idx_clab_desired
% M100 is the matrix such that  feature = norm(M100*ica_pattern(idx_clab_desired), 'fro')
%
% (c) Stefan Haufe

lambda = 100;

load inv_matrix_icbm152; %L (forward matrix 115 x 2124 x 3), clab (channel labels)

[cl_ ia idx_clab_desired] = intersect(clab, clab_desired);
F = L(ia, :, :); %forward matrix for desired channels labels
[n_channels m foo] = size(F);  %m = 2124, number of dipole locations
F = reshape(F, n_channels, 3*m);

%H - matrix that centralizes the pattern, i.e. mean(H*pattern) = 0
H = eye(n_channels) -  ones(n_channels, n_channels)./ n_channels;
%W - inverse of the depth compensation matrix Lambda
W = sloreta_invweights(L);

L = H*F*W;

%We have inv(L'L +lambda eye(size(L'*L))* L' = L'*inv(L*L' + lambda
%eye(size(L*L')), which is easier to calculate as number of dimensions is
%much smaller

%calulate the inverse of L*L' + lambda * eye(size(L*L')
[U D] = eig(L*L');
d = diag(D);
di = d+lambda;
di = 1./di;
di(d < 1e-10) = 0;
inv1 = U*diag(di)*U';  %inv1 = inv(L*L' + lambda *eye(size(L*L'))

%get M100
M100 = L'*inv1*H;



function W = sloreta_invweights(LL)
% inverse sLORETA-based weighting
%
% Synopsis:
%   W = sloreta_invweights(LL);
%
% Arguments:
%   LL: [M N 3] leadfield tensor
%
% Returns:
%   W: [3*N 3*N] block-diagonal matrix of weights
%
% Stefan Haufe, 2007, 2008
%
% License
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see http://www.gnu.org/licenses/.

[M N NDUM]=size(LL);
L=reshape(permute(LL, [1 3 2]), M, N*NDUM);

L = L - repmat(mean(L, 1), M, 1);

T = L'*pinv(L*L');

W = spalloc(N*NDUM, N*NDUM, N*NDUM*NDUM);
for ivox = 1:N
    W(NDUM*(ivox-1)+(1:NDUM), NDUM*(ivox-1)+(1:NDUM)) = (T(NDUM*(ivox-1)+(1:NDUM), :)*L(:, NDUM*(ivox-1)+(1:NDUM)))^-.5;
end

ind = [];
for idum = 1:NDUM
    ind = [ind idum:NDUM:N*NDUM];
end
W = W(ind, ind);




function [i_te, i_tr] = findconvertedlabels(pos_3d, chanlocs)
% IN  pos_3d  - 3d-positions of training channel labels
%     chanlocs - EEG.chanlocs structure of data to be classified

%compute spherical coordinates theta and phi for the training channel
%label
[theta, phi, r] = cart2sph(pos_3d(1,:),pos_3d(2,:), pos_3d(3,:));
theta = theta - pi/2;
theta(theta < -pi) = theta(theta < -pi) + 2*pi;
theta = theta*180/pi;
phi = phi * 180/pi;
theta(find(pos_3d(1,:) == 0 & pos_3d(2,:) == 0)) = 0; %exception for Cz


clab_common = {};
i_te = [];
i_tr = [];

%For each channel in EEG.chanlocs, try to find matching channel in
%training data
for chan = 1:length(chanlocs)
    if not(isempty(chanlocs(chan).sph_phi))
        idx = find((theta <= chanlocs(chan).sph_theta + 6) ...
            & (theta >= chanlocs(chan).sph_theta - 6) ...
            & (phi <= chanlocs(chan).sph_phi + 6) ...
            & (phi >= chanlocs(chan).sph_phi - 6));
        if not(isempty(idx))
            i_tr = [i_tr, idx(1)];
            i_te = [i_te, chan];
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%   END MARA CODE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function textprogressbar(c)
% This function creates a text progress bar. It should be called with a
% STRING argument to initialize and terminate. Otherwise the number correspoding
% to progress in % should be supplied.
% INPUTS:   C   Either: Text string to initialize or terminate
%                       Percentage number to show progress
% OUTPUTS:  N/A
% Example:  Please refer to demo_textprogressbar.m

% Author: Paul Proteus (e-mail: proteus.paul (at) yahoo (dot) com)
% Version: 1.0
% Changes tracker:  29.06.2010  - First version

% Inspired by: http://blogs.mathworks.com/loren/2007/08/01/monitoring-progress-of-a-calculation/

%% Initialization
persistent strCR prevc strCRtitle;           %   Carriage return pesistent variable

% Vizualization parameters
strPercentageLength = 10;   %   Length of percentage string (must be >5)
strDotsMaximum      = 10;   %   The total number of dots in a progress bar

%% Main
if nargin == 0
    % Progress bar  - force termination/initialization
    fprintf('\n');
    strCR = [];
    strCRtitle = [];
    prevc = [];
elseif ischar(c)
    % Progress bar - set/reset title
    if not(isempty(strCR)) && all(strCR ~= -1)
        fprintf(strCR);
    end
    if not(isempty(strCRtitle))
        fprintf(strCRtitle);
    end
    % add trailing space if not one already
    if isempty(regexp(c,'\s$', 'once'))
        c = [c ' '];
    end
    fprintf('%s',c);
    strCR = -1;strCRtitle = repmat('\b',1,numel(c));
elseif isnumeric(c)
    % Progress bar - normal progress
    if isempty(prevc)
        prevc = 0;
    end
    c = floor(c);
    if c == prevc
        return
    else
        prevc = c;
    end
    percentageOut = [num2str(c) '%%'];
    percentageOut = [percentageOut repmat(' ',1,strPercentageLength-length(percentageOut)-1)];
    nDots = floor(c/100*strDotsMaximum);
    dotOut = ['[' repmat('.',1,nDots) repmat(' ',1,strDotsMaximum-nDots) ']'];
    strOut = [percentageOut dotOut];
    
    % Print it on the screen
    if strCR == -1,
        % Don't do carriage return during first run
        fprintf(strOut);
    else
        % Do it during all the other runs
        fprintf([strCR strOut]);
    end
    
    % Update carriage return
    strCR = repmat('\b',1,length(strOut)-1);
    
else
    % Any other unexpected input
    error('Unsupported argument type');
end

% eeg_getica() - get ICA component activation. Recompute if necessary.
%
% >> mergelocs = eeg_getica(EEG, comp);
%
% Inputs:
%     EEG     - EEGLAB dataset structure
%     comp    - component index
%
% Output:
%     icaact  - ICA component activity
%
% Author: Arnaud Delorme, 2006

% Copyright (C) Arnaud Delorme, CERCO, 2006, arno@salk.edu
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

function icaact = eeg_getica(EEG, comp)

if nargin < 1
    help eeg_getica;
    return;
end;
if nargin < 2
    comp = 1:size(EEG.icaact,1);
end;

if ~isempty(EEG.icaact)
    icaact = EEG.icaact(comp,:,:);
else
    disp('Recomputing ICA activations');
    if isempty(EEG.icachansind)
        EEG.icachansind = 1:EEG.nbchan;
        disp('Channels indices are assumed to be in regular order and arranged accordingly');
    end
    icaact = (EEG.icaweights(comp,:)*EEG.icasphere)*reshape(EEG.data(EEG.icachansind,:,:), length(EEG.icachansind), EEG.trials*EEG.pnts);
    icaact = reshape( icaact, size(icaact,1), EEG.pnts, EEG.trials);
end;
% kurt() - return kurtosis of input data distribution
%
% Usage:
%   >> k=kurt(data)
%
% Algorithm:
%   Calculates kurtosis or normalized 4th moment of an input data vector
%   Given a matrix, returns a row vector giving the kurtosis' of the columns
%   (Ref: "Numerical Recipes," p. 612)
%
% Author: Martin Mckeown, CNL / Salk Institute, La Jolla, 10/2/96

% Copyright (C) Martin Mckeown, CNL / Salk Institute, La Jolla, 7/1996
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

% 2/28/97 - made to return separate kurtosis estimates of columns -Scott Makeig
% 01-25-02 reformated help & license, added links -ad

function [k] = kurt(data)

[r,c]=size(data);
if r==1,
    kdata = data';  % if a row vector, make into a column vector
    r = c;
else
    kdata = data;
end
%fprintf('size of kdata = [%d,%d]\n',size(kdata,1),size(kdata,2));

mn = mean(kdata);              % find the column means
diff = kdata-ones(r,1)*mn;     % remove the column means
dsq = diff.*diff;              % square the data

k =  (sum(dsq.*dsq)./std(kdata).^4)./r - 3;

% topoplot() - plot a topographic map of a scalp data field in a 2-D circular view
%              (looking down at the top of the head) using interpolation on a fine
%              cartesian grid. Can also show specified channnel location(s), or return
%              an interpolated value at an arbitrary scalp location (see 'noplot').
%              By default, channel locations below head center (arc_length 0.5) are
%              shown in a 'skirt' outside the cartoon head (see 'plotrad' and 'headrad'
%              options below). Nose is at top of plot; left is left; right is right.
%              Using option 'plotgrid', the plot may be one or more rectangular grids.
% Usage:
%        >>  topoplot(datavector, EEG.chanlocs);   % plot a map using an EEG chanlocs structure
%        >>  topoplot(datavector, 'my_chan.locs'); % read a channel locations file and plot a map
%        >>  topoplot('example');                  % give an example of an electrode location file
%        >>  [h grid_or_val plotrad_or_grid, xmesh, ymesh]= ...
%                           topoplot(datavector, chan_locs, 'Input1','Value1', ...);
% Required Inputs:
%   datavector        - single vector of channel values. Else, if a vector of selected subset
%                       (int) channel numbers -> mark their location(s) using 'style' 'blank'.
%   chan_locs         - name of an EEG electrode position file (>> topoplot example).
%                       Else, an EEG.chanlocs structure (>> help readlocs or >> topoplot example)
% Optional inputs:
%   'maplimits'       - 'absmax'   -> scale map colors to +/- the absolute-max (makes green 0);
%                       'maxmin'   -> scale colors to the data range (makes green mid-range);
%                       [lo.hi]    -> use user-definined lo/hi limits
%                       {default: 'absmax'}
%   'style'           - 'map'      -> plot colored map only
%                       'contour'  -> plot contour lines only
%                       'both'     -> plot both colored map and contour lines
%                       'fill'     -> plot constant color between contour lines
%                       'blank'    -> plot electrode locations only {default: 'both'}
%   'electrodes'      - 'on','off','labels','numbers','ptslabels','ptsnumbers'. To set the 'pts'
%                       marker,,see 'Plot detail options' below. {default: 'on' -> mark electrode
%                       locations with points ('.') unless more than 64 channels, then 'off'}.
%   'plotchans'       - [vector] channel numbers (indices) to use in making the head plot.
%                       {default: [] -> plot all chans}
%   'chantype'        - cell array of channel type(s) to plot. Will also accept a single quoted
%                       string type. Channel type for channel k is field EEG.chanlocs(k).type.
%                       If present, overrides 'plotchans' and also 'chaninfo' with field
%                       'chantype'. Ex. 'EEG' or {'EEG','EOG'} {default: all, or 'plotchans' arg}
%   'plotgrid'        - [channels] Plot channel data in one or more rectangular grids, as
%                       specified by [channels],  a position matrix of channel numbers defining
%                       the topographic locations of the channels in the grid. Zero values are
%                       given the figure background color; negative integers, the color of the
%                       polarity-reversed channel values.  Ex: >> figure; ...
%                        >> topoplot(values,'chanlocs','plotgrid',[11 12 0; 13 14 15]);
%                       % Plot a (2,3) grid of data values from channels 11-15 with one empty
%                       grid cell (top right) {default: no grid plot}
%   'nosedir'         - ['+X'|'-X'|'+Y'|'-Y'] direction of nose {default: '+X'}
%   'chaninfo'        - [struct] optional structure containing fields 'nosedir', 'plotrad'
%                       and/or 'chantype'. See these (separate) field definitions above, below.
%                       {default: nosedir +X, plotrad 0.5, all channels}
%   'plotrad'         - [0.15<=float<=1.0] plotting radius = max channel arc_length to plot.
%                       See >> topoplot example. If plotrad > 0.5, chans with arc_length > 0.5
%                       (i.e. below ears-eyes) are plotted in a circular 'skirt' outside the
%                       cartoon head. See 'intrad' below. {default: max(max(chanlocs.radius),0.5);
%                       If the chanlocs structure includes a field chanlocs.plotrad, its value
%                       is used by default}.
%   'headrad'         - [0.15<=float<=1.0] drawing radius (arc_length) for the cartoon head.
%                       NOTE: Only headrad = 0.5 is anatomically correct! 0 -> don't draw head;
%                       'rim' -> show cartoon head at outer edge of the plot {default: 0.5}
%   'intrad'          - [0.15<=float<=1.0] radius of the scalp map interpolation area (square or
%                       disk, see 'intsquare' below). Interpolate electrodes in this area and use
%                       this limit to define boundaries of the scalp map interpolated data matrix
%                       {default: max channel location radius}
%   'intsquare'       - ['on'|'off'] 'on' -> Interpolate values at electrodes located in the whole
%                       square containing the (radius intrad) interpolation disk; 'off' -> Interpolate
%                       values from electrodes shown in the interpolation disk only {default: 'on'}.
%   'conv'            - ['on'|'off'] Show map interpolation only out to the convext hull of
%                       the electrode locations to minimize extrapolation.  {default: 'off'}
%   'noplot'          - ['on'|'off'|[rad theta]] do not plot (but return interpolated data).
%                       Else, if [rad theta] are coordinates of a (possibly missing) channel,
%                       returns interpolated value for channel location.  For more info,
%                       see >> topoplot 'example' {default: 'off'}
%   'verbose'         - ['on'|'off'] comment on operations on command line {default: 'on'}.
%
% Plot detail options:
%   'drawaxis'        - ['on'|'off'] draw axis on the top left corner.
%   'emarker'         - Matlab marker char | {markerchar color size linewidth} char, else cell array
%                       specifying the electrode 'pts' marker. Ex: {'s','r',32,1} -> 32-point solid
%                       red square. {default: {'.','k',[],1} where marker size ([]) depends on the number
%                       of channels plotted}.
%   'emarker2'        - {markchans}|{markchans marker color size linewidth} cell array specifying
%                       an alternate marker for specified 'plotchans'. Ex: {[3 17],'s','g'}
%                       {default: none, or if {markchans} only are specified, then {markchans,'o','r',10,1}}
%   'hcolor'          - color of the cartoon head. Use 'hcolor','none' to plot no head. {default: 'k' = black}
%   'shading'         - 'flat','interp'  {default: 'flat'}
%   'numcontour'      - number of contour lines {default: 6}
%   'contourvals'     - values for contour {default: same as input values}
%   'pmask'           - values for masking topoplot. Array of zeros and 1 of the same size as the input
%                       value array {default: []}
%   'color'           - color of the contours {default: dark grey}
%   'whitebk '        -  ('on'|'off') make the background color white (e.g., to print empty plotgrid channels)
%                       {default: 'off'}
%   'gridscale'       - [int > 32] size (nrows) of interpolated scalp map data matrix {default: 67}
%   'colormap'        -  (n,3) any size colormap {default: existing colormap}
%   'circgrid'        - [int > 100] number of elements (angles) in head and border circles {201}
%
% Dipole plotting options:
%   'dipole'          - [xi yi xe ye ze] plot dipole on the top of the scalp map
%                       from coordinate (xi,yi) to coordinates (xe,ye,ze) (dipole head
%                       model has radius 1). If several rows, plot one dipole per row.
%                       Coordinates returned by dipplot() may be used. Can accept
%                       an EEG.dipfit.model structure (See >> help dipplot).
%                       Ex: ,'dipole',EEG.dipfit.model(17) % Plot dipole(s) for comp. 17.
%   'dipnorm'         - ['on'|'off'] normalize dipole length {default: 'on'}.
%   'diporient'       - [-1|1] invert dipole orientation {default: 1}.
%   'diplen'          - [real] scale dipole length {default: 1}.
%   'dipscale'        - [real] scale dipole size {default: 1}.
%   'dipsphere'       - [real] size of the dipole sphere. {default: 85 mm}.
%   'dipcolor'        - [color] dipole color as Matlab code code or [r g b] vector
%                       {default: 'k' = black}.
% Outputs:
%                   h - handle of the colored surface. If no surface is plotted,
%                       return "gca", the handle of the current plot.
%         grid_or_val - [matrix] the interpolated data image (with off-head points = NaN).
%                       Else, single interpolated value at the specified 'noplot' arg channel
%                       location ([rad theta]), if any.
%     plotrad_or_grid - IF grid image returned above, then the 'plotrad' radius of the grid.
%                       Else, the grid image
%     xmesh, ymesh    - x and y values of the returned grid (above)
%
% Chan_locs format:
%    See >> topoplot 'example'
%
% Examples:
%
%    To plot channel locations only:
%    >> figure; topoplot([],EEG.chanlocs,'style','blank','electrodes','labelpoint','chaninfo',EEG.chaninfo);
%
% Notes: - To change the plot map masking ring to a new figure background color,
%            >> set(findobj(gca,'type','patch'),'facecolor',get(gcf,'color'))
%        - Topoplots may be rotated. From the commandline >> view([deg 90]) {default: [0 90])
%
% Authors: Andy Spydell, Colin Humphries, Arnaud Delorme & Scott Makeig
%          CNL / Salk Institute, 8/1996-/10/2001; SCCN/INC/UCSD, Nov. 2001 -
%
% See also: timtopo(), envtopo()

% Deprecated options:
%           'shrink' - ['on'|'off'|'force'|factor] Deprecated. 'on' -> If max channel arc_length
%                       > 0.5, shrink electrode coordinates towards vertex to plot all channels
%                       by making max arc_length 0.5. 'force' -> Normalize arc_length
%                       so the channel max is 0.5. factor -> Apply a specified shrink
%                       factor (range (0,1) = shrink fraction). {default: 'off'}
%   'electcolor' {'k'}  ... electrode marking details and their {defaults}.
%   'emarker' {'.'}|'emarkersize' {14}|'emarkersizemark' {40}|'efontsize' {var} -
%                       electrode marking details and their {defaults}.
%   'ecolor'          - color of the electrode markers {default: 'k' = black}
%   'interplimits'    - ['electrodes'|'head'] 'electrodes'-> interpolate the electrode grid;
%                       'head'-> interpolate the whole disk {default: 'head'}.

% Unimplemented future options:

% Copyright (C) Colin Humphries & Scott Makeig, CNL / Salk Institute, Aug, 1996
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

% Topoplot Version 2.1
% Early development history:
% Begun by Andy Spydell and Scott Makeig, NHRC,  7-23-96
% 8-96 Revised by Colin Humphries, CNL / Salk Institute, La Jolla CA
%   -changed surf command to imagesc (faster)
%   -can now handle arbitrary scaling of electrode distances
%   -can now handle non integer angles in chan_locs
% 4-4-97 Revised again by Colin Humphries, reformatted by SM
%   -added parameters
%   -changed chan_locs format
% 2-26-98 Revised by Colin
%   -changed image back to surface command
%   -added fill and blank styles
%   -removed extra background colormap entry (now use any colormap)
%   -added parameters for electrode colors and labels
%   -now each topoplot axes use the caxis command again.
%   -removed OUTPUT parameter
% 3-11-98 changed default emarkersize, improve help msg -sm
% 5-24-01 made default emarkersize vary with number of channels -sm
% 01-25-02 reformated help & license, added link -ad
% 03-15-02 added readlocs and the use of eloc input structure -ad
% 03-25-02 added 'labelpoint' options and allow Values=[] -ad &sm
% 03-25-02 added details to "Unknown parameter" warning -sm & ad

function [handle,Zi,grid,Xi,Yi] = topoplot(Values,loc_file,varargin)

%
%%%%%%%%%%%%%%%%%%%%%%%% Set defaults %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
icadefs                 % read defaults MAXTOPOPLOTCHANS and DEFAULT_ELOC and BACKCOLOR
if ~exist('BACKCOLOR')  % if icadefs.m does not define BACKCOLOR
    BACKCOLOR = [.93 .96 1];  % EEGLAB standard
end
whitebk = 'off';  % by default, make gridplot background color = EEGLAB screen background color

plotgrid = 'off';
plotchans = [];
noplot  = 'off';
handle = [];
Zi = [];
chanval = NaN;
rmax = 0.5;             % actual head radius - Don't change this!
INTERPLIMITS = 'head';  % head, electrodes
INTSQUARE = 'on';       % default, interpolate electrodes located though the whole square containing
% the plotting disk
default_intrad = 1;     % indicator for (no) specified intrad
MAPLIMITS = 'absmax';   % absmax, maxmin, [values]
GRID_SCALE = 67;        % plot map on a 67X67 grid
CIRCGRID   = 201;       % number of angles to use in drawing circles
AXHEADFAC = 1.3;        % head to axes scaling factor
CONTOURNUM = 6;         % number of contour levels to plot
STYLE = 'both';         % default 'style': both,straight,fill,contour,blank
HEADCOLOR = [0 0 0];    % default head color (black)
CCOLOR = [0.2 0.2 0.2]; % default contour color
ELECTRODES = [];        % default 'electrodes': on|off|label - set below
MAXDEFAULTSHOWLOCS = 64;% if more channels than this, don't show electrode locations by default
EMARKER = '.';          % mark electrode locations with small disks
ECOLOR = [0 0 0];       % default electrode color = black
EMARKERSIZE = [];       % default depends on number of electrodes, set in code
EMARKERLINEWIDTH = 1;   % default edge linewidth for emarkers
EMARKERSIZE1CHAN = 20;  % default selected channel location marker size
EMARKERCOLOR1CHAN = 'red'; % selected channel location marker color
EMARKER2CHANS = [];      % mark subset of electrode locations with small disks
EMARKER2 = 'o';          % mark subset of electrode locations with small disks
EMARKER2COLOR = 'r';     % mark subset of electrode locations with small disks
EMARKERSIZE2 = 10;      % default selected channel location marker size
EMARKER2LINEWIDTH = 1;
EFSIZE = get(0,'DefaultAxesFontSize'); % use current default fontsize for electrode labels
HLINEWIDTH = 1.7;         % default linewidth for head, nose, ears
BLANKINGRINGWIDTH = .035;% width of the blanking ring
HEADRINGWIDTH    = .007;% width of the cartoon head ring
SHADING = 'flat';       % default 'shading': flat|interp
shrinkfactor = [];      % shrink mode (dprecated)
intrad       = [];      % default interpolation square is to outermost electrode (<=1.0)
plotrad      = [];      % plotting radius ([] = auto, based on outermost channel location)
headrad      = [];      % default plotting radius for cartoon head is 0.5
squeezefac = 1.0;
MINPLOTRAD = 0.15;      % can't make a topoplot with smaller plotrad (contours fail)
VERBOSE = 'off';
MASKSURF = 'off';
CONVHULL = 'off';       % dont mask outside the electrodes convex hull
DRAWAXIS = 'off';
CHOOSECHANTYPE = 0;
ContourVals = Values;
PMASKFLAG   = 0;

%%%%%% Dipole defaults %%%%%%%%%%%%
DIPOLE  = [];
DIPNORM   = 'on';
DIPSPHERE = 85;
DIPLEN    = 1;
DIPSCALE  = 1;
DIPORIENT  = 1;
DIPCOLOR  = [0 0 0];
NOSEDIR   = '+X';
CHANINFO  = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%%%%%%%%%%%%%%%%%%%%%%% Handle arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if nargin< 1
    help topoplot;
    return
end

% calling topoplot from Fieldtrip
% -------------------------------
fieldtrip = 0;
if nargin < 2, loc_file = []; end;
if isstruct(Values) | ~isstruct(loc_file), fieldtrip == 1; end;
if isstr(loc_file), if exist(loc_file) ~= 2, fieldtrip == 1; end; end;
if fieldtrip
    disp('Calling topoplot from Fieldtrip');
    dir1 = which('topoplot');           dir1 = fileparts(dir1);
    dir2 = which('electrodenormalize'); dir2 = fileparts(dir2);
    addpath(dir2);
    try,
        topoplot(Values, loc_file, varargin{:});
    catch,
    end;
    addpath(dir1);
    return;
end;

nargs = nargin;
if nargs == 1
    if isstr(Values)
        if any(strcmp(lower(Values),{'example','demo'}))
            fprintf(['This is an example of an electrode location file,\n',...
                'an ascii file consisting of the following four columns:\n',...
                ' channel_number degrees arc_length channel_name\n\n',...
                'Example:\n',...
                ' 1               -18    .352       Fp1 \n',...
                ' 2                18    .352       Fp2 \n',...
                ' 5               -90    .181       C3  \n',...
                ' 6                90    .181       C4  \n',...
                ' 7               -90    .500       A1  \n',...
                ' 8                90    .500       A2  \n',...
                ' 9              -142    .231       P3  \n',...
                '10               142    .231       P4  \n',...
                '11                 0    .181       Fz  \n',...
                '12                 0    0          Cz  \n',...
                '13               180    .181       Pz  \n\n',...
                ...
                'In topoplot() coordinates, 0 deg. points to the nose, positive\n',...
                'angles point to the right hemisphere, and negative to the left.\n',...
                'The model head sphere has a circumference of 2; the vertex\n',...
                '(Cz) has arc_length 0. Locations with arc_length > 0.5 are below\n',...
                'head center and are plotted outside the head cartoon.\n',...
                'Option plotrad controls how much of this lower-head "skirt" is shown.\n',...
                'Option headrad controls if and where the cartoon head will be drawn.\n',...
                'Option intrad controls how many channels will be included in the interpolation.\n',...
                ])
            return
        end
    end
end
if nargs < 2
    loc_file = DEFAULT_ELOC;
    if ~exist(loc_file)
        fprintf('default locations file "%s" not found - specify chan_locs in topoplot() call.\n',loc_file)
        error(' ')
    end
end
if isempty(loc_file)
    loc_file = 0;
end
if isnumeric(loc_file) & loc_file == 0
    loc_file = DEFAULT_ELOC;
end

if nargs > 2
    if ~(round(nargs/2) == nargs/2)
        error('Odd number of input arguments??')
    end
    for i = 1:2:length(varargin)
        Param = varargin{i};
        Value = varargin{i+1};
        if ~isstr(Param)
            error('Flag arguments must be strings')
        end
        Param = lower(Param);
        switch Param
            case 'conv'
                CONVHULL = lower(Value);
                if ~strcmp(CONVHULL,'on') & ~strcmp(CONVHULL,'off')
                    error('Value of ''conv'' must be ''on'' or ''off''.');
                end
            case 'colormap'
                if size(Value,2)~=3
                    error('Colormap must be a n x 3 matrix')
                end
                colormap(Value)
            case 'intsquare'
                INTSQUARE = lower(Value);
                if ~strcmp(INTSQUARE,'on') & ~strcmp(INTSQUARE,'off')
                    error('Value of ''intsquare'' must be ''on'' or ''off''.');
                end
            case {'interplimits','headlimits'}
                if ~isstr(Value)
                    error('''interplimits'' value must be a string')
                end
                Value = lower(Value);
                if ~strcmp(Value,'electrodes') & ~strcmp(Value,'head')
                    error('Incorrect value for interplimits')
                end
                INTERPLIMITS = Value;
            case 'verbose'
                VERBOSE = Value;
            case 'nosedir'
                NOSEDIR = Value;
                if isempty(strmatch(lower(NOSEDIR), { '+x', '-x', '+y', '-y' }))
                    error('Invalid nose direction');
                end;
            case 'chaninfo'
                CHANINFO = Value;
                if isfield(CHANINFO, 'nosedir'), NOSEDIR      = CHANINFO.nosedir; end;
                if isfield(CHANINFO, 'shrink' ), shrinkfactor = CHANINFO.shrink;  end;
                if isfield(CHANINFO, 'plotrad') & isempty(plotrad), plotrad = CHANINFO.plotrad; end;
                if isfield(CHANINFO, 'chantype')
                    chantype = CHANINFO.chantype;
                    if ischar(chantype), chantype = cellstr(chantype); end
                    CHOOSECHANTYPE = 1;
                end
            case 'chantype'
                chantype = Value;
                CHOOSECHANTYPE = 1;
                if ischar(chantype), chantype = cellstr(chantype); end
                if ~iscell(chantype), error('chantype must be cell array. e.g. {''EEG'', ''EOG''}'); end
            case 'drawaxis'
                DRAWAXIS = Value;
            case 'maplimits'
                MAPLIMITS = Value;
            case 'masksurf'
                MASKSURF = Value;
            case 'circgrid'
                CIRCGRID = Value;
                if isstr(CIRCGRID) | CIRCGRID<100
                    error('''circgrid'' value must be an int > 100');
                end
            case 'style'
                STYLE = lower(Value);
            case 'numcontour'
                CONTOURNUM = Value;
            case 'electrodes'
                ELECTRODES = lower(Value);
                if strcmpi(ELECTRODES,'pointlabels') | strcmpi(ELECTRODES,'ptslabels') ...
                        | strcmpi(ELECTRODES,'labelspts') | strcmpi(ELECTRODES,'ptlabels') ...
                        | strcmpi(ELECTRODES,'labelpts')
                    ELECTRODES = 'labelpoint'; % backwards compatability
                elseif strcmpi(ELECTRODES,'pointnumbers') | strcmpi(ELECTRODES,'ptsnumbers') ...
                        | strcmpi(ELECTRODES,'numberspts') | strcmpi(ELECTRODES,'ptnumbers') ...
                        | strcmpi(ELECTRODES,'numberpts')  | strcmpi(ELECTRODES,'ptsnums')  ...
                        | strcmpi(ELECTRODES,'numspts')
                    ELECTRODES = 'numpoint'; % backwards compatability
                elseif strcmpi(ELECTRODES,'nums')
                    ELECTRODES = 'numbers'; % backwards compatability
                elseif strcmpi(ELECTRODES,'pts')
                    ELECTRODES = 'on'; % backwards compatability
                elseif ~strcmp(ELECTRODES,'off') ...
                        & ~strcmpi(ELECTRODES,'on') ...
                        & ~strcmp(ELECTRODES,'labels') ...
                        & ~strcmpi(ELECTRODES,'numbers') ...
                        & ~strcmpi(ELECTRODES,'labelpoint') ...
                        & ~strcmpi(ELECTRODES,'numpoint')
                    error('Unknown value for keyword ''electrodes''');
                end
            case 'dipole'
                DIPOLE = Value;
            case 'dipsphere'
                DIPSPHERE = Value;
            case 'dipnorm'
                DIPNORM = Value;
            case 'diplen'
                DIPLEN = Value;
            case 'dipscale'
                DIPSCALE = Value;
            case 'contourvals'
                ContourVals = Value;
            case 'pmask'
                ContourVals = Value;
                PMASKFLAG   = 1;
            case 'diporient'
                DIPORIENT = Value;
            case 'dipcolor'
                DIPCOLOR = Value;
            case 'emarker'
                if ischar(Value)
                    EMARKER = Value;
                elseif ~iscell(Value) | length(Value) > 4
                    error('''emarker'' argument must be a cell array {marker color size linewidth}')
                else
                    EMARKER = Value{1};
                end
                if length(Value) > 1
                    ECOLOR = Value{2};
                end
                if length(Value) > 2
                    EMARKERSIZE = Value{3};
                end
                if length(Value) > 3
                    EMARKERLINEWIDTH = Value{4};
                end
            case 'emarker2'
                if ~iscell(Value) | length(Value) > 5
                    error('''emarker2'' argument must be a cell array {chans marker color size linewidth}')
                end
                EMARKER2CHANS = abs(Value{1}); % ignore channels < 0
                if length(Value) > 1
                    EMARKER2 = Value{2};
                end
                if length(Value) > 2
                    EMARKER2COLOR = Value{3};
                end
                if length(Value) > 3
                    EMARKERSIZE2 = Value{4};
                end
                if length(Value) > 4
                    EMARKER2LINEWIDTH = Value{5};
                end
            case 'shrink'
                shrinkfactor = Value;
            case 'intrad'
                intrad = Value;
                if isstr(intrad) | (intrad < MINPLOTRAD | intrad > 1)
                    error('intrad argument should be a number between 0.15 and 1.0');
                end
            case 'plotrad'
                plotrad = Value;
                if isstr(plotrad) | (plotrad < MINPLOTRAD | plotrad > 1)
                    error('plotrad argument should be a number between 0.15 and 1.0');
                end
            case 'headrad'
                headrad = Value;
                if isstr(headrad) & ( strcmpi(headrad,'off') | strcmpi(headrad,'none') )
                    headrad = 0;       % undocumented 'no head' alternatives
                end
                if isempty(headrad) % [] -> none also
                    headrad = 0;
                end
                if ~isstr(headrad)
                    if ~(headrad==0) & (headrad < MINPLOTRAD | headrad>1)
                        error('bad value for headrad');
                    end
                elseif  ~strcmpi(headrad,'rim')
                    error('bad value for headrad');
                end
            case {'headcolor','hcolor'}
                HEADCOLOR = Value;
            case {'contourcolor','ccolor'}
                CCOLOR = Value;
            case {'electcolor','ecolor'}
                ECOLOR = Value;
            case {'emarkersize','emsize'}
                EMARKERSIZE = Value;
            case {'emarkersize1chan','emarkersizemark'}
                EMARKERSIZE1CHAN= Value;
            case {'efontsize','efsize'}
                EFSIZE = Value;
            case 'shading'
                SHADING = lower(Value);
                if ~any(strcmp(SHADING,{'flat','interp'}))
                    error('Invalid shading parameter')
                end
            case 'noplot'
                noplot = Value;
                if ~isstr(noplot)
                    if length(noplot) ~= 2
                        error('''noplot'' location should be [radius, angle]')
                    else
                        chanrad = noplot(1);
                        chantheta = noplot(2);
                        noplot = 'on';
                    end
                end
            case 'gridscale'
                GRID_SCALE = Value;
                if isstr(GRID_SCALE) | GRID_SCALE ~= round(GRID_SCALE) | GRID_SCALE < 32
                    error('''gridscale'' value must be integer > 32.');
                end
            case {'plotgrid','gridplot'}
                plotgrid = 'on';
                gridchans = Value;
            case 'plotchans'
                plotchans = Value(:);
                if find(plotchans<=0)
                    error('''plotchans'' values must be > 0');
                end
                % if max(abs(plotchans))>max(Values) | max(abs(plotchans))>length(Values) -sm ???
            case {'whitebk','whiteback','forprint'}
                whitebk = Value;
            otherwise
                error(['Unknown input parameter ''' Param ''' ???'])
        end
    end
end
if strcmpi(whitebk, 'on')
    BACKCOLOR = [ 1 1 1 ];
end;

cmap = colormap;
cmaplen = size(cmap,1);

if strcmp(STYLE,'blank')    % else if Values holds numbers of channels to mark
    if length(Values) < length(loc_file)
        ContourVals = zeros(1,length(loc_file));
        ContourVals(Values) = 1;
        Values = ContourVals;
    end;
end;

%
%%%%%%%%%%%%%%%%%%%%%%%%%%% test args for plotting an electrode grid %%%%%%%%%%%%%%%%%%%%%%
%
if strcmp(plotgrid,'on')
    STYLE = 'grid';
    gchans = sort(find(abs(gridchans(:))>0));
    
    % if setdiff(gchans,unique(gchans))
    %      fprintf('topoplot() warning: ''plotgrid'' channel matrix has duplicate channels\n');
    % end
    
    if ~isempty(plotchans)
        if intersect(gchans,abs(plotchans))
            fprintf('topoplot() warning: ''plotgrid'' and ''plotchans'' have channels in common\n');
        end
    end
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%% misc arg tests %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if isempty(ELECTRODES)                     % if electrode labeling not specified
    if length(Values) > MAXDEFAULTSHOWLOCS   % if more channels than default max
        ELECTRODES = 'off';                    % don't show electrodes
    else                                     % else if fewer chans,
        ELECTRODES = 'on';                     % do
    end
end

if isempty(Values)
    STYLE = 'blank';
end
[r,c] = size(Values);
if r>1 & c>1,
    error('input data must be a single vector');
end
Values = Values(:); % make Values a column vector
ContourVals = ContourVals(:); % values for contour

if ~isempty(intrad) & ~isempty(plotrad) & intrad < plotrad
    error('intrad must be >= plotrad');
end

if ~strcmpi(STYLE,'grid')                     % if not plot grid only
    
    %
    %%%%%%%%%%%%%%%%%%%% Read the channel location information %%%%%%%%%%%%%%%%%%%%%%%%
    %
    if isstr(loc_file)
        [tmpeloc labels Th Rd indices] = readlocs( loc_file);
    elseif isstruct(loc_file) % a locs struct
        [tmpeloc labels Th Rd indices] = readlocs( loc_file );
        % Note: Th and Rd correspond to indices channels-with-coordinates only
    else
        error('loc_file must be a EEG.locs struct or locs filename');
    end
    Th = pi/180*Th;                              % convert degrees to radians
    allchansind = 1:length(Th);
    
    
    if ~isempty(plotchans)
        if max(plotchans) > length(Th)
            error('''plotchans'' values must be <= max channel index');
        end
    end
    
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% channels to plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    if ~isempty(plotchans)
        plotchans = intersect(plotchans, indices);
    end;
    if ~isempty(Values) & ~strcmpi( STYLE, 'blank') & isempty(plotchans)
        plotchans = indices;
    end
    if isempty(plotchans) & strcmpi( STYLE, 'blank')
        plotchans = indices;
    end
    
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%% filter for channel type(s), if specified %%%%%%%%%%%%%%%%%%%%%
    %
    
    if CHOOSECHANTYPE,
        newplotchans = eeg_chantype(loc_file,chantype);
        plotchans = intersect(newplotchans, plotchans);
    end
    
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%% filter channels used for components %%%%%%%%%%%%%%%%%%%%%
    %
    if isfield(CHANINFO, 'icachansind') & ~isempty(Values) & length(Values) ~= length(tmpeloc)
        
        % test if ICA component
        % ---------------------
        if length(CHANINFO.icachansind) == length(Values)
            
            % if only a subset of channels are to be plotted
            % and ICA components also use a subject of channel
            % we must find the new indices for these channels
            
            plotchans = intersect(CHANINFO.icachansind, plotchans);
            tmpvals   = zeros(1, length(tmpeloc));
            tmpvals(CHANINFO.icachansind) = Values;
            Values    = tmpvals;
            tmpvals   = zeros(1, length(tmpeloc));
            tmpvals(CHANINFO.icachansind) = ContourVals;
            ContourVals = tmpvals;
            
        end;
    end;
    
    %
    %%%%%%%%%%%%%%%%%%% last channel is reference? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    if length(tmpeloc) == length(Values) + 1 % remove last channel if necessary
        % (common reference channel)
        if plotchans(end) == length(tmpeloc)
            plotchans(end) = [];
        end;
        
    end;
    
    %
    %%%%%%%%%%%%%%%%%%% remove infinite and NaN values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    if length(Values) > 1
        inds          = union(find(isnan(Values)), find(isinf(Values))); % NaN and Inf values
        plotchans     = setdiff(plotchans, inds);
    end;
    if strcmp(plotgrid,'on')
        plotchans = setxor(plotchans,gchans);   % remove grid chans from head plotchans
    end
    
    [x,y]     = pol2cart(Th,Rd);  % transform electrode locations from polar to cartesian coordinates
    plotchans = abs(plotchans);   % reverse indicated channel polarities
    allchansind = allchansind(plotchans);
    Th        = Th(plotchans);
    Rd        = Rd(plotchans);
    x         = x(plotchans);
    y         = y(plotchans);
    labels    = labels(plotchans); % remove labels for electrodes without locations
    labels    = strvcat(labels); % make a label string matrix
    if ~isempty(Values) & length(Values) > 1
        Values      = Values(plotchans);
        ContourVals = ContourVals(plotchans);
    end;
    
    %
    %%%%%%%%%%%%%%%%%% Read plotting radius from chanlocs  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    if isempty(plotrad) & isfield(tmpeloc, 'plotrad'),
        plotrad = tmpeloc(1).plotrad;
        if isstr(plotrad)                        % plotrad shouldn't be a string
            plotrad = str2num(plotrad)           % just checking
        end
        if plotrad < MINPLOTRAD | plotrad > 1.0
            fprintf('Bad value (%g) for plotrad.\n',plotrad);
            error(' ');
        end
        if strcmpi(VERBOSE,'on') & ~isempty(plotrad)
            fprintf('Plotting radius plotrad (%g) set from EEG.chanlocs.\n',plotrad);
        end
    end;
    if isempty(plotrad)
        plotrad = min(1.0,max(Rd)*1.02);            % default: just outside the outermost electrode location
        plotrad = max(plotrad,0.5);                 % default: plot out to the 0.5 head boundary
    end                                           % don't plot channels with Rd > 1 (below head)
    
    if isempty(intrad)
        default_intrad = 1;     % indicator for (no) specified intrad
        intrad = min(1.0,max(Rd)*1.02);             % default: just outside the outermost electrode location
    else
        default_intrad = 0;                         % indicator for (no) specified intrad
        if plotrad > intrad
            plotrad = intrad;
        end
    end                                           % don't interpolate channels with Rd > 1 (below head)
    if isstr(plotrad) | plotrad < MINPLOTRAD | plotrad > 1.0
        error('plotrad must be between 0.15 and 1.0');
    end
    
    %
    %%%%%%%%%%%%%%%%%%%%%%% Set radius of head cartoon %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    if isempty(headrad)  % never set -> defaults
        if plotrad >= rmax
            headrad = rmax;  % (anatomically correct)
        else % if plotrad < rmax
            headrad = 0;    % don't plot head
            if strcmpi(VERBOSE, 'on')
                fprintf('topoplot(): not plotting cartoon head since plotrad (%5.4g) < 0.5\n',...
                    plotrad);
            end
        end
    elseif strcmpi(headrad,'rim') % force plotting at rim of map
        headrad = plotrad;
    end
    
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Shrink mode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    if ~isempty(shrinkfactor) | isfield(tmpeloc, 'shrink'),
        if isempty(shrinkfactor) & isfield(tmpeloc, 'shrink'),
            shrinkfactor = tmpeloc(1).shrink;
            if strcmpi(VERBOSE,'on')
                if isstr(shrinkfactor)
                    fprintf('Automatically shrinking coordinates to lie above the head perimter.\n');
                else
                    fprintf('Automatically shrinking coordinates by %3.2f\n', shrinkfactor);
                end;
            end
        end;
        
        if isstr(shrinkfactor)
            if strcmpi(shrinkfactor, 'on') | strcmpi(shrinkfactor, 'force') | strcmpi(shrinkfactor, 'auto')
                if abs(headrad-rmax) > 1e-2
                    fprintf('     NOTE -> the head cartoon will NOT accurately indicate the actual electrode locations\n');
                end
                if strcmpi(VERBOSE,'on')
                    fprintf('     Shrink flag -> plotting cartoon head at plotrad\n');
                end
                headrad = plotrad; % plot head around outer electrodes, no matter if 0.5 or not
            end
        else % apply shrinkfactor
            plotrad = rmax/(1-shrinkfactor);
            headrad = plotrad;  % make deprecated 'shrink' mode plot
            if strcmpi(VERBOSE,'on')
                fprintf('    %g%% shrink  applied.');
                if abs(headrad-rmax) > 1e-2
                    fprintf(' Warning: With this "shrink" setting, the cartoon head will NOT be anatomically correct.\n');
                else
                    fprintf('\n');
                end
            end
        end
    end; % if shrink
    
    %
    %%%%%%%%%%%%%%%%% Issue warning if headrad ~= rmax  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    
    if headrad ~= 0.5 & strcmpi(VERBOSE, 'on')
        fprintf('     NB: Plotting map using ''plotrad'' %-4.3g,',plotrad);
        fprintf(    ' ''headrad'' %-4.3g\n',headrad);
        fprintf('Warning: The plotting radius of the cartoon head is NOT anatomically correct (0.5).\n')
    end
    %
    %%%%%%%%%%%%%%%%%%%%% Find plotting channels  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    
    pltchans = find(Rd <= plotrad); % plot channels inside plotting circle
    
    if strcmpi(INTSQUARE,'on') % interpolate channels in the radius intrad square
        intchans = find(x <= intrad & y <= intrad); % interpolate and plot channels inside interpolation square
    else
        intchans = find(Rd <= intrad); % interpolate channels in the radius intrad circle only
    end
    
    %
    %%%%%%%%%%%%%%%%%%%%% Eliminate channels not plotted  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    
    allx      = x;
    ally      = y;
    intchans; % interpolate using only the 'intchans' channels
    pltchans; % plot using only indicated 'plotchans' channels
    
    if length(pltchans) < length(Rd) & strcmpi(VERBOSE, 'on')
        fprintf('Interpolating %d and plotting %d of the %d scalp electrodes.\n', ...
            length(intchans),length(pltchans),length(Rd));
    end;
    
    
    % fprintf('topoplot(): plotting %d channels\n',length(pltchans));
    if ~isempty(EMARKER2CHANS)
        if strcmpi(STYLE,'blank')
            error('emarker2 not defined for style ''blank'' - use marking channel numbers in place of data');
        else % mark1chans and mark2chans are subsets of pltchans for markers 1 and 2
            [tmp1 mark1chans tmp2] = setxor(pltchans,EMARKER2CHANS);
            [tmp3 tmp4 mark2chans] = intersect(EMARKER2CHANS,pltchans);
        end
    end
    
    if ~isempty(Values)
        if length(Values) == length(Th)  % if as many map Values as channel locs
            intValues      = Values(intchans);
            intContourVals = ContourVals(intchans);
            Values         = Values(pltchans);
            ContourVals    = ContourVals(pltchans);
        end;
    end;   % now channel parameters and values all refer to plotting channels only
    
    allchansind = allchansind(pltchans);
    intTh = Th(intchans);           % eliminate channels outside the interpolation area
    intRd = Rd(intchans);
    intx  = x(intchans);
    inty  = y(intchans);
    Th    = Th(pltchans);              % eliminate channels outside the plotting area
    Rd    = Rd(pltchans);
    x     = x(pltchans);
    y     = y(pltchans);
    
    labels= labels(pltchans,:);
    %
    %%%%%%%%%%%%%%% Squeeze channel locations to <= rmax %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    
    squeezefac = rmax/plotrad;
    intRd = intRd*squeezefac; % squeeze electrode arc_lengths towards the vertex
    Rd = Rd*squeezefac;       % squeeze electrode arc_lengths towards the vertex
    % to plot all inside the head cartoon
    intx = intx*squeezefac;
    inty = inty*squeezefac;
    x    = x*squeezefac;
    y    = y*squeezefac;
    allx    = allx*squeezefac;
    ally    = ally*squeezefac;
    % Note: Now outermost channel will be plotted just inside rmax
    
else % if strcmpi(STYLE,'grid')
    intx = rmax; inty=rmax;
end % if ~strcmpi(STYLE,'grid')

%
%%%%%%%%%%%%%%%% rotate channels based on chaninfo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if strcmpi(lower(NOSEDIR), '+x')
    rotate = 0;
else
    if strcmpi(lower(NOSEDIR), '+y')
        rotate = 3*pi/2;
    elseif strcmpi(lower(NOSEDIR), '-x')
        rotate = pi;
    else rotate = pi/2;
    end;
    allcoords = (inty + intx*sqrt(-1))*exp(sqrt(-1)*rotate);
    intx = imag(allcoords);
    inty = real(allcoords);
    allcoords = (ally + allx*sqrt(-1))*exp(sqrt(-1)*rotate);
    allx = imag(allcoords);
    ally = real(allcoords);
    allcoords = (y + x*sqrt(-1))*exp(sqrt(-1)*rotate);
    x = imag(allcoords);
    y = real(allcoords);
end;

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Make the plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if ~strcmpi(STYLE,'blank') % if draw interpolated scalp map
    if ~strcmpi(STYLE,'grid') %  not a rectangular channel grid
        %
        %%%%%%%%%%%%%%%% Find limits for interpolation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        if default_intrad % if no specified intrad
            if strcmpi(INTERPLIMITS,'head') % intrad is 'head'
                xmin = min(-rmax,min(intx)); xmax = max(rmax,max(intx));
                ymin = min(-rmax,min(inty)); ymax = max(rmax,max(inty));
                
            else % INTERPLIMITS = rectangle containing electrodes -- DEPRECATED OPTION!
                xmin = max(-rmax,min(intx)); xmax = min(rmax,max(intx));
                ymin = max(-rmax,min(inty)); ymax = min(rmax,max(inty));
            end
        else % some other intrad specified
            xmin = -intrad*squeezefac; xmax = intrad*squeezefac;   % use the specified intrad value
            ymin = -intrad*squeezefac; ymax = intrad*squeezefac;
        end
        %
        %%%%%%%%%%%%%%%%%%%%%%% Interpolate scalp map data %%%%%%%%%%%%%%%%%%%%%%%%
        %
        xi = linspace(xmin,xmax,GRID_SCALE);   % x-axis description (row vector)
        yi = linspace(ymin,ymax,GRID_SCALE);   % y-axis description (row vector)
        ws = warning('off','MATLAB:griddata:DuplicateDataPoints');
        try
            [Xi,Yi,Zi] = griddata(inty,intx,double(intValues),yi',xi,'v4'); % interpolate data
            [Xi,Yi,ZiC] = griddata(inty,intx,double(intContourVals),yi',xi,'v4'); % interpolate data
        catch,
            [Xi,Yi,Zi] = griddata(inty,intx,intValues',yi,xi'); % interpolate data (Octave)
            [Xi,Yi,ZiC] = griddata(inty,intx,intContourVals',yi,xi'); % interpolate data
        end;
        warning(ws)
        %
        %%%%%%%%%%%%%%%%%%%%%%% Mask out data outside the head %%%%%%%%%%%%%%%%%%%%%
        %
        mask = (sqrt(Xi.^2 + Yi.^2) <= rmax); % mask outside the plotting circle
        ii = find(mask == 0);
        Zi(ii)  = NaN;                         % mask non-plotting voxels with NaNs
        ZiC(ii) = NaN;                         % mask non-plotting voxels with NaNs
        grid = plotrad;                       % unless 'noplot', then 3rd output arg is plotrad
        %
        %%%%%%%%%% Return interpolated value at designated scalp location %%%%%%%%%%
        %
        if exist('chanrad')   % optional first argument to 'noplot'
            chantheta = (chantheta/360)*2*pi;
            chancoords = round(ceil(GRID_SCALE/2)+GRID_SCALE/2*2*chanrad*[cos(-chantheta),...
                -sin(-chantheta)]);
            if chancoords(1)<1 ...
                    | chancoords(1) > GRID_SCALE ...
                    | chancoords(2)<1 ...
                    | chancoords(2)>GRID_SCALE
                error('designated ''noplot'' channel out of bounds')
            else
                chanval = Zi(chancoords(1),chancoords(2));
                grid = Zi;
                Zi = chanval;  % return interpolated value instead of Zi
            end
        end
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%% Return interpolated image only  %%%%%%%%%%%%%%%%%
        %
        if strcmpi(noplot, 'on')
            if strcmpi(VERBOSE,'on')
                fprintf('topoplot(): no plot requested.\n')
            end
            return;
        end
        %
        %%%%%%%%%%%%%%%%%%%%%%% Calculate colormap limits %%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        if isstr(MAPLIMITS)
            if strcmp(MAPLIMITS,'absmax')
                amax = max(max(abs(Zi)));
                amin = -amax;
            elseif strcmp(MAPLIMITS,'maxmin') | strcmp(MAPLIMITS,'minmax')
                amin = min(min(Zi));
                amax = max(max(Zi));
            else
                error('unknown ''maplimits'' value.');
            end
        elseif length(MAPLIMITS) == 2
            amin = MAPLIMITS(1);
            amax = MAPLIMITS(2);
        else
            error('unknown ''maplimits'' value');
        end
        delta = xi(2)-xi(1); % length of grid entry
        
    end % if ~strcmpi(STYLE,'grid')
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%% Scale the axes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %cla  % clear current axis
    hold on
    h = gca; % uses current axes
    
    % instead of default larger AXHEADFAC
    if squeezefac<0.92 & plotrad-headrad > 0.05  % (size of head in axes)
        AXHEADFAC = 1.05;     % do not leave room for external ears if head cartoon
        % shrunk enough by the 'skirt' option
    end
    
    set(gca,'Xlim',[-rmax rmax]*AXHEADFAC,'Ylim',[-rmax rmax]*AXHEADFAC);
    % specify size of head axes in gca
    
    unsh = (GRID_SCALE+1)/GRID_SCALE; % un-shrink the effects of 'interp' SHADING
    
    %
    %%%%%%%%%%%%%%%%%%%%%%%% Plot grid only %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    if strcmpi(STYLE,'grid')                     % plot grid only
        
        %
        % The goal below is to make the grid cells square - not yet achieved in all cases? -sm
        %
        g1 = size(gridchans,1);
        g2 = size(gridchans,2);
        gmax = max([g1 g2]);
        Xi = linspace(-rmax*g2/gmax,rmax*g2/gmax,g1+1);
        Xi = Xi+rmax/g1; Xi = Xi(1:end-1);
        Yi = linspace(-rmax*g1/gmax,rmax*g1/gmax,g2+1);
        Yi = Yi+rmax/g2; Yi = Yi(1:end-1); Yi = Yi(end:-1:1); % by trial and error!
        %
        %%%%%%%%%%% collect the gridchans values %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        gridvalues = zeros(size(gridchans));
        for j=1:size(gridchans,1)
            for k=1:size(gridchans,2)
                gc = gridchans(j,k);
                if gc > 0
                    gridvalues(j,k) = Values(gc);
                elseif gc < 0
                    gridvalues(j,k) = -Values(gc);
                else
                    gridvalues(j,k) = nan; % not-a-number = no value
                end
            end
        end
        %
        %%%%%%%%%%% reset color limits for grid plot %%%%%%%%%%%%%%%%%%%%%%%%%
        %
        if isstr(MAPLIMITS)
            if strcmp(MAPLIMITS,'maxmin') | strcmp(MAPLIMITS,'minmax')
                amin = min(min(gridvalues(~isnan(gridvalues))));
                amax = max(max(gridvalues(~isnan(gridvalues))));
            elseif strcmp(MAPLIMITS,'absmax')
                % 11/21/2005 Toby edit
                % This should now work as specified. Before it only crashed (using
                % "plotgrid" and "maplimits>absmax" options).
                amax = max(max(abs(gridvalues(~isnan(gridvalues)))));
                amin = -amax;
                %amin = -max(max(abs([amin amax])));
                %amax = max(max(abs([amin amax])));
            else
                error('unknown ''maplimits'' value');
            end
        elseif length(MAPLIMITS) == 2
            amin = MAPLIMITS(1);
            amax = MAPLIMITS(2);
        else
            error('unknown ''maplimits'' value');
        end
        %
        %%%%%%%%%% explicitly compute grid colors, allowing BACKCOLOR  %%%%%%
        %
        gridvalues = 1+floor(cmaplen*(gridvalues-amin)/(amax-amin));
        gridvalues(find(gridvalues == cmaplen+1)) = cmaplen;
        gridcolors = zeros([size(gridvalues),3]);
        for j=1:size(gridchans,1)
            for k=1:size(gridchans,2)
                if ~isnan(gridvalues(j,k))
                    gridcolors(j,k,:) = cmap(gridvalues(j,k),:);
                else
                    if strcmpi(whitebk,'off')
                        gridcolors(j,k,:) = BACKCOLOR; % gridchans == 0 -> background color
                        % This allows the plot to show 'space' between separate sub-grids or strips
                    else % 'on'
                        gridcolors(j,k,:) = [1 1 1]; BACKCOLOR; % gridchans == 0 -> white for printing
                    end
                end
            end
        end
        
        %
        %%%%%%%%%% draw the gridplot image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        handle=imagesc(Xi,Yi,gridcolors); % plot grid with explicit colors
        axis square
        
        %
        %%%%%%%%%%%%%%%%%%%%%%%% Plot map contours only %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
    elseif strcmp(STYLE,'contour')                     % plot surface contours only
        [cls chs] = contour(Xi,Yi,ZiC,CONTOURNUM,'k');
        % for h=chs, set(h,'color',CCOLOR); end
        %
        %%%%%%%%%%%%%%%%%%%%%%%% Else plot map and contours %%%%%%%%%%%%%%%%%%%%%%%%%
        %
    elseif strcmp(STYLE,'both')  % plot interpolated surface and surface contours
        if strcmp(SHADING,'interp')
            tmph = surface(Xi*unsh,Yi*unsh,zeros(size(Zi))-0.1,Zi,...
                'EdgeColor','none','FaceColor',SHADING);
        else % SHADING == 'flat'
            tmph = surface(Xi-delta/2,Yi-delta/2,zeros(size(Zi))-0.1,Zi,...
                'EdgeColor','none','FaceColor',SHADING);
        end
        if strcmpi(MASKSURF, 'on')
            set(tmph, 'visible', 'off');
            handle = tmph;
        end;
        
        warning off;
        if ~PMASKFLAG
            [cls chs] = contour(Xi,Yi,ZiC,CONTOURNUM,'k');
        else
            ZiC(find(ZiC > 0.5 )) = NaN;
            [cls chs] = contourf(Xi,Yi,ZiC,0,'k');
            subh = get(chs, 'children');
            for indsubh = 1:length(subh)
                numfaces = size(get(subh(indsubh), 'XData'),1);
                set(subh(indsubh), 'FaceVertexCData', ones(numfaces,3), 'Cdatamapping', 'direct', 'facealpha', 0.5, 'linewidth', 2);
            end;
        end;
        for h=chs, set(h,'color',CCOLOR); end
        warning on;
        %
        %%%%%%%%%%%%%%%%%%%%%%%% Else plot map only %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
    elseif strcmp(STYLE,'straight') | strcmp(STYLE,'map') % 'straight' was former arg
        
        if strcmp(SHADING,'interp') % 'interp' mode is shifted somehow... but how?
            tmph = surface(Xi*unsh,Yi*unsh,zeros(size(Zi)),Zi,'EdgeColor','none',...
                'FaceColor',SHADING);
        else
            tmph = surface(Xi-delta/2,Yi-delta/2,zeros(size(Zi)),Zi,'EdgeColor','none',...
                'FaceColor',SHADING);
        end
        if strcmpi(MASKSURF, 'on')
            set(tmph, 'visible', 'off');
            handle = tmph;
        end;
        %
        %%%%%%%%%%%%%%%%%% Else fill contours with uniform colors  %%%%%%%%%%%%%%%%%%
        %
    elseif strcmp(STYLE,'fill')
        [cls chs] = contourf(Xi,Yi,Zi,CONTOURNUM,'k');
        
        % for h=chs, set(h,'color',CCOLOR); end
        %     <- 'not line objects.' Why does 'both' work above???
        
    else
        error('Invalid style')
    end
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set color axis  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    caxis([amin amax]); % set coloraxis
    
else % if STYLE 'blank'
    %
    %%%%%%%%%%%%%%%%%%%%%%% Draw blank head %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    if strcmpi(noplot, 'on')
        if strcmpi(VERBOSE,'on')
            fprintf('topoplot(): no plot requested.\n')
        end
        return;
    end
    %cla
    hold on
    
    set(gca,'Xlim',[-rmax rmax]*AXHEADFAC,'Ylim',[-rmax rmax]*AXHEADFAC)
    % pos = get(gca,'position');
    % fprintf('Current axes size %g,%g\n',pos(3),pos(4));
    
    if strcmp(ELECTRODES,'labelpoint') |  strcmp(ELECTRODES,'numpoint')
        text(-0.6,-0.6, ...
            [ int2str(length(Rd)) ' of ' int2str(length(tmpeloc)) ' electrode locations shown']);
        text(-0.6,-0.7, [ 'Click on electrodes to toggle name/number']);
        tl = title('Channel locations');
        set(tl, 'fontweight', 'bold');
    end;
end % STYLE 'blank'

if exist('handle') ~= 1
    handle = gca;
end;

if ~strcmpi(STYLE,'grid')                     % if not plot grid only
    
    %
    %%%%%%%%%%%%%%%%%%% Plot filled ring to mask jagged grid boundary %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    hwidth = HEADRINGWIDTH;                   % width of head ring
    hin  = squeezefac*headrad*(1- hwidth/2);  % inner head ring radius
    
    if strcmp(SHADING,'interp')
        rwidth = BLANKINGRINGWIDTH*1.3;             % width of blanking outer ring
    else
        rwidth = BLANKINGRINGWIDTH;         % width of blanking outer ring
    end
    rin    =  rmax*(1-rwidth/2);              % inner ring radius
    if hin>rin
        rin = hin;                              % dont blank inside the head ring
    end
    
    if strcmp(CONVHULL,'on') %%%%%%%%% mask outside the convex hull of the electrodes %%%%%%%%%
        cnv = convhull(allx,ally);
        cnvfac = round(CIRCGRID/length(cnv)); % spline interpolate the convex hull
        if cnvfac < 1, cnvfac=1; end;
        CIRCGRID = cnvfac*length(cnv);
        
        startangle = atan2(allx(cnv(1)),ally(cnv(1)));
        circ = linspace(0+startangle,2*pi+startangle,CIRCGRID);
        rx = sin(circ);
        ry = cos(circ);
        
        allx = allx(:)';  % make x (elec locations; + to nose) a row vector
        ally = ally(:)';  % make y (elec locations, + to r? ear) a row vector
        erad = sqrt(allx(cnv).^2+ally(cnv).^2);  % convert to polar coordinates
        eang = atan2(allx(cnv),ally(cnv));
        eang = unwrap(eang);
        eradi =spline(linspace(0,1,3*length(cnv)), [erad erad erad], ...
            linspace(0,1,3*length(cnv)*cnvfac));
        eangi =spline(linspace(0,1,3*length(cnv)), [eang+2*pi eang eang-2*pi], ...
            linspace(0,1,3*length(cnv)*cnvfac));
        xx = eradi.*sin(eangi);           % convert back to rect coordinates
        yy = eradi.*cos(eangi);
        yy = yy(CIRCGRID+1:2*CIRCGRID);
        xx = xx(CIRCGRID+1:2*CIRCGRID);
        eangi = eangi(CIRCGRID+1:2*CIRCGRID);
        eradi = eradi(CIRCGRID+1:2*CIRCGRID);
        xx = xx*1.02; yy = yy*1.02;           % extend spline outside electrode marks
        
        splrad = sqrt(xx.^2+yy.^2);           % arc radius of spline points (yy,xx)
        oob = find(splrad >= rin);            %  enforce an upper bound on xx,yy
        xx(oob) = rin*xx(oob)./splrad(oob);   % max radius = rin
        yy(oob) = rin*yy(oob)./splrad(oob);   % max radius = rin
        
        splrad = sqrt(xx.^2+yy.^2);           % arc radius of spline points (yy,xx)
        oob = find(splrad < hin);             % don't let splrad be inside the head cartoon
        xx(oob) = hin*xx(oob)./splrad(oob);   % min radius = hin
        yy(oob) = hin*yy(oob)./splrad(oob);   % min radius = hin
        
        ringy = [[ry(:)' ry(1) ]*(rin+rwidth) yy yy(1)];
        ringx = [[rx(:)' rx(1) ]*(rin+rwidth) xx xx(1)];
        
        ringh2= patch(ringy,ringx,ones(size(ringy)),get(gcf,'color'),'edgecolor','none','tag','toporingmask'); hold on
        
        % plot(ry*rmax,rx*rmax,'b') % debugging line
        
    else %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mask the jagged border around rmax %%%%%%%%%%%%%%%5%%%%%%
        
        circ = linspace(0,2*pi,CIRCGRID);
        rx = sin(circ);
        ry = cos(circ);
        ringx = [[rx(:)' rx(1) ]*(rin+rwidth)  [rx(:)' rx(1)]*rin];
        ringy = [[ry(:)' ry(1) ]*(rin+rwidth)  [ry(:)' ry(1)]*rin];
        
        if ~strcmpi(STYLE,'blank')
            ringh= patch(ringx,ringy,0.01*ones(size(ringx)),BACKCOLOR,'edgecolor','none','tag','toporingmask'); hold on
        end
        % plot(ry*rmax,rx*rmax,'b') % debugging line
    end
    
    %f1= fill(rin*[rx rX],rin*[ry rY],BACKCOLOR,'edgecolor',BACKCOLOR); hold on
    %f2= fill(rin*[rx rX*(1+rwidth)],rin*[ry rY*(1+rwidth)],BACKCOLOR,'edgecolor',BACKCOLOR);
    
    % Former line-style border smoothing - width did not scale with plot
    %  brdr=plot(1.015*cos(circ).*rmax,1.015*sin(circ).*rmax,...      % old line-based method
    %      'color',HEADCOLOR,'Linestyle','-','LineWidth',HLINEWIDTH);    % plot skirt outline
    %  set(brdr,'color',BACKCOLOR,'linewidth',HLINEWIDTH + 4);        % hide the disk edge jaggies
    
    %
    %%%%%%%%%%%%%%%%%%%%%%%%% Plot cartoon head, ears, nose %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    if headrad > 0                         % if cartoon head to be plotted
        %
        %%%%%%%%%%%%%%%%%%% Plot head outline %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        headx = [[rx(:)' rx(1) ]*(hin+hwidth)  [rx(:)' rx(1)]*hin];
        heady = [[ry(:)' ry(1) ]*(hin+hwidth)  [ry(:)' ry(1)]*hin];
        
        if ~isstr(HEADCOLOR) | ~strcmpi(HEADCOLOR,'none')
            ringh= patch(headx,heady,ones(size(headx)),HEADCOLOR,'edgecolor',HEADCOLOR); hold on
        end
        
        % rx = sin(circ); rX = rx(end:-1:1);
        % ry = cos(circ); rY = ry(end:-1:1);
        % for k=2:2:CIRCGRID
        %   rx(k) = rx(k)*(1+hwidth);
        %   ry(k) = ry(k)*(1+hwidth);
        % end
        % f3= fill(hin*[rx rX],hin*[ry rY],HEADCOLOR,'edgecolor',HEADCOLOR); hold on
        % f4= fill(hin*[rx rX*(1+hwidth)],hin*[ry rY*(1+hwidth)],HEADCOLOR,'edgecolor',HEADCOLOR);
        
        % Former line-style head
        %  plot(cos(circ).*squeezefac*headrad,sin(circ).*squeezefac*headrad,...
        %      'color',HEADCOLOR,'Linestyle','-','LineWidth',HLINEWIDTH);    % plot head outline
        
        %
        %%%%%%%%%%%%%%%%%%% Plot ears and nose %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        base  = rmax-.0046;
        basex = 0.18*rmax;                   % nose width
        tip   = 1.15*rmax;
        tiphw = .04*rmax;                    % nose tip half width
        tipr  = .01*rmax;                    % nose tip rounding
        q = .04; % ear lengthening
        EarX  = [.497-.005  .510  .518  .5299 .5419  .54    .547   .532   .510   .489-.005]; % rmax = 0.5
        EarY  = [q+.0555 q+.0775 q+.0783 q+.0746 q+.0555 -.0055 -.0932 -.1313 -.1384 -.1199];
        sf    = headrad/plotrad;                                          % squeeze the model ears and nose
        % by this factor
        if ~isstr(HEADCOLOR) | ~strcmpi(HEADCOLOR,'none')
            plot3([basex;tiphw;0;-tiphw;-basex]*sf,[base;tip-tipr;tip;tip-tipr;base]*sf,...
                2*ones(size([basex;tiphw;0;-tiphw;-basex])),...
                'Color',HEADCOLOR,'LineWidth',HLINEWIDTH);                 % plot nose
            plot3(EarX*sf,EarY*sf,2*ones(size(EarX)),'color',HEADCOLOR,'LineWidth',HLINEWIDTH)    % plot left ear
            plot3(-EarX*sf,EarY*sf,2*ones(size(EarY)),'color',HEADCOLOR,'LineWidth',HLINEWIDTH)   % plot right ear
        end
    end
    
    %
    % %%%%%%%%%%%%%%%%%%% Show electrode information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    plotax = gca;
    axis square                                           % make plotax square
    axis off
    
    pos = get(gca,'position');
    xlm = get(gca,'xlim');
    ylm = get(gca,'ylim');
    % textax = axes('position',pos,'xlim',xlm,'ylim',ylm);  % make new axes so clicking numbers <-> labels
    % will work inside head cartoon patch
    % axes(textax);
    axis square                                           % make textax square
    
    pos = get(gca,'position');
    set(plotax,'position',pos);
    
    xlm = get(gca,'xlim');
    set(plotax,'xlim',xlm);
    
    ylm = get(gca,'ylim');
    set(plotax,'ylim',ylm);                               % copy position and axis limits again
    
    axis equal;
    set(gca, 'xlim', [-0.525 0.525]); set(plotax, 'xlim', [-0.525 0.525]);
    set(gca, 'ylim', [-0.525 0.525]); set(plotax, 'ylim', [-0.525 0.525]);
    
    %get(textax,'pos')    % test if equal!
    %get(plotax,'pos')
    %get(textax,'xlim')
    %get(plotax,'xlim')
    %get(textax,'ylim')
    %get(plotax,'ylim')
    
    if isempty(EMARKERSIZE)
        EMARKERSIZE = 10;
        if length(y)>=160
            EMARKERSIZE = 3;
        elseif length(y)>=128
            EMARKERSIZE = 3;
        elseif length(y)>=100
            EMARKERSIZE = 3;
        elseif length(y)>=80
            EMARKERSIZE = 4;
        elseif length(y)>=64
            EMARKERSIZE = 5;
        elseif length(y)>=48
            EMARKERSIZE = 6;
        elseif length(y)>=32
            EMARKERSIZE = 8;
        end
    end
    %
    %%%%%%%%%%%%%%%%%%%%%%%% Mark electrode locations only %%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    ELECTRODE_HEIGHT = 2.1;  % z value for plotting electrode information (above the surf)
    
    if strcmp(ELECTRODES,'on')   % plot electrodes as spots
        if isempty(EMARKER2CHANS)
            hp2 = plot3(y,x,ones(size(x))*ELECTRODE_HEIGHT,...
                EMARKER,'Color',ECOLOR,'markersize',EMARKERSIZE,'linewidth',EMARKERLINEWIDTH);
        else % plot markers for normal chans and EMARKER2CHANS separately
            hp2 = plot3(y(mark1chans),x(mark1chans),ones(size((mark1chans)))*ELECTRODE_HEIGHT,...
                EMARKER,'Color',ECOLOR,'markersize',EMARKERSIZE,'linewidth',EMARKERLINEWIDTH);
            hp2b = plot3(y(mark2chans),x(mark2chans),ones(size((mark2chans)))*ELECTRODE_HEIGHT,...
                EMARKER2,'Color',EMARKER2COLOR,'markerfacecolor',EMARKER2COLOR,'linewidth',EMARKER2LINEWIDTH,'markersize',EMARKERSIZE2);
        end
        %
        %%%%%%%%%%%%%%%%%%%%%%%% Print electrode labels only %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
    elseif strcmp(ELECTRODES,'labels')  % print electrode names (labels)
        for i = 1:size(labels,1)
            text(double(y(i)),double(x(i)),...
                ELECTRODE_HEIGHT,labels(i,:),'HorizontalAlignment','center',...
                'VerticalAlignment','middle','Color',ECOLOR,...
                'FontSize',EFSIZE)
        end
        %
        %%%%%%%%%%%%%%%%%%%%%%%% Mark electrode locations plus labels %%%%%%%%%%%%%%%%%%%
        %
    elseif strcmp(ELECTRODES,'labelpoint')
        if isempty(EMARKER2CHANS)
            hp2 = plot3(y,x,ones(size(x))*ELECTRODE_HEIGHT,...
                EMARKER,'Color',ECOLOR,'markersize',EMARKERSIZE,'linewidth',EMARKERLINEWIDTH);
        else
            hp2 = plot3(y(mark1chans),x(mark1chans),ones(size((mark1chans)))*ELECTRODE_HEIGHT,...
                EMARKER,'Color',ECOLOR,'markersize',EMARKERSIZE,'linewidth',EMARKERLINEWIDTH);
            hp2b = plot3(y(mark2chans),x(mark2chans),ones(size((mark2chans)))*ELECTRODE_HEIGHT,...
                EMARKER2,'Color',EMARKER2COLOR,'markerfacecolor',EMARKER2COLOR,'linewidth',EMARKER2LINEWIDTH,'markersize',EMARKERSIZE2);
        end
        for i = 1:size(labels,1)
            hh(i) = text(double(y(i)+0.01),double(x(i)),...
                ELECTRODE_HEIGHT,labels(i,:),'HorizontalAlignment','left',...
                'VerticalAlignment','middle','Color', ECOLOR,'userdata', num2str(allchansind(i)), ...
                'FontSize',EFSIZE, 'buttondownfcn', ...
                ['tmpstr = get(gco, ''userdata'');'...
                'set(gco, ''userdata'', get(gco, ''string''));' ...
                'set(gco, ''string'', tmpstr); clear tmpstr;'] );
        end
        %
        %%%%%%%%%%%%%%%%%%%%%%% Mark electrode locations plus numbers %%%%%%%%%%%%%%%%%%%
        %
    elseif strcmp(ELECTRODES,'numpoint')
        if isempty(EMARKER2CHANS)
            hp2 = plot3(y,x,ones(size(x))*ELECTRODE_HEIGHT,...
                EMARKER,'Color',ECOLOR,'markersize',EMARKERSIZE,'linewidth',EMARKERLINEWIDTH);
        else
            hp2 = plot3(y(mark1chans),x(mark1chans),ones(size((mark1chans)))*ELECTRODE_HEIGHT,...
                EMARKER,'Color',ECOLOR,'markersize',EMARKERSIZE,'linewidth',EMARKERLINEWIDTH);
            hp2b = plot3(y(mark2chans),x(mark2chans),ones(size((mark2chans)))*ELECTRODE_HEIGHT,...
                EMARKER2,'Color',EMARKER2COLOR,'markerfacecolor',EMARKER2COLOR,'linewidth',EMARKER2LINEWIDTH,'markersize',EMARKERSIZE2);
        end
        for i = 1:size(labels,1)
            hh(i) = text(double(y(i)+0.01),double(x(i)),...
                ELECTRODE_HEIGHT,num2str(allchansind(i)),'HorizontalAlignment','left',...
                'VerticalAlignment','middle','Color', ECOLOR,'userdata', labels(i,:) , ...
                'FontSize',EFSIZE, 'buttondownfcn', ...
                ['tmpstr = get(gco, ''userdata'');'...
                'set(gco, ''userdata'', get(gco, ''string''));' ...
                'set(gco, ''string'', tmpstr); clear tmpstr;'] );
        end
        %
        %%%%%%%%%%%%%%%%%%%%%% Print electrode numbers only %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
    elseif strcmp(ELECTRODES,'numbers')
        for i = 1:size(labels,1)
            text(double(y(i)),double(x(i)),...
                ELECTRODE_HEIGHT,int2str(allchansind(i)),'HorizontalAlignment','center',...
                'VerticalAlignment','middle','Color',ECOLOR,...
                'FontSize',EFSIZE)
        end
        %
        %%%%%%%%%%%%%%%%%%%%%% Mark emarker2 electrodes only  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
    elseif strcmp(ELECTRODES,'off') & ~isempty(EMARKER2CHANS)
        hp2b = plot3(y(mark2chans),x(mark2chans),ones(size((mark2chans)))*ELECTRODE_HEIGHT,...
            EMARKER2,'Color',EMARKER2COLOR,'markerfacecolor',EMARKER2COLOR,'linewidth',EMARKER2LINEWIDTH,'markersize',EMARKERSIZE2);
    end
    %
    %%%%%%%% Mark specified electrode locations with red filled disks  %%%%%%%%%%%%%%%%%%%%%%
    %
    try,
        if strcmpi(STYLE,'blank') % if mark-selected-channel-locations mode
            for kk = 1:length(1:length(x))
                if Values(kk) == 3
                    hp2 = plot3(y(kk),x(kk),ELECTRODE_HEIGHT,EMARKER,'Color', [0 0 0], 'markersize', EMARKERSIZE1CHAN);
                elseif Values(kk) == 2
                    hp2 = plot3(y(kk),x(kk),ELECTRODE_HEIGHT,EMARKER,'Color', [0.5 0 0], 'markersize', EMARKERSIZE1CHAN);
                elseif Values(kk) == 1
                    hp2 = plot3(y(kk),x(kk),ELECTRODE_HEIGHT,EMARKER,'Color', [1 0 0], 'markersize', EMARKERSIZE1CHAN);
                elseif strcmpi(ELECTRODES,'on')
                    hp2 = plot3(y(kk),x(kk),ELECTRODE_HEIGHT,EMARKER,'Color', ECOLOR, 'markersize', EMARKERSIZE);
                end
            end
        end
    catch, end;
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%% Plot dipole(s) on the scalp map  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    if ~isempty(DIPOLE)
        hold on;
        tmp = DIPOLE;
        if isstruct(DIPOLE)
            if ~isfield(tmp,'posxyz')
                error('dipole structure is not an EEG.dipfit.model')
            end
            DIPOLE = [];  % Note: invert x and y from dipplot usage
            DIPOLE(:,1) = -tmp.posxyz(:,2)/DIPSPHERE; % -y -> x
            DIPOLE(:,2) =  tmp.posxyz(:,1)/DIPSPHERE; %  x -> y
            DIPOLE(:,3) = -tmp.momxyz(:,2);
            DIPOLE(:,4) =  tmp.momxyz(:,1);
        else
            DIPOLE(:,1) = -tmp(:,2);                    % same for vector input
            DIPOLE(:,2) =  tmp(:,1);
            DIPOLE(:,3) = -tmp(:,4);
            DIPOLE(:,4) =  tmp(:,3);
        end;
        for index = 1:size(DIPOLE,1)
            if ~any(DIPOLE(index,:))
                DIPOLE(index,:) = [];
            end
        end;
        DIPOLE(:,1:4)   = DIPOLE(:,1:4)*rmax*(rmax/plotrad); % scale radius from 1 -> rmax (0.5)
        DIPOLE(:,3:end) = (DIPOLE(:,3:end))*rmax/100000*(rmax/plotrad);
        if strcmpi(DIPNORM, 'on')
            for index = 1:size(DIPOLE,1)
                DIPOLE(index,3:4) = DIPOLE(index,3:4)/norm(DIPOLE(index,3:end))*0.2;
            end;
        end;
        DIPOLE(:, 3:4) =  DIPORIENT*DIPOLE(:, 3:4)*DIPLEN;
        
        PLOT_DIPOLE=1;
        if sum(DIPOLE(1,3:4).^2) <= 0.00001
            if strcmpi(VERBOSE,'on')
                fprintf('Note: dipole is length 0 - not plotted\n')
            end
            PLOT_DIPOLE = 0;
        end
        if 0 % sum(DIPOLE(1,1:2).^2) > plotrad
            if strcmpi(VERBOSE,'on')
                fprintf('Note: dipole is outside plotting area - not plotted\n')
            end
            PLOT_DIPOLE = 0;
        end
        if PLOT_DIPOLE
            for index = 1:size(DIPOLE,1)
                hh = plot( DIPOLE(index, 1), DIPOLE(index, 2), '.');
                set(hh, 'color', DIPCOLOR, 'markersize', DIPSCALE*30);
                hh = line( [DIPOLE(index, 1) DIPOLE(index, 1)+DIPOLE(index, 3)]', ...
                    [DIPOLE(index, 2) DIPOLE(index, 2)+DIPOLE(index, 4)]',[10 10]);
                set(hh, 'color', DIPCOLOR, 'linewidth', DIPSCALE*30/7);
            end;
        end;
    end;
    
end % if ~ 'gridplot'

%
%%%%%%%%%%%%% Plot axis orientation %%%%%%%%%%%%%%%%%%%%
%
if strcmpi(DRAWAXIS, 'on')
    axes('position', [0 0.85 0.08 0.1]);
    axis off;
    coordend1 = sqrt(-1)*3;
    coordend2 = -3;
    coordend1 = coordend1*exp(sqrt(-1)*rotate);
    coordend2 = coordend2*exp(sqrt(-1)*rotate);
    
    line([5 5+round(real(coordend1))]', [5 5+round(imag(coordend1))]', 'color', 'k');
    line([5 5+round(real(coordend2))]', [5 5+round(imag(coordend2))]', 'color', 'k');
    if round(real(coordend2))<0
        text( 5+round(real(coordend2))*1.2, 5+round(imag(coordend2))*1.2-2, '+Y');
    else text( 5+round(real(coordend2))*1.2, 5+round(imag(coordend2))*1.2, '+Y');
    end;
    if round(real(coordend1))<0
        text( 5+round(real(coordend1))*1.2, 5+round(imag(coordend1))*1.2+1.5, '+X');
    else text( 5+round(real(coordend1))*1.2, 5+round(imag(coordend1))*1.2, '+X');
    end;
    set(gca, 'xlim', [0 10], 'ylim', [0 10]);
end;

%
%%%%%%%%%%%%% Set EEGLAB background color to match head border %%%%%%%%%%%%%%%%%%%%%%%%
%
try,
    icadefs;
    set(gcf, 'color', BACKCOLOR);
catch,
end;

hold off
axis off
return
% readlocs() - read electrode location coordinates and other information from a file.
%              Several standard file formats are supported. Users may also specify
%              a custom column format. Defined format examples are given below
%              (see File Formats).
% Usage:
%   >>  eloc = readlocs( filename );
%   >>  EEG.chanlocs = readlocs( filename, 'key', 'val', ... );
%   >>  [eloc, labels, theta, radius, indices] = ...
%                                               readlocs( filename, 'key', 'val', ... );
% Inputs:
%   filename   - Name of the file containing the electrode locations
%                {default: 2-D polar coordinates} (see >> help topoplot )
%
% Optional inputs:
%   'filetype'  - ['loc'|'sph'|'sfp'|'xyz'|'asc'|'polhemus'|'besa'|'chanedit'|'custom']
%                 Type of the file to read. By default the file type is determined
%                 using the file extension (see below under File Formats),
%                  'loc'   an EEGLAB 2-D polar coordinates channel locations file
%                          Coordinates are theta and radius (see definitions below).
%                  'sph'   Matlab spherical coordinates (Note: spherical
%                          coordinates used by Matlab functions are different
%                          from spherical coordinates used by BESA - see below).
%                  'sfp'   EGI Cartesian coordinates (NOT Matlab Cartesian - see below).
%                  'xyz'   Matlab/EEGLAB Cartesian coordinates (NOT EGI Cartesian).
%                          z is toward nose; y is toward left ear; z is toward vertex
%                  'asc'   Neuroscan polar coordinates.
%                  'polhemus' or 'polhemusx' - Polhemus electrode location file recorded
%                          with 'X' on sensor pointing to subject (see below and readelp()).
%                  'polhemusy' - Polhemus electrode location file recorded with
%                          'Y' on sensor pointing to subject (see below and readelp()).
%                  'besa' BESA-'.elp' spherical coordinates. (Not MATLAB spherical -
%                           see below).
%                  'chanedit' - EEGLAB channel location file created by pop_chanedit().
%                  'custom' - Ascii file with columns in user-defined 'format' (see below).
%   'importmode' - ['eeglab'|'native'] for location files containing 3-D cartesian electrode
%                  coordinates, import either in EEGLAB format (nose pointing toward +X).
%                  This may not always be possible since EEGLAB might not be able to
%                  determine the nose direction for scanned electrode files. 'native' import
%                  original carthesian coordinates (user can then specify the position of
%                  the nose when calling the topoplot() function; in EEGLAB the position
%                  of the nose is stored in the EEG.chaninfo structure). {default 'eeglab'}
%   'format'    -  [cell array] Format of a 'custom' channel location file (see above).
%                  {default: if no file type is defined. The cell array contains
%                  labels defining the meaning of each column of the input file.
%                           'channum'   [positive integer] channel number.
%                           'labels'    [string] channel name (no spaces).
%                           'theta'     [real degrees] 2-D angle in polar coordinates.
%                                       positive => rotating from nose (0) toward left ear
%                           'radius'    [real] radius for 2-D polar coords; 0.5 is the head
%                                       disk radius and limit for topoplot() plotting).
%                           'X'         [real] Matlab-Cartesian X coordinate (to nose).
%                           'Y'         [real] Matlab-Cartesian Y coordinate (to left ear).
%                           'Z'         [real] Matlab-Cartesian Z coordinate (to vertex).
%                           '-X','-Y','-Z' Matlab-Cartesian coordinates pointing opposite
%                                       to the above.
%                           'sph_theta' [real degrees] Matlab spherical horizontal angle.
%                                       positive => rotating from nose (0) toward left ear.
%                           'sph_phi'   [real degrees] Matlab spherical elevation angle.
%                                       positive => rotating from horizontal (0) upwards.
%                           'sph_radius' [real] distance from head center (unused).
%                           'sph_phi_besa' [real degrees] BESA phi angle from vertical.
%                                       positive => rotating from vertex (0) towards right ear.
%                           'sph_theta_besa' [real degrees] BESA theta horiz/azimuthal angle.
%                                       positive => rotating from right ear (0) toward nose.
%                           'ignore'    ignore column}.
%     The input file may also contain other channel information fields.
%                           'type'      channel type: 'EEG', 'MEG', 'EMG', 'ECG', others ...
%                           'calib'     [real near 1.0] channel calibration value.
%                           'gain'      [real > 1] channel gain.
%                           'custom1'   custom field #1.
%                           'custom2', 'custom3', 'custom4', etc.    more custom fields
%   'skiplines' - [integer] Number of header lines to skip (in 'custom' file types only).
%                 Note: Characters on a line following '%' will be treated as comments.
%   'readchans' - [integer array] indices of electrodes to read. {default: all}
%   'center'    - [(1,3) real array or 'auto'] center of xyz coordinates for conversion
%                 to spherical or polar, Specify the center of the sphere here, or 'auto'.
%                 This uses the center of the sphere that best fits all the electrode
%                 locations read. {default: [0 0 0]}
% Outputs:
%   eloc        - structure containing the channel names and locations (if present).
%                 It has three fields: 'eloc.labels', 'eloc.theta' and 'eloc.radius'
%                 identical in meaning to the EEGLAB struct 'EEG.chanlocs'.
%   labels      - cell array of strings giving the names of the electrodes. NOTE: Unlike the
%                 three outputs below, includes labels of channels *without* location info.
%   theta       - vector (in degrees) of polar angles of the electrode locations.
%   radius      - vector of polar-coordinate radii (arc_lengths) of the electrode locations
%   indices     - indices, k, of channels with non-empty 'locs(k).theta' coordinate
%
% File formats:
%   If 'filetype' is unspecified, the file extension determines its type.
%
%   '.loc' or '.locs' or '.eloc':
%               polar coordinates. Notes: angles in degrees:
%               right ear is 90; left ear -90; head disk radius is 0.5.
%               Fields:   N    angle  radius    label
%               Sample:   1    -18    .511       Fp1
%                         2     18    .511       Fp2
%                         3    -90    .256       C3
%                         4     90    .256       C4
%                           ...
%               Note: In previous releases, channel labels had to contain exactly
%               four characters (spaces replaced by '.'). This format still works,
%               though dots are no longer required.
%   '.sph':
%               Matlab spherical coordinates. Notes: theta is the azimuthal/horizontal angle
%               in deg.: 0 is toward nose, 90 rotated to left ear. Following this, performs
%               the elevation (phi). Angles in degrees.
%               Fields:   N    theta    phi    label
%               Sample:   1      18     -2      Fp1
%                         2     -18     -2      Fp2
%                         3      90     44      C3
%                         4     -90     44      C4
%                           ...
%   '.elc':
%               Cartesian 3-D electrode coordinates scanned using the EETrak software.
%               See readeetraklocs().
%   '.elp':
%               Polhemus-.'elp' Cartesian coordinates. By default, an .elp extension is read
%               as PolhemusX-elp in which 'X' on the Polhemus sensor is pointed toward the
%               subject. Polhemus files are not in columnar format (see readelp()).
%   '.elp':
%               BESA-'.elp' spherical coordinates: Need to specify 'filetype','besa'.
%               The elevation angle (phi) is measured from the vertical axis. Positive
%               rotation is toward right ear. Next, perform azimuthal/horizontal rotation
%               (theta): 0 is toward right ear; 90 is toward nose, -90 toward occiput.
%               Angles are in degrees.  If labels are absent or weights are given in
%               a last column, readlocs() adjusts for this. Default labels are E1, E2, ...
%               Fields:   Type  label      phi  theta
%               Sample:   EEG   Fp1        -92   -72
%                         EEG   Fp2         92    72
%                         EEG   C3         -46    0
%                         EEG   C4          46    0
%                           ...
%   '.xyz':
%               Matlab/EEGLAB Cartesian coordinates. Here. x is towards the nose,
%               y is towards the left ear, and z towards the vertex. Note that the first
%               column (x) is -Y in a Matlab 3-D plot, the second column (y) is X in a
%               matlab 3-D plot, and the third column (z) is Z.
%               Fields:   channum   x           y         z     label
%               Sample:   1       .950        .308     -.035     Fp1
%                         2       .950       -.308     -.035     Fp2
%                         3        0           .719      .695    C3
%                         4        0          -.719      .695    C4
%                           ...
%   '.asc', '.dat':
%               Neuroscan-.'asc' or '.dat' Cartesian polar coordinates text file.
%   '.sfp':
%               BESA/EGI-xyz Cartesian coordinates. Notes: For EGI, x is toward right ear,
%               y is toward the nose, z is toward the vertex. EEGLAB converts EGI
%               Cartesian coordinates to Matlab/EEGLAB xyz coordinates.
%               Fields:   label   x           y          z
%               Sample:   Fp1    -.308        .950      -.035
%                         Fp2     .308        .950      -.035
%                         C3     -.719        0          .695
%                         C4      .719        0          .695
%                           ...
%   '.ced':
%               ASCII file saved by pop_chanedit(). Contains multiple MATLAB/EEGLAB formats.
%               Cartesian coordinates are as in the 'xyz' format (above).
%               Fields:   channum  label  theta  radius   x      y      z    sph_theta   sph_phi  ...
%               Sample:   1        Fp1     -18    .511   .950   .308  -.035   18         -2       ...
%                         2        Fp2      18    .511   .950  -.308  -.035  -18         -2       ...
%                         3        C3      -90    .256   0      .719   .695   90         44       ...
%                         4        C4       90    .256   0     -.719   .695  -90         44       ...
%                           ...
%               The last columns of the file may contain any other defined fields (gain,
%               calib, type, custom).
%
%    Fieldtrip structure:
%               If a Fieltrip structure is given as input, an EEGLAB
%               chanlocs structure is returned
%
% Author: Arnaud Delorme, Salk Institute, 8 Dec 2002
%
% See also: readelp(), writelocs(), topo2sph(), sph2topo(), sph2cart()

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 28 Feb 2002
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


function [eloc, labels, theta, radius, indices] = readlocs( filename, varargin );

if nargin < 1
    help readlocs;
    return;
end;

% NOTE: To add a new channel format:
% ----------------------------------
% 1) Add a new element to the structure 'chanformat' (see 'ADD NEW FORMATS HERE' below):
% 2)  Enter a format 'type' for the new file format,
% 3)  Enter a (short) 'typestring' description of the format
% 4)  Enter a longer format 'description' (possibly multiline, see ex. (1) below)
% 5)  Enter format file column labels in the 'importformat' field (see ex. (2) below)
% 6)  Enter the number of header lines to skip (if any) in the 'skipline' field
% 7)  Document the new channel format in the help message above.
% 8)  After testing, please send the new version of readloca.m to us
%       at eeglab@sccn.ucsd.edu with a sample locs file.
% The 'chanformat' structure is also used (automatically) by the writelocs()
% and pop_readlocs() functions. You do not need to edit these functions.

chanformat(1).type         = 'polhemus';
chanformat(1).typestring   = 'Polhemus native .elp file';
chanformat(1).description  = [ 'Polhemus native coordinate file containing scanned electrode positions. ' ...
    'User must select the direction ' ...
    'for the nose after importing the data file.' ];
chanformat(1).importformat = 'readelp() function';
% ---------------------------------------------------------------------------------------------------
chanformat(2).type         = 'besa';
chanformat(2).typestring   = 'BESA spherical .elp file';
chanformat(2).description  = [ 'BESA spherical coordinate file. Note that BESA spherical coordinates ' ...
    'are different from Matlab spherical coordinates' ];
chanformat(2).skipline     = 0; % some BESA files do not have headers
chanformat(2).importformat = { 'type' 'labels' 'sph_theta_besa' 'sph_phi_besa' 'sph_radius' };
% ---------------------------------------------------------------------------------------------------
chanformat(3).type         = 'xyz';
chanformat(3).typestring   = 'Matlab .xyz file';
chanformat(3).description  = [ 'Standard 3-D cartesian coordinate files with electrode labels in ' ...
    'the first column and X, Y, and Z coordinates in columns 2, 3, and 4' ];
chanformat(3).importformat = { 'channum' '-Y' 'X' 'Z' 'labels'};
% ---------------------------------------------------------------------------------------------------
chanformat(4).type         = 'sfp';
chanformat(4).typestring   = 'BESA or EGI 3-D cartesian .sfp file';
chanformat(4).description  = [ 'Standard BESA 3-D cartesian coordinate files with electrode labels in ' ...
    'the first column and X, Y, and Z coordinates in columns 2, 3, and 4.' ...
    'Coordinates are re-oriented to fit the EEGLAB standard of having the ' ...
    'nose along the +X axis.' ];
chanformat(4).importformat = { 'labels' '-Y' 'X' 'Z' };
chanformat(4).skipline     = 0;
% ---------------------------------------------------------------------------------------------------
chanformat(5).type         = 'loc';
chanformat(5).typestring   = 'EEGLAB polar .loc file';
chanformat(5).description  = [ 'EEGLAB polar .loc file' ];
chanformat(5).importformat = { 'channum' 'theta' 'radius' 'labels' };
% ---------------------------------------------------------------------------------------------------
chanformat(6).type         = 'sph';
chanformat(6).typestring   = 'Matlab .sph spherical file';
chanformat(6).description  = [ 'Standard 3-D spherical coordinate files in Matlab format' ];
chanformat(6).importformat = { 'channum' 'sph_theta' 'sph_phi' 'labels' };
% ---------------------------------------------------------------------------------------------------
chanformat(7).type         = 'asc';
chanformat(7).typestring   = 'Neuroscan polar .asc file';
chanformat(7).description  = [ 'Neuroscan polar .asc file, automatically recentered to fit EEGLAB standard' ...
    'of having ''Cz'' at (0,0).' ];
chanformat(7).importformat = 'readneurolocs';
% ---------------------------------------------------------------------------------------------------
chanformat(8).type         = 'dat';
chanformat(8).typestring   = 'Neuroscan 3-D .dat file';
chanformat(8).description  = [ 'Neuroscan 3-D cartesian .dat file. Coordinates are re-oriented to fit ' ...
    'the EEGLAB standard of having the nose along the +X axis.' ];
chanformat(8).importformat = 'readneurolocs';
% ---------------------------------------------------------------------------------------------------
chanformat(9).type         = 'elc';
chanformat(9).typestring   = 'ASA .elc 3-D file';
chanformat(9).description  = [ 'ASA .elc 3-D coordinate file containing scanned electrode positions. ' ...
    'User must select the direction ' ...
    'for the nose after importing the data file.' ];
chanformat(9).importformat = 'readeetraklocs';
% ---------------------------------------------------------------------------------------------------
chanformat(10).type         = 'chanedit';
chanformat(10).typestring   = 'EEGLAB complete 3-D file';
chanformat(10).description  = [ 'EEGLAB file containing polar, cartesian 3-D, and spherical 3-D ' ...
    'electrode locations.' ];
chanformat(10).importformat = { 'channum' 'labels'  'theta' 'radius' 'X' 'Y' 'Z' 'sph_theta' 'sph_phi' ...
    'sph_radius' 'type' };
chanformat(10).skipline     = 1;
% ---------------------------------------------------------------------------------------------------
chanformat(11).type         = 'custom';
chanformat(11).typestring   = 'Custom file format';
chanformat(11).description  = 'Custom ASCII file format where user can define content for each file columns.';
chanformat(11).importformat = '';
% ---------------------------------------------------------------------------------------------------
% ----- ADD MORE FORMATS HERE -----------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------

listcolformat = { 'labels' 'channum' 'theta' 'radius' 'sph_theta' 'sph_phi' ...
    'sph_radius' 'sph_theta_besa' 'sph_phi_besa' 'gain' 'calib' 'type' ...
    'X' 'Y' 'Z' '-X' '-Y' '-Z' 'custom1' 'custom2' 'custom3' 'custom4' 'ignore' 'not def' };

% ----------------------------------
% special mode for getting the info
% ----------------------------------
if isstr(filename) & strcmp(filename, 'getinfos')
    eloc = chanformat;
    labels = listcolformat;
    return;
end;

g = finputcheck( varargin, ...
    { 'filetype'	   'string'  {}                 '';
    'importmode'  'string'  { 'eeglab','native' } 'eeglab';
    'defaultelp'  'string'  { 'besa','polhemus' } 'polhemus';
    'skiplines'   'integer' [0 Inf] 			[];
    'elecind'     'integer' [1 Inf]	    	[];
    'format'	   'cell'	 []					{} }, 'readlocs');
if isstr(g), error(g); end;
if ~isempty(g.format), g.filetype = 'custom'; end;

if isstr(filename)
    
    % format auto detection
    % --------------------
    if strcmpi(g.filetype, 'autodetect'), g.filetype = ''; end;
    g.filetype = strtok(g.filetype);
    periods = find(filename == '.');
    fileextension = filename(periods(end)+1:end);
    g.filetype = lower(g.filetype);
    if isempty(g.filetype)
        switch lower(fileextension),
            case {'loc' 'locs' 'eloc'}, g.filetype = 'loc'; % 5/27/2014 Ramon: 'eloc' option introduced.
            case 'xyz', g.filetype = 'xyz';
                fprintf( [ 'WARNING: Matlab Cartesian coord. file extension (".xyz") detected.\n' ...
                    'If importing EGI Cartesian coords, force type "sfp" instead.\n'] );
            case 'sph', g.filetype = 'sph';
            case 'ced', g.filetype = 'chanedit';
            case 'elp', g.filetype = g.defaultelp;
            case 'asc', g.filetype = 'asc';
            case 'dat', g.filetype = 'dat';
            case 'elc', g.filetype = 'elc';
            case 'eps', g.filetype = 'besa';
            case 'sfp', g.filetype = 'sfp';
            otherwise, g.filetype =  '';
        end;
        fprintf('readlocs(): ''%s'' format assumed from file extension\n', g.filetype);
    else
        if strcmpi(g.filetype, 'locs'),  g.filetype = 'loc'; end
        if strcmpi(g.filetype, 'eloc'),  g.filetype = 'loc'; end
    end;
    
    % assign format from filetype
    % ---------------------------
    if ~isempty(g.filetype) & ~strcmpi(g.filetype, 'custom') ...
            & ~strcmpi(g.filetype, 'asc') & ~strcmpi(g.filetype, 'elc') & ~strcmpi(g.filetype, 'dat')
        indexformat = strmatch(lower(g.filetype), { chanformat.type }, 'exact');
        g.format = chanformat(indexformat).importformat;
        if isempty(g.skiplines)
            g.skiplines = chanformat(indexformat).skipline;
        end;
        if isempty(g.filetype)
            error( ['readlocs() error: The filetype cannot be detected from the \n' ...
                '                  file extension, and custom format not specified']);
        end;
    end;
    
    % import file
    % -----------
    if strcmp(g.filetype, 'asc') | strcmp(g.filetype, 'dat')
        eloc = readneurolocs( filename );
        eloc = rmfield(eloc, 'sph_theta'); % for the conversion below
        eloc = rmfield(eloc, 'sph_theta_besa'); % for the conversion below
        if isfield(eloc, 'type')
            for index = 1:length(eloc)
                type = eloc(index).type;
                if type == 69,     eloc(index).type = 'EEG';
                elseif type == 88, eloc(index).type = 'REF';
                elseif type >= 76 & type <= 82, eloc(index).type = 'FID';
                else eloc(index).type = num2str(eloc(index).type);
                end;
            end;
        end;
    elseif strcmp(g.filetype, 'elc')
        eloc = readeetraklocs( filename );
        %eloc = read_asa_elc( filename ); % from fieldtrip
        %eloc = struct('labels', eloc.label, 'X', mattocell(eloc.pnt(:,1)'), 'Y', ...
        %                        mattocell(eloc.pnt(:,2)'), 'Z', mattocell(eloc.pnt(:,3)'));
        eloc = convertlocs(eloc, 'cart2all');
        eloc = rmfield(eloc, 'sph_theta'); % for the conversion below
        eloc = rmfield(eloc, 'sph_theta_besa'); % for the conversion below
    elseif strcmp(lower(g.filetype(1:end-1)), 'polhemus') | ...
            strcmp(g.filetype, 'polhemus')
        try,
            [eloc labels X Y Z]= readelp( filename );
            if strcmp(g.filetype, 'polhemusy')
                tmp = X; X = Y; Y = tmp;
            end;
            for index = 1:length( eloc )
                eloc(index).X = X(index);
                eloc(index).Y = Y(index);
                eloc(index).Z = Z(index);
            end;
        catch,
            disp('readlocs(): Could not read Polhemus coords. Trying to read BESA .elp file.');
            [eloc, labels, theta, radius, indices] = readlocs( filename, 'defaultelp', 'besa', varargin{:} );
        end;
    else
        % importing file
        % --------------
        if isempty(g.skiplines), g.skiplines = 0; end;
        if strcmpi(g.filetype, 'chanedit')
            array = loadtxt( filename, 'delim', 9, 'skipline', g.skiplines, 'blankcell', 'off');
        else
            array = load_file_or_array( filename, g.skiplines);
        end;
        if size(array,2) < length(g.format)
            fprintf(['readlocs() warning: Fewer columns in the input than expected.\n' ...
                '                    See >> help readlocs\n']);
        elseif size(array,2) > length(g.format)
            fprintf(['readlocs() warning: More columns in the input than expected.\n' ...
                '                    See >> help readlocs\n']);
        end;
        
        % removing lines BESA
        % -------------------
        if isempty(array{1,2})
            disp('BESA header detected, skipping three lines...');
            array = load_file_or_array( filename, g.skiplines-1);
            if isempty(array{1,2})
                array = load_file_or_array( filename, g.skiplines-1);
            end;
        end;
        
        % xyz format, is the first col absent
        % -----------------------------------
        if strcmp(g.filetype, 'xyz')
            if size(array, 2) == 4
                array(:, 2:5) = array(:, 1:4);
            end;
        end;
        
        % removing comments and empty lines
        % ---------------------------------
        indexbeg = 1;
        while isempty(array{indexbeg,1}) | ...
                (isstr(array{indexbeg,1}) & array{indexbeg,1}(1) == '%' )
            indexbeg = indexbeg+1;
        end;
        array = array(indexbeg:end,:);
        
        % converting file
        % ---------------
        for indexcol = 1:min(size(array,2), length(g.format))
            [str mult] = checkformat(g.format{indexcol});
            for indexrow = 1:size( array, 1)
                if mult ~= 1
                    eval ( [ 'eloc(indexrow).'  str '= -array{indexrow, indexcol};' ]);
                else
                    eval ( [ 'eloc(indexrow).'  str '= array{indexrow, indexcol};' ]);
                end;
            end;
        end;
    end;
    
    % handling BESA coordinates
    % -------------------------
    if isfield(eloc, 'sph_theta_besa')
        if isfield(eloc, 'type')
            if isnumeric(eloc(1).type)
                disp('BESA format detected ( Theta | Phi )');
                for index = 1:length(eloc)
                    eloc(index).sph_phi_besa   = eloc(index).labels;
                    eloc(index).sph_theta_besa = eloc(index).type;
                    eloc(index).labels         = '';
                    eloc(index).type           = '';
                end;
                eloc = rmfield(eloc, 'labels');
            end;
        end;
        if isfield(eloc, 'labels')
            if isnumeric(eloc(1).labels)
                disp('BESA format detected ( Elec | Theta | Phi )');
                for index = 1:length(eloc)
                    eloc(index).sph_phi_besa   = eloc(index).sph_theta_besa;
                    eloc(index).sph_theta_besa = eloc(index).labels;
                    eloc(index).labels         = eloc(index).type;
                    eloc(index).type           = '';
                    eloc(index).radius         = 1;
                end;
            end;
        end;
        
        try
            eloc = convertlocs(eloc, 'sphbesa2all');
            eloc = convertlocs(eloc, 'topo2all'); % problem with some EGI files (not BESA files)
        catch, disp('Warning: coordinate conversion failed'); end;
        fprintf('Readlocs: BESA spherical coords. converted, now deleting BESA fields\n');
        fprintf('          to avoid confusion (these fields can be exported, though)\n');
        eloc = rmfield(eloc, 'sph_phi_besa');
        eloc = rmfield(eloc, 'sph_theta_besa');
        
        % converting XYZ coordinates to polar
        % -----------------------------------
    elseif isfield(eloc, 'sph_theta')
        try
            eloc = convertlocs(eloc, 'sph2all');
        catch, disp('Warning: coordinate conversion failed'); end;
    elseif isfield(eloc, 'X')
        try
            eloc = convertlocs(eloc, 'cart2all');
        catch, disp('Warning: coordinate conversion failed'); end;
    else
        try
            eloc = convertlocs(eloc, 'topo2all');
        catch, disp('Warning: coordinate conversion failed'); end;
    end;
    
    % inserting labels if no labels
    % -----------------------------
    if ~isfield(eloc, 'labels')
        fprintf('readlocs(): Inserting electrode labels automatically.\n');
        for index = 1:length(eloc)
            eloc(index).labels = [ 'E' int2str(index) ];
        end;
    else
        % remove trailing '.'
        for index = 1:length(eloc)
            if isstr(eloc(index).labels)
                tmpdots = find( eloc(index).labels == '.' );
                eloc(index).labels(tmpdots) = [];
            end;
        end;
    end;
    
    % resorting electrodes if number not-sorted
    % -----------------------------------------
    if isfield(eloc, 'channum')
        if ~isnumeric(eloc(1).channum)
            error('Channel numbers must be numeric');
        end;
        allchannum = [ eloc.channum ];
        if any( sort(allchannum) ~= allchannum )
            fprintf('readlocs(): Re-sorting channel numbers based on ''channum'' column indices\n');
            [tmp newindices] = sort(allchannum);
            eloc = eloc(newindices);
        end;
        eloc = rmfield(eloc, 'channum');
    end;
else
    if isstruct(filename)
        % detect Fieldtrip structure and convert it
        % -----------------------------------------
        if isfield(filename, 'pnt')
            neweloc = [];
            for index = 1:length(filename.label)
                neweloc(index).labels = filename.label{index};
                neweloc(index).X      = filename.pnt(index,1);
                neweloc(index).Y      = filename.pnt(index,2);
                neweloc(index).Z      = filename.pnt(index,3);
            end;
            eloc = neweloc;
            eloc = convertlocs(eloc, 'cart2all');
        else
            eloc = filename;
        end;
    else
        disp('readlocs(): input variable must be a string or a structure');
    end;
end;
if ~isempty(g.elecind)
    eloc = eloc(g.elecind);
end;
if nargout > 2
    if isfield(eloc, 'theta')
        tmptheta = { eloc.theta }; % check which channels have (polar) coordinates set
    else tmptheta = cell(1,length(eloc));
    end;
    if isfield(eloc, 'theta')
        tmpx = { eloc.X }; % check which channels have (polar) coordinates set
    else tmpx = cell(1,length(eloc));
    end;
    
    indices           = find(~cellfun('isempty', tmptheta));
    indices           = intersect_bc(find(~cellfun('isempty', tmpx)), indices);
    indices           = sort(indices);
    
    indbad            = setdiff_bc(1:length(eloc), indices);
    tmptheta(indbad)  = { NaN };
    theta             = [ tmptheta{:} ];
end;
if nargout > 3
    if isfield(eloc, 'theta')
        tmprad = { eloc.radius }; % check which channels have (polar) coordinates set
    else tmprad = cell(1,length(eloc));
    end;
    tmprad(indbad)    = { NaN };
    radius            = [ tmprad{:} ];
end;

%tmpnum = find(~cellfun('isclass', { eloc.labels }, 'char'));
%disp('Converting channel labels to string');
for index = 1:length(eloc)
    if ~isstr(eloc(index).labels)
        eloc(index).labels = int2str(eloc(index).labels);
    end;
end;
labels = { eloc.labels };
if isfield(eloc, 'ignore')
    eloc = rmfield(eloc, 'ignore');
end;

% process fiducials if any
% ------------------------
fidnames = { 'nz' 'lpa' 'rpa' 'nasion' 'left' 'right' 'nazion' 'fidnz' 'fidt9' 'fidt10' 'cms' 'drl' };
for index = 1:length(fidnames)
    ind = strmatch(fidnames{index}, lower(labels), 'exact');
    if ~isempty(ind), eloc(ind).type = 'FID'; end;
end;

return;

% interpret the variable name
% ---------------------------
function array = load_file_or_array( varname, skiplines );
if isempty(skiplines),
    skiplines = 0;
end;
if exist( varname ) == 2
    array = loadtxt(varname,'verbose','off','skipline',skiplines,'blankcell','off');
else % variable in the global workspace
    % --------------------------
    try, array = evalin('base', varname);
    catch, error('readlocs(): cannot find the named file or variable, check syntax');
    end;
end;
return;

% check field format
% ------------------
function [str, mult] = checkformat(str)
mult = 1;
if strcmpi(str, 'labels'),         str = lower(str); return; end;
if strcmpi(str, 'channum'),        str = lower(str); return; end;
if strcmpi(str, 'theta'),          str = lower(str); return; end;
if strcmpi(str, 'radius'),         str = lower(str); return; end;
if strcmpi(str, 'ignore'),         str = lower(str); return; end;
if strcmpi(str, 'sph_theta'),      str = lower(str); return; end;
if strcmpi(str, 'sph_phi'),        str = lower(str); return; end;
if strcmpi(str, 'sph_radius'),     str = lower(str); return; end;
if strcmpi(str, 'sph_theta_besa'), str = lower(str); return; end;
if strcmpi(str, 'sph_phi_besa'),   str = lower(str); return; end;
if strcmpi(str, 'gain'),           str = lower(str); return; end;
if strcmpi(str, 'calib'),          str = lower(str); return; end;
if strcmpi(str, 'type') ,          str = lower(str); return; end;
if strcmpi(str, 'X'),              str = upper(str); return; end;
if strcmpi(str, 'Y'),              str = upper(str); return; end;
if strcmpi(str, 'Z'),              str = upper(str); return; end;
if strcmpi(str, '-X'),             str = upper(str(2:end)); mult = -1; return; end;
if strcmpi(str, '-Y'),             str = upper(str(2:end)); mult = -1; return; end;
if strcmpi(str, '-Z'),             str = upper(str(2:end)); mult = -1; return; end;
if strcmpi(str, 'custom1'), return; end;
if strcmpi(str, 'custom2'), return; end;
if strcmpi(str, 'custom3'), return; end;
if strcmpi(str, 'custom4'), return; end;
error(['readlocs(): undefined field ''' str '''']);

% finputcheck() - check Matlab function {'key','value'} input argument pairs
%
% Usage: >> result = finputcheck( varargin, fieldlist );
%        >> [result varargin] = finputcheck( varargin, fieldlist, ...
%                                              callingfunc, mode, verbose );
% Input:
%   varargin  - Cell array 'varargin' argument from a function call using 'key',
%               'value' argument pairs. See Matlab function 'varargin'.
%               May also be a structure such as struct(varargin{:})
%   fieldlist - A 4-column cell array, one row per 'key'. The first
%               column contains the key string, the second its type(s),
%               the third the accepted value range, and the fourth the
%               default value.  Allowed types are 'boolean', 'integer',
%               'real', 'string', 'cell' or 'struct'.  For example,
%                       {'key1' 'string' { 'string1' 'string2' } 'defaultval_key1'}
%                       {'key2' {'real' 'integer'} { minint maxint } 'defaultval_key2'}
%  callingfunc - Calling function name for error messages. {default: none}.
%  mode        - ['ignore'|'error'] ignore keywords that are either not specified
%                in the fieldlist cell array or generate an error.
%                {default: 'error'}.
%  verbose     - ['verbose', 'quiet'] print information. Default: 'verbose'.
%
% Outputs:
%   result     - If no error, structure with 'key' as fields and 'value' as
%                content. If error this output contain the string error.
%   varargin   - residual varagin containing unrecognized input arguments.
%                Requires mode 'ignore' above.
%
% Note: In case of error, a string is returned containing the error message
%       instead of a structure.
%
% Example (insert the following at the beginning of your function):
%	result = finputcheck(varargin, ...
%               { 'title'         'string'   []       ''; ...
%                 'percent'       'real'     [0 1]    1 ; ...
%                 'elecamp'       'integer'  [1:10]   [] });
%   if isstr(result)
%       error(result);
%   end
%
% Note:
%   The 'title' argument should be a string. {no default value}
%   The 'percent' argument should be a real number between 0 and 1. {default: 1}
%   The 'elecamp' argument should be an integer between 1 and 10 (inclusive).
%
%   Now 'g.title' will contain the title arg (if any, else the default ''), etc.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 10 July 2002

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 10 July 2002, arno@salk.edu
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

function [g, varargnew] = finputcheck( vararg, fieldlist, callfunc, mode, verbose )

if nargin < 2
    help finputcheck;
    return;
end;
if nargin < 3
    callfunc = '';
else
    callfunc = [callfunc ' ' ];
end;
if nargin < 4
    mode = 'do not ignore';
end;
if nargin < 5
    verbose = 'verbose';
end;
NAME = 1;
TYPE = 2;
VALS = 3;
DEF  = 4;
SIZE = 5;

varargnew = {};
% create structure
% ----------------
if ~isempty(vararg)
    if isstruct(vararg)
        g = vararg;
    else
        for index=1:length(vararg)
            if iscell(vararg{index})
                vararg{index} = {vararg{index}};
            end;
        end;
        try
            g = struct(vararg{:});
        catch
            vararg = removedup(vararg, verbose);
            try
                g = struct(vararg{:});
            catch
                g = [ callfunc 'error: bad ''key'', ''val'' sequence' ]; return;
            end;
        end;
    end;
else
    g = [];
end;

for index = 1:size(fieldlist,NAME)
    % check if present
    % ----------------
    if ~isfield(g, fieldlist{index, NAME})
        g = setfield( g, fieldlist{index, NAME}, fieldlist{index, DEF});
    end;
    tmpval = getfield( g, {1}, fieldlist{index, NAME});
    
    % check type
    % ----------
    if ~iscell( fieldlist{index, TYPE} )
        res = fieldtest( fieldlist{index, NAME},  fieldlist{index, TYPE}, ...
            fieldlist{index, VALS}, tmpval, callfunc );
        if isstr(res), g = res; return; end;
    else
        testres = 0;
        tmplist = fieldlist;
        for it = 1:length( fieldlist{index, TYPE} )
            if ~iscell(fieldlist{index, VALS})
                res{it} = fieldtest(  fieldlist{index, NAME},  fieldlist{index, TYPE}{it}, ...
                    fieldlist{index, VALS}, tmpval, callfunc );
            else res{it} = fieldtest(  fieldlist{index, NAME},  fieldlist{index, TYPE}{it}, ...
                    fieldlist{index, VALS}{it}, tmpval, callfunc );
            end;
            if ~isstr(res{it}), testres = 1; end;
        end;
        if testres == 0,
            g = res{1};
            for tmpi = 2:length(res)
                g = [ g 10 'or ' res{tmpi} ];
            end;
            return;
        end;
    end;
end;

% check if fields are defined
% ---------------------------
allfields = fieldnames(g);
for index=1:length(allfields)
    if isempty(strmatch(allfields{index}, fieldlist(:, 1)', 'exact'))
        if ~strcmpi(mode, 'ignore')
            g = [ callfunc 'error: undefined argument ''' allfields{index} '''']; return;
        end;
        varargnew{end+1} = allfields{index};
        varargnew{end+1} = getfield(g, {1}, allfields{index});
    end;
end;


function g = fieldtest( fieldname, fieldtype, fieldval, tmpval, callfunc );
NAME = 1;
TYPE = 2;
VALS = 3;
DEF  = 4;
SIZE = 5;
g = [];

switch fieldtype
    case { 'integer' 'real' 'boolean' 'float' },
        if ~isnumeric(tmpval) && ~islogical(tmpval)
            g = [ callfunc 'error: argument ''' fieldname ''' must be numeric' ]; return;
        end;
        if strcmpi(fieldtype, 'boolean')
            if tmpval ~=0 && tmpval ~= 1
                g = [ callfunc 'error: argument ''' fieldname ''' must be 0 or 1' ]; return;
            end;
        else
            if strcmpi(fieldtype, 'integer')
                if ~isempty(fieldval)
                    if (any(isnan(tmpval(:))) && ~any(isnan(fieldval))) ...
                            && (~ismember(tmpval, fieldval))
                        g = [ callfunc 'error: wrong value for argument ''' fieldname '''' ]; return;
                    end;
                end;
            else % real or float
                if ~isempty(fieldval) && ~isempty(tmpval)
                    if any(tmpval < fieldval(1)) || any(tmpval > fieldval(2))
                        g = [ callfunc 'error: value out of range for argument ''' fieldname '''' ]; return;
                    end;
                end;
            end;
        end;
        
        
    case 'string'
        if ~isstr(tmpval)
            g = [ callfunc 'error: argument ''' fieldname ''' must be a string' ]; return;
        end;
        if ~isempty(fieldval)
            if isempty(strmatch(lower(tmpval), lower(fieldval), 'exact'))
                g = [ callfunc 'error: wrong value for argument ''' fieldname '''' ]; return;
            end;
        end;
        
        
    case 'cell'
        if ~iscell(tmpval)
            g = [ callfunc 'error: argument ''' fieldname ''' must be a cell array' ]; return;
        end;
        
        
    case 'struct'
        if ~isstruct(tmpval)
            g = [ callfunc 'error: argument ''' fieldname ''' must be a structure' ]; return;
        end;
        
    case 'function_handle'
        if ~isa(tmpval, 'function_handle')
            g = [ callfunc 'error: argument ''' fieldname ''' must be a function handle' ]; return;
        end;
        
    case '';
    otherwise, error([ 'finputcheck error: unrecognized type ''' fieldname '''' ]);
end;

% remove duplicates in the list of parameters
% -------------------------------------------
function cella = removedup(cella, verbose)
% make sure if all the values passed to unique() are strings, if not, exist
%try
[tmp indices] = unique_bc(cella(1:2:end));
if length(tmp) ~= length(cella)/2
    myfprintf(verbose,'Note: duplicate ''key'', ''val'' parameter(s), keeping the last one(s)\n');
end;
cella = cella(sort(union(indices*2-1, indices*2)));
%catch
% some elements of cella were not string
%    error('some ''key'' values are not string.');
%end;

function myfprintf(verbose, varargin)

if strcmpi(verbose, 'verbose')
    fprintf(varargin{:});
end;
% intersect_bc - intersect backward compatible with Matlab versions prior to 2013a

function [C,IA,IB] = intersect_bc(A,B,varargin);

errorFlag = false;

v = version;
indp = find(v == '.');
v = str2num(v(1:indp(2)-1));
if v > 7.19, v = floor(v) + rem(v,1)/10; end;

if nargin > 2
    ind = strmatch('legacy', varargin);
    if ~isempty(ind)
        varargin(ind) = [];
    end;
end;

if v >= 7.14
    [C,IA,IB] = intersect(A,B,varargin{:},'legacy');
    if errorFlag
        [C2,IA2,IB2] = intersect(A,B,varargin{:});
        if (~isequal(C, C2) || ~isequal(IA, IA2) || ~isequal(IB, IB2))
            warning('backward compatibility issue with call to intersect function');
        end;
    end;
else
    [C,IA,IB] = intersect(A,B,varargin{:});
end;
% setdiff_bc - setdiff backward compatible with Matlab versions prior to 2013a

function [C,IA] = setdiff_bc(A,B,varargin);

errorFlag = false;

v = version;
indp = find(v == '.');
v = str2num(v(1:indp(2)-1));
if v > 7.19, v = floor(v) + rem(v,1)/10; end;

if nargin > 2
    ind = strmatch('legacy', varargin);
    if ~isempty(ind)
        varargin(ind) = [];
    end;
end;

if v >= 7.14
    [C,IA] = setdiff(A,B,varargin{:},'legacy');
    if errorFlag
        [C2,IA2] = setdiff(A,B,varargin{:});
        if (~isequal(C, C2) || ~isequal(IA, IA2))
            warning('backward compatibility issue with call to setdiff function');
        end;
    end;
else
    [C,IA] = setdiff(A,B,varargin{:});
end;
% fastif() - fast if function.
%
% Usage:
%  >> res = fastif(test, s1, s2);
%
% Input:
%   test   - logical test with result 0 or 1
%   s1     - result if 1
%   s2     - result if 0
%
% Output:
%   res    - s1 or s2 depending on the value of the test
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001

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

function res = fastif(s1, s2, s3);

if s1
    res = s2;
else
    res = s3;
end;
return;
% vararg2str() - transform arguments into string for evaluation
%                using the eval() command
%
% Usage:
%   >> strout = vararg2str( allargs );
%   >> strout = vararg2str( allargs, inputnames, inputnum, nostrconv );
%
% Inputs:
%   allargs    - Cell array containing all arguments
%   inputnames - Cell array of input names for these arguments, if any.
%   inputnum   - Vector of indices for all inputs. If present, the
%                string output may by replaced by varargin{num}.
%                Include NaN in the vector to avoid specific parameters
%                being converted in this way.
%   nostrconv  - Vector of 0s and 1s indicating where the string
%                should be not be converted.
%
% Outputs:
%   strout     - output string
%
% Author: Arnaud Delorme, CNL / Salk Institute, 9 April 2002

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 9 April 2002
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

function strout = vararg2str(allargs, inputnam, inputnum, int2str );

if nargin < 1
    help vararg2str;
    return;
end;
if isempty(allargs)
    strout = '';
    return;
end;

% default arguments
% -----------------
if nargin < 2
    inputnam(1:length(allargs)) = {''};
else
    if length(inputnam) < length(allargs)
        inputnam(end+1:length(allargs)) = {''};
    end;
end;
if nargin < 3
    inputnum(1:length(allargs)) = NaN;
else
    if length(inputnum) < length(allargs)
        inputnum(end+1:length(allargs)) = NaN;
    end;
end;
if nargin < 4
    int2str(1:length(allargs)) = 0;
else
    if length(int2str) < length(allargs)
        int2str(end+1:length(allargs)) = 0;
    end;
end;
if ~iscell( allargs )
    allargs = { allargs };
end;

% actual conversion
% -----------------
strout = '';
for index = 1:length(allargs)
    tmpvar = allargs{index};
    if ~isempty(inputnam{index})
        strout = [ strout ',' inputnam{index} ];
    else
        if isstr( tmpvar )
            if int2str(index)
                strout = [ strout ',' tmpvar ];
            else
                strout = [ strout ',' str2str( tmpvar ) ];
            end;
        elseif isnumeric( tmpvar ) | islogical( tmpvar )
            strout = [ strout ',' array2str( tmpvar ) ];
        elseif iscell( tmpvar )
            tmpres = vararg2str( tmpvar );
            comas  = find( tmpres == ',' );
            tmpres(comas) = ' ';
            strout = [ strout ',{' tmpres '}' ];
        elseif isstruct(tmpvar)
            strout = [ strout ',' struct2str( tmpvar ) ];
        else
            error('Unrecognized input');
        end;
    end;
    
end;
if ~isempty(strout)
    strout = strout(2:end);
end;

% convert string to string
% ------------------------
function str = str2str( array )
if isempty( array), str = ''''''; return; end;
str = '';
for index = 1:size(array,1)
    tmparray = deblank(array(index,:));
    if isempty(tmparray)
        str = [ str ','' ''' ];
    else
        str = [ str ',''' doublequotes(tmparray) '''' ];
    end;
end;
if size(array,1) > 1
    str = [ 'strvcat(' str(2:end) ')'];
else
    str = str(2:end);
end;
return;

% convert array to string
% -----------------------
function str = array2str( array )
if isempty( array), str = '[]'; return; end;
if prod(size(array)) == 1, str = num2str(array); return; end;
if size(array,1) == 1, str = [ '[' contarray(array) '] ' ]; return; end;
if size(array,2) == 1, str = [ '[' contarray(array') ']'' ' ]; return; end;
str = '';
for index = 1:size(array,1)
    str = [ str ';' contarray(array(index,:)) ];
end;
str = [ '[' str(2:end) ']' ];
return;

% convert struct to string
% ------------------------
function str = struct2str( structure )
if isempty( structure )
    str = 'struct([])'; return;
end;
str = '';
allfields = fieldnames( structure );
for index = 1:length( allfields )
    strtmp = '';
    eval( [ 'allcontent = { structure.' allfields{index} ' };' ] ); % getfield generates a bug
    str = [ str, '''' allfields{index} ''',{' vararg2str( allcontent ) '},' ];
end;
str = [ 'struct(' str(1:end-1) ')' ];
return;

% double the quotes in strings
% ----------------------------
function str = doublequotes( str )
quoteloc = union_bc(findstr( str, ''''), union(findstr(str, '%'), findstr(str, '\')));
if ~isempty(quoteloc)
    for index = length(quoteloc):-1:1
        str = [ str(1:quoteloc(index)) str(quoteloc(index):end) ];
    end;
end;
return;

% test continuous arrays
% ----------------------
function str = contarray( array )
array = double(array);
tmpind = find( round(array) ~= array );
if prod(size(array)) == 1
    str =  num2str(array);
    return;
end;
if size(array,1) == 1 & size(array,2) == 2
    str =  [num2str(array(1)) ' ' num2str(array(2))];
    return;
end;
if isempty(tmpind) | all(isnan(array(tmpind)))
    str = num2str(array(1));
    skip = 0;
    indent = array(2) - array(1);
    for index = 2:length(array)
        if array(index) ~= array(index-1)+indent | indent == 0
            if skip <= 1
                if skip == 0
                    str = [str ' ' num2str(array(index))];
                else
                    str = [str ' ' num2str(array(index-1)) ' ' num2str(array(index))];
                end;
            else
                if indent == 1
                    str = [str ':' num2str(array(index-1)) ' ' num2str(array(index))];
                else
                    str = [str ':' num2str(indent) ':' num2str(array(index-1)) ' ' num2str(array(index))];
                end;
            end;
            skip = 0;
            indent = array(index) - array(index-1);
        else
            skip = skip + 1;
        end;
    end;
    if array(index) == array(index-1)+indent
        if skip ~= 0
            if indent == 1
                str = [str ':' num2str(array(index)) ];
            elseif indent == 0
                str = [str ' ' num2str(array(index)) ];
            else
                str = [str ':' num2str(indent) ':' num2str(array(index)) ];
            end;
        end;
    end;
else
    if length(array) < 10
        str = num2str(array(1));
        for index = 2:length(array)
            str = [str ' ' num2str(array(index)) ];
        end;
    else
        str = num2str(double(array));
    end;
end;

