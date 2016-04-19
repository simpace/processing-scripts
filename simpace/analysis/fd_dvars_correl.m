
clear all;

pipes = {'NR_0';'NR_1';'NR_2';'NR_3';'NR_4'};
pipe_nms = {'BF';'Mn';'Mn+GS';'C';'C+GS'};
sub = 'sub01';
analysis_dir = '/home/despo/simpace/analyses/';

% load([analysis_dir 'pipeline_list_hpf.mat'])
% pipes = pipes(2:end);
% pipe_nms = pipe_nms(2:end);
% 
% pipes(end:end+2) = {'localWMreg+GS';'localWMreg_smooth';'localWMreg+GS_smooth'};
% pipe_nms(end:end+2) = {'LWM+GS';'LWM(smo)';'LWM(smo)+GS'};

Npipes = length(pipes);
workflow_idx = ones(Npipes,1);
workflow_nm = {'preproc'};

% pipes(end:end+2) = {'localWMreg_smooth';'HPF';'NR_4'};
% pipe_nms(end:end+2) = {'LWM(smo)';'HP-noMC';'TF-noMC'};
% workflow_idx = cat(1,ones(8,1),2*ones(2,1));
% workflow_oth = 'preproc_noMC';
% Npipes = 8;%length(pipes);

load([analysis_dir 'mot_levels.mat'])

data_dir = [analysis_dir 'fd-dvars/' sub '/'];
save_dir = [analysis_dir 'fd-dvars-figs/' sub '/'];

data_types = {'FD';'DVARS_LFF'};
Ndata_types = length(data_types);

Nsessions = 13;
Ntrs = 195-1;

Meas = NaN(Nsessions,Npipes,Nruns,2);
    
data = NaN(Ntrs,Nruns,Ndata_types,Npipes,Nsessions);

reliability = NaN(Npipes,Npipes,Nruns,Nsessions,Ndata_types);

for isess = 1:Nsessions

    if isess < 10; sess_str = ['0' num2str(isess)]; else sess_str = num2str(isess); end
    
    for ipipe = 1:Npipes
        
        pipe = pipes{ipipe};

        for itype = 1:Ndata_types
            
            type = data_types{itype};

            file_nm = [sess_str '_' workflow_nm{workflow_idx(ipipe)} '_' pipe '_' type '.txt'];

            file = deal_names_full(file_nm,data_dir);
            
            if itype==1 && ipipe==1; % read in FD file just once for each session
                fd_curr = textread(file);
                temp = fd_curr;
                
            elseif itype==1; % FD file
                temp = fd_curr;
            
            elseif itype==2; % read in DVARS file
                temp = textread(file);
            end
            
            data(:,:,itype,ipipe,isess) = temp;
            clear temp;
             
            a = 1;

        end

        for irun = 1:Nruns
            curr = squeeze(data(:,irun,:,ipipe,isess));
            
            % correlation
            Meas(isess,ipipe,irun,1) = rtoz(corr(curr(:,1),curr(:,2)));
            
            % dot product
            Meas(isess,ipipe,irun,2) = curr(:,1)'*curr(:,2);
        end
    end
    
    % assess similarity of FD / Dvars across pipelines
    for ipipe = 1:Npipes-1
        for ipipe_p = ipipe+1:Npipes
            for itype = 1:Ndata_types
                for irun = 1:Nruns
                    reliability(ipipe,ipipe_p,irun,isess,itype) = corr(data(:,irun,itype,ipipe,isess),...
                        data(:,irun,itype,ipipe_p,isess));
                end

            end
        end
    end
    
end

%%

for irun = 1:Nruns
    figure;imagesc(squeeze(nanmean(reliability(:,:,irun,:,2),4)),[0 1]);
    colorbar
    title([mot_levels{irun} ' ; Dvars reliability across pipelines'])
end

%% FD-DVARS correlation by pipeline separately for each motion condition

ylims = [-.08 .4];

imeas = 1;
meas_name = 'FD-DVARS correlation (z)';
meas_name_save = 'fd-dvars-corr_12_nox';

to_save = 0;
jitter = .1;

% clrs = [0 0 0; 0 0 .8; .2 .2 .6; 0 .8 0; .2 .6 .2; .8 0 0; .6 .2 .2; .8 0 .4; .6 .2 .4];
% % clrs = [0 0 0; .2 .2 .6; 0 0 .8; .2 .6 .2; 0 .8 0; .6 .2 .2; .8 0 0; .6 .2 .6;.8 0 .4;];
% x_ax = (1:Npipes) + jitter*[0 1 -1 1 -1 1 -1 1 -1];

clrs = [0 0 0; 0 0 .8; .2 .2 .6; 0 .8 0; .2 .6 .2];
x_ax = (1:Npipes) + jitter*[0 1 -1 1 -1];

x_tick = x_ax;

%separate plots by condition
for irun = 1:Nruns+1
        
    set_figsize([1100 500],1)
    
    if irun <= Nruns;
        dat = squeeze(Meas(:,:,irun,imeas));
        name = [save_dir mot_levels{irun} '_cond_' meas_name_save];

    else
        dat = squeeze(mean(Meas(:,:,:,imeas),3));
        name = [save_dir 'AvgMOTION_cond_' meas_name_save];

    end
    
    for ipipe = 1:Npipes
       
        abar_dat(dat(:,ipipe),x_ax(ipipe),clrs(ipipe,:),.75);
        hold on;
        
    end
    
    errorbar(x_ax,mean(dat),std(dat)/sqrt(size(dat,1)),'k','LineStyle','none','LineWidth',3)

    xlim([.3 Npipes+.7-jitter])
    
    ylabel(meas_name)
    
%     set(gca,'FontSize',24,'XTick',1:Npipes,'XTickLabel',pipe_nms,...
%         'YTick',0:.1:.3,'box','off')
    set(gca,'FontSize',18,'XTick',x_tick,'XTickLabel',pipe_nms,...
        'YTick',0:.1:.3,'box','off')
    
    if irun <= Nruns;
        title([mot_levels{irun} ' motion'])
    else title(['Avg across runs'])
    end
    ylim(ylims)
    
    if to_save;
        saveas(gcf,[name '.eps'],'epsc')
        saveas(gcf,[name '.fig'],'fig')
        close gcf
    end
    
end




%% FD-DVARS correlations by motion condition

to_save = 0;

imeas = 1;
meas_name = 'FD-DVARS corr (z)';
meas_name_save = 'fd-dvars-corr';

ylims = [0 .5];

for ipipe = 1
    
    pipe_name = 'Temp Filt only';
   
    set_figsize([600 500],1)
    
    dat = squeeze(Meas(:,ipipe,:,imeas));

    abar_dat(dat);
    
    set(gca,'FontSize',24,'XTick',1:Nruns,'XTickLabel',mot_levels,'box','off','YTick',0:.1:.3)
    ylim(ylims)
    
    xlim([.5 Nruns+.5])
    ylabel(meas_name)
    title([pipe_name ' pipeline'])
    
    if to_save;
        saveas(gcf,[save_dir pipe_name '_motconds.eps'],'epsc')
    end

    r_mot = corr(dat',[1:Nruns]');
    
end

pipe_name = 'Mean over pipelines';

set_figsize([600 500],1)

dat = squeeze(mean(Meas(:,:,:,imeas),2));

abar_dat(dat);

set(gca,'FontSize',24,'XTick',1:Nruns,'XTickLabel',mot_levels,'box','off','YTick',0:.1:.3)
ylim(ylims)

xlim([.5 Nruns+.5])
ylabel(meas_name)
title(pipe_name)
if to_save;
    saveas(gcf,[save_dir pipe_name '_motconds.eps'],'epsc')
end

%% mean DVARS plots

to_save = 0;

y = squeeze(mean(data(:,:,2,:,:)));
meas_name = 'mean DVARS';
meas_name_save = 'dvars-mn';
ylims = [0 8];

for ipipe = 1:2%:Npipes
    
    pipe_name = cat(2,pipe_nms{ipipe},' only');
   
    set_figsize([600 500],1)
    
    dat = squeeze(y(:,ipipe,:))';

    abar_dat(dat);
    
    set(gca,'FontSize',24,'XTick',1:Nruns,'XTickLabel',mot_levels,'box','off')
%         ,'YTick',0:.1:.3)
    ylim(ylims)
    
    xlim([.5 Nruns+.5])
    ylabel(meas_name)
    title([pipe_name ' pipeline'])
    
    if to_save;
        saveas(gcf,[save_dir pipe_name '_meanFD_motconds.eps'],'epsc')
    end

    r_mot = corr(dat',[1:Nruns]');
    
end

pipe_name = 'Mean over pipelines';

set_figsize([600 500],1)

dat = squeeze(mean(y,2))';

abar_dat(dat);

set(gca,'FontSize',24,'XTick',1:Nruns,'XTickLabel',mot_levels,'box','off')
% ,'YTick',0:.1:.3)
ylim(ylims)

xlim([.5 Nruns+.5])
ylabel(meas_name)
title(pipe_name)
if to_save;
    saveas(gcf,[save_dir pipe_name '_meanFD_motconds.eps'],'epsc')
end




%% individual run plots

save_dir_individual = [save_dir 'individual_run_plots/'];

to_save = 0;

pipe_idx = 1:5;
% pipe_idx = [1 2 3 5 7 8];
% pipe_idx = 1:Npipes;

clrs = [0*ones(1,3); .3*ones(1,3); 0 .8 0; 0 0 .8; 1 .7 .9; .9 0 0];
% clrs = [0*ones(1,3); 0 .8 0; 0 .8 0; 0 0 .8; 0 0 .8; 1 .7 .9; .9 0 0];
save_name_append = '_setof5__';
line_wd = [3;3;2*ones(length(pipe_idx)-2,1)];

% pipe_idx = [1 2 9 10];
% clrs = [0*ones(1,3); .3*ones(1,3); .2 .4 .2; .2 .6 .2];
% save_name_append = '_MCcomparison';
% line_wd = [3;3;2;2];


straight = '-';
dotted = '--';
line_style = cell(length(pipe_idx),1);
for ipipe = 1:length(pipe_idx);
    if sum(pipe_idx(ipipe) == [4 6])>0;
        line_style{ipipe} = dotted;
    else line_style{ipipe} = straight;
    end
end

for isess = 1:Nsessions
    
    for irun = Nruns
    
        dat = squeeze(data(:,irun,:,pipe_idx,isess));
        
        FD = dat(:,1,1);
        DVARS = squeeze(dat(:,2,:));
        
        set_figsize([1000 600],1)
        
        subplot(2,1,1)
        plot(FD,'Color',[.6 .6 .6],'LineWidth',3)
        xlim([0 Ntrs+1])
        ylim([0 1.03*max(FD)])
        ylabel('FD (mm)')
        title(['Sess ' num2str(isess) ' ; ' mot_levels{irun}])
        
        subplot(2,1,2)
        for ipipe = 1:length(pipe_idx)
            plot(DVARS(:,ipipe),line_style{ipipe},'LineWidth',line_wd(ipipe),'Color',clrs(ipipe,:))
            hold on;
        end
        xlim([0 Ntrs+1])
%         legend(pipe_nms(pipe_idx))
        ylim([.98*min(DVARS(:)) 1.01*max(DVARS(:))])
        xlabel('TR')
        ylabel('DVARS')
        
        if to_save==1;
            saveas(gcf,[save_dir_individual 'sess_' num2str(isess) ...
                '_' mot_levels{irun} save_name_append '.eps'],'epsc')
            close gcf
        end
    
    end
    
    
    
end



%%

dat_an = [];
for irun = 1:Nruns
    dat_an = [dat_an squeeze(Meas(:,:,irun,imeas))];
end

dat_an = squeeze(Meas(:,1,:,imeas));

% idx = 2
dat_an = [];
for irun = 1:Nruns
    dat_an = [dat_an squeeze(Meas(:,2:7,irun,imeas))];
end

idx = 3:6;
dat_an = [];
for irun = 1:Nruns
    dat_an = [dat_an squeeze(Meas(:,idx,irun,imeas))];
end
save('anova_MnC_GS.mat','dat_an')

idx = [3 5];
idx = [3 5 7 8]; %or 8
dat_an = [];
for irun = 1:Nruns
    dat_an = [dat_an squeeze(Meas(:,idx,irun,imeas))];
end
save('anova_MnCLocal.mat','dat_an')


idx = 3:10; %or 8
dat_an = [];
for irun = 1:Nruns
    dat_an = [dat_an squeeze(Meas(:,idx,irun,imeas))];
end
save('anova_4_GS.mat','dat_an')


