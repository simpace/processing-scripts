% 
clear all;
sub = 'sub01';

Ntrs_min = 60;

atlas = 'power';
scrub_type = 'SCRUB_v3';
scrub_label = 'Scrub_v_3 (.15 FD)';

% scrub_type = 'SCRUB_v2';
% scrub_label = 'Scrub_v_2 (4 min.)';

% scrub_type = 'SCRUB_v1';
% scrub_label = 'Scrub_v_1 (5 min.)';

corr_txt = 'CorrMtx_LFF';

analysis_dir = '/home/despo/simpace/analyses/';

data_dir = [analysis_dir 'roi-pair-distance/' sub '/'];
save_dir = [analysis_dir 'roi-pair-distance-figs/' sub '/scrub/' atlas '/'];

% load([analysis_dir 'pipeline_list_wo56.mat'])
load([analysis_dir 'mot_levels.mat'])

ipipe = 1;
pipes = {'NR_0'};
pipe_nms = {'BP'};

% pipes{ismember(pipes,'localWMreg')} = 'WMr';
pipeline = pipes{ipipe}; 
pipe_nm = pipe_nms{ipipe};

to_equate_ROIs = 1;

scrub_append{1} = ['*' scrub_type '*'];
scrub_append{2} = [];
Nscrub = 2;

Nsess = 13;
Dist_vec = cell(Nsess,1);
Z_mtx = Dist_vec;
Z_vec = Dist_vec;

Ntrs_full = 195;

if to_equate_ROIs==1;

    for isess = 1:Nsess
        
        if isess < 10; extra = '0'; else extra = ''; end
        sess_str = cat(2,'sess',extra,num2str(isess)); 
        
        roi_name_file = [pipeline '*' sess_str '*' atlas '*ROI_Name*.txt'];
        name_file = deal_names_full(roi_name_file, data_dir);
        roi_names = textread(name_file,'%s');

        if isess==1;
            roi_list = roi_names;
            roi_tally = ones(size(roi_list,1),1);
        else
            idx = ismember(roi_list,roi_names);
            roi_tally(idx) = roi_tally(idx)+1;
        end
        
        clear roi_names name_file idx;
    end
    
    roi_list = roi_list(roi_tally==Nsess);
    loc = cell(Nsess,1);

end



for isess = 1:Nsess

    if isess < 10; extra = '0'; else extra = ''; end
    sess_str = cat(2,'sess',extra,num2str(isess)); 
        
    if to_equate_ROIs==1;
        
        roi_name_file = [pipeline '*' sess_str '*' atlas '*ROI_Name*.txt'];
        name_file = deal_names_full(roi_name_file, data_dir);
        roi_names = textread(name_file,'%s');

        [~,loc{isess}] = ismember(roi_names,roi_list);
        loc{isess}(loc{isess}==0) = NaN;
        idx = loc{isess};
        idx = idx(~isnan(idx));
        
        clear name_file roi_names;
    end
    
    
    dist_file = deal_names_full([pipeline '*' sess_str '*' atlas '*ROI_Distance*.txt'], data_dir);
    dist = textread(dist_file);
    
    if to_equate_ROIs==0; idx = 1:size(dist,1); end       
    dist = dist(idx,idx);

    temp = triu(dist,1);
    Dist_vec{isess} = temp(temp~=0);
    
    for irun = 1:Nruns

        for iscrub = 1:Nscrub
            
            corr_file = deal_names_full([pipeline '*' sess_str '*' mot_levels{irun} '*'...
                atlas '*' corr_txt scrub_append{iscrub} '.txt'],data_dir);

            if iscrub==1;
                Ntrs(isess,irun) = str2double(corr_file(max(findstr(corr_file,'_'))+1:end-4));
            end
            if isess==1;
                corr_file
            end
            
            r = textread(corr_file);
            r = r(idx,idx);

            if irun==1 && iscrub==1;
                Z_mtx{isess} = NaN(size(r,1),size(r,2),Nruns,Nscrub);
                Z_vec{isess} = NaN(nchoosek(size(r,1),2),Nruns,Nscrub);
            end
            temp = rtoz(r); clear r;
            temp(isinf(temp)) = NaN;

            if Ntrs(isess,irun) > Ntrs_min;

                Z_mtx{isess}(:,:,irun,iscrub) = temp; 
                temp = triu(temp,1);
                Z_vec{isess}(:,irun,iscrub) = temp(temp~=0);

                clear temp;

            end
        end
        

    end


end

xcr_mean = cell2mat(cellfun(@mean,Z_vec,'UniformOutput',false));

%%

for isess = 1
    dat = Z_mtx{isess}(:,:,:,1);
    set_figsize([1400 400],1)

    ylims = [min(dat(~isnan(dat(:)))) 2];

    for irun = 1:Nruns

        subplot(1,Nruns,irun);imagesc(dat(:,:,irun),ylims);
        if irun==Nruns; colorbar; end
        title([mot_levels{irun}])
        set(gca,'XTickLabel',[],'YTickLabel',[])
%         title([mot_levels{irun} ', mean z = ' num2str(around(xcr_mean(irun), 3))])

    end
end


%%

clrs = [.5 .5 .5; .8 0 0; .6 0 .6; 0 0 .8];

x_string = 'ROI-pair distance (mm)';
to_save = 1;

if to_equate_ROIs==1;

    Npairs_per = length(Dist_vec{1});
    
    Dist = NaN(Npairs_per,Nsess);
    Z = NaN(Npairs_per,Nsess,Nruns);
    for isess = 1:Nsess
        Dist(:,isess) = Dist_vec{isess};
        Z(:,isess,:) = Z_vec{isess}(:,:,1)-Z_vec{isess}(:,:,2);
    end
    Dist_mn = mean(Dist,2);

    
    % Z_sort: ROI-pair x Nsess x Nruns
    % Z_block: ROI-pair x ROI-bin x Nsess x Nruns
    
    Nbins = 24;
    Ntot = length(Dist_mn);
    [Dist_sort,idx] = sort(Dist_mn);
    Z_sort = Z(idx,:,:);

    Nper = floor(Ntot/Nbins);

    Z_block = reshape(Z_sort(1:Nper*Nbins,:,:),[Nper,Nbins,Nsess,Nruns]);
    Dist_block = reshape(Dist_sort(1:Nper*Nbins),[Nper,Nbins]);
        
    x_ax = mean(Dist_block);

    Z_sort_mnsess = squeeze(nanmean(Z_sort,2)); % ROI-pair x Nruns
    Z_block_mn = squeeze(nanmean(mean(Z_block,1),3)); % ROI-bin x Nruns
    
    xl = [min(Dist_sort)-2 max(Dist_sort)+2];
    yl = [min(Z_sort_mnsess(:))-.05 max(Z_sort_mnsess(:))+.05];
    y_string = 'Mean corr. change (z)';
    save_string = [save_dir pipe_nm '_' scrub_type '_Mean_ROI_r_change_vs_distance_'];
    
    for irun = 1:Nruns
        set_figsize([900 400],1)
        y = Z_sort_mnsess(:,irun);
        y_mn = Z_block_mn(:,irun);
        
        plot(Dist_sort,y,'.','MarkerSize',15,'MarkerEdgeColor',clrs(irun,:))
        hold on;
        plot(x_ax,y_mn,'k','LineWidth',3)
        plot(x_ax,y_mn,'k.','MarkerSize',40)
        plot(xl,[0 0],'k--')
        xlim(xl);ylim(yl)
        xlabel(x_string)
        ylabel(y_string)
        title([pipe_nm ' ; ' scrub_label ' - NotScrub ; ' mot_levels{irun}]) 
        
        if to_save==1;
            save_nm = [save_string mot_levels{irun}];
            saveas(gcf,[save_nm '.eps'],'epsc')
            close gcf;
        end
    end
    
    set_figsize([900 400],1)
    y_mn = Z_block_mn;
    xl = [min(x_ax)*.85 max(x_ax)*1.05];
    dat = squeeze(mean(Z_block));    
    ntot = squeeze(sum(~isnan(dat),2));
    title_str = [pipe_nm ' ; ' scrub_label ' - NotScrub ; (N = ' num2str(ntot(1,:)) ')'];
    
    for irun = 1:Nruns
        plot(x_ax,y_mn(:,irun),'Color',clrs(irun,:),'LineWidth',3); hold on
    end
    for irun = 1:Nruns
        plot(x_ax,y_mn(:,irun),'.','MarkerEdgeColor',clrs(irun,:),'MarkerSize',40)
    end
%     legend(mot_levels)
    xlabel(x_string); ylabel(y_string)
    title(title_str) 
    plot(xl,[0 0],'k--')
    set(gca,'box','off','XLim',xl)
    ylim([-.08 .04])
    if to_save==1;
        save_nm = [save_string 'Summary'];
        saveas(gcf,[save_nm '.eps'],'epsc')
        close gcf;
    end
    

    set_figsize([900 400],1)
    y_mn = squeeze(nanmean(dat,2));
    y_sem = squeeze(nanstd(dat,0,2))./sqrt(ntot);
    for irun = 1:Nruns
        errorbar(x_ax,y_mn(:,irun),y_sem(:,irun),'Color',clrs(irun,:),'LineWidth',2); hold on
    end
    for irun = 1:Nruns
        plot(x_ax,y_mn(:,irun),'.','MarkerEdgeColor',clrs(irun,:),'MarkerSize',25)
    end
    plot(xl,[0 0],'k--')
%     legend(mot_levels)
    xlabel(x_string); ylabel(y_string)
    title(title_str) 
    set(gca,'box','off','XLim',xl,'YLim',1.5*[min(y_mn(:)-y_sem(:)) max(y_mn(:)+y_sem(:))])
    ylim([-.08 .04])

    if to_save==1;
        save_nm = [save_string 'Summary_errbars'];
        saveas(gcf,[save_nm '.eps'],'epsc')
        close gcf;
    end
    

end

%%
% Dist_all = [];
% Z_all = [];
% 
% for isess = 1:Nsess
%     
%     Dist_all = cat(1,Dist_all,Dist_vec{isess});
%     Z_all = cat(1,Z_all,Z_vec{isess});
%     
% end
% 
% Nbins = 13;
% Ntot = length(Z_all);
% [Dist_sort,idx] = sort(Dist_all);
% Z_sort = Z_all(idx,:);
% 
% Nper = floor(size(Z_sort,1)/Nbins);
% 
% Z_block = reshape(Z_sort(1:Nper*Nbins,:),[Nper,Nbins,Nruns]);





%%

% xlims = [min(Dist_all)-10 max(Dist_all)+10];
% 
% set_figsize([1200 600],1)
% for irun = 1:Nruns-1:Nruns
% 
%     plot(Dist_all,Z_all(:,irun),[clrs(irun) '.'],'MarkerSize',10)
%     hold on;
%     
% end
% xlim(xlims)
% 
% set_figsize([1200 600],1)
% plot(Dist_all,Z_all(:,end)-Z_all(:,1),'.','MarkerSize',10);
% hold on;
% plot(xlims,[0 0],'k','LineWidth',2);
% xlim(xlims)
