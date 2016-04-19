
clear all;
sub = 'sub01';

atlas = 'aal';
corr_txt = 'CorrMtx_LFF';

to_save_file=0;

analysis_dir = '/home/despo/simpace/analyses/';

data_dir = [analysis_dir 'roi-pair-distance/' sub '/'];
save_dir = [analysis_dir 'roi-pair-distance-figs/' sub '/' atlas '/'];

% load([analysis_dir 'pipeline_list_wo56.mat'])


load([analysis_dir 'mot_levels.mat'])
Nmotion = Nruns;

ipipe = 2;

% pipes{ismember(pipes,'localWMreg')} = 'WMr';
pipeline = 'NR_0';
pipe_nm = 'BP';

to_equate_ROIs = 1;

Nsess = 13;
Dist_vec = cell(Nsess,1);
Z_mtx = Dist_vec;
Z_vec = Dist_vec;


if to_equate_ROIs==1;

    for isess = 1:Nsess
        
        if isess < 10; extra = '0'; else extra = ''; end
        sess_str = cat(2,'sess',extra,num2str(isess)); 
        
        roi_name_file = [pipeline '*' sess_str '*' atlas '*ROI_Name*.txt'];
        name_file = deal_names_full(roi_name_file, data_dir);
        
        if length(name_file)>1; 
            roi_name_file = [pipeline '_' sess_str '*' atlas '*ROI_Name*.txt'];
            name_file = deal_names_full(roi_name_file, data_dir);        
        end
        
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
        
         if length(name_file)>1; 
            roi_name_file = [pipeline '_' sess_str '*' atlas '*ROI_Name*.txt'];
            name_file = deal_names_full(roi_name_file, data_dir);        
         end
        
        roi_names = textread(name_file,'%s');

        [~,loc{isess}] = ismember(roi_names,roi_list);
        loc{isess}(loc{isess}==0) = NaN;
        idx = loc{isess};
        idx = idx(~isnan(idx));
        
        clear name_file roi_names;
    end
    
    
    dist_file = deal_names_full([pipeline '*' sess_str '*' atlas '*ROI_Distance*.txt'], data_dir);
    if length(dist_file)>1; 
        dist_file = deal_names_full([pipeline '_' sess_str '*' atlas '*ROI_Distance*.txt'], data_dir);
    end
    dist = textread(dist_file);        

    if to_equate_ROIs==0; idx = 1:size(dist,1); end       
    dist = dist(idx,idx);

    temp = triu(dist,1);
    Dist_vec{isess} = temp(temp~=0);
    
    for imotion = 1:Nmotion

        corr_file = deal_names_full([pipeline '_' sess_str '*' mot_levels{imotion} '*' atlas '*' corr_txt '.txt'],...
            data_dir);
        r = textread(corr_file);
        r = r(idx,idx);
        if imotion==1 && isess==1; 
            figure;imagesc(r);colorbar
        end

        if imotion==1;
            Z_mtx{isess} = NaN(size(r,1),size(r,2),Nmotion);
            Z_vec{isess} = NaN(nchoosek(size(r,1),2),Nmotion);
        end
        temp = rtoz(r); clear r;
        temp(isinf(temp)) = NaN;

        Z_mtx{isess}(:,:,imotion) = temp; 
        temp = triu(temp,1);
        Z_vec{isess}(:,imotion) = temp(temp~=0);

        clear temp;

    end


end

xcr_mean = cell2mat(cellfun(@mean,Z_vec,'UniformOutput',false));

%%

for isess = 1%Nsess
    dat = Z_mtx{isess};
    set_figsize([1400 400],1)

    ylims = [min(dat(:)) 2];

    for imotion = 1:Nmotion

        subplot(1,Nmotion,imotion);imagesc(dat(:,:,imotion),ylims);
        if imotion==Nmotion; colorbar; end
        title([mot_levels{imotion}])
        set(gca,'XTickLabel',[],'YTickLabel',[])
%         title([mot_levels{imotion} ', mean z = ' num2str(around(xcr_mean(imotion), 3))])

    end
end

%% for an individual session

% clrs = ['k';'g';'b';'r'];
% xlims = [min(Dist_vec)-10 max(Dist_vec)+10];
% 
% set_figsize([1200 600],1)
% for imotion = 1:Nmotion-1:Nmotion
% 
%     plot(Dist_vec,Z_vec(:,imotion),[clrs(imotion) '.'],'MarkerSize',10)
%     hold on;
%     
% end
% xlim(xlims)
% 
% set_figsize([1200 600],1)
% plot(Dist_vec,Z_vec(:,end)-Z_vec(:,1),'.','MarkerSize',10);
% hold on;
% plot(xlims,[0 0],'k','LineWidth',2);
% xlim(xlims)

%%

clrs = [.5 .5 .5; .8 0 0; .6 0 .6; 0 0 .8];

x_string = 'ROI-pair distance (mm)';
to_save = 0;
to_disp = 1;

Nbins = 20;

if to_equate_ROIs==1;

    Npairs_per = length(Dist_vec{1});
    
    Dist = NaN(Npairs_per,Nsess);
    Z = NaN(Npairs_per,Nsess,Nmotion);
    for isess = 1:Nsess
        Dist(:,isess) = Dist_vec{isess};
        Z(:,isess,:) = Z_vec{isess};
    end
    Dist_mn = mean(Dist,2);
    
    % Z_sort: ROI-pair x Nsess x Nmotion
    % Z_block: ROI-pair x ROI-bin x Nsess x Nmotion
    
    Ntot = length(Dist_mn);
    [Dist_sort,idx] = sort(Dist_mn);
    Z_sort = Z(idx,:,:);

    Nper = floor(Ntot/Nbins);

    Z_block = reshape(Z_sort(1:Nper*Nbins,:,:),[Nper,Nbins,Nsess,Nmotion]);
    Dist_block = reshape(Dist_sort(1:Nper*Nbins),[Nper,Nbins]);
        
    x_ax = mean(Dist_block);

    Z_sort_mnsess = squeeze(nanmean(Z_sort,2)); % ROI-pair x Nmotion
    Z_block_mn = squeeze(mean(nanmean(Z_block,1),3)); % ROI-bin x Nmotion
    
    xl = [min(Dist_sort)-2 max(Dist_sort)+2];
    yl = [min(Z_sort_mnsess(:))-.05 max(Z_sort_mnsess(:))+.05];
    
    y_string = 'Mean correlation (z)';
    save_string = [save_dir pipe_nm '_Mean_ROI_r_vs_distance_'];
    
    % connectivity x ROI distance for each motion condition
    for imotion = 1:Nmotion
        set_figsize([900 400],to_disp)
        y = Z_sort_mnsess(:,imotion);
        y_mn = Z_block_mn(:,imotion);
        
        plot(Dist_sort,y,'.','MarkerSize',15,'MarkerEdgeColor',clrs(imotion,:))
        hold on;
        plot(x_ax,y_mn,'k','LineWidth',3)
        plot(x_ax,y_mn,'k.','MarkerSize',40)
        plot(xl,[0 0],'k--')
        xlim(xl);ylim(yl)
        xlabel(x_string)
        ylabel(y_string)
        title([pipe_nm ' ; ' mot_levels{imotion}]) 
        
        if to_save==1;
            save_nm = [save_string mot_levels{imotion}];
            saveas(gcf,[save_nm '.eps'],'epsc')
            close gcf;
        end
    end
    
    % mean connectivity x ROI pair distance for all motion conditions
    set_figsize([900 400],1)
    y_mn = Z_block_mn;
    for imotion = 1:Nmotion
        plot(x_ax,y_mn(:,imotion),'Color',clrs(imotion,:),'LineWidth',3); hold on
    end
    for imotion = 1:Nmotion
        plot(x_ax,y_mn(:,imotion),'.','MarkerEdgeColor',clrs(imotion,:),'MarkerSize',40)
    end
%     legend(mot_levels)
    xlim([min(x_ax)*.85 max(x_ax)*1.05])
    xlabel(x_string)
    ylabel(y_string)
    title([pipe_nm]) 
    set(gca,'box','off')
    ylim([.2 .8])
    if to_save==1;
        save_nm = [save_string 'Summary'];
        saveas(gcf,[save_nm '.eps'],'epsc')
        close gcf;
    end
    
    % average connectivity difference from NONE for each motion condition
    Z_diff = repmat(Z_block(:,:,:,1),[1 1 1 Nmotion-1]) - Z_block(:,:,:,2:Nmotion);
    
    Z_diff_dat = squeeze(nanmean(Z_diff,3));
    Z_diff_mn = squeeze(nanmean(Z_diff_dat,1));
    
    Z_out = squeeze(nanmean(Z_diff));
    if to_save_file==1;
        save([data_dir 'r_change_summary_' pipeline '_' atlas '.mat'],'Z_out','x_ax')
        disp([pipeline ' mean  = ' around(squeeze(mean(mean(Z_out),2))')])
    end
    
    yl = [min(Z_diff_dat(:))-.05 max(Z_diff_dat(:))+.05];
%     yl = [-.8 max(Z_diff_dat(:))+.05];
    y_string = 'Mean corr. change (z)';
    save_string = [save_dir pipe_nm '_Mean_ROI_r_change_vs_distance_'];
    
    for imotion = 2:Nmotion
        set_figsize([900 400],to_disp)
        y = Z_diff_dat(:,:,imotion-1);
        y_mn = Z_diff_mn(:,imotion-1);
        plot(Dist_block(:),y(:),'.','MarkerSize',15,'MarkerEdgeColor',clrs(imotion,:))
        hold on;
        plot(x_ax,y_mn,'k','LineWidth',3)
        plot(x_ax,y_mn,'k.','MarkerSize',40)
        plot(xl,[0 0],'k--')
        xlim(xl);ylim(yl)
        xlabel(x_string)
        ylabel(y_string)
        title([pipe_nm ' ; NONE - ' mot_levels{imotion} ' corr. difference']) 
        if to_save==1;
            save_nm = [save_string mot_levels{imotion}];
            saveas(gcf,[save_nm '.eps'],'epsc')
            close gcf;
        end
    end
    
    % average connectivity difference from NONE for all motion conditions
    set_figsize([900 400],to_disp)
    y_mn = Z_diff_mn;
    for imotion = 2:Nmotion
        plot(x_ax,y_mn(:,imotion-1),'Color',clrs(imotion,:),'LineWidth',3); hold on
    end
    for imotion = 2:Nmotion
        plot(x_ax,y_mn(:,imotion-1),'.','MarkerEdgeColor',clrs(imotion,:),'MarkerSize',40)
    end
%     legend(mot_levels(2:Nmotion),'Location','Best')
    xlim([min(x_ax)*.85 max(x_ax)*1.05])
    xlabel(x_string)
    ylabel(y_string)
    title([pipe_nm ' ; NONE - {mot} corr. difference']) 
    set(gca,'box','off')

    if to_save==1;
        save_nm = [save_string 'Summary'];
        saveas(gcf,[save_nm '.eps'],'epsc')
        close gcf;
    end
    
    % connectivity difference from NONE per motion condition for each
    % session
    for imotion = 2:Nmotion
        dat = squeeze(mean(Z_diff(:,:,:,imotion-1)));
        
        set_figsize([900 400],to_disp)
        plot(x_ax, dat, 'Color', clrs(imotion,:), 'LineWidth', 2)
        hold on;
        xlim([-5 5]+minmax(x_ax))
        xlabel(x_string)
        ylabel(y_string)
        title([pipe_nm ' ; NONE - ' mot_levels{imotion} ' corr. diff. (PER session)']) 
    
    end
end
