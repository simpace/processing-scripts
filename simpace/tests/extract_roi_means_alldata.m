
main_dir = '/home/despo/simpace/rename_files/';

sub_dir = [main_dir 'sub01/'];

data_types = {'smooth';'NR_0';'NR_1';'NR_2';'NR_3';'NR_4';'HPF';'NR_5'};

Ndata_types = length(data_types);
Nsess = 13;
Nruns = 4;
Ntrs = 195;
% xmax = 20;
% save_dir = [pwd '/mean_figs/'];
% mkdir(save_dir)

for isess = 1:Nsess

    if isess < 10; sess_str = ['0' num2str(isess)];
    else sess_str = num2str(isess);
    end

    preproc_dir = [sub_dir 'sess' sess_str '/preproc/'];

    roi_dir = [preproc_dir 'coreg/'];
    roi_files = cell2mat( deal_names_full('rraal*.nii', roi_dir) );
    Nrois = length(roi_files);
    
    sess_mask = deal_names_full('sess_mask.nii',[preproc_dir 'sess_mask']);
    
    for itype = 1:Ndata_types

        data_type = data_types{itype};
        
        data_dir = [preproc_dir data_type '/'];

        for irun = 1:Nruns
            
            TS = NaN(Ntrs, Nrois);

            files = cell2mat( deal_names_full(['*ra*run0' num2str(irun) '*.nii'], data_dir) );

            for iroi = 1:Nrois
                try
                    TS(:,iroi) = rex(files,roi_files(iroi,:) , 'conjunction_mask' , sess_mask);
                catch
                    disp(['roi not good: ' roi_files(iroi,:)])
                end
            end
            
            out_filename = [data_type '_sess' num2str(isess) '_run' num2str(irun)];
            
            
            disp(['sess ' num2str(isess) ', run ' num2str(irun) ' : ' files(1,:)])
            disp(['    ' files(end,:)])

%             for iplot = 1:2
%                 subplot(Nruns,2, iplot+ (irun-1)*2)
%                 set(gca,'FontSize',18)
% 
%                 if iplot==1; 
%                     plot(means)
%                     ylabel(['run ' num2str(irun)])
%                     if irun==1; title(['session ' num2str(isess) ', raw']); end
% 
%                 elseif iplot==2;                 
%                     plot(zscore(means))
%                     if irun==1; title(['session ' num2str(isess) ', z-scored']); end
% 
%                 end
%                 xlim([1 xmax])
% 
%             end

        end
%         saveas(gcf,[save_dir 'sub01_sess' sess_str '.eps'],'epsc')
%         close gcf
        a = 1;


    end

end
