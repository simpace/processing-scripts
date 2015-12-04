
main_dir = '/home/despo/simpace/rename_files/';

sub_dir = [main_dir 'sub01/'];

Nsess = 13;
Nruns = 4;
xmax = 20;
save_dir = [pwd '/mean_figs/'];
mkdir(save_dir)

for isess = 8:Nsess
    
    if isess < 10; sess_str = ['0' num2str(isess)];
    else sess_str = num2str(isess);
    end
    
    preproc_dir = [sub_dir 'sess' sess_str '/preproc/'];
    data_dir = [preproc_dir 'smooth/'];
    
    roi_dir = [preproc_dir 'coreg/'];
    roi_files = cell2mat( deal_names_full('rraal*.nii', roi_dir) );
    
    set_figsize([1000 750],1)
    
    for irun = 1:Nruns
        
        files = cell2mat( deal_names_full(['*ra*run0' num2str(irun) '*.nii'], data_dir) );
        means = rex(files,roi_files );
        disp(['sess ' num2str(isess) ', run ' num2str(irun) ' : ' files(1,:)])
        disp(['    ' files(end,:)])
    
        for iplot = 1:2
            subplot(Nruns,2, iplot+ (irun-1)*2)
            set(gca,'FontSize',18)
            
            if iplot==1; 
                plot(means)
                ylabel(['run ' num2str(irun)])
                if irun==1; title(['session ' num2str(isess) ', raw']); end
             
            elseif iplot==2;                 
                plot(zscore(means))
                if irun==1; title(['session ' num2str(isess) ', z-scored']); end

            end
            xlim([1 xmax])
        
        end
        
    end
    saveas(gcf,[save_dir 'sub01_sess' sess_str '.eps'],'epsc')
    close gcf
    a = 1;
        
    
end

