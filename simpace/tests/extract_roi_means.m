function extract_roi_means( data_path, apply_mask, save_dir)
%data_path, ex: /home/despo/simpace/rename_files/sub01/sess01/preproc/NR_0
%apply_mask 0 or 1
%save_dir: where to save, defaults to data_path/../extracted_signals

Nruns = 4;

if ~exist('deal_names','file'); addpath('helpers'); end

end_path = max(strfind(data_path,'/'));
preproc_dir = [data_path(1:end_path-1) '/'];

pipeline_name = data_path(end_path+1:end);

s = strfind(preproc_dir, 'sub');
sub_str = preproc_dir( s:s+4 );

s = strfind(preproc_dir, 'sess');
sess_str = preproc_dir( s:s+5 );

roi_dir = [preproc_dir 'coreg/'];

if nargin < 3;
    save_dir = [preproc_dir 'extracted_signals/'];
elseif isempty(strmatch(save_dir(end),'/'));
    save_dir = [save_dir '/'];
end

if ~isdir(save_dir); mkdir(save_dir); end

roi_files = cell2mat( deal_names_full('rraal*.nii', roi_dir) );
roi_names = deal_names('rraal*.nii', roi_dir);
Nrois = length(roi_files);
    
if apply_mask==1;
    sess_mask = deal_names_full('sess_mask.nii',[preproc_dir 'sess_mask']);
end


for irun = 1:Nruns

%     files = cell2mat( deal_names_full(['*ra*run0' num2str(irun) '*.nii'], data_dir) );
    files = cell2mat( deal_names_full(['*run0' num2str(irun) '*.nii'], data_path) );
    
    for iroi = 1:Nrois
        try
            if apply_mask==1;
                ts_roi = rex(files,roi_files(iroi,:),'conjunction_mask',sess_mask);
            else
                ts_roi = rex(files,roi_files(iroi,:));
            end        
            
            if iroi==1; TS = NaN(length(ts_roi),Nrois); end
            
            TS(:,iroi) = ts_roi; clear ts_roi
            
        catch
            disp(['roi not good: ' roi_files(iroi,:) ', run ' num2str(irun)])
        end
    end

    out_filename = [pipeline_name '_signal_matlab_' sub_str '_' sess_str '_run0' num2str(irun) '_mask' num2str(apply_mask)];

    save([save_dir out_filename '.mat'], 'TS', 'roi_names', 'roi_files')
    csvwrite([save_dir out_filename '.txt'], TS)
    dlmwrite([save_dir out_filename '_rois.txt'], roi_names ,'')

    disp([pipeline_name ' ' sess_str ', run ' num2str(irun) ' : ' files(1,:)])
    disp(['    ' files(end,:)])
    
end
