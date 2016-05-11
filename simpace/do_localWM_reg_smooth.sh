# script to run AFNI's ANATICOR (local WM regression) on simpace data


for sub in sub-01 sub-02; do

	WD="/home/despo/simpace/rename_files/derivatives/${sub}/"

	cd ${WD}
	for session in $(/bin/ls -d ses-t*); do #$(ls -d sess*)
		
		cd ${WD}/${session}/func_proc

		rm -rf ${WD}/${session}/func_proc/localWMreg_smooth/
		# where to put stuff
		if [ ! -d ${WD}/${session}/func_proc/localWMreg_smooth ]; then
			mkdir ${WD}/${session}/func_proc/localWMreg_smooth/
		fi
		

		#split motion regressors
		#cd ${WD}/${session}/func_proc/motion_params/
		ln -s ${WD}/${session}/func_proc/motion_params/MP_NONE.txt ${WD}/${session}/func_proc/localWMreg_smooth/${session}_MP_NONE.1D
		ln -s ${WD}/${session}/func_proc/motion_params/MP_LOW.txt ${WD}/${session}/func_proc/localWMreg_smooth/${session}_MP_LOW.1D
		ln -s ${WD}/${session}/func_proc/motion_params/MP_MED.txt ${WD}/${session}/func_proc/localWMreg_smooth/${session}_MP_MED.1D
		ln -s ${WD}/${session}/func_proc/motion_params/MP_HIGH.txt ${WD}/${session}/func_proc/localWMreg_smooth/${session}_MP_HIGH.1D

		#create WM mask
		if [ ! -e ${WD}/${session}/func_proc/localWMreg_smooth/${session}_WM_mask_erode1x.nii.gz ]; then
			3dmask_tool -input ${WD}/${session}/func_proc/wm_in_func_res/wm_func_res_0.85.nii.gz -dilate_inputs -1 \
			-prefix ${WD}/${session}/func_proc/localWMreg_smooth/${session}_WM_mask_erode1x.nii.gz
		fi	


		for run in NONE LOW MED HIGH; do

			cd ${WD}/${session}/func_proc/localWMreg_smooth/

			#merge motion corrected preproc files into one 4D nifit
			#gzip ${WD}/${session}/func_proc/realign/asub01_${session}_${run}.nii 
			if [ ! -e ${WD}/${session}/func_proc/localWMreg_smooth/${session}_preproc_${run}.nii.gz ]; then
				#fslmerge -t ${WD}/${session}/func_proc/localWMreg_smooth/${session}_preproc_${run}.nii.gz $(ls ${WD}/${session}/func_proc/smooth/*${run}*.nii)
				# fix TR info
				#3drefit -TR 2 ${WD}/${session}/func_proc/localWMreg_smooth/${session}_preproc_${run}.nii.gz 
				ln -s ${WD}/${session}/func_proc/smooth/sra${sub}_${session}_task-rest_acq-${run}_bold.nii ${WD}/${session}/func_proc/localWMreg_smooth/${session}_preproc_${run}.nii
			fi

			
			#get csf regressors
			3dmaskave -mask ${WD}/${session}/func_proc/csf_in_func_res/csf_func_res_final.nii.gz -quiet \
			${WD}/${session}/func_proc/localWMreg_smooth/${session}_preproc_${run}.nii > ${WD}/${session}/func_proc/localWMreg_smooth/${session}_csf_${run}.1D

			#get global signal
			3dmaskave -mask ${WD}/${session}/func_proc/sess_mask/sess_mask.nii -quiet \
			${WD}/${session}/func_proc/localWMreg_smooth/${session}_preproc_${run}.nii > ${WD}/${session}/func_proc/localWMreg_smooth/${session}_gs_${run}.1D

			#get local WM regressors	
			if [ ! -e ${WD}/${session}/func_proc/localWMreg_smooth/${session}_WMeLOCAL_${run}.nii.gz ]; then
				3dLocalstat -prefix ${WD}/${session}/func_proc/localWMreg_smooth/${session}_WMeLOCAL_${run}.nii.gz \
				-nbhd 'SPHERE(50)' \
				-stat mean \
				-mask ${WD}/${session}/func_proc/localWMreg_smooth/${session}_WM_mask_erode1x.nii.gz \
				-use_nonmask ${WD}/${session}/func_proc/localWMreg_smooth/${session}_preproc_${run}.nii
			fi	

			##do regression
			#the model at this point is Y = [ Motion + localWMreg_smooth + Ventricles ] X + Residual 
			# bandpass filtering at 0.01 0.08, simultaneously with regression
			# spatial smooth after regression
			
			if [ ! -e ${WD}/${session}/func_proc/localWMreg_smooth/${session}_preproc_localWMreg_${run}.nii.gz ]; then
				3dTproject \
				-input ${WD}/${session}/func_proc/localWMreg_smooth/${session}_preproc_${run}.nii \
				-prefix ${WD}/${session}/func_proc/localWMreg_smooth/${session}_preproc_localWMreg_${run}.nii.gz \
				-ort ${WD}/${session}/func_proc/localWMreg_smooth/${session}_MP_${run}.1D \
				-ort ${WD}/${session}/func_proc/localWMreg_smooth/${session}_csf_${run}.1D \
				-dsort ${WD}/${session}/func_proc/localWMreg_smooth/${session}_WMeLOCAL_${run}.nii.gz \
				-automask \
				-passband 0.01 0.08 \
				-polort 2
			fi

		done
	done
done