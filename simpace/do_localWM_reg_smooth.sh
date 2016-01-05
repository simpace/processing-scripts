# script to run AFNI's ANATICOR (local WM regression) on simpace data


WD='/home/despo/simpace/rename_files/sub01/'
cd ${WD}
for session in $(ls -d sess*); do #
	
	cd ${WD}/${session}/preproc

	# where to put stuff
	if [ ! -d ${WD}/${session}/preproc/localWMreg+GS_smooth ]; then
		mkdir ${WD}/${session}/preproc/localWMreg+GS_smooth/
	fi
	

	#split motion regressors
	cd ${WD}/${session}/preproc/motion_params/
	1d_tool.py -infile ${WD}/${session}/preproc/motion_params/MP.txt -select_rows '0..194' \
	-write ${WD}/${session}/preproc/localWMreg+GS_smooth/${session}_MP_run01.1D
	
	1d_tool.py -infile ${WD}/${session}/preproc/motion_params/MP.txt -select_rows '195..389' \
	-write ${WD}/${session}/preproc/localWMreg+GS_smooth/${session}_MP_run02.1D
	
	1d_tool.py -infile ${WD}/${session}/preproc/motion_params/MP.txt -select_rows '390..584' \
	-write ${WD}/${session}/preproc/localWMreg+GS_smooth/${session}_MP_run03.1D
	
	1d_tool.py -infile ${WD}/${session}/preproc/motion_params/MP.txt -select_rows '585..779' \
	-write ${WD}/${session}/preproc/localWMreg+GS_smooth/${session}_MP_run04.1D

	#create WM mask
	if [ ! -e ${WD}/${session}/preproc/localWMreg+GS_smooth/${session}_WM_mask_erode1x.nii.gz ]; then
		3dmask_tool -input ${WD}/${session}/preproc/wm_in_func_res/wm_func_res_0.85.nii.gz -dilate_inputs -1 \
		-prefix ${WD}/${session}/preproc/localWMreg+GS_smooth/${session}_WM_mask_erode1x.nii.gz
	fi	


	for run in run01 run02 run03 run04; do

		cd ${WD}/${session}/preproc/localWMreg+GS_smooth/

		#merge motion corrected preproc files into one 4D nifit
		#gzip ${WD}/${session}/preproc/realign/asub01_${session}_${run}.nii 
		if [ ! -e ${WD}/${session}/preproc/localWMreg+GS_smooth/${session}_preproc_${run}.nii.gz ]; then
			fslmerge -t ${WD}/${session}/preproc/localWMreg+GS_smooth/${session}_preproc_${run}.nii.gz $(ls ${WD}/${session}/preproc/smooth/*${run}*.nii)
			# fix TR info
			3drefit -TR 2 ${WD}/${session}/preproc/localWMreg+GS_smooth/${session}_preproc_${run}.nii.gz 
		fi

		
		#get csf regressors
		3dmaskave -mask ${WD}/${session}/preproc/csf_in_func_res/csf_func_res_final.nii.gz -quiet \
		${WD}/${session}/preproc/localWMreg+GS_smooth/${session}_preproc_${run}.nii.gz > ${WD}/${session}/preproc/localWMreg+GS_smooth/${session}_csf_${run}.1D

		#get global signal
		3dmaskave -mask ${WD}/${session}/preproc/sess_mask/sess_mask.nii -quiet \
		${WD}/${session}/preproc/localWMreg+GS_smooth/${session}_preproc_${run}.nii.gz > ${WD}/${session}/preproc/localWMreg+GS_smooth/${session}_GS_${run}.1D

		#get local WM regressors	
		if [ ! -e ${WD}/${session}/preproc/localWMreg+GS_smooth/${session}_WMeLOCAL_${run}.nii.gz ]; then
			3dLocalstat -prefix ${WD}/${session}/preproc/localWMreg+GS_smooth/${session}_WMeLOCAL_${run}.nii.gz \
			-nbhd 'SPHERE(50)' \
			-stat mean \
			-mask ${WD}/${session}/preproc/localWMreg+GS_smooth/${session}_WM_mask_erode1x.nii.gz \
			-use_nonmask ${WD}/${session}/preproc/localWMreg+GS_smooth/${session}_preproc_${run}.nii.gz
		fi	

		##do regression
		#the model at this point is Y = [ Motion + localWMreg+GS_smooth + Ventricles ] X + Residual 
		# bandpass filtering at 0.01 0.08, simultaneously with regression
		# spatial smooth after regression
		
		if [ ! -e ${WD}/${session}/preproc/localWMreg+GS_smooth/${session}_preproc_localWMreg+GS_${run}.nii.gz ]; then
			3dTproject \
			-input ${WD}/${session}/preproc/localWMreg+GS_smooth/${session}_preproc_${run}.nii.gz \
			-prefix ${WD}/${session}/preproc/localWMreg+GS_smooth/${session}_preproc_localWMreg+GS_${run}.nii.gz \
			-ort ${WD}/${session}/preproc/localWMreg+GS_smooth/${session}_MP_${run}.1D \
			-ort ${WD}/${session}/preproc/localWMreg+GS_smooth/${session}_csf_${run}.1D \
			-ort ${WD}/${session}/preproc/localWMreg+GS_smooth/${session}_GS_${run}.1D \
			-dsort ${WD}/${session}/preproc/localWMreg+GS_smooth/${session}_WMeLOCAL_${run}.nii.gz \
			-automask \
			-passband 0.01 0.08 \
			-polort 2
		fi

	done
done