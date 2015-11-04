# script to run AFNI's ANATICOR (local WM regression) on simpace data


WD='/home/despo/simpace/rename_files/sub01/'

for session in sess01; do #$(ls -d sess*)
	
	cd ${WD}/${session}/preproc

	# where to put stuff
	if [ ! -d ${WD}/${session}/preproc/localWMreg ]; then
		mkdir ${WD}/${session}/preproc/localWMreg/
	fi
	

	#split motion regressors
	cd ${WD}/${session}/preproc/motion_params/
	1d_tool.py -infile ${WD}/${session}/preproc/motion_params/MP.txt -select_rows '0..195' \
	-write ${WD}/${session}/preproc/localWMreg/${session}_MP_run01.1D
	
	1d_tool.py -infile ${WD}/${session}/preproc/motion_params/MP.txt -select_rows '196..391' \
	-write ${WD}/${session}/preproc/localWMreg/${session}_MP_run02.1D
	
	1d_tool.py -infile ${WD}/${session}/preproc/motion_params/MP.txt -select_rows '392..587' \
	-write ${WD}/${session}/preproc/localWMreg/${session}_MP_run03.1D
	
	1d_tool.py -infile ${WD}/${session}/preproc/motion_params/MP.txt -select_rows '588..783' \
	-write ${WD}/${session}/preproc/localWMreg/${session}_MP_run04.1D

	#create WM mask
	3dmask_tool -input ${WD}/${session}/preproc/wm_in_func_res/wm_func_res_0.85.nii.gz -dilate_inputs -1 -prefix ${WD}/${session}/preproc/localWMreg/${session}_WM_mask_erode1x.nii.gz
		


	for run in run01 run02 run03 run04; do

		cd ${WD}/${session}/preproc/localWMreg/

		#merge motion corrected preproc files into one 4D nifit
		#gzip ${WD}/${session}/preproc/realign/asub01_${session}_${run}.nii 
		fslmerge -t ${WD}/${session}/preproc/localWMreg/${session}_preproc_${run}.nii.gz $(ls ${WD}/${session}/preproc/realign/*${run}*.nii.gz)
		
		# fix TR info
		3drefit -TR 2 ${WD}/${session}/preproc/localWMreg/${session}_preproc_${run}.nii.gz 

		#get csf regressors
		3dmaskave -mask ${WD}/${session}/preproc/csf_in_func_res/csf_func_res_final.nii.gz -quiet \
		${WD}/${session}/preproc/localWMreg/${session}_preproc_${run}.nii.gz > ${WD}/${session}/preproc/localWMreg/${session}_csf_${run}.1D

		#get local WM regressors	
		3dLocalstat -prefix ${WD}/${session}/preproc/localWMreg/${session}_WMeLOCAL_${run}.nii.gz \
		-nbhd 'SPHERE(50)' \
		-stat mean \
		-mask ${WD}/${session}/preproc/localWMreg/${session}_WM_mask_erode1x.nii.gz \
		-use_nonmask ${WD}/${session}/preproc/localWMreg/${session}_preproc_${run}.nii.gz


		##do regression
		#the model at this point is Y = [ Motion + localWMreg + Ventricles ] X + Residual 
		# bandpass filtering at 0.01 0.08, simultaneously with regression
		# spatial smooth after regression
		3dTproject \
		-input ${WD}/${session}/preproc/localWMreg/${session}_preproc_${run}.nii.gz \
		-prefix ${WD}/${session}/preproc/localWMreg/${session}_preproc_localWMreg_${run}.nii.gz \
		-ort ${WD}/${session}/preproc/localWMreg/${session}_MP_${run}.1D \
		-ort ${WD}/${session}/preproc/localWMreg/${session}_csf_${run}.1D \
		-dsort ${WD}/${session}/preproc/localWMreg/${session}_WMeLOCAL_${run}.nii.gz \
		-automask \
		-passband 0.01 0.08 \
		-blur 6 \
		-polort 2

	done
done