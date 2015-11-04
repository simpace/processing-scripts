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

	for run in run01 run02 run03 run04; do

		cd ${WD}/${session}/preproc/localWMreg/${session}

		#merge smoothed preproc files into one 4D nifit
		fslmerge -t ${WD}/${session}/preproc/localWMreg/${session}_preproc_smooth_${run}.nii.gz $(ls ${WD}/${session}/preproc/smooth/*${run}*.nii)
		
		# fix TR info
		3drefit -TR 2 ${WD}/${session}/preproc/localWMreg/${session}_preproc_smooth_${run}.nii.gz 

		#get csf regressors
		3dmaskave -mask ${WD}/${session}/preproc/csf_in_func_res/csf_func_res_final.nii.gz -quiet \
		${WD}/${session}/preproc/localWMreg/${session}_preproc_smooth_${run}.nii.gz > ${WD}/${session}/preproc/localWMreg/${session}_csf_${run}.1D

		#get local WM regressors
		3dLocalstat -prefix ${WD}/${session}/preproc/localWMreg/${session}_WMeLOCAL_${run}.nii.gz \
		-nbhd 'SPHERE(30)' \
		-stat mean \
		-mask ${WD}/${session}/preproc/wm_in_func_res/wm_func_res_final.nii.gz \
		-use_nonmask ${WD}/${session}/preproc/localWMreg/${session}_preproc_smooth_${run}.nii.gz


		##do regression
		#the model at this point is Y = [ Motion + localWMreg + Ventricles ] X + Residual 
		# did not remove linear drift
		# bandpass filtering at 0.008 0.09
		3dTproject \
		-input ${WD}/${session}/preproc/localWMreg/${session}_preproc_smooth_${run}.nii.gz \
		-prefix ${WD}/${session}/preproc/localWMreg/${session}_preproc_smooth_localWMreg_${run}.nii.gz \
		-ort ${WD}/${session}/preproc/localWMreg/${session}_MP_${run}.1D \
		-ort ${WD}/${session}/preproc/localWMreg/${session}_csf_${run}.1D \
		-dsort ${WD}/${session}/preproc/localWMreg/${session}_WMeLOCAL_${run}.nii.gz \
		-automask \
		-passband 0.008 0.09

	done
done