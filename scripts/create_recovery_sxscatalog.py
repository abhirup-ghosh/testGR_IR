import numpy as np
import os

data = np.genfromtxt("./SXS_campaign.dat", names=True, dtype=None)
sub_dir_list = ['IMR', 'inspiral', 'post-inspiral']

for idx in range(len(data)):
	print idx, data[idx]['tag']

	tag = data[idx]['tag']
	geocentric_end_time = data[idx]['geocentric_end_time']

	out_dir = '/home/abhirup/Documents/Work/testGR_IR/runs/systematics_error_characterisation/sxs_injection_o1o2_noise/%s'%tag
	for sub_dir in sub_dir_list:

		# creating output directory and sub-directories: IMR, inspiral, post-inspiral
		os.system('mkdir -p %s/IMR %s/inspiral %s/post-inspiral'%(out_dir, out_dir, out_dir))

		# copying gwf, psd and injection files to output sub-directories
	        os.system('cp %s* %s/%s'%(tag, out_dir, sub_dir))

		# creating a gpstime.txt file in each output sub-directory
	        os.system('echo %f > %s/%s/gpstime.txt'%(geocentric_end_time, out_dir, sub_dir))

		# changing nomenclature of gwf, psd and injection xml files
		os.system("cp %s/%s/%s_H-H1HWINJ.gwf %s/%s/H-H1HWINJ.gwf"%(out_dir, sub_dir, tag, out_dir, sub_dir))
		os.system("cp %s/%s/%s_L-L1HWINJ.gwf %s/%s/L-L1HWINJ.gwf"%(out_dir, sub_dir, tag, out_dir, sub_dir))
		os.system("cp %s/%s/%s_V-V1HWINJ.gwf %s/%s/V-V1HWINJ.gwf"%(out_dir, sub_dir, tag, out_dir, sub_dir))

		os.system("cp %s/%s/%s_H1_psd.txt %s/%s/H1_psd.txt"%(out_dir, sub_dir, tag, out_dir, sub_dir))
		os.system("cp %s/%s/%s_L1_psd.txt %s/%s/L1_psd.txt"%(out_dir, sub_dir, tag, out_dir, sub_dir))
		os.system("cp %s/%s/%s_V1_psd.txt %s/%s/V1_psd.txt"%(out_dir, sub_dir, tag, out_dir, sub_dir))

		os.system("cp %s/%s/%s_*.xml.gz %s/%s/injection.xml.gz"%(out_dir, sub_dir, tag, out_dir, sub_dir))

		# copying copying example_config.ini to each sub-directory
		os.system("cp /home/abhirup/Documents/Work/testGR_IR/runs/systematics_error_characterisation/sxs_injection_o1o2_noise/example_config.ini %s/%s"%(out_dir, sub_dir))

		# making appropriate changes inside example_config.ini
		os.system("sed -i 's/SXS_BBH_0150o2HL_1170730903\/IMR/%s\/%s/g' %s/%s/example_config.ini"%(tag, sub_dir, out_dir, sub_dir))

		# change detector configuration inside config.ini
		if "HLV" in tag:
			os.system("""sed -i "s/'H1','L1'/'H1','L1','V1'/" %s/%s/example_config.ini"""%(out_dir, sub_dir))

		# insert f_isco
		f_isco = data[idx]['f_isco']
		if sub_dir == 'inspiral':
			os.system("""sed -i "/^flow/a fhigh=\{'H1':%.2f,'L1':%.2f,'V1':%.2f\}" %s/%s/example_config.ini"""%(f_isco, f_isco, f_isco, out_dir, sub_dir))
		elif sub_dir == 'post-inspiral':
			os.system("""sed -i "s/^flow.*/flow=\{'H1':%.2f,'L1':%.2f,'V1':%.2f\}/g" %s/%s/example_config.ini"""%(f_isco, f_isco, f_isco, out_dir, sub_dir))

		# launch dag
		os.system("lalinference_pipe -g %s/%s/gpstime.txt -r %s/%s -p %s/%s %s/%s/example_config.ini --condor-submit"%(out_dir, sub_dir, out_dir, sub_dir, out_dir, sub_dir, out_dir, sub_dir))

