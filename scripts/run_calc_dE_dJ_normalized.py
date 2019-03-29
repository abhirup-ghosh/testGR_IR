import numpy as np
import os

M_inj_list = [30., 50., 100., 150., 200.]
q_inj_list = [1., 2., 4.]
spin_fit_formula = 'nospin_Pan2011'
#spin_fit_formula = 'nonprecspin_Healy2014'
out_dir_root = '../runs/2015-06-05/imrtestgr_nonspinprior_%s' %spin_fit_formula
run_nohup = True 
comp_spin_min = 0. 
comp_spin_max = 0. 
post_loc = '/home/abhirup/Documents/Work/testGR_IR/runs/2015-06-05'

for M_inj in M_inj_list:
  for q_inj in q_inj_list:
		m1_inj = M_inj/(1.+q_inj)
		m2_inj = M_inj*q_inj/(1.+q_inj)
		
		print '=== M_inj = %2.1f q = %2.1f m1 = %2.1f m2 = %2.1f ===' %(M_inj, q_inj, m1_inj, m2_inj)

		post_samples_i = '%s/inspiral/%s_%s/nisco_1.0/cbc_bayes/posterior_samples.dat' %(post_loc, M_inj, q_inj)
		post_samples_r = '%s/ringdown_more_samples/%s_%s/nisco_1.0/cbc_bayes/posterior_samples.dat' %(post_loc, M_inj, q_inj)
		post_samples_imr = '%s/IMR/%s_%s/nisco_0.0/cbc_bayes/posterior_samples.dat' %(post_loc, M_inj, q_inj)
		comp_mass_min = 0.25*m1_inj 
		comp_mass_max = 2*m2_inj 
		waveform = 'IMRPhenomB'
		
		insp_fhigh = (1./6.)**(3./2)/(LAL_PI*M_inj*LAL_MTSUN_SI)
		if q_inj == 1.:
        		f_qnm = 0.0879
		elif q_inj == 2.:
        		f_qnm = 0.0829
		elif q_inj == 4.:
        		f_qnm = 0.0744
		ring_flow = f_qnm/(M_inj*LAL_MTSUN_SI)


		out_dir = '%s/M_inj%3.2f_q_inj%3.2f' %(out_dir_root, M_inj, q_inj) 
		#os.system('mkdir -p %s' %out_dir)
		log_file = '%s/log.txt' %out_dir
		err_file = '%s/err.txt' %out_dir

		os.system('mkdir -p %s' %out_dir)

		if run_nohup == True: 
			run_cmd = 'nohup python imrconsistencytest.py --insp-post=%s --ring-post=%s --imr-post=%s --fit-formula=%s --out-dir=%s --insp_fhigh=%.2f --ring_flow=%.2f --mtot_inj=%f --mratio_inj=%f --comp-mass-min=%f --comp-mass-max=%f --comp-spin-min=%f --comp-spin-max=%f --waveform=%s > %s 2> %s &' %(post_samples_i, post_samples_r, post_samples_imr, spin_fit_formula, out_dir, insp_fhigh, ring_flow, M_inj, q_inj, comp_mass_min, comp_mass_max, comp_spin_min, comp_spin_max, waveform, log_file, err_file)
		else: 
			run_cmd = 'python imrconsistencytest.py --insp-post=%s --ring-post=%s --imr-post=%s --fit-formula=%s --out-dir=%s --insp_fhigh=%.2f --ring_flow=%.2f --mtot_inj=%f --mratio_inj=%f --comp-mass-min=%f --comp-mass-max=%f --comp-spin-min=%f --comp-spin-max=%f --waveform=%s' %(post_samples_i, post_samples_r, post_samples_imr, spin_fit_formula, out_dir, insp_fhigh, ring_flow, M_inj, q_inj, comp_mass_min, comp_mass_max, comp_spin_min, comp_spin_max, waveform)
		
		#print run_cmd 
		os.system(run_cmd)  
