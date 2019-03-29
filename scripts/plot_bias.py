import numpy as np
import matplotlib.pyplot as plt
"""
# Change according to the total mass and mass ratio required 

i = 1 # index for total mass 
j = 0 # index for assymetric mass ratio


n_isco = [0.125, 0.25, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 4.0]

if i == 0:
	nisco = n_isco
elif i == 1:
	nisco = n_isco[1:]
elif i == 2:
	nisco = n_isco[2:]

date = '2015-01-12'

m = [30., 50., 100.]
q = [1., 2., 4.]

m_1 = np.zeros((len(m), len(q)))
m_2 = np.zeros((len(m), len(q)))
eta_0 = np.zeros((len(m), len(q)))
m_c = np.zeros((len(m), len(q)))

for k in range(len(m)):
	for l in range(len(q)):
		m_1[k,l] = m[k]/(1. + q[l])
		m_2[k,l] = (m[k]*q[l])/(1. + q[l])
		eta_0[k,l] = (m_1[k,l]*m_2[k,l])/((m_1[k,l]+m_2[k,l])**2.)
		m_c[k,l] = ((m_1[k,l]*m_2[k,l])**(3./5.))/((m_1[k,l]+m_2[k,l])**(1./5.))


delta = np.zeros((len(nisco),5))
sigma = np.zeros((len(nisco),5))

for k in range(len(nisco)):
	address = '../runs/%s/%s_%s/nisco_%s/cbc_bayes/posterior_samples.dat'%(date, str(m[i]), str(q[j]), str(nisco[k]))
	M, m1, m2, mc, eta = np.loadtxt(address, usecols = (5, 13, 16, 23, 32), skiprows = 1, unpack = True)
	M_mean = np.mean(M)
	delta[k,0] = m[i] - M_mean
	sigma[k,0] = np.std(M)
	m1_mean = np.mean(m1)
	delta[k,1] = m_1[i,j] - m1_mean
	sigma[k,1] = np.std(m1)
	m2_mean = np.mean(m2)
	delta[k,2] = m_2[i,j] - m2_mean
	sigma[k,2] = np.std(m2)
	mc_mean = np.mean(mc)
	delta[k,3] = m_c[i,j] - mc_mean
	sigma[k,3] = np.std(mc)
	eta_mean = np.mean(eta)	
	delta[k,4] = eta_0[i,j] - eta_mean
	sigma[k,4] = np.std(eta)

xlabel = ['M', '$m_1$', '$m_2$', '$M_c$', '$\eta$']
	
plt.figure(figsize=(16,9))
for k in range(5):
	plt.subplot(2,3,k+1)
	plt.semilogx(nisco, delta[:,k], 'r', marker = '.', label='$\Delta$: systematic error')
	plt.semilogx(nisco, sigma[:,k], 'g', marker = '.', label='$\sigma$: statistical error')
	plt.legend(loc='best')
	plt.xlabel('$f_{high}/f_{isco}$')
	plt.ylabel('%s'%(xlabel[k]))
	#plt.title('$M$='+str(M_0)+'; q='+str(q))
	#plt.title('$m_1$=15.; $m_2$=15.; \n Injection: EOBNRv2; Recovery: TaylorF2(no noise); \n $f_{isco}$=146.57')
	plt.grid(True, which='both', ls='-')
plt.suptitle('M='+str(m[i])+'; q='+ str(q[j]))
#plt.savefig('../runs/%s/%s_%s/bias.png'%(date, str(M_0), str(q)))
plt.tight_layout()
plt.show()
"""

date = '2015-01-28'

m = [30., 50., 100.]
q = [1., 2., 4.]
nisco = [0.125, 0.25, 0.5, 1.0]

delta = np.zeros(len(nisco))
sigma = np.zeros(len(nisco))


plt.figure(figsize=(16,8))
for i in range(1):
	for j in range(1):
		for k in range(len(nisco)):	
			address = '../runs/%s/%s_%s/nisco_%s/cbc_bayes/posterior_samples.dat'%(date, str(m[i]), str(q[j]), str(nisco[k]))
			M, m1, m2, mc, eta = np.loadtxt(address, usecols = (5, 13, 16, 23, 32), skiprows = 1, unpack = True)
			m_1 = m[i]/(1. + q[j])
			m_2 = (m[i]*q[j])/(1. + q[j])
			m_c = ((m_1*m_2)**(3./5.))/((m_1+m_2)**(1./5.))
			M_mean = np.mean(M)
			delta[k] = abs(m[i] - M_mean)
			sigma[k] = np.std(M)
		plt.subplot(3,3,3*i+j+1)
		plt.plot(nisco, delta, 'r', marker = '.', label='$\Delta$: systematic error')
		plt.plot(nisco, sigma, 'g', marker = '.', label='$\sigma$: statistical error')
		plt.legend(loc='best')
		plt.xlabel('$f_{high}/f_{isco}$')
		plt.ylabel('$M$ ($M_\odot$)')
		#plt.title('$M$='+str(M_0)+'; q='+str(q))
		#plt.title('$m_1$=15.; $m_2$=15.; \n Injection: EOBNRv2; Recovery: TaylorF2(no noise); \n $f_{isco}$=146.57')
		plt.grid(True, which='both', ls='-')
plt.tight_layout()
plt.show()
