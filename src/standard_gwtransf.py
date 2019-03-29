"""
Module for basic GW parameter transformations.

CONVENTIONS:
===========

1. 0 < q <= 1.
2. m1 > m2 in the tuple (m1, m2). m2 = q*m1.

(c) Archisman Ghosh, 2014-10-22
"""


import numpy as np

def eta_from_q(q):
  return q/(1.+q)**2

def q_from_eta(eta):
  return ((1.-2.*eta)+np.sqrt(1.-4.*eta))/(2.*eta)

def q_from_comp(m1, m2):
  if m1<m2:
    return m1/float(m2)
  return m2/float(m1)

def eta_from_comp(m1, m2):
  return m1*m2/(m1+m2)**2.

def tot_from_comp(m1, m2):
  return m1+m2

def mc_from_comp(m1, m2):
  return ((m1*m2)**(3./5.))/((m1+m2)**(1./5.))

def mc_from_toteta(m, eta):
  return m*eta**(3./5.)

def tot_from_mceta(mc, eta):
  return mc*eta**(-3./5.)

def tot_from_mcq(mc, q):
  return tot_from_mceta(mc, eta_from_q(q))

def mcq_from_comp(m1, m2):
  return (mc_from_comp(m1, m2), q_from_comp(m1, m2))

def mceta_from_comp(m1, m2):
  return (mc_from_comp(m1, m2), eta_from_comp(m1, m2))

def comp_from_totq(m, q):
  m1 = m/(1.+q)
  m2 = q*m1
  return (m1, m2)

def comp_from_toteta(m, eta):
  q = q_from_eta(eta)
  return comp_from_totq(m, q)

def comp_from_mcq(mc, q):
  eta = eta_from_q(q)
  m = tot_from_mceta(mc, eta)
  return comp_from_totq(m, q)

def comp_from_mceta(mc, eta):
  m = tot_from_mceta(mc, eta)
  return comp_from_toteta(m, eta)

def mf_from_toteta(m, eta):
  return m*(1. + (np.sqrt(8./9.)-1.)*eta - 0.4333*(eta**2.) - (0.4392*(eta**3)))

def mf_from_comp(m1, m2):
  return mf_from_toteta(tot_from_comp(m1, m2), eta_from_comp(m1, m2))

def mf_from_totq(m, q):
  return mf_from_toteta(m, eta_from_q(q))

def sf_from_eta(eta):
  return np.sqrt(12.)*eta - 3.871*(eta**2.) + 4.028*(eta**3)

def sf_from_q(q):
  return sf_from_eta(eta_from_q(q))
