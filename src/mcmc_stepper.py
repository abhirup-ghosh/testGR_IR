"""
A general purpose MCMC stepper

(C) Archisman Ghosh, Siddharth Mohite
"""

import copy
import numpy as np
from scipy import random

def inrange(num, left, right):
  return left < num < right

class parameters(object):
  def __init__(self, *args):
    self._dimension = np.size(args)
    self._names = list(args)
    for name in self._names:
      self.__dict__[name] = None
      self.__dict__['_'+name+'_sigma'] = None
      self.__dict__['_'+name+'_bounds'] = None
  def get_bounds(self, name):
    return self.__dict__['_'+name+'_bounds']
  def get_sigma(self, name):
    return self.__dict__['_'+name+'_sigma']
  def get_value(self, name):
    return self.__dict__[name]
  def set_bounds(self, name, bounds):
    self.__dict__['_'+name+'_bounds'] = bounds
  def set_sigma(self, name, sigma):
    self.__dict__['_'+name+'_sigma'] = sigma
  def set_value(self, name, value):
    self.__dict__[name] = value
  def init_params(self):
    for name in self._names:
      self.set_value(name, random.uniform(*self.get_bounds(name)))
  def in_bounds(self):
    for name in self._names:
      if not inrange(self.get_value(name), *self.get_bounds(name)):
        return False
    return True
  def try_params(self):
    params_trial = copy.copy(self)
    for name in self._names:
      params_trial.set_value(name, random.normal(self.get_value(name), self.get_sigma(name)))
    return params_trial
  def walk_params(self):
    while True:
      params_trial = self.try_params()
      if params_trial.in_bounds():
        return params_trial
  def values(self):
    return tuple(self.get_value(name) for name in self._names)
  def value_string(self, delimiter='\t', formatter='%f'):
    return delimiter.join([(formatter)%(self.get_value(name)) for name in self._names])
  def write_string(self, delimiter='\t', formatter='%.3f', assigner='='):
    return delimiter.join([('%s'+assigner+formatter)%(name, self.get_value(name)) for name in self._names])

class mcmc_stepper(object):
  def __init__(self, params, loglikelihood, data, beta=1.):
    self.params = params
    self.loglikelihood = loglikelihood
    self.data = data
    self.beta = beta
  def mc_step(self):
    #if self.params.logl is None:
      #self.params.logl = self.loglikelihood(self.params, self.data)
    params_new = self.params.walk_params()
    params_new.logl = self.loglikelihood(params_new, self.data)
    rr = random.random()
    if np.log(rr)<self.beta*(params_new.logl-self.params.logl):
      self.params = params_new
      return True
    return False
  def mc_run(self, N_cycle, N_step, ofile, verbose=True):
    self.params.logl = self.loglikelihood(self.params, self.data)
    for cycle in range(N_cycle):
      acc = 0.
      for step in range(N_step):
        acc += self.mc_step()/float(N_step)
      if verbose:
        print('Cycle: %d, %s, logl=%.3f, acc=%.3f'%(cycle, self.params.write_string(delimiter=', ', formatter='%f'), self.params.logl, acc))
      ofile.write('%d\t%s\t%f\n'%(cycle, self.params.value_string(), self.params.logl))
  