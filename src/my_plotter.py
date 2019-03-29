"""
This header defines some pretty plotting routines.
"""

import numpy as np
import bayesian as ba
import matplotlib.pyplot as plt
from matplotlib.pyplot import *

def plot1d(data, centers=None, binned=False, bins=20, range=None, injection=None, instance=None, limits=True, y_lim=None, xlabel=None, significance=False, sbins=None, n_smooth=40, histtype='kde', orientation='vertical', color='g', icolor='r', scolors=['m', 'c', 'y'], linestyle='-', ilinestyle='-', slinestyle='--', linewidth=2, ilinewidth=1, slinewidth=1, **kwargs):
  hn = [0., 0.]
  if sbins is None:
    sbins = bins
  if type(data) is ba.prob_distr:
    pdf = data
    binned = True
  else:
    pdf = ba.prob_distr(data, centers, binned=binned, bins=sbins, range=range)
  (pos, cen) = pdf.get_prob()
  if significance or histtype in ['smooth', 'fsmooth', 'kde', 'fkde']:
    spdf = pdf.smooth(n_smooth=n_smooth)
    (spos, scen) = pdf.smooth(n_smooth=n_smooth).get_prob()
    if histtype in ['kde', 'fkde']:
      spos = pdf.kde(scen)
  if significance:
    spline1 = pdf.spline()
    r1 = spdf.nsigma_range(1); print '1-sigma range =', r1
    r2 = spdf.nsigma_range(2); print '2-sigma range =', r2
    r3 = spdf.nsigma_range(3); print '3-sigma range =', r3
  if instance is None:
    instance = plt
    instance.figure()
    instance.clf()
  if histtype in ['smooth', 'kde']:
    if orientation=='vertical':
      instance.plot(scen, spos, color=color, linestyle=linestyle, linewidth=linewidth, **kwargs)
    else:
      instance.plot(spos, scen, color=color, linestyle=linestyle, linewidth=linewidth, **kwargs)
  elif histtype in ['fsmooth', 'fkde']:
    if orientation=='vertical':
      #instance.plot(scen, spos, color=color, linestyle=linestyle, linewidth=linewidth, **kwargs)
      instance.fill_between(scen, 0, spos, color=color, linestyle=linestyle, linewidth=linewidth, **kwargs)
    else:
      #instance.plot(spos, scen, color=color, linestyle=linestyle, linewidth=linewidth, **kwargs)
      instance.fill_betweenx(scen, 0, spos, color=color, linestyle=linestyle, linewidth=linewidth, **kwargs)
  else:
    (hn, hb, hp) = instance.hist(data, range=range, bins=bins, normed=True, histtype=histtype, orientation=orientation, color=color, **kwargs)
  #print 'hn =', hn
  #print 'max(hn) =', max(hn)
  axoline = eval('plt.ax%sline'%(orientation[0]))
  oxlim = eval('plt.%slim'%((orientation=='vertical')*'x'+(orientation=='horizontal')*'y'))
  oylim = eval('plt.%slim'%((orientation=='vertical')*'y'+(orientation=='horizontal')*'x'))
  oxlabel = eval('plt.%slabel'%((orientation=='vertical')*'x'+(orientation=='horizontal')*'y'))
  if y_lim is None:
    y_lim = 1.1*max(max(pos), max(hn));
  if limits:
    oxlim(min(cen),max(cen))
    oylim(0,y_lim*1.1)
  if range is not None:
    oxlim(*range)
  if significance:
    axoline(r3[0], 0, 1, linestyle=slinestyle, linewidth=slinewidth, color=scolors[2])
    axoline(r3[1], 0, 1, linestyle=slinestyle, linewidth=slinewidth, color=scolors[2])
    axoline(r2[0], 0, 1, linestyle=slinestyle, linewidth=slinewidth, color=scolors[1])
    axoline(r2[1], 0, 1, linestyle=slinestyle, linewidth=slinewidth, color=scolors[1])
    axoline(r1[0], 0, 1, linestyle=slinestyle, linewidth=slinewidth, color=scolors[0])
    axoline(r1[1], 0, 1, linestyle=slinestyle, linewidth=slinewidth, color=scolors[0])
  if injection is not None:
    axoline(injection, linestyle=ilinestyle, linewidth=ilinewidth, color=icolor)
  if xlabel is not None:
    instance.oxlabel(xlabel)
  #if xlim is not None:
    #oxlim(xlim[0],xlim[1])
  if significance:
    return (np.average(cen, weights=pos), r1[0], r1[1])
  else:
    return (spos, scen)

def plot2d(data_x, data_y=None, ycenters=None, binned=False, bins=100, range=None, sbins=20, injection=None, instance=None, significance=False, limits=True, labels=None, histtype='hist2d', color='g', cmap='BuGn', icolor='r', scolors=['m', 'c', 'y'], ilinestyle='-', slinestyle='--', ilinewidth=1, slinewidth=1, **kwargs):
  if type(data_x) is ba.prob_distr_2d:
    pdf2 = data_x
    binned = True
  else:
    pdf2 = ba.prob_distr_2d(data_x, data_y, ycenters, binned=binned, bins=sbins, range=range)
  (pos2, xcen, ycen) = pdf2.get_prob()
  if range is None:
    xmin = min(xcen)
    xmax = max(xcen)
    ymin = min(ycen)
    ymax = max(ycen)
  else:
    xmin = range[0][0]
    xmax = range[0][1]
    ymin = range[1][0]
    ymax = range[1][1]
  if significance:
    spline2 = pdf2.spline()
    s1 = pdf2.nsigma_value(1); print 's1 =', s1
    s2 = pdf2.nsigma_value(2); print 's2 =', s2
    s3 = pdf2.nsigma_value(3); print 's3 =', s3
  if instance is None:
    instance = plt
    instance.figure(figsize=(6,6))
    instance.clf()
  if binned or histtype=='imshow':
    instance.imshow(pos2.T, origin='lower', aspect='auto', extent=(xmin, xmax, ymin, ymax), cmap=cmap, **kwargs)
  elif histtype=='scatter':
      plt.scatter(data_x, data_y, color=color, **kwargs)
  elif histtype=='hist2d':
      plt.hist2d(data_x, data_y, range=range, bins=bins, cmap=cmap, **kwargs)
  if limits:
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
  if significance:
    (XX, YY) = np.meshgrid(xcen, ycen, indexing='ij')
    plt.contour(XX, YY, pos2, levels=[s1,s2,s3], colors=scolors, linewidths=slinewidth, linestyles=slinestyle)
  if injection is not None:
    plt.axvline(x=injection[0], linestyle=ilinestyle, linewidth=ilinewidth, color=icolor)
    plt.axhline(y=injection[1], linestyle=ilinestyle, linewidth=ilinewidth, color=icolor)
  if labels:
    plt.xlabel(legend[0])
    plt.ylabel(legend[1])
  return (pos2, xcen, ycen)

def plttext(x, y, text, instance=None, **kwargs):
  if instance is None:
    instance = plt
  instance.text(plt.axis()[0]*(1-x)+plt.axis()[1]*x, plt.axis()[-2]*(1-y)+plt.axis()[-1]*y, text, **kwargs)
