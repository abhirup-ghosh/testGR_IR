""" 
Module for computation of confidence intervals, sky areas, etc.

(c) Archisman Ghosh, 2015-11-17
"""

import numpy as np
from scipy import interpolate

# Simple functions to translate between bin edges, bin centers, bin widths
def width(edg_or_cen):
  return np.average(np.diff(edg_or_cen))

def cen_from_edg(edges):
  return (edges[1:]+edges[:-1])/2.

def edg_from_cen(centers):
  bin_width = width(centers)
  return  np.array([centers[0]-bin_width/2.]+list((centers[1:]+centers[:-1])/2.)+[centers[-1]+bin_width/2.])

# Basic module for confidence calculations
class confidence(object):
  def __init__(self, counts):
    self.counts = counts
    # Sort in descending order in frequency
    self.counts_sorted = np.sort(counts.flatten())[::-1]
    # Get a normalized cumulative distribution from the mode
    self.norm_cumsum_counts_sorted = np.cumsum(self.counts_sorted) / np.float(np.sum(counts))
    # Set interpolations between heights, bins and levels
    self._set_interp()
  def _set_interp(self):
    self._length = len(self.counts_sorted)
    # height from index
    self._height_from_idx = interpolate.interp1d(np.arange(self._length), self.counts_sorted, bounds_error=False, fill_value=0.)
    # index from height
    self._idx_from_height = lambda ht: interpolate.interp1d(self.counts_sorted[::-1], np.arange(self._length)[::-1], bounds_error=False, fill_value=0.)(ht) if ht>0. else self._length
    # level from index
    self._level_from_idx = interpolate.interp1d(np.arange(self._length), self.norm_cumsum_counts_sorted, bounds_error=False, fill_value=1.)
    # index from level
    self._idx_from_level = lambda lev: interpolate.interp1d(self.norm_cumsum_counts_sorted, np.arange(self._length), bounds_error=False, fill_value=0.)(lev) if lev<1. else self._length 
  def level_from_height(self, height):
    """
    Returns value for a given confidence level.
    """
    return self._level_from_idx(self._idx_from_height(height))
  def height_from_level(self, level):
    """
    Returns confidence level for a given value
    """
    return self._height_from_idx(self._idx_from_level(level))
  def _idx_above(self, level):
    """
    Returns indices of values above a given confidence level.
    """
    return np.where(self.counts>self.height_from_level(level))

# 1d confidence
class confidence_1d(confidence):
  def __init__(self, counts, edges):
    confidence.__init__(self, counts)
    self.centers = cen_from_edg(edges)
    self.width = width(edges)
  def level_at(self, X):
    if not '_height_interp' in dir(self):
      self._height_interp = interpolate.interp1d(self.centers, self.counts, bounds_error=False, fill_value=0.)
    return self.level_from_height(self._height_interp(X))
  def edges(self, level=0.68):
    """
    Returns edges of the confidence interval (assumes connected).
    """
    valid_x = self.centers[self._idx_above(level)]
    return (min(valid_x)-self.width/2., max(valid_x)+self.width/2.)
  def range(self, level=0.68):
    """
    Returns size of 1d confidence region.
    """
    return np.size(self._idx_above(level))*self.width

# 2d confidence
class confidence_2d(confidence):
  def __init__(self, counts, xedges, yedges):
    confidence.__init__(self, counts)
    self.xcenters = cen_from_edg(xedges)
    self.ycenters = cen_from_edg(yedges)
  def level_at(self, X, Y):
    if not '_height_interp' in dir(self):
      self._height_interp = interpolate.interp2d(self.xcenters, self.ycenters, self.counts, bounds_error=False, fill_value=0.)
    return self.level_from_height(self._height_interp(X, Y)[0])
  def area(self, level=0.68, xweights=None, yweights=None):
    """
    Returns area of 2d confidence region.
    """
    idx = self._idx_above(level)
    if xweights is None:
      xweights = self.xcenters*0.+width(self.xcenters)
    elif callable(xweights):
      xweights = np.vectorize(xweights)(self.xcenters)*width(self.xcenters)
    if yweights is None:
      yweights = self.ycenters*0.+width(self.ycenters)
    elif callable(yweights):
      yweights = np.vectorize(yweights)(self.ycenters)*width(self.ycenters)
    (XX, YY) = np.meshgrid(self.xcenters, self.ycenters, indexing='ij')
    (Xw, Yw) = np.meshgrid(xweights, yweights, indexing='ij')
    return sum(Xw[idx]*Yw[idx])
  def sky_area(self, level=0.68):
    """
    Returns sky area. Arguments have to be in the order (ra, dec).
    """
    return self.area(level=level, yweights=np.cos)
  
# 3d confidence
class confidence_3d(confidence):
  def __init__(self, counts, xedges, yedges, zedges):
    confidence.__init__(self, counts)
    self.xcenters = cen_from_edg(xedges)
    self.ycenters = cen_from_edg(yedges)
    self.zcenters = cen_from_edg(zedges)
  def volume(self, level=0.68, xweights=None, yweights=None, zweights=None):
    """
    Returns volume of 3d confidence region.
    """
    idx = self._idx_above(level)
    if xweights is None:
      xweights = self.xcenters*0.+width(self.xcenters)
    elif callable(xweights):
      xweights = np.vectorize(xweights)(self.xcenters)*width(self.xcenters)
    if yweights is None:
      yweights = self.ycenters*0.+width(self.ycenters)
    elif callable(yweights):
      yweights = np.vectorize(yweights)(self.ycenters)*width(self.ycenters)
    if zweights is None:
      zweights = self.zcenters*0.+width(self.zcenters)
    elif callable(zweights):
      zweights = np.vectorize(zweights)(self.zcenters)*width(self.zcenters)
    (XX, YY, ZZ) = np.meshgrid(self.xcenters, self.ycenters, self.zcenters, indexing='ij')
    (Xw, Yw, Zw) = np.meshgrid(xweights, yweights, zweights, indexing='ij')
    return sum(Xw[idx]*Yw[idx]*Zw[idx])
  def sky_volume(self, level=0.68):
    """
    Returns sky volume. Arguments have to be in the order (ra, dec, dist).
    """
    return self.volume(level=level, yweights=np.cos, zweights=(lambda z: z**2))
