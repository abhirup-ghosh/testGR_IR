{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "metadata": {},
   "outputs": [],
   "source": [
    "import os, sys, commands, numpy as np"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "metadata": {},
   "outputs": [],
   "source": [
    "import matplotlib.pyplot as plt"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "metadata": {},
   "outputs": [],
   "source": [
    "import standard_gwtransf as gw"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 20,
   "metadata": {},
   "outputs": [],
   "source": [
    "m1, m2 = np.random.uniform(2.43,182.59,100000), np.random.uniform(2.43,182.59,100000)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 21,
   "metadata": {},
   "outputs": [],
   "source": [
    "mc = gw.mc_from_comp(m1,m2)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 25,
   "metadata": {
    "scrolled": false
   },
   "outputs": [
    {
     "data": {
      "image/png": "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\n",
      "text/plain": [
       "<matplotlib.figure.Figure at 0x7fc9cd68e1d0>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "plt.hist(mc, bins=50, histtype='step', normed=True); plt.xlabel('$M_c$'); plt.ylabel('P($M_c$)');"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 2",
   "language": "python",
   "name": "python2"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 2
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython2",
   "version": "2.7.9"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
