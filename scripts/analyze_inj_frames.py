#!/usr/bin/env python

import sys, os, commands, glob, re, numpy as np
from pylal import Fr
import matplotlib.pyplot as plt
from mynewwaveforms import fd_from_td
import LIGOpsd as psds

start_time = 501
stop_time = 503

frame_folder = '/home/archis/Work/WaveformSystematics/frames'
out_folder = '/home/archis/Work/WaveformSystematics/plot_frames'
os.system('mkdir -p %s'%(out_folder))

#frame_dir = 'NR_0002_HM'
#frame_loc = '/home/archis/Work/WaveformSystematics/frames/%s'%(frame_dir)

for frame_loc in glob.glob(frame_folder+'/*'):
  if os.path.isdir(frame_loc):
    frame_dir = os.path.basename(frame_loc)
    print frame_dir
    
    try:

      if frame_dir[-3:] == '_HM':
        H1_nonoise_frame = glob.glob(os.path.join(frame_loc, 'H1*_HM.gwf'))[0]
        L1_nonoise_frame = glob.glob(os.path.join(frame_loc, 'L1*_HM.gwf'))[0]
      else:
        H1_nonoise_frame = glob.glob(os.path.join(frame_loc, 'H1*_l2m2.gwf'))[0]
        L1_nonoise_frame = glob.glob(os.path.join(frame_loc, 'L1*_l2m2.gwf'))[0]

      H1_C00_frame = glob.glob(os.path.join(frame_loc, 'H1*_C00.gwf'))[0]
      L1_C00_frame = glob.glob(os.path.join(frame_loc, 'L1*_C00.gwf'))[0]

      # Load data
      for prefix in ['H1_nonoise', 'L1_nonoise', 'H1_C00', 'L1_C00']:
        frame = vars()[prefix+'_frame']
        frame_dump = commands.getoutput('FrDump -i %s'%(frame))
        channel = re.findall('ProcData:.*', frame_dump)[0][10:]
        (vars()[prefix+'_frdata'], vars()[prefix+'_gps_start'], vars()[prefix+'_xoffset'], vars()[prefix+'_xspacing'], vars()[prefix+'_xunit'], vars()[prefix+'_yunit']) = Fr.frgetvect(frame, channel)
        vars()[prefix+'_time'] = np.arange(0, len(vars()[prefix+'_frdata']))*(vars()[prefix+'_xspacing'])

      # Do FFTs
      for prefix in ['H1_nonoise', 'L1_nonoise', 'H1_C00', 'L1_C00']:
        start_idx = np.abs(vars()[prefix+'_time']-start_time).argmin()
        stop_idx = np.abs(vars()[prefix+'_time']-stop_time).argmin()
        (vars()[prefix+'_freq'], vars()[prefix+'_frdata_fd']) = fd_from_td(vars()[prefix+'_time'], vars()[prefix+'_frdata'])
        (vars()[prefix+'_freq_trunc'], vars()[prefix+'_frdata_fd_trunc']) = fd_from_td(vars()[prefix+'_time'][start_idx:stop_idx], vars()[prefix+'_frdata'][start_idx:stop_idx])

      # Calculate PSDs
      f_low = 20.
      noise_model = 'aLIGO_EARLY_HIGH'
      for prefix in ['H1_nonoise', 'L1_nonoise']:
        vars()[prefix+'_snr'] = psds.SNR(vars()[prefix+'_freq'], vars()[prefix+'_frdata_fd'], noise_model, f_low)
        vars()[prefix+'_snr_trunc'] = psds.SNR(vars()[prefix+'_freq_trunc'], vars()[prefix+'_frdata_fd_trunc'], noise_model, f_low)
        vars()[prefix+'_snr_C00'] = psds.SNR(vars()[prefix+'_freq'], vars()[prefix+'_frdata_fd'], [vars()[prefix.replace('nonoise', 'C00')+'_freq'], np.abs(vars()[prefix.replace('nonoise', 'C00')+'_frdata_fd'])], f_low)
        vars()[prefix+'_snr_trunc_C00'] = psds.SNR(vars()[prefix+'_freq_trunc'], vars()[prefix+'_frdata_fd_trunc'], [vars()[prefix.replace('nonoise', 'C00')+'_freq'], np.abs(vars()[prefix.replace('nonoise', 'C00')+'_frdata_fd'])], f_low)
        print("With %s f_low = %d Hz, %s SNR = %.2f (full chunk), %.2f (2s of data)."%(noise_model, f_low, prefix[:2], vars()[prefix+'_snr'], vars()[prefix+'_snr_trunc']))
        print("With C00 f_low = %d Hz, %s SNR = %.2f (full chunk), %.2f (2s of data)."%(f_low, prefix[:2], vars()[prefix+'_snr_C00'], vars()[prefix+'_snr_trunc_C00']))
      ofile = open('%s/%s_SNR.txt'%(out_folder, frame_dir), 'w')
      ofile.write('#Noise\tH1_SNR\tL1_SNR\tNetwork_SNR\n')
      ofile.write('%s\t%.2f\t%.2f\t%.2f\n'%(noise_model, H1_nonoise_snr, L1_nonoise_snr, np.sqrt(H1_nonoise_snr**2 + L1_nonoise_snr**2)))
      ofile.write('%s_C00\t%.2f\t%.2f\t%.2f\n'%(frame_dir, H1_nonoise_snr_C00, L1_nonoise_snr_C00, np.sqrt(H1_nonoise_snr_C00**2 + L1_nonoise_snr_C00**2)))
      ofile.close()

      ff = np.linspace(1e1, 1e4, 1e3)
      aleh_psd = psds.LIGOPSD(ff, noise_model)

      # TD plots
      plt.figure(figsize=(9,12))
      ax = plt.subplot2grid((2,1),(0,0))
      ax.xaxis.tick_top()
      ax.xaxis.set_label_position('top') 
      plt.xlabel(r'Time (s) after start of frame')
      for prefix in ['H1_nonoise', 'L1_nonoise']:
        plt.plot(vars()[prefix+'_time'], vars()[prefix+'_frdata'], label=prefix)
      plt.xlim(start_time, stop_time)
      plt.legend()
      plt.ylabel(r'Strain')
      plt.subplot2grid((2,1),(1,0))
      for prefix in ['H1_nonoise', 'L1_nonoise', 'H1_C00', 'L1_C00']:
        plt.plot(vars()[prefix+'_time'], vars()[prefix+'_frdata'], label=prefix)
      plt.xlim(start_time, stop_time)
      plt.legend()
      plt.xlabel(r'Time (s) after start of frame')
      plt.ylabel(r'Strain')
      plt.tight_layout()
      plt.savefig('%s/%s_TD.png'%(out_folder, frame_dir))

      # FD plots
      plt.figure(figsize=(9,12))
      ax = plt.subplot2grid((2,1),(0,0))
      ax.xaxis.tick_top()
      ax.xaxis.set_label_position('top') 
      plt.xlabel(r'Freq (Hz)')
      for prefix in ['H1_nonoise', 'L1_nonoise', 'H1_C00', 'L1_C00']:
        plt.loglog(vars()[prefix+'_freq'], np.abs(vars()[prefix+'_frdata_fd']), label=prefix)
      plt.loglog(ff, np.sqrt(aleh_psd), 'k--', linewidth=2, label='AdvLIGO Early High')
      plt.xlim(1e1, 1e4)
      plt.legend(loc='lower left', title='full 512s')
      plt.ylabel(r'Strain (Hz$^{-1}$)')
      plt.subplot2grid((2,1),(1,0))
      for prefix in ['H1_nonoise', 'L1_nonoise', 'H1_C00', 'L1_C00']:
        plt.loglog(vars()[prefix+'_freq_trunc'], np.abs(vars()[prefix+'_frdata_fd_trunc']), label=prefix)
      plt.loglog(ff, np.sqrt(aleh_psd), 'k--', linewidth=2, label='AdvLIGO Early High')
      plt.xlim(1e1, 1e4)
      plt.legend(loc='lower left', title='2s of data')
      plt.xlabel(r'Freq (Hz)')
      plt.ylabel(r'Strain (Hz$^{-1}$)')
      plt.tight_layout()
      plt.savefig('%s/%s_FD.png'%(out_folder, frame_dir))

      #plt.show()
    
    except Exception,e:
      print 'Failed for reason: %s'%(str(e))
