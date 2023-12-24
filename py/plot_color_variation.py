import sys
import warnings
from scipy.signal import argrelextrema
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import binned_statistic
from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as ticker
from scipy.optimize import curve_fit
from scipy import interpolate
from scipy import stats

import h5py
warnings.filterwarnings("ignore")
# matplotlib.use('agg')
plt.rcParams.update({
    'font.size': 12, 
    'font.family': 'Times',#sans-serif
    'axes.linewidth': 0.5,
    'axes.spines.bottom': True,
    'axes.spines.left': True,
    'axes.spines.right': True,
    'axes.spines.top': True,
    'xtick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
    'xtick.bottom': True,
    'ytick.left': True,
    'ytick.direction': 'in',
    'xtick.major.size': 3,
    'xtick.major.width': 0.8,
    'ytick.major.size': 3,
    'ytick.major.width': 0.8,
    'xtick.minor.size': 1.5,
    'xtick.minor.width': 0.6,
    'ytick.minor.size': 1.5,
    'ytick.minor.width': 0.6,
    "text.usetex": True,
    "axes.titlesize":"medium",
    "figure.dpi" : 500
})
def geterror(x,perclim=84.1):
    """
    x: tau list.
    return:
        tau, tau_lerr, tau_uperr
    """
    tau = stats.scoreatpercentile(x, 50)
    centau_uperr = (stats.scoreatpercentile(x, perclim))-tau
    centau_loerr = tau-(stats.scoreatpercentile(x, (100.-perclim)))
    return [round(tau,3), round(centau_loerr,3), round(centau_uperr,3)]
    
# f = h5py.File('./run_result/run_bin_103.68_montano22_n1_nsim_1000_mode_both_nsigma_3_band_i_z_sun14.h5', 'r') 
f = h5py.File('./run_i_z.h5', 'r')

fit = f['fit']
band1, band2 = f['band'].value.flatten()

num_all = f['num_all'].value
num_cut = f['num_cut'].value
num_pos = f['num_pos'].value

err1, err2 = f['par'].value[0][6], f['par'].value[1][6]

t1_min, t2_min = f['t_min']
t1_min, t2_min = 10 ** t1_min, 10 ** t2_min

sf1_min, sf2_min = f['sf_min']

t1_max, t2_max = f['t_max']
t1_max, t2_max = 10 ** t1_max, 10 ** t2_max

t_used_min = max(t1_min, t2_min)
t_used_max = min(t1_max, t2_max)

flux_ratio = f['flux_ratio'].value.flatten()


sf = f['sf'].value
bin_sf_1_t = 10 ** sf[0][0]
bin_sf_1_t_err = 10 ** sf[0][1]
bin_sf_1_sf = sf[0][2]
bin_sf_1_sf_err = sf[0][3]
# bin_sf_1_sf_err_hig = sf[0][4]

idx_br_1 = np.where(bin_sf_1_t <= t1_max, True, False)

bin_sf_2_t = 10 ** sf[1][0]
bin_sf_2_t_err = 10 ** sf[1][1]
bin_sf_2_sf = sf[1][2]
bin_sf_2_sf_err = sf[1][3]
# bin_sf_2_sf_err_hig = sf[1][4]

idx_br_2 = np.where(bin_sf_2_t <= t2_max, True, False)


# flux-flux

cv_1_t = 10 ** f['cv'][0][0]

cv_1_t_err = 10 ** f['cv'][0][1]
# cv_1_t_err = [(cv_1_t[i+1] - cv_1_t[i]) / 2 for i in range(len(cv_1_t)- 1)]
cv_1_y = f['cv'][0][2]
cv_1_y_err = f['cv'][0][3]

idx_cv_1 = np.where((cv_1_t>=t_used_min) & (cv_1_t<=t_used_max), True, False)


# mag-mag

cv_2_t = 10 ** f['cv'][1][0]
cv_2_t_err = 10 ** f['cv'][1][1]
cv_2_y = f['cv'][1][2]
cv_2_y_err = f['cv'][1][3]

idx_cv_2 = np.where((cv_2_t>=t_used_min) & (cv_2_t<=t_used_max), True, False)



cv_use_bin= 10 ** np.arange(0., 5.1, 0.1)
# cv_use_bin= 10 ** np.arange(0., 5.2, 0.2)
bin_width_all = [(cv_use_bin[i+1] - cv_use_bin[i]) / 2 for i in range(len(cv_use_bin)- 1)]
xerr_all = np.vstack((bin_width_all, bin_width_all))
bin_width_all = np.array(bin_width_all)
bin_centers_all = cv_use_bin[1:] - bin_width_all


"""
TODO ERROR.
"""

# For plotting -------------------------------------------------------------->
mosaic = ['A', 'B']
figure_mosaic = """
AC
BD
"""
colors = ['xkcd:mid blue', 'darkred']

fig, axes = plt.subplot_mosaic(mosaic = figure_mosaic, figsize=(4.2, 3.5))


axes['A'].axhline(np.sqrt(2) * err1, lw=0.2, ls='-', color=colors[0])
axes['A'].axhline(np.sqrt(2) * err2, lw=0.2, ls='--', color=colors[1])

axes['A'].scatter(t1_min, sf1_min, marker='x', zorder=3, s=8, linewidth=1.2, color='blue')
axes['A'].scatter(t2_min, sf2_min, marker='x', zorder=3, s=8, linewidth=1.2, color='red')

axes['A'].scatter(bin_sf_1_t[idx_br_1], bin_sf_1_sf[idx_br_1], 
                  marker='^', s=5, color='black', linewidth=0.2,facecolors='none',)
axes['A'].scatter(bin_sf_1_t[~idx_br_1], bin_sf_1_sf[~idx_br_1], 
                  marker='^', s=5, linewidth=0.2, facecolors='none', color='gray')

axes['A'].errorbar(bin_sf_1_t[idx_br_1], bin_sf_1_sf[idx_br_1],
                   yerr=np.vstack((bin_sf_1_sf_err[idx_br_1], bin_sf_1_sf_err[idx_br_1])),
                   ms=1, ls='none',mew=0.2, elinewidth=0.1, mec='black',ecolor='black', capsize=0.,mfc='none')

axes['A'].errorbar(bin_sf_1_t[~idx_br_1], bin_sf_1_sf[~idx_br_1],
                   yerr=np.vstack((bin_sf_1_sf_err[~idx_br_1], bin_sf_1_sf_err[~idx_br_1])),
                   ms=2, ls='none', mew=0.2, elinewidth=0.1, mec='gray', ecolor='gray',mfc='none', capsize=0.)

axes['A'].scatter(bin_sf_2_t[idx_br_2], bin_sf_2_sf[idx_br_2], 
                  marker='s', s=5, color='black', linewidth=0.2,facecolors='none',)
axes['A'].scatter(bin_sf_2_t[~idx_br_2], bin_sf_2_sf[~idx_br_2], 
                  marker='s', s=5, linewidth=0.2, facecolors='none', color='gray')

axes['A'].errorbar(bin_sf_2_t[idx_br_2], bin_sf_2_sf[idx_br_2],
                   yerr=np.vstack((bin_sf_2_sf_err[idx_br_2], bin_sf_2_sf_err[idx_br_2])),
                   ms=1, ls='none',mew=0.2, elinewidth=0.1, mec='black',ecolor='black', capsize=0.,mfc='none')
axes['A'].errorbar(bin_sf_2_t[~idx_br_2], bin_sf_2_sf[~idx_br_2],
                   yerr=np.vstack((bin_sf_2_sf_err[~idx_br_2], bin_sf_2_sf_err[~idx_br_2])),
                   ms=2, ls='none', mew=0.2, elinewidth=0.1, mec='gray', ecolor='gray',mfc='none', capsize=0.)

axes['A'].plot(fit[0][0], fit[0][1], lw=0.6, color=colors[0])#, label=r'$%s$-band'%str(band_pair[0]))
axes['A'].plot(fit[1][0], fit[1][1], lw=0.6, color=colors[1], ls='--')#, label=r'$%s$-band'%str(band_pair[0]))


axes['A'].set_yscale('log')
axes['A'].set_ylabel('SF [mag]', fontsize=10)

axes['A'].set_ylim(6e-3, 1e-1)

axes['A'].yaxis.set_major_formatter(ScalarFormatter()) 
axes['A'].set_yticklabels(['','', 0.01, 0.1])


axes['B'].errorbar(cv_1_t[~idx_cv_1], cv_1_y[~idx_cv_1],
                   xerr= bin_width_all[~idx_cv_1], 
                   yerr=cv_1_y_err[~idx_cv_1], fmt='o', mfc='none',mew=0.2,elinewidth=0.2, ms=2, lw=0.5,ls='none',
                   color='black', zorder=2, alpha=0.5,capsize=0.5)

axes['B'].scatter(cv_1_t[idx_cv_1], cv_1_y[idx_cv_1], 
                  marker='o', s=5,  linewidth=0.2, color='magenta')#facecolors='none', 

axes['B'].errorbar(cv_1_t[idx_cv_1], cv_1_y[idx_cv_1],
                   xerr= np.vstack((bin_width_all[idx_cv_1], bin_width_all[idx_cv_1])), 
                   yerr=np.vstack((cv_1_y_err[idx_cv_1], cv_1_y_err[idx_cv_1])), color='black',
                   mew=0.2, elinewidth=0.2, lw=0.3,ls='-', capsize=0.5)#, label='%s vs. %s'%(band_pair[0], band_pair[1]), color='black', ecolor='black', zorder=2)

axes['B'].axhline(y=1 / flux_ratio[0], xmin=0, xmax=1e4, color='gray', ls='--', zorder=0, lw=0.5) # ewidth


axes['B'].set_ylabel(r'$C_f(\tau)$', fontsize=10)

axes['B'].set_yticks(np.arange(0.30, 1.5, 0.1), major=True)
axes['B'].set_yticks(np.arange(0.30, 1.52, 0.02), minor=True)
axes['B'].set_ylim(0.5, 1.1)


axes['B'].text(0.05, 0.05, 'BWB', transform=axes['B'].transAxes, color='blue', fontsize=7)
axes['B'].text(0.05, 0.9, 'RWB', transform=axes['B'].transAxes, color='red', fontsize=7)


axes['B'].set_xlabel(r'$\tau$ [sec]', fontsize=10)
axes['A'].text(0.9, 0.9, '(a)', transform=axes['A'].transAxes, color='black', fontsize=7)
axes['B'].text(0.9, 0.9, '(c)', transform=axes['B'].transAxes, color='black', fontsize=7)

axes['A'].text(0.8, 0.1, r'$\sqrt{2}\sigma_n^i$', transform=axes['A'].transAxes, color=colors[0], fontsize=5)
axes['A'].text(0.8, 0.2, r'$\sqrt{2}\sigma_n^z$', transform=axes['A'].transAxes, color=colors[1], fontsize=5)

axes['A'].scatter(150, 0.047, marker='s', s=8, color='black',
                  facecolors='none', lw=0.3)
axes['A'].plot([130, 170], [0.047, 0.047], ls='--', color=colors[1], lw=0.5)

axes['A'].scatter(150, 0.063, marker='^', s=8, color='black',
                  facecolors='none', lw=0.3)

axes['A'].plot([130, 170], [0.063, 0.063], ls='-', color=colors[0], lw=0.5)
axes['A'].text(200, 0.06, r'$%s$-band'%str(band1), color=colors[0], fontsize=6)
axes['A'].text(195, 0.045, r'$%s$-band'%str(band2),color=colors[1], fontsize=6)

plt.setp(axes['A'].get_xticklabels(), visible=False)


use_bin = 10 ** np.arange(1, 5.2, 0.1)
axes['C'].hist(num_all, bins=use_bin, histtype='step', color='black', ls=(5, (11, 3)),lw=0.5, label='All pairs')
axes['C'].hist(num_cut, bins=use_bin, histtype='step',ls=(0, (1, 0.5)), lw=0.5, color='blue',  label=r'$3\sigma$-cut', zorder=4)
axes['C'].hist(num_pos, bins=use_bin, histtype='step',lw=0.7, color='magenta', label=r'$3\sigma$-${\rm cut}~\&~C_{f/m}>0$', zorder=3)
# axes['A'].legend(loc='upper left', frameon=False, fontsize=5)
axes['C'].set_ylabel('number of pairs', fontsize=10)
axes['C'].set_yscale('log')
axes['C'].set_ylim(2, 3e4)

axes['D'].set_xlabel(r'$\tau$ [sec]', fontsize=10)
axes['C'].text(0.9, 0.9, '(b)', transform=axes['C'].transAxes, color='black', fontsize=7)
axes['D'].text(0.9, 0.9, '(d)', transform=axes['D'].transAxes, color='black', fontsize=7)

axes['D'].axhline(y=1.0, xmin=0, xmax=1e4, color='gray', ls='--', zorder=0, lw=0.5) # ewidth


axes['D'].errorbar(cv_2_t[~idx_cv_2], cv_2_y[~idx_cv_2],
                   xerr= bin_width_all[~idx_cv_2], 
                   yerr=cv_2_y_err[~idx_cv_2], fmt='o', mfc='none',
                   mew=0.2,elinewidth=0.2, ms=2, lw=0.5,ls='none',
                   color='black', zorder=2, alpha=0.5,capsize=0.5)

axes['D'].scatter(cv_2_t[idx_cv_2], cv_2_y[idx_cv_2], 
                  marker='o', s=5,  linewidth=0.2, color='blue')#facecolors='none', 

axes['D'].errorbar(cv_2_t[idx_cv_2], cv_2_y[idx_cv_2],
                   xerr= np.vstack((bin_width_all[idx_cv_2], bin_width_all[idx_cv_2])), 
                   yerr=np.vstack((cv_2_y_err[idx_cv_2], cv_2_y_err[idx_cv_2])), color='black',
                   mew=0.2, elinewidth=0.2, lw=0.3,ls='-', capsize=0.5)#, label='%s vs. %s'%(band_pair[0], band_pair[1]), color='black', ecolor='black', zorder=2)


axes['D'].text(0.05, 0.05, 'BWB', transform=axes['D'].transAxes, color='blue', fontsize=7)
axes['D'].text(0.05, 0.9, 'RWB', transform=axes['D'].transAxes, color='red', fontsize=7)

axes['D'].set_yticks(np.arange(0.5, 1.5, 0.1), major=True)
axes['D'].set_yticks(np.arange(0.5, 1.5, 0.02), minor=True)
axes['D'].set_ylim(0.7,1.3)
# axes['B'].set_ylabel(r'$\bar{\vartheta}(\tau)$ [deg]', fontsize=10)
axes['D'].set_ylabel(r'$C_m(\tau)$', fontsize=10)

plt.setp(axes['A'].get_xticklabels(), visible=False)
plt.setp(axes['C'].get_xticklabels(), visible=False)
# plt.setp(axes['B'].get_xticklabels(), visible=False) 
for ax in ['A', 'B', 'C', 'D']:#, 'C'
    axes[ax].set_xscale('log')
    axes[ax].set_xlim(80, 3e4)
    axes[ax].tick_params(axis='x', labelsize=8)
    axes[ax].tick_params(axis='y', labelsize=8)
    
axes['A'].fill_between(np.arange(t_used_min, t_used_max, 0.1), 2e-3, 1e-1,
                       facecolor='palevioletred', alpha=0.15,zorder=0)
axes['B'].fill_between(np.arange(t_used_min, t_used_max, 0.1), 0., 2,
                       facecolor='palevioletred', alpha=0.15, zorder=0)
axes['C'].fill_between(np.arange(t_used_min, t_used_max, 0.1), 1, 3e4,
                       facecolor='palevioletred', alpha=0.15, zorder=0)
axes['D'].fill_between(np.arange(t_used_min, t_used_max, 0.1), 0, 2,
                       facecolor='palevioletred', alpha=0.15, zorder=0)
                       
# ==========================================================================


axes['C'].legend(loc='upper left', frameon=False, fontsize=4)


# ==========================================================================
plt.tight_layout()
plt.subplots_adjust(wspace=0.3, hspace=0.03)
# plt.savefig('./figure/mag_flux_new_sun14_test.png', dpi=500, bbox_inches='tight')
# plt.savefig('./figure/mag_flux_new_sun14_sf_0.1_cv_0.5.png', dpi=500, bbox_inches='tight')
plt.show()