#%%
import variation as var
import sys
import warnings
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as ticker
from scipy.optimize import curve_fit
from scipy import stats
import matplotlib
# print(matplotlib.rcParams['mathtext.fontset'])
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
    # "figure.dpi" : 500
})

band = ["g", "r", "i", "z"]

band_pair = []

for i in range(len(band)):
   for j in range(1, len(band)):
       if i<j:
           band_pair.append([band[i], band[j]])
       

used_band = ['i', 'z']

band1, band2 = used_band[0], used_band[1]

t_cad = 103.68 # bin size of lightcurve
t_start = 0 - (t_cad / 2) # make sure that the bin starts centered at 0 second
t_end = 25000 + (t_cad /2 ) # the bin end 
lc_edges = np.arange(t_start, t_end + t_cad, t_cad)

sf_bin = 10 ** np.arange(0, 5.05, 0.05) # to bin the structure function
cv_bin = 10 ** np.arange(0, 5.1, 0.1) # to bin the color variation


# --------------------------------------------------------------------------->
# Loading data
# --------------------------------------------------------------------------->

lc1_dat = np.loadtxt("../test/ngc4395/Montano22_n1/%s_4395.txt"%str(used_band[0]), usecols=[0, 1, 2])
lc2_dat = np.loadtxt("../test/ngc4395/Montano22_n1/%s_4395.txt"%str(used_band[1]), usecols=[0, 1, 2])

t1, y1, e1 = lc1_dat[:, 0], lc1_dat[:, 1], lc1_dat[:, 2]
t1 = (t1 - t1[0]) * 24 * 3600 - 1e-9 # for some reasons, we shift the time by 1e-9 to match the results from `julia` version, and it do not influence our final analysis.
lc1 = var.lightcurve(t1, y1, e1)

t2, y2, e2 = lc2_dat[:, 0], lc2_dat[:, 1], lc2_dat[:, 2]
t2 = (t2 - t2[0]) * 24 * 3600 - 1e-9 # for some reasons, we shift the time by 1e-9 to match the results from `julia` version, and it do not influence our final analysis.
lc2 = var.lightcurve(t2, y2, e2)

# --------------------------------------------------------------------------->
# Binning lightcurves
# --------------------------------------------------------------------------->
lc1_bin = var.bin_lightcurve(lc1, lc_edges)
lc2_bin = var.bin_lightcurve(lc2, lc_edges)

# --------------------------------------------------------------------------->
# Calculating structure function
# --------------------------------------------------------------------------->
binsf1 = var.structure_function_mcmc(lc1_bin, sf_bin, mode="both", nsim=1000)
binsf2 = var.structure_function_mcmc(lc2_bin, sf_bin, mode="both", nsim=1000)

# --------------------------------------------------------------------------->
# Finding the proper time range
# --------------------------------------------------------------------------->
fit1_res = var.find_proper_time(binsf1, p0 = [1, 1e3, 1, 0.005],  lower_bounds = [0, 0, 0, 0.001], upper_bounds = [10, 2e4, 2, 0.01 ],)
fit2_res = var.find_proper_time(binsf2, p0 = [1, 1e3, 1, 0.005],lower_bounds = [0, 0, 0, 0.001], upper_bounds = [10, 2e4, 2, 0.01 ], )


t1_max, sf1_max, t1_min, sf1_min, t1_fit, sf1_fit, par1, par1_err = fit1_res.t_max, fit1_res.sf_max, fit1_res.t_min, fit1_res.sf_min, fit1_res.t_fit, fit1_res.sf_fit, fit1_res.par, fit1_res.par_err

t2_max, sf2_max, t2_min, sf2_min, t2_fit, sf2_fit, par2, par2_err = fit2_res.t_max, fit2_res.sf_max, fit2_res.t_min, fit2_res.sf_min, fit2_res.t_fit, fit2_res.sf_fit, fit2_res.par, fit2_res.par_err

t_used_min = max(t1_min, t2_min)
t_used_max = min(t1_max, t2_max)

bin_sf_1_t = binsf1.x
bin_sf_1_t_err = binsf1.xerr
bin_sf_1_sf = binsf1.y
bin_sf_1_sf_err = binsf1.yerr

idx_br_1 = np.where(bin_sf_1_t <= t1_max, True, False)

bin_sf_2_t = binsf2.x
bin_sf_2_t_err = binsf2.xerr
bin_sf_2_sf = binsf2.y
bin_sf_2_sf_err = binsf2.yerr

idx_br_2 = np.where(bin_sf_2_t <= t2_max, True, False)

# --------------------------------------------------------------------------->
# getting the new lightcurves for calculating color variation in flux-flux and mag-mag
# --------------------------------------------------------------------------->
lc1_new, lc2_new = var.get_common_lc(lc1_bin, lc2_bin)
cv_flux = var.color_variation(lc1_new, lc2_new, erron=True, nsigma=3, used_bins=10 ** np.arange(0, 5.1, 0.1), showhist=True)
cv_mag = var.color_variation(lc1_new, lc2_new, erron=True, nsigma=3, used_bins=10 ** np.arange(0, 5.1, 0.1), mode='mag')

flux_ratio = np.mean(lc1_new.flux) / np.mean(lc2_new.flux) #f['flux_ratio'].value.flatten()

# flux-flux

num_all = cv_flux.num_all
num_cut = cv_flux.num_cut
num_pos = cv_flux.num_pos

cv_1_t = cv_flux.x

cv_1_t_err = cv_flux.xerr
# cv_1_t_err = [(cv_1_t[i+1] - cv_1_t[i]) / 2 for i in range(len(cv_1_t)- 1)]
cv_1_y = cv_flux.y
cv_1_y_err = cv_flux.yerr

idx_cv_1 = np.where((cv_1_t>=t_used_min) & (cv_1_t<=t_used_max), True, False)

# mag-mag

cv_2_t = cv_mag.x
cv_2_t_err = cv_mag.xerr
cv_2_y = cv_mag.y
cv_2_y_err = cv_mag.yerr

idx_cv_2 = np.where((cv_2_t>=t_used_min) & (cv_2_t<=t_used_max), True, False)

# cv_use_bin= 10 ** np.arange(0., 5.1, 0.1)
# bin_width_all = [(cv_use_bin[i+1] - cv_use_bin[i]) / 2 for i in range(len(cv_use_bin)- 1)]
# xerr_all = np.vstack((bin_width_all, bin_width_all))
# bin_width_all = np.array(bin_width_all)
# bin_centers_all = cv_use_bin[1:] - bin_width_all

# -------------------------------------------------------------->
# For plotting 
# -------------------------------------------------------------->
figure_mosaic = """
AC
BD
"""

colors = [(46/255, 89/255, 167/255),'xkcd:orange red']

fig, axes = plt.subplot_mosaic(mosaic = figure_mosaic, figsize=(4.2, 3.5))

axes['A'].axhline(np.sqrt(2) * par1[-1], lw=0.2, ls='-', color=colors[0])
axes['A'].axhline(np.sqrt(2) * par2[-1], lw=0.2, ls=(0, (15, 3)), color=colors[1])
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
axes['A'].plot(t1_fit, sf1_fit, lw=0.6, color=colors[0])
axes['A'].plot(t2_fit, sf2_fit, lw=0.6, color=colors[1], ls=(0, (5, 3)))#, 
axes['A'].set_yscale('log')
axes['A'].set_ylabel('SF [mag]', fontsize=10)
axes['A'].set_ylim(6e-3, 1e-1)
axes['A'].yaxis.set_major_formatter(ScalarFormatter()) 
axes['A'].set_yticklabels(['','', 0.01, 0.1])
axes['A'].text(0.9, 0.9, '(a)', transform=axes['A'].transAxes, color='black', fontsize=7)
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

axes['B'].errorbar(cv_1_t[~idx_cv_1], cv_1_y[~idx_cv_1],
                   xerr= cv_1_t_err[~idx_cv_1], 
                   yerr=cv_1_y_err[~idx_cv_1], fmt='o', mfc='none',mew=0.2,elinewidth=0.2, ms=2, lw=0.5,ls='none',
                   color='black', zorder=2, alpha=0.5,capsize=0.5)
axes['B'].scatter(cv_1_t[idx_cv_1], cv_1_y[idx_cv_1], 
                  marker='o', s=5,  linewidth=0.2, color='magenta')#facecolors='none', 
axes['B'].errorbar(cv_1_t[idx_cv_1], cv_1_y[idx_cv_1],
                   xerr= np.vstack((cv_1_t_err[idx_cv_1], cv_1_t_err[idx_cv_1])), 
                   yerr=np.vstack((cv_1_y_err[idx_cv_1], cv_1_y_err[idx_cv_1])), color='black',
                   mew=0.2, elinewidth=0.2, lw=0.3,ls='-', capsize=0.5)#, label='%s vs. %s'%(band_pair[0], band_pair[1]), color='black', ecolor='black', zorder=2)
axes['B'].axhline(y=1 / flux_ratio, xmin=0, xmax=1e4, color='gray', ls='--', zorder=0, lw=0.5) # ewidth
axes['B'].set_ylabel(r'$C_f(\tau)$', fontsize=10)
axes['B'].set_yticks(np.arange(0.30, 1.5, 0.1), major=True)
axes['B'].set_yticks(np.arange(0.30, 1.52, 0.02), minor=True)
axes['B'].set_ylim(0.5, 1.1)
axes['B'].text(0.05, 0.05, 'BWB', transform=axes['B'].transAxes, color='blue', fontsize=7)
axes['B'].text(0.05, 0.9, 'RWB', transform=axes['B'].transAxes, color='red', fontsize=7)
axes['B'].set_xlabel(r'$\tau$ [sec]', fontsize=10)
axes['B'].text(0.9, 0.9, '(c)', transform=axes['B'].transAxes, color='black', fontsize=7)


use_bin = 10 ** np.arange(1, 5.2, 0.1)
axes['C'].hist(num_all, bins=use_bin, histtype='step', color='black', ls=(5, (11, 3)),lw=0.5, label='All pairs')
axes['C'].hist(num_cut, bins=use_bin, histtype='step',ls=(0, (1, 0.5)), lw=0.5, color='blue',  label=r'$3\sigma$-cut', zorder=4)
axes['C'].hist(num_pos, bins=use_bin, histtype='step',lw=0.7, color='magenta', label=r'$3\sigma$-${\rm cut}~\&~C_{f/m}>0$', zorder=3)
axes['C'].legend(loc='upper left', frameon=False, fontsize=5)
axes['C'].set_ylabel('number of pairs', fontsize=10)
axes['C'].set_yscale('log')
axes['C'].set_ylim(2, 3e4)

axes['D'].set_xlabel(r'$\tau$ [sec]', fontsize=10)
axes['C'].text(0.9, 0.9, '(b)', transform=axes['C'].transAxes, color='black', fontsize=7)
axes['D'].text(0.9, 0.9, '(d)', transform=axes['D'].transAxes, color='black', fontsize=7)

axes['D'].axhline(y=1.0, xmin=0, xmax=1e4, color='gray', ls='--', zorder=0, lw=0.5) # ewidth


axes['D'].errorbar(cv_2_t[~idx_cv_2], cv_2_y[~idx_cv_2],
                   xerr= cv_2_t_err[~idx_cv_2], 
                   yerr=cv_2_y_err[~idx_cv_2], fmt='o', mfc='none',
                   mew=0.2,elinewidth=0.2, ms=2, lw=0.5,ls='none',
                   color='black', zorder=2, alpha=0.5,capsize=0.5)

axes['D'].scatter(cv_2_t[idx_cv_2], cv_2_y[idx_cv_2], 
                  marker='o', s=5,  linewidth=0.2, color='blue')#facecolors='none', 

axes['D'].errorbar(cv_2_t[idx_cv_2], cv_2_y[idx_cv_2],
                   xerr= np.vstack((cv_2_t_err[idx_cv_2], cv_2_t_err[idx_cv_2])), 
                   yerr=np.vstack((cv_2_y_err[idx_cv_2], cv_2_y_err[idx_cv_2])), color='black',
                   mew=0.2, elinewidth=0.2, lw=0.3,ls='-', capsize=0.5)

axes['D'].text(0.05, 0.05, 'BWB', transform=axes['D'].transAxes, color='blue', fontsize=7)
axes['D'].text(0.05, 0.9, 'RWB', transform=axes['D'].transAxes, color='red', fontsize=7)

axes['D'].set_yticks(np.arange(0.5, 1.5, 0.1), major=True)
axes['D'].set_yticks(np.arange(0.5, 1.5, 0.02), minor=True)
axes['D'].set_ylim(0.7,1.3)
axes['D'].set_ylabel(r'$C_m(\tau)$', fontsize=10)

plt.setp(axes['A'].get_xticklabels(), visible=False)
plt.setp(axes['C'].get_xticklabels(), visible=False)

for ax in ['A', 'B', 'C', 'D']:#, 'C'
    axes[ax].set_xscale('log')
    axes[ax].set_xlim(80, 3e4)
    axes[ax].tick_params(axis='x', labelsize=8)
    axes[ax].tick_params(axis='y', labelsize=8)

color_filled = "palevioletred"

axes['A'].fill_between(np.arange(t_used_min, t_used_max, 0.1), 2e-3, 1e-1,
                       facecolor=color_filled, alpha=0.15, zorder=0)
axes['B'].fill_between(np.arange(t_used_min, t_used_max, 0.1), 0., 2,
                       facecolor=color_filled, alpha=0.15, zorder=0)
axes['C'].fill_between(np.arange(t_used_min, t_used_max, 0.1), 1, 3e4,
                       facecolor=color_filled, alpha=0.15, zorder=0)
axes['D'].fill_between(np.arange(t_used_min, t_used_max, 0.1), 0, 2,
                       facecolor=color_filled, alpha=0.15, zorder=0)
           
plt.tight_layout()
plt.subplots_adjust(wspace=0.3, hspace=0.03)
plt.savefig('./mag_flux_new_sun14_%s_%s.png'%(str(used_band[0]), str(used_band[1])), dpi=500, bbox_inches='tight')
