import numpy as np
from scipy.stats import binned_statistic
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from collections import namedtuple
import matplotlib.pyplot as plt


lightcurve = namedtuple('lc', ('time', 'flux', 'err'))
binned_data = namedtuple('bin', ('x', 'xerr', 'y', 'yerr'))
binned_data_cv = namedtuple('bin', ('x', 'xerr', 'y', 'yerr', 'num_all', 'num_cut', 'num_pos'))
fit_result = namedtuple('fit', ('t_max', 'sf_max', 't_min', 'sf_min', 't_fit', 'sf_fit', 'par', 'par_err'))

def bin_lightcurve(lc, used_bins):
    """
    binning the light curves

    Inputs:

        lc -- input lightcurves, should be in type of namedtuple-like
        used_bins -- the bins for binning lightcurve
    
    Outputs:

        lightcurve(bin_dt_center[~idx], bin_flux[~idx], bin_err[~idx]) -- the resultant lightcurves after binning
    """
    if not isinstance(lc, lightcurve):
        print("light curves should be input with namedtuple-like type!!!")
        
    bin_dt_width = np.array([(used_bins[i+1] - used_bins[i]) / 2 for i in range(len(used_bins) - 1)])
    bin_dt_center = used_bins[1:] - bin_dt_width

    bin_flux = binned_statistic(lc.time, lc.flux, statistic='mean', bins=used_bins)[0]

    bin_err = binned_statistic(lc.time, lc.err, statistic=lambda x: np.sqrt(np.sum(np.array(x) ** 2)) / len(x), bins=used_bins)[0]

    idx = np.isnan(bin_flux)

    return lightcurve(bin_dt_center[~idx], bin_flux[~idx], bin_err[~idx])


def get_common_lc(lc1, lc2):
    """
    calculating color variation requires quasi-simultaneou observations, thus this function help get the overlap observed time bewteen 2 lightcurves.

    Inputs: 
        lc1 -- first lightcurve
        lc2 -- second lightcurve

    Outputs:
        lc1 -- new first lightcurve
        lc2 -- new first lightcurve
    """


    common_t = np.intersect1d(lc1.time, lc2.time, return_indices=True)
    
    return lightcurve(lc1.time[common_t[1]], lc1.flux[common_t[1]], lc1.err[common_t[1]]), lightcurve(lc2.time[common_t[2]], lc2.flux[common_t[2]], lc2.err[common_t[2]])

def lc_mcmc(lc, mode='both'):
    """
    fr/rss methods, flux randomness and random subset selection for lightcurves.
    reference: Fausnaugh+16

    Inputs: 
        lc -- lightcurve
        mode -- "both": fr/rss; "fr": fr only; "rss": rss only.
    
    Outputs: lc -- lightcurve after sampling
    """    
    
    if not isinstance(lc, lightcurve):
        print("light curves should be input with namedtuple-like type!!!")


    numt = len(lc.time)
    indx = np.random.randint(0, numt, numt)
    unique, counts = np.unique(indx, return_counts=True) # sorted unique value
    t_rss = lc.time[unique]
    y_rss = lc.flux[unique]
    e_rss = lc.err[unique]/np.sqrt(counts)
    
    y_fr_rss = np.random.normal(y_rss, e_rss)
    idx_sort = np.argsort(t_rss, kind='mergesort')
    
    if mode == "rss":
        return lightcurve(t_rss[idx_sort], y_rss[idx_sort], e_rss[idx_sort])
    elif mode == "both":
        return lightcurve(t_rss[idx_sort], y_fr_rss[idx_sort], e_rss[idx_sort])
    else:
        return lightcurve(lc.time, np.random.normal(lc.flux, lc.err), lc.err)

def structure_function(lc, used_bins):
    """
    calculate the structure function without removing photometric uncertainties.

    reference: kozlowski+16
    
    Inputs: 
        lc -- lightcurve
        used_bins -- used for bin structure function    

    Outputs: the binned structure function
    """

    if not isinstance(lc, lightcurve):
        print("light curves should be input with namedtuple-like type!!!")
    
    time = lc.time    
    flux = -2.5 * np.log10(lc.flux)
    
    time_2d = np.tile(time, (len(time), 1))
    dt_all = time_2d - np.transpose(time_2d)    
    x_flux_2d = np.tile(flux, (len(flux), 1))
    dx_flux_all = x_flux_2d - np.transpose(x_flux_2d)
    
    idx = np.where(dt_all>0, True, False)

    dt_all = dt_all[idx]
    dx_flux_all = dx_flux_all[idx]

    int0 = np.argsort(dt_all, kind='mergesort')
    dt_all_sort = dt_all[int0]
    flux_sort = np.sqrt(np.pi / 2 * dx_flux_all[int0] ** 2)
    
    bin_dt_width = np.array([(used_bins[i+1] - used_bins[i]) / 2 for i in range(len(used_bins) - 1)])
    bin_dt_center = used_bins[1:] - bin_dt_width

    bin_sf = binned_statistic(dt_all_sort, flux_sort, statistic=lambda x: np.nanmean(x), bins=used_bins)[0]

    bin_sf_err = binned_statistic(dt_all_sort, flux_sort, statistic=lambda x : np.nanstd(x) / np.sqrt(np.count_nonzero(~np.isnan(x))), bins=used_bins)[0]

    return binned_data(bin_dt_center, bin_dt_width, bin_sf, bin_sf_err)

def structure_function_mcmc(lc,  used_bins, mode='both', nsim=1000):
    """
    same as `structure_function`, but use the fr/rss to sampling the lightcurves.

    """


    if not isinstance(lc, lightcurve):
        print("light curves should be input with namedtuple-like type!!!")
    
    _tmp_sf = np.zeros((nsim, len(used_bins) - 1))

    for i in range(nsim):
        tmp_lc = lc_mcmc(lc, mode=mode)
        tmp_sf = structure_function(tmp_lc, used_bins)
        _tmp_sf[i] = tmp_sf.y

    sf_err_mcmc = np.zeros(len(used_bins) - 1)

    for i in range(len(used_bins) - 1):
        sf_err_mcmc[i] = np.nanstd(_tmp_sf[:, i])

    mysf = structure_function(lc, used_bins)
    
    idx = np.isnan(mysf.y)
    
    # Updating the structure function errors!

    return binned_data(mysf.x[~idx], mysf.xerr[~idx], mysf.y[~idx], sf_err_mcmc.flatten()[~idx])


def fit_sf(t, sf_inf, tau, beta, sigma):
    return np.sqrt(sf_inf ** 2 * (1-np.exp(-(t / tau) ** beta)) + 2 * sigma ** 2)

def find_t_max(binsf, prominence=1e-3, width=1e-3, check_plot=False):
    """
    find the local maximum structure function, and you can check if the maximum time is the proper one by set `check_plot=True`.

    for more details, see: https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html
    """
    peaks, properties = find_peaks(binsf.y, prominence=prominence, width=width)

    if check_plot:
        fig = plt.figure(figsize=(6, 4))
        plt.errorbar(binsf.x, binsf.y, binsf.yerr, label='SF')
        plt.scatter(binsf.x[peaks], binsf.y[peaks], marker='*', color='red', label='local maximum')
        plt.xlabel(r'$\Delta t$')
        plt.ylabel('SF [mag]')
        plt.xscale('log')
        plt.yscale('log')
        plt.legend()
        plt.show()

    return binsf.x[peaks], binsf.y[peaks]


def find_proper_time(binsf, prominence=1e-3, width=1e-3,  sf_noise_sigma = 2, p0 = [1, 1e3, 1, 0.005],lower_bounds = [], upper_bounds = [],  t_fit = 10 ** np.arange(0, 5.1, 0.1)):
    
    """
    
    find the proper time range of a given SF, see Su+24 for more details.

    """
    peaks, properties = find_peaks(binsf.y, prominence=prominence, width=width)

    t_max, sf_max = binsf.x[peaks][0], binsf.y[peaks][0]
    
    fit_idx = np.where(binsf.x <= t_max, True, False)

    if upper_bounds[1] >= t_max:
        upper_bounds[1] = t_max

    popt, pcov = curve_fit(fit_sf, binsf.x[fit_idx], binsf.y[fit_idx], sigma= binsf.yerr[fit_idx], p0=p0, bounds=(lower_bounds, upper_bounds))
    
    par, par_err = popt, np.sqrt(np.diag(pcov))
    
    sf_fit = fit_sf(t_fit, *popt)

    noise_cut = sf_noise_sigma * np.sqrt(2) * popt[-1] 

    t_min = np.interp(noise_cut, sf_fit, t_fit)

    sf_min = noise_cut

    return fit_result(t_max, sf_max, t_min, sf_min, t_fit, sf_fit, par, par_err)

def err_std(x, num=500):
    res = []
    x = np.array(x)
    for i in range(num):
        idx = np.random.randint(0, len(x), len(x))
        unique, counts = np.unique(idx, return_counts=True)
        res.append(np.median(x[unique]))
    
    return np.std(res)

def flux2mag(flux, err):
    mag = -2.5 * np.log10(flux) - 21.175
    mag_err = 2.5 / np.log(10) * (err / flux) # approximately!!!
    return mag, mag_err

def color_variation(lc1, lc2, erron, nsigma, used_bins, mode='flux', showhist=False):
    
    """
    
    calculte the timescale-dependent color variation, see Su+24 for more details.

    """

    if not (isinstance(lc1, lightcurve) & isinstance(lc2, lightcurve)):
        print("light curves should be input with namedtuple-like type!!!")
    
    if not (lc1.time == lc2.time).all():
        print("The Light curves should be quasi-simultaneou observed!!!")

    t1 = lc1.time

    if mode == 'mag':
        y1, e1 = flux2mag(lc1.flux, lc1.err)
        y2, e2 = flux2mag(lc2.flux, lc2.err)
    else:
        y1, e1 = lc1.flux, lc1.err
        y2, e2 = lc2.flux, lc2.err

    time_2d = np.tile(t1, (len(t1), 1))
    dt_all = time_2d - np.transpose(time_2d)     

    x_wmin_2d = np.tile(y1, (len(y1), 1))
    dx_wmin_all = x_wmin_2d - np.transpose(x_wmin_2d)

    x_wmax_2d = np.tile(y2, (len(y2), 1))
    dx_wmax_all = x_wmax_2d - np.transpose(x_wmax_2d)

    idx = np.where(dt_all > 0, True, False)

    dt_all = dt_all[idx]
    dx_wmin_all = dx_wmin_all[idx]
    dx_wmax_all = dx_wmax_all[idx]

    theta = dx_wmax_all / dx_wmin_all

    num_all = dt_all

    if erron:
        ex_wmin_1e = np.tile(e1, (len(e1), 1))
        ex_wmax_1e = np.tile(e2, (len(e2), 1))

        ex_wmin_2e, ex_wmax_2e = np.transpose(ex_wmin_1e), np.transpose(ex_wmax_1e)
        
        ex_wmin_1e, ex_wmax_1e = ex_wmin_1e[idx], ex_wmax_1e[idx]
        ex_wmin_2e, ex_wmax_2e = ex_wmin_2e[idx], ex_wmax_2e[idx]

        dx_ratio2 = (dx_wmax_all / dx_wmin_all) ** 2
        mmd = np.sqrt(
            dx_wmin_all ** 2 + dx_wmax_all ** 2
        )

        err_mmd1 = ex_wmin_1e * ex_wmax_1e * np.sqrt((1 + dx_ratio2) / (ex_wmax_1e ** 2 + ex_wmin_1e ** 2 * dx_ratio2))
        err_mmd2 = ex_wmin_2e * ex_wmax_2e * np.sqrt((1 + dx_ratio2) / (ex_wmax_2e ** 2 + ex_wmin_2e ** 2 * dx_ratio2)) 

        err_mmd = np.sqrt(err_mmd1 ** 2 + err_mmd2 ** 2)

        err_mmd_nsigma = nsigma * err_mmd

        criterion = np.where(mmd >= err_mmd_nsigma, True, False)

        dt_all = dt_all[criterion]
        theta = theta[criterion]
        num_cut = dt_all

    int0 = np.argsort(dt_all, kind='mergesort')
    dt_all_sort = dt_all[int0]
    theta_sort = theta[int0]

    idx_pos = np.where(theta_sort > 0, True, False)

    dt_all_pos, theta_sort_pos = dt_all_sort[idx_pos], theta_sort[idx_pos]

    num_pos = dt_all_pos

    bin_dt_width = np.array([(used_bins[i+1] - used_bins[i]) / 2 for i in range(len(used_bins) - 1)])
    bin_dt_center = used_bins[1:] - bin_dt_width

    bin_cv = binned_statistic(dt_all_pos, theta_sort_pos, statistic='median', bins=used_bins)[0]

    bin_cv_err = binned_statistic(dt_all_pos, theta_sort_pos, statistic=lambda x : err_std(x), bins=used_bins)[0]

    if showhist:
        return binned_data_cv(bin_dt_center, bin_dt_width, bin_cv, bin_cv_err, num_all, num_cut, num_pos)
    else:
        return binned_data(bin_dt_center, bin_dt_width, bin_cv, bin_cv_err)
    

def main():
    return 0

if __name__ == '__main__':
    main()