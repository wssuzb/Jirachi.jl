import numpy as np
from scipy.stats import binned_statistic

def color_variation(t, r1, r2, e1, e2, mode='flux', erron=True, nsigma=3):

    time_2d = np.tile(t, (len(t), 1))
    dt_all = time_2d - np.transpose(time_2d)     

    x_wmin_2d = np.tile(r1, (len(r1), 1))
    dx_wmin_all = x_wmin_2d - np.transpose(x_wmin_2d)

    x_wmax_2d = np.tile(r2, (len(r2), 1))
    dx_wmax_all = x_wmax_2d - np.transpose(x_wmax_2d)
    idx = np.where(dt_all>0, True, False)

    dt_all = dt_all[idx]
    dx_wmin_all = dx_wmin_all[idx]
    dx_wmax_all = dx_wmax_all[idx]

    theta = dx_wmax_all / dx_wmin_all

    if erron:
        ex_wmin_1e = np.tile(e1, (len(e1), 1))
        ex_wmax_1e = np.tile(e2, (len(e2), 1))

        ex_wmin_2e, ex_wmax_2e = np.transpose(ex_wmin_1e), np.transpose(ex_wmax_1e)

        ex_wmin_1e, ex_wmax_1e = ex_wmin_1e[idx], ex_wmax_1e[idx]
        ex_wmin_2e, ex_wmax_2e = ex_wmin_2e[idx], ex_wmax_2e[idx]
        
        dx_ratio2 = (dx_wmax_all / dx_wmin_all) ** 2

        mmd = np.sqrt(dx_wmin_all ** 2 + dx_wmax_all ** 2)

        err_mmd1 = ex_wmin_1e * ex_wmax_1e * np.sqrt((1 + dx_ratio2) / (ex_wmax_1e ** 2 + ex_wmin_1e ** 2 * dx_ratio2))
        err_mmd2 = ex_wmin_2e * ex_wmax_2e * np.sqrt((1 + dx_ratio2) / (ex_wmax_2e ** 2 + ex_wmin_2e ** 2 * dx_ratio2)) 

        err_mmd = np.sqrt(err_mmd1 ** 2 + err_mmd2 ** 2)

        err_mmd_nsigma = nsigma * err_mmd

        criterion = np.where(mmd >= err_mmd_nsigma, True, False)
        
        dt_all = dt_all[criterion]
        theta = theta[criterion]
        # print(len(theta))
    else:
        pass
    
    int0 = np.argsort(dt_all, kind='mergesort')
    dt_all_sort = dt_all[int0]
    theta_sort = theta[int0]

    idx_pos = np.where(theta_sort > 0, True, False)
    return dt_all_sort[idx_pos], theta_sort[idx_pos]


def get_bin(nbins, dt, theta, method='mean'):

    def err_bootstrapped(x, num):
        res = []
        for i in range(num):
            idx = np.random.randint(0, len(x), len(x))
            unique, counts = np.unique(idx, return_counts=True)
            res.append(np.median(x[unique]))
        return np.std(res)

    bin_num = binned_statistic(dt, dt, statistic='count', bins=nbins)[0]
    bin_edge = binned_statistic(dt, dt, statistic=lambda x: len(x), bins=nbins)[1]

    bin_dt = binned_statistic(dt, dt, statistic='median', bins=nbins)[0]
    bin_theta = binned_statistic(dt, theta, statistic='median', bins=nbins)[0]
    bin_rbst = binned_statistic(dt, theta, statistic=lambda x: err_bootstrapped(np.array(x), 100), bins = nbins)[0]
    
    yerr = np.vstack((bin_rbst, bin_rbst))

    return [bin_dt, bin_theta, yerr, bin_edge, bin_num]

if __name__ == "__main__":

    t1 = np.loadtxt("../test/data/montano22_n1_i_binned.txt", usecols = [0, 2, 3])
    t2 = np.loadtxt("../test/data/montano22_n1_z_binned.txt", usecols = [0, 2, 3])

    it, iy, ie = t1[:, 0], t1[:, 1], t1[:, 2]
    zt, zy, ze = t2[:, 0], t2[:, 1], t2[:, 2]

    res = color_variation(it, iy, zy, ie, ze, mode="flux", erron=True, nsigma=3)
