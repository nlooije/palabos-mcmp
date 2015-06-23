
# run with: python analysis.py

import os, glob
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

# Windows only: set working directory
#os.chdir('C:/Users/niels/Desktop/large_density_droplet')

def calculate_radius_width(fname, plots_on=True):
    """ Calculate radius & width of interface from density profile according to:
            rho = alpha + beta * tanh[-(z-z_int)/width]
        where alpha = 0.5*(rho_l + rho_v), beta = 0.5*(rho_l - rho_v)
        The code performs a linear fit on:
            atanh[(rho-alpha)/beta] = -z/width + radius/width
        where the slope provides a measure for the width and the intercept gives
        the radius.
    """
    # load files
    rho = np.loadtxt(fname)
    nx, z = (len(rho)-1)/2-1, np.arange(len(rho))
    z, rho = z[0:nx], rho[0:nx]

    # calculate parameters
    rho_l, rho_v = rho[0], rho[-1]
    alpha, beta = 0.5*(rho_l+rho_v), 0.5*(rho_l-rho_v)

    # linearize & linear fit
    theta = np.arctanh((rho-alpha)/beta)
    mask = ~np.isnan(z) & ~np.isinf(z) & ~np.isnan(theta) & ~np.isinf(theta) # remove nans & infs
    slope, intercept, r_value, p_value, std_err = stats.linregress(z[mask],theta[mask])
    if plots_on:
        print slope, intercept, r_value
        plt.plot(z, theta, 'o', z, slope*z + intercept, '-')
        plt.show()

    # calculate radius, interface location & thickness
    width = -1./slope
    z_int = intercept*width
    radius = nx-z_int

    # plot numerical and analytical profiles
    rho_ana = alpha + beta * np.tanh(-(z-z_int)/width)
    if plots_on:
        plt.plot(z, rho, 'o', z, rho_ana, '-')
        plt.show()
        
    # print 'fname: ', fname, ', radius: ', radius, ', width: ', width
    return (radius, width, rho_l, rho_v)

def calculate_pressure_vdW(rho, Tr):
    """van der Waals equation of state"""
    a, b = 9./49., 2./21.
    RTc = 8./27.*a/b
    return rho*RTc*Tr/(1.-b*rho)-a*rho**2
    
if __name__ == "__main__":
    # Initialization
    reducedT = 0.85
    plots_on = False

    # loop over all density profile files
    press_diff, inv_radii = [], []
    files = glob.glob('rho_profile*')
    for fname in files:
        r, w, rho_l, rho_v = calculate_radius_width(fname, plots_on)
        press_l = calculate_pressure_vdW(rho_l, reducedT)
        press_v = calculate_pressure_vdW(rho_v, reducedT) 
        press_diff.append(press_v - press_l)
        inv_radii.append(1./r)
    
    # create numpy arrays from python lists
    press_diff = np.asarray(press_diff)
    inv_radii = np.asarray(inv_radii)
    # print pressures, inv_radii
    
    # calculate surface tension as slope of pressure difference 
    # vs radius of curvature (and plot/fit, calculate residual sum of squares error)
    slope = np.dot(press_diff,inv_radii)/np.dot(inv_radii, inv_radii)
    rssq = np.dot(press_diff - slope*inv_radii, press_diff - slope*inv_radii)

    # plot and save laplace law
    plt.plot(inv_radii, press_diff, 'o', inv_radii, slope*inv_radii, '-')
    plt.savefig('laplace_law.png')
    plt.show()

    print 'surface_tension: ', slope
    print 'error: ', rssq
