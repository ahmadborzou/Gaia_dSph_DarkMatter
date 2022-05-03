import pandas as pd
import numpy as np
import uncertainties.unumpy as unp
import matplotlib.pyplot as plt
from matplotlib import cm
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='large')
plt.rc('ytick', labelsize='large')
labelsize = 20
font = {
	'family': 'serif',
	'size': 20,
	'style': 'normal',
	'weight': 'medium',
	'fontname': 'Arial'
}

def Preprocess(dSph):
##########################################################
## distance of the dwarf spheroidal from the earth       #
  Distance = float(input("Insert the distance of the Galaxy from the earth in kilo-parsecs\
  (for example use table 1 of arXiv:0908.2995 [astro-ph.CO]) For example: use 147 for Fornax\n\n "))
  Dist_err = float(input("Insert the error of the distance of the Galaxy from the earth in kilo-parsecs\
  (for example use table 1 of arXiv:0908.2995 [astro-ph.CO]). For example: use 3 for Fornax\n\n "))

## D = Distance +/- Dist_err (kpc)                       #
  D = unp.uarray((Distance*1000.,Dist_err*1000.)) ## in pc #
##########################################################



## r hat unit vector
  def r_hat(th,ph):
    ## return the i_hat, j_hat, k_hat components
    return np.array([unp.sin(th)*unp.cos(ph), unp.sin(th)*unp.sin(ph), unp.cos(th)])

## theta hat 
  def th_hat(th,ph):
    return np.array([unp.cos(th)*unp.cos(ph), unp.cos(th)*unp.sin(ph), -unp.sin(th)])

## phi hat
  def ph_hat(th,ph):
    ## the third component is zero. The trick is to vectorize the computation
    return np.array([-unp.sin(ph), unp.cos(ph), unp.sin(th-th)])

  txt = "%s_dSph.csv"%(dSph)
  df = pd.read_csv(txt)

## choose the members
  df = df[df['member']==True]

## theta and phi in radian
  th_s = np.pi/2. - unp.uarray((df['dec'],df['dec_error']))*np.pi/180.
  ph_s   = unp.uarray((df['ra'],df['ra_error']))*np.pi/180.

## time derivative of angles in 
#th_ds = 

## we define our coordinate system such that its z is toward the dwarf
## the x hat is in theta hat direction at the mean of the angles theta, phi
## the y hat is in phi hat direction at the mean of the angels theta, phi
## first find the mean of the angles
  th_m = th_s.mean()
  ph_m = ph_s.mean()

  print("Plotting the histograms of the co-moving frame ...")
## make sure v_x ~ v_theta, and v_y ~ v_phi, i.e. the other components are negligible (here assumption is v_r, v_theta, and v_phi are comparable in values)
  temp_arr = [
  np.dot(th_hat(th_m,ph_m), ph_hat(th_s,ph_s))/np.dot(th_hat(th_m,ph_m), th_hat(th_s,ph_s)),
  np.dot(th_hat(th_m,ph_m), r_hat(th_s,ph_s))/np.dot(th_hat(th_m,ph_m), th_hat(th_s,ph_s)),
  np.dot(ph_hat(th_m,ph_m), th_hat(th_s,ph_s))/np.dot(ph_hat(th_m,ph_m), ph_hat(th_s,ph_s)),
  np.dot(ph_hat(th_m,ph_m), r_hat(th_s,ph_s))/np.dot(ph_hat(th_m,ph_m), ph_hat(th_s,ph_s))
  ]

  FigLabel_arr = [
  'PhiThm_ThThm.pdf',
  'RThm_ThThm.pdf',
  'ThPhim_PhiPhim.pdf',
  'RPhim_PhiPhim.pdf'
  ]
  xlabel_arr = [
  r"$|(\hat{\phi}\cdot\hat{\theta}_m)/(\hat{\theta}\cdot\hat{\theta}_m)|$",
  r"$|(\hat{r}\cdot\hat{\theta}_m)/(\hat{\theta}\cdot\hat{\theta}_m)|$",
  r"$|(\hat{\theta}\cdot\hat{\phi}_m)/(\hat{\phi}\cdot\hat{\phi}_m)|$",
  r"$|(\hat{r}\cdot\hat{\phi}_m)/(\hat{\phi}\cdot\hat{\phi}_m)|$"
  ]
  for i, temp in enumerate(temp_arr):
    fig, ax = plt.subplots()
    plt.locator_params(axis="x", nbins=4)
    ax.hist(np.abs([a.n for a in temp]), bins=5, color="grey")
    ax.hist(np.abs([a.n for a in temp]), bins=5, histtype="step",color="k")
    ax.set_yscale("log")
    ax.set_xlabel(xlabel_arr[i], fontdict=font)
    ax.set_ylabel(r"Number of stars", fontdict=font)
    ax.tick_params(labelsize=20)
    plt.tight_layout()
    plt.show()
    fig.savefig(FigLabel_arr[i])

  print("Saved the plots in the local directory.") 

  print("Computing x and y of stars and their errors ...")

## x = \vec{r}\cdot\hat{x}. Here assumption is that D is very large
## unit is pc
  x_arr = D*np.dot(th_hat(th_m,ph_m), r_hat(th_s,ph_s))
  df['x']     = [a.n for a in x_arr]
  df['x_err'] = [a.s for a in x_arr]

## y = \vec{r}\cdot\hat{y}
  y_arr = D*np.dot(ph_hat(th_m,ph_m), r_hat(th_s,ph_s))
  df['y']     = [a.n for a in y_arr]
  df['y_err'] = [a.s for a in y_arr]

  print("Computing the velocites of the stars and their errors ...")
## vx in units of pc/yr 
## Each arc-second=4.8481e-6 radian => mas = 4.8481e-9 rad
## unit of df['pmdec'] is mas/yr
  vx_arr = -D*(unp.uarray((df['pmdec'],df['pmdec_error']))*(4.8481e-9))
  df['vx']     = [a.n for a in vx_arr-vx_arr.mean()]
  df['vx_err'] = [a.s for a in vx_arr-vx_arr.mean()]

## vy in units of pc/yr
  vy_arr = D*(unp.uarray((df['pmra'],df['pmra_error']))*(4.8481e-9))
  df['vy']     = [a.n for a in vy_arr-vy_arr.mean()] 
  df['vy_err'] = [a.s for a in vy_arr-vy_arr.mean()]

## store the results
  df_result = pd.DataFrame({
    'x (pc)':df['x'], 'x_err':df['x_err'], 
    'y (pc)':df['y'], 'y_err':df['y_err'], 
    'vx (pc/yr)':df['vx'], 'vx_err':df['vx_err'],
    'vy (pc/yr)':df['vy'], 'vy_err':df['vy_err']
    })
  df_result.to_csv("Data.csv",index=None)


  print("Output file Data.csv stored in the local directory!")
