from sklearn.model_selection import GridSearchCV
from sklearn.neighbors import KernelDensity
import uncertainties.unumpy as unp
import time
import Constants as co
import pandas as pd
import numpy as np
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


def EstimateDM(dSph):
  doPlot = False


## import data  
## (x, y, vx, vy) and their errors are the columns & each row should be one star
  txt = "%s_Data.csv"%(dSph)
  df = pd.read_csv(txt)
## number of stars and the dimension of phase-space
  Npar, dim = df.shape
  dim/=2

  if dim != 4:
    raise(Exception("dimension of the data should be 4"))


## a summary of data 
  print("#############\ndata summary:\n%s\n#############"%(df.describe()))

  if doPlot:
    alpha=1
    xlim=1000
    fig, axs = plt.subplots(1,2)
    ax = axs[0]
    cb=ax.scatter(df['x (pc)'],df['y (pc)'], c=df['vx (pc/yr)']*10000., marker='.',alpha=alpha, cmap=plt.cm.seismic)
    fig.colorbar(cb,ax=ax)
    ax.set_xlim(-xlim,xlim)
    ax.set_ylim(-xlim,xlim)
    ax = axs[1]
    cb=ax.scatter(df['x (pc)'],df['y (pc)'], c=df['vy (pc/yr)']*10000., marker='.',alpha=alpha, cmap=plt.cm.seismic)
    fig.colorbar(cb,ax=ax)
    ax.set_xlim(-xlim,xlim)
    ax.set_ylim(-xlim,xlim)
    plt.tight_layout()
    plt.show()




  if doPlot:
    for i in range(dim):
      for j in range(dim):
        col1 = df.columns[i]
        col2 = df.columns[j]
        fig, ax = plt.subplots()
        ax.scatter(df[col1]-np.mean(df[col1]),df[col2]-np.mean(df[col2]),marker='.',color='k',alpha=0.3)
        ax.set_xlabel(df.columns[i])
        ax.set_ylabel(df.columns[j])
        plt.tight_layout()
        plt.show()



## rotate the axes around the z direction such that the correlation between vx-vy is the least. 
## this means that later on the xy component of the presure will be negligible
  RotAng = 0*np.pi/180.
  Corr   = 1.
  for ang in np.arange(0, np.pi, np.pi/180*0.1):
    ## rotate the axes so that the system is symmetric around x axis
    x  =  np.cos(ang)*df['x (pc)'] + np.sin(ang)*df['y (pc)']
    y  = -np.sin(ang)*df['x (pc)'] + np.cos(ang)*df['y (pc)']
    vx =  np.cos(ang)*df['vx (pc/yr)'] + np.sin(ang)*df['vy (pc/yr)']
    vy = -np.sin(ang)*df['vx (pc/yr)'] + np.cos(ang)*df['vy (pc/yr)']
    df_rot = pd.DataFrame({'x':x, 'y':y, 'vx':vx, 'vy':vy})
    cor = np.corrcoef(df_rot['vx'],df_rot['vy'])
    if abs(cor[0,1]) < Corr:
      Corr = abs(cor[0,1])
      RotAng = ang
  print("ang: %g\nCorr: %s\n##"%(RotAng*180/np.pi,Corr))

## rotate the axes so that the system is symmetric around x axis
  x_ =  np.cos(RotAng)*df['x (pc)'] + np.sin(RotAng)*df['y (pc)']
  y_ = -np.sin(RotAng)*df['x (pc)'] + np.cos(RotAng)*df['y (pc)']
  vx =  np.cos(RotAng)*df['vx (pc/yr)'] + np.sin(RotAng)*df['vy (pc/yr)']
  vy = -np.sin(RotAng)*df['vx (pc/yr)'] + np.cos(RotAng)*df['vy (pc/yr)']
  df_rot = pd.DataFrame({'x':x_, 'y':y_, 'vx':vx, 'vy':vy})

## for error propagation 
  x__ =  np.cos(RotAng)*unp.uarray((df['x (pc)'],df['x_err'])) + np.sin(RotAng)*unp.uarray((df['y (pc)'],df['y_err']))
  y__ = -np.sin(RotAng)*unp.uarray((df['x (pc)'],df['x_err'])) + np.cos(RotAng)*unp.uarray((df['y (pc)'],df['y_err']))
  vx_ =  np.cos(RotAng)*unp.uarray((df['vx (pc/yr)'],df['vx_err'])) + np.sin(RotAng)*unp.uarray((df['vy (pc/yr)'],df['vy_err']))
  vy_ = -np.sin(RotAng)*unp.uarray((df['vx (pc/yr)'],df['vx_err'])) + np.cos(RotAng)*unp.uarray((df['vy (pc/yr)'],df['vy_err']))



## the mean of each column
  mean = np.mean(df_rot,axis=0)

## scale (change units)
  xscale = np.max([abs(df_rot['x']),abs(df_rot['y'])])
  vscale = np.max([abs(df_rot['vx']),abs(df_rot['vy'])])
  scale = np.array([xscale,xscale,vscale,vscale])

## scale the data (Nxdim) 
  data = np.array(df_rot/scale)

  print("\n\nPlotting the position and velocity of stars:")
  i, j = 0, 1
  fig, ax = plt.subplots()
  ax.scatter(data[:,i]*scale[i],data[:,j]*scale[j],marker='*',color='k',alpha=0.3)
  xup = 1000. ## (pc)
  ax.set_xlim(-xup,xup)
  ax.set_ylim(-xup,xup)
  ax.set_xlabel(r"x (pc)", fontdict=font)
  ax.set_ylabel(r"y (pc)", fontdict=font)
  ax.tick_params(labelsize=20)
  plt.tight_layout()
  plt.show()
  txt = "%s_StarPositions.pdf"%(dSph)
  fig.savefig(txt)

## one year
  yr = 365*24*60*60.
  i, j = 2, 3
  fig, ax = plt.subplots()
  ax.scatter(data[:,i]*scale[i]*(co.parsec/yr)/1000.,data[:,j]*scale[j]*(co.parsec/yr)/1000.,marker='o',color='k',alpha=0.3)
  vup = 200 # km/s
  ax.set_xlim(-vup-100,vup+100)
  ax.set_ylim(-vup,vup)
  ax.set_xlabel(r"$v_x$ (km/s)", fontdict=font)
  ax.set_ylabel(r"$v_y$ (km/s)", fontdict=font)
  ax.tick_params(labelsize=20)
  plt.tight_layout()
  plt.show()
  txt = "%s_StarVelocities.pdf"%(dSph)
  fig.savefig(txt)


  print("\n\nPlotting the histograms of the positions and velocities of the stars:")
  xlabel_arr   = ['x (pc)', 'y (pc)', r'v$_x$ (km/s)', r'v$_y$ (km/s)']
  FigLabel_arr = ['Star_x_hist.pdf', 'Star_y_hist.pdf', 'Star_vx_hist.pdf', 'Star_vy_hist.pdf']
  for i in range(4):
    fig, ax = plt.subplots()
    plt.locator_params(axis="x", nbins=4)
    if i < 2:
      sc = scale[i]
    else:
      sc = scale[i]*(co.parsec/yr)/1000.
    ax.hist(data[:,i]*sc, bins=15, color="grey")
    ax.hist(data[:,i]*sc, bins=15, histtype="step",color="k")
    if i < 2:
      ax.set_yscale("log")
    if i == 3:
      ax.set_xlim(-203,250)
    ax.set_xlabel(xlabel_arr[i], fontdict=font)
    ax.set_ylabel(r"Number of stars", fontdict=font)
    ax.tick_params(labelsize=20)
    plt.tight_layout()
    plt.show()
    #fig.savefig(FigLabel_arr[i])




  if doPlot:
    for i in range(dim):
      for j in range(dim):
        if j < i: continue
        print("(%s,%s): %g"%(df_rot.columns[i],df_rot.columns[j],np.mean(data[:,i]*data[:,j])))
        fig, ax = plt.subplots()
        ax.scatter(data[:,i],data[:,j],marker='.',color='k',alpha=0.3)
        ax.set_xlabel(df_rot.columns[i])
        ax.set_ylabel(df_rot.columns[j])
        plt.show()

##analyze part of data only
#condition = (abs(data[:,0])<0.1)&(abs(data[:,1])<0.1)


  print("\n\nEstimating the phase-space density:")

## determine the best bandwidth
## use grid search cross-validation to optimize the bandwidth
  params = {'bandwidth': np.linspace(0.01, 0.2, 50)}
  grid = GridSearchCV(KernelDensity(), params, cv=5)
  grid.fit(data)
  print(grid.best_params_) ## result is 0.14
  print("resolution: %g (pc) %g (km/s)"%(grid.best_params_['bandwidth']*scale[0],grid.best_params_['bandwidth']*scale[2]*(co.parsec/yr)/1000.))


## use a gaussian kernel density esstimator to find single-body phase-space density
  kde   = KernelDensity(kernel='gaussian', bandwidth=grid.best_params_['bandwidth']).fit(data)



## lower and upper coordinates
  qmin = np.min(data,axis=0)
  qmax = np.max(data,axis=0)


#import pickle
#with open("needed.pickle","wb") as fi:  
#  pickle.dump([0, 0, Npar, scale, qmin, qmax, kde],fi) 




  Nbin = 500 ## the number of bins
  vx_arr = np.linspace(qmin[2], qmax[2], Nbin)
  vy_arr = np.linspace(qmin[3], qmax[3], Nbin)

## one year in seconds
  yr = 365*24*60*60
## the mass of each particle (star). The value is not important since it cancels in P/rho below
  m = 1.


  x = float(input("\n\nEstimate dark matter mass density at x = ? (pc)\n"))/scale[0]
  y = float(input("\n\nEstimate dark matter mass density at y = ? (pc)\n"))/scale[1]
  dx = 0.1/scale[0]#np.max(abs(data[:,0]))/1.e6
  dy = 0.1/scale[1]#np.max(abs(data[:,1]))/1.e6


  print("\n\nWorking on DM mass density estimation:  ")


## mesh of the 4D space
  qqs = np.meshgrid(x, y, vx_arr, vy_arr)
  vx  = qqs[2]
  vy  = qqs[3]
  D4space = np.vstack([qqs[0].ravel(), qqs[1].ravel(), qqs[2].ravel(), qqs[3].ravel()]).T

  start  = time.time()
  Log_f = kde.score_samples(D4space)  
  print("time: %g"%(time.time()-start))
  f     = np.exp(Log_f).reshape(qqs[0].shape)


## density and stress tensor
  rho = m*np.sum(f,axis=(2,3))

## assuming <vx> = <vy> = 0
  Pxx = m*np.sum(f*vx**2,axis=(2,3))
  Pyy = m*np.sum(f*vy**2,axis=(2,3))
  Pxy = m*np.sum(f*vx*vy,axis=(2,3))


  print("Finished computations at (x,y). Working on their errors ...")
#############################
## evaluate the error of f ##
#############################
## error due to measurement uncertainties
  rho_ = 0
  Pxx_ = 0
  Pyy_ = 0
  h_   = grid.best_params_['bandwidth']
  vx_arr_ = np.linspace(qmin[2], qmax[2], 10)
  vy_arr_ = np.linspace(qmin[3], qmax[3], 10)


  for q2 in vx_arr_*vscale:
    for q3 in vy_arr_*vscale:
      dist2_h2 = ((x__ - x)/(h_*xscale))**2 + ((y__ - y)/(h_*xscale))**2 +((vx_ - q2)/(h_*vscale))**2 +((vy_ - q3)/(h_*vscale))**2
      rho_ += (unp.exp(-dist2_h2)).sum()
      Pxx_ += (unp.exp(-dist2_h2)*q2**2).sum() 
      Pyy_ += (unp.exp(-dist2_h2)*q3**2).sum() 


###############################
## statistical error estimation
###############################
  NEnsemb = 10
  rho_Ens = np.ones(NEnsemb)*(-99)
  Pxx_Ens = np.ones(NEnsemb)*(-99)
  Pyy_Ens = np.ones(NEnsemb)*(-99) 
  for iEns in range(10):

    if iEns % 10 == 0:
      print(iEns)

## draw a sample of the ensemble
    data2 = kde.sample(Npar)

## use a gaussian kernel density esstimator to find single-body phase-space density
    kde2   = KernelDensity(kernel='gaussian', bandwidth=grid.best_params_['bandwidth']).fit(data2)

    Log_f = kde2.score_samples(D4space)  
    f     = np.exp(Log_f).reshape(qqs[0].shape)

    rho_Ens[iEns] = m*np.sum(f,axis=(2,3))
    Pxx_Ens[iEns] = m*np.sum(f*vx**2,axis=(2,3))
    Pyy_Ens[iEns] = m*np.sum(f*vy**2,axis=(2,3))

    print("iEns: %d"%(iEns))
    print(np.std(rho_Ens[:iEns+1]))
    print(np.std(Pxx_Ens[:iEns+1]))
    print(np.std(Pyy_Ens[:iEns+1]))

## plot stat error
  fig, ax = plt.subplots()
  ax.plot([np.std(rho_Ens[:iEns+1])/rho[0,0] for iEns in range(0,NEnsemb)],label=r"$\delta\rho/\rho$")
  ax.plot([np.std(Pxx_Ens[:iEns+1])/Pxx[0,0] for iEns in range(0,NEnsemb)],label=r"$\delta$P$_{xx}$/P$_{xx}$")
  ax.plot([np.std(Pyy_Ens[:iEns+1])/Pyy[0,0] for iEns in range(0,NEnsemb)],label=r"$\delta$P$_{yy}$/P$_{yy}$")
  ax.set_xlabel("Ensemble size", fontdict=font)
  ax.set_ylabel("std", fontdict=font)
  ax.tick_params(labelsize=20)
  ax.legend(prop={'size':15})
  plt.show()

## Statistical error 
  rho_stat_err = np.std(rho_Ens)
  Pxx_stat_err = np.std(Pxx_Ens)
  Pyy_stat_err = np.std(Pyy_Ens)


  print("drho/rho=> systematic: %g statistic: %g"%(rho_.s/rho_.n, rho_stat_err/rho))
  print("dPxx/Pxx=> systematic: %g statistic: %g"%(Pxx_.s/Pxx_.n, Pxx_stat_err/Pxx))
  print("dPyy/Pyy=> systematic: %g statistic: %g"%(Pyy_.s/Pyy_.n, Pyy_stat_err/Pyy))
#############################
#############################
## x+dx

  print("Working on estimations at x+dx ...")

## mesh of the 4D space
  qqs = np.meshgrid(x+dx, y, vx_arr, vy_arr)
  vx  = qqs[2]
  vy  = qqs[3]
  D4space = np.vstack([qqs[0].ravel(), qqs[1].ravel(), qqs[2].ravel(), qqs[3].ravel()]).T

  start  = time.time()
  Log_f = kde.score_samples(D4space)
  print("time: %g"%(time.time()-start))
  f_dx     = np.exp(Log_f).reshape(qqs[0].shape)


## density and stress tensor
  rho_dx = m*np.sum(f_dx,axis=(2,3))

## assuming <vx> = <vy> = 0
  Pxx_dx = m*np.sum(f_dx*vx**2,axis=(2,3))
  Pyy_dx = m*np.sum(f_dx*vy**2,axis=(2,3))
  Pxy_dx = m*np.sum(f_dx*vx*vy,axis=(2,3))


## y+dy
  print("Working on estimations at y+dy ...")

## mesh of the 4D space
  qqs = np.meshgrid(x, y+dy, vx_arr, vy_arr)
  vx  = qqs[2]
  vy  = qqs[3]
  D4space = np.vstack([qqs[0].ravel(), qqs[1].ravel(), qqs[2].ravel(), qqs[3].ravel()]).T

  start  = time.time()
  Log_f = kde.score_samples(D4space)
  print("time: %g"%(time.time()-start))
  f_dy     = np.exp(Log_f).reshape(qqs[0].shape)


## density and stress tensor
  rho_dy = m*np.sum(f_dy,axis=(2,3))

## assuming <vx> = <vy> = 0
  Pxx_dy = m*np.sum(f_dy*vx**2,axis=(2,3))
  Pyy_dy = m*np.sum(f_dy*vy**2,axis=(2,3))
  Pxy_dy = m*np.sum(f_dy*vx*vy,axis=(2,3))



## x+2dx
  print("Working on estimations at x+2dx ...")

## mesh of the 4D space
  qqs = np.meshgrid(x+2.*dx, y, vx_arr, vy_arr)
  vx  = qqs[2]
  vy  = qqs[3]
  D4space = np.vstack([qqs[0].ravel(), qqs[1].ravel(), qqs[2].ravel(), qqs[3].ravel()]).T

  start  = time.time()
  Log_f = kde.score_samples(D4space)
  print("time: %g"%(time.time()-start))
  f_2dx     = np.exp(Log_f).reshape(qqs[0].shape)


## density and stress tensor
  rho_2dx = m*np.sum(f_2dx,axis=(2,3))

## assuming <vx> = <vy> = 0
  Pxx_2dx = m*np.sum(f_2dx*vx**2,axis=(2,3))
  Pyy_2dx = m*np.sum(f_2dx*vy**2,axis=(2,3))
  Pxy_2dx = m*np.sum(f_2dx*vx*vy,axis=(2,3))


## y+2dy
  print("Working on estimations at y+2dy ...")

## mesh of the 4D space
  qqs = np.meshgrid(x, y+2.*dy, vx_arr, vy_arr)
  vx  = qqs[2]
  vy  = qqs[3]
  D4space = np.vstack([qqs[0].ravel(), qqs[1].ravel(), qqs[2].ravel(), qqs[3].ravel()]).T

  start  = time.time()
  Log_f = kde.score_samples(D4space)
  print("time: %g"%(time.time()-start))
  f_2dy     = np.exp(Log_f).reshape(qqs[0].shape)


## density and stress tensor
  rho_2dy = m*np.sum(f_2dy,axis=(2,3))

## assuming <vx> = <vy> = 0
  Pxx_2dy = m*np.sum(f_2dy*vx**2,axis=(2,3))
  Pyy_2dy = m*np.sum(f_2dy*vy**2,axis=(2,3))
  Pxy_2dy = m*np.sum(f_2dy*vx*vy,axis=(2,3))


  print("Using the finite difference method to estimate the derivatives ...")
############################
## calculate the derivatives
## and their errors
############################
## combine the syst. and stat. errors
  err    = np.sqrt( (rho*rho_.s/rho_.n)**2 + rho_stat_err**2 )
  err_dx = np.sqrt( (rho_dx*rho_.s/rho_.n)**2 +  rho_stat_err**2)
  corr   = 1. ## for variables aprat by dx the correlation is ~1. 
  diff_err = np.sqrt(err**2 + err_dx**2 - 2*corr*err*err_dx)
  dx_rho = unp.uarray((rho_dx-rho, diff_err))/dx
##
  err      = np.sqrt( (Pxx*Pxx_.s/Pxx_.n)**2 + Pxx_stat_err**2 )
  err_dx   = np.sqrt( (Pxx_dx*Pxx_.s/Pxx_.n)**2 + Pxx_stat_err**2 )
  diff_err = np.sqrt(err**2 + err_dx**2 - 2*corr*err*err_dx)
  dx_Pxx   = unp.uarray((Pxx_dx - Pxx, diff_err))/dx
##
  err      = np.sqrt( (rho*rho_.s/rho_.n)**2 + rho_stat_err**2 )
  err_dy   = np.sqrt( (rho_dy*rho_.s/rho_.n)**2 + rho_stat_err**2 )
  corr     = 1. ## for variables aprat by dx the correlation is ~1.
  diff_err = np.sqrt(err**2 + err_dy**2 - 2*corr*err*err_dy)
  dy_rho   = unp.uarray((rho_dy-rho, diff_err))/dy
##
  err      = np.sqrt( (Pyy*Pyy_.s/Pyy_.n)**2 + Pyy_stat_err**2 )
  err_dy   = np.sqrt( (Pyy_dy*Pyy_.s/Pyy_.n)**2 + Pyy_stat_err**2 )
  diff_err = np.sqrt(err**2 + err_dy**2 - 2*corr*err*err_dy)
  dy_Pyy   = unp.uarray((Pyy_dy - Pyy, diff_err))/dy
##
  err      = np.sqrt( (Pxx_dx*Pxx_.s/Pxx_.n)**2 + Pxx_stat_err**2 )
  err_dx   = np.sqrt( (Pxx_2dx*Pxx_.s/Pxx_.n)**2 + Pxx_stat_err**2 )
  diff_err = np.sqrt(err**2 + err_dx**2 - 2*corr*err*err_dx) ## error on difference (P_xx(x+2dx) - P_xx(x+dx))
  err1 = diff_err/dx ## error on dP_xx/dx at x+dx
  err2 = dx_Pxx[0,0].s ## erro on dP_xx/dx at x
  diff_err2 = np.sqrt(err1**2 + err2**2 - 2*corr*err1*err2)
  dxx_Pxx = unp.uarray(( (Pxx_2dx-Pxx_dx)/dx-dx_Pxx[0,0].n, diff_err2))/dx
##
  err      = np.sqrt( (Pyy_dy*Pyy_.s/Pyy_.n)**2 + Pyy_stat_err**2 )
  err_dy   = np.sqrt( (Pyy_2dy*Pyy_.s/Pyy_.n)**2 + Pyy_stat_err**2 )
  diff_err = np.sqrt(err**2 + err_dy**2 - 2*corr*err*err_dy) ## error on difference (P_yy(y+2dy) - P_yy(y+dy))
  err1 = diff_err/dy ## error on dP_yy/dy at y+dy
  err2 = dy_Pyy[0,0].s ## erro on dP_yy/dy at y
  diff_err2 = np.sqrt(err1**2 + err2**2 - 2*corr*err1*err2)
  dyy_Pyy = unp.uarray(( (Pyy_2dy-Pyy_dy)/dy-dy_Pyy[0,0].n, diff_err2))/dy


  print("dx_rho: %s \ndy_rho: %s"%(dx_rho[0,0],dy_rho[0,0]))
  print("dx_Pxx: %s \ndy_Pyy: %s"%(dx_Pxx[0,0],dy_Pyy[0,0]))
  print("dxx_Pxx: %s \ndyy_Pyy: %s"%(dxx_Pxx[0,0],dyy_Pyy[0,0]))


## change the units
## in SI
  y_SI = y*scale[1]*co.parsec
## Due to the later cancelation, we don't change \int f dvx dvy. Only the v*v is converted to SI
  Pxx_SI     = unp.uarray((Pxx, np.sqrt( (Pxx*Pxx_.s/Pxx_.n)**2 + Pxx_stat_err**2 ) )) *scale[2]**2 *(co.parsec/yr)**2
  Pyy_SI     = unp.uarray((Pyy, np.sqrt( (Pyy*Pyy_.s/Pyy_.n)**2 + Pyy_stat_err**2 ) )) *scale[3]**2 *(co.parsec/yr)**2
  dx_rho_SI  = dx_rho /(scale[0]*co.parsec)
  dy_rho_SI  = dy_rho /(scale[1]*co.parsec)
  dx_Pxx_SI  = dx_Pxx /(scale[0]*co.parsec) *scale[2]**2 *(co.parsec/yr)**2 
  dy_Pyy_SI  = dy_Pyy /(scale[1]*co.parsec) *scale[3]**2 *(co.parsec/yr)**2 
  dxx_Pxx_SI = dxx_Pxx / (scale[0]*co.parsec)**2 *scale[2]**2 *(co.parsec/yr)**2
  dyy_Pyy_SI = dyy_Pyy / (scale[1]*co.parsec)**2 *scale[3]**2 *(co.parsec/yr)**2

#
  DoubleDiv = 1./(y_SI*(unp.uarray((rho, np.sqrt( (rho*rho_.s/rho_.n)**2 + rho_stat_err**2) )))**2)*(
  -y_SI*dx_Pxx_SI*dx_rho_SI
  -y_SI*dy_Pyy_SI*dy_rho_SI
  +rho*(y_SI*dxx_Pxx_SI+dy_Pyy_SI+y_SI*dyy_Pyy_SI)
  )


  rhoDM_SI = -1./(4.*np.pi*co.G)*DoubleDiv # -rho*UnitConversion_normalization  





  print("DM rho: %s (kg/m^3)"%(rhoDM_SI[0,0]))


#### discritize dimensions 
##Nbin = 20
##qs = []
##for i in range(int(dim)):
##  if i < 2:
##    qs.append(np.linspace(-0.106, 0.106,Nbin))
##  else:
##    qs.append(np.linspace(qmin[i], qmax[i], Nbin))
##
##
#### mesh of the 4D space
##qqs = np.meshgrid(qs[0], qs[1], qs[2], qs[3])
##D4space = np.vstack([qqs[0].ravel(), qqs[1].ravel(), qqs[2].ravel(), qqs[3].ravel()]).T
##
##start  = time.time()
##Log_f = kde.score_samples(D4space)  
##print("time: %g"%(time.time()-start))
##f     = np.exp(Log_f)
##
##
##from Plot_f import *
##Plot4D(f.reshape(qqs[0].shape),qqs,Nbin) 
##
##
##import pickle
##with open("f_qs_ForPlot.pickle","wb") as fi:
##  pickle.dump([scale, Nbin, qs, qqs, f], fi)




