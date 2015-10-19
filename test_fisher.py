
import scipy.linalg as spla
import scipy as sp
import testvis as tv

if do_constSurfBright:
    print "SURFACE BRIGHTNESS of components held fixed"
else:
    print "TOTAL FLUX density of components held fixed"

# need to set do_constSurfBright={True,False} externally...

# better pardelta - numerical derivative scheme

# ALMA B3 numbers
nu=100.0
lam=3e8/(nu*1e9)

cfg_files = ['../alma.C36-1n.cfg','../alma.C36-2n.cfg', \
             '../alma.C36-3n.cfg','../alma.C36-4n.cfg','../alma.C36-5n.cfg', \
             '../alma.C36-6n.cfg','../alma.C36-7n.cfg','../alma.C36-8n.cfg']
n_configs=len(cfg_files)
# resolution and LAS in arcsec from C3 THB - convert to radians - aca is 15", 42.8"
cfg_resolution=sp.array([3.4,1.8,1.2,0.7,0.5,0.3,0.1,0.08])*3.1415/180.0/3600.0
min_resolution = min(cfg_resolution)
cfg_las = sp.array([25.0,25.0,25.0,10.0,8.0,5.0,1.5,1.1])*3.1415/180.0/3600.0
# beam fwhm in rad- (B3 alma) - 1' fwhm = 2.91e-4 radians fwhm 
beamfwhm= sp.array([2.91e-4,2.91e-4,2.91e-4,2.91e-4,2.91e-4,2.91e-4,2.91e-4,])*(100.0/nu)
# TBD - the ACA one crashed for some reason :/

# get some typical #s...
# master normalization
master_norm=1.0
# if doing SB case, this will be the flux at 1" fwhm...

# if doing SB case, I believe the changes needed are
#  a) invoke fish w/flux_norm=False 
#  b) figure out appropriate typical_err. can keep master_norm=1.0 and just apply the size
#  c) interpret norm as the SB not the integrated flux.

# assuming 1 mJy typical error per visibility
typical_err = 1.0

# number of x-axis samples at which to evaluate the parameter uncertainties of interest-
n_samps=1000
# n_configs=1
# loop over configurations
#for j in range(n_configs):
j=0
parvec=sp.array([master_norm,0.0,0.0,cfg_resolution[j],cfg_resolution[j],0.0])
# param order is - norm,l0,m0,fwhm_1,fwhm_2,axis_angle
# initialize variables
norm_snr=sp.zeros(n_samps)
fwhm_snr=sp.zeros(n_samps)
pos_err=sp.zeros(n_samps)
# go from resolution/5 to LAS*2
#min_scale=cfg_resolution[j]*0.1
min_scale=min_resolution*0.05
max_scale=cfg_las[j]*2.5
# compute range of models to consider-
step_size=(max_scale-min_scale)/n_samps
signal_fwhm= (sp.arange(n_samps)*step_size+step_size) 
# this will return the uv baselines in inverse radians-
bl=tv.getbaselines(cfg_files[j],lam=lam)
# loop over gaussian component sizes-
print '***',j,step_size,max_scale,min_scale,cfg_files[j]
#for i in range(n_samps):
i=500
parvec[3] = signal_fwhm[i] 
parvec[4] = signal_fwhm[i]
if do_constSurfBright:
    # normalize total flux so that master_norm is the flux in mJy
    #  of a component with 1" FWHM - note signal_fwhm here
    #  is in radians
    parvec[0] = master_norm * (signal_fwhm[i] * 206264.8)**2 
# set default deltas for calculating the numerical derivative
#  default to 1% for nonzero params; 0.1xsynth beam for positions;
#  and 0.5 deg for the axis angle-
default_par_delta=sp.copy(parvec*0.01)
default_par_delta[1]=cfg_resolution[j]*0.1
default_par_delta[2]=cfg_resolution[j]*0.1
default_par_delta[5]=0.5
f=tv.make_fisher_mx(bl,typical_err,default_par_delta,beamfwhm[j],parvec,brute_force=False,flux_norm=True)
# use SVD pseudo inverse instead of direct
#  inverse for stability
finv0=spla.inv(f)
u,s,v=spla.svd(f)
s2=sp.zeros([s.size,s.size])
for i in range(s.size):
    if (s[i] > 0):
        s2[i,i] = 1.0/s[i]
finv = v.dot( s2.dot( sp.transpose(u) ) )
norm_snr[i]= parvec[0] / (finv[0,0])**0.5
# save position error = average 1D error-
pos_err[i]= 0.5 * (finv[1,1]**0.5 + finv[2,2]**0.5)
fwhm_snr[i]= parvec[3] / (finv[3,3])**0.5
fwhm_snr[i]= parvec[3] / (finv[3,3])**0.5
# save fisher mx i here

print '*****'
print f
print '*****'
print finv0
print u
print s
print v
print '*****'
print finv
print '*****'
print norm_snr[i],fwhm_snr[i],pos_err[i]

