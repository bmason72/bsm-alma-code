
import scipy.linalg as spla
import scipy as sp
import testvis as tv

nu=100.0
lam=3e8/(nu*1e9)

# this will return the uv baselines in inverse radians-
bl=tv.getbaselines('Alma_cycle1_2.cfg.txt',lam=lam)

# beam fwhm in rad- (B3 alma)
beamfwhm=2.91e-4
parvec=sp.array([0.01,0.0,0.0,2e-5,2e-5,0.0])
# param order is - norm,l0,m0,fwhm_1,fwhm_2,axis_angle
parvec += 1e-9

print tv.visval(beamfwhm,0,0,parvec[0],parvec[1],parvec[2],parvec[3],parvec[4],parvec[5],brute_force=False)
# this is a typical visibility value. set the error in make_fisher_mx to something of this order

signal_fwhm= (sp.arange(100)*0.02+0.02) * beamfwhm
norm_snr=sp.zeros(100)
fwhm_snr=sp.zeros(100)
for i in range(100):
	parvec[3] = signal_fwhm[i] 
	parvec[4] = signal_fwhm[i]
	f=tv.make_fisher_mx(bl,1e-12,0.01,beamfwhm,parvec,brute_force=False)
	finv=spla.inv(f)
	norm_snr[i]= parvec[0] / (finv[0,0])**0.5
	fwhm_snr[i]= parvec[3] / (finv[3,3])**0.5

pl.plot(signal_fwhm/beamfwhm,fwhm_snr)
pl.plot(signal_fwhm/beamfwhm,norm_snr,'-r')

##save comptued array 2 values for comparison-
#norm_snr_2=sp.copy(norm_snr)
#fwhm_snr_2=sp.copy(fwhm_snr)

# check the visibility values
#  integral of this should be pi*sigma_beam^2 = pi*(0.1/2.3548)^2
u=0
v=0
l0=0
m0=0
beamfwhm=0.1
fwhm_1=beamfwhm
fwhm_2=beamfwhm
axis_angle=0.0
norm=1.0
print tv.visval(beamfwhm,u,v,norm,l0,m0,fwhm_1,fwhm_2,axis_angle,brute_force=True)
print tv.visval(beamfwhm,u,v,norm,l0,m0,fwhm_1,fwhm_2,axis_angle,brute_force=False)

# with v=0,u=0.5/sig_a the right answer is exp(-pi^2/2)
u=0.5/(beamfwhm/2.3548)
print tv.visval(beamfwhm,u,v,norm,l0,m0,fwhm_1,fwhm_2,axis_angle,brute_force=True)
print tv.visval(beamfwhm,u,v,norm,l0,m0,fwhm_1,fwhm_2,axis_angle,brute_force=False)
# with v=u, the right answer is exp(-pi^2)
v=0.5/(beamfwhm/2.3548)
print tv.visval(beamfwhm,u,v,norm,l0,m0,fwhm_1,fwhm_2,axis_angle,brute_force=True)
print tv.visval(beamfwhm,u,v,norm,l0,m0,fwhm_1,fwhm_2,axis_angle,brute_force=False)
# the analytic form gives these results but the brute_force integral is having
#   problems. maybe b/c it is oscillatory. but it is only for checking anyway! ...
