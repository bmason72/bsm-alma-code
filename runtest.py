
import scipy as sp
import testvis as tv

nu=100.0
lam=3e8/(nu*1e9)

# this will return the uv baselines in inverse radians-
bl=tv.getbaselines('Alma_cycle1_2.cfg.txt',lam=lam)

# beam fwhm in rad- (B3 alma)
beamfwhm=2.91e-4
parvec=sp.array([0.01,0.0,0.0,2e-5,2e-5,0.0])

f=tv.make_fisher_mx(bl,1e-4,1e-4,beamfwhm,parvec)

