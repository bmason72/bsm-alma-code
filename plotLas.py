
import pylab as pl
import scipy as sp
import testvis as tv

# hi there this is a test --

i_i=0
i_theta=1
i_pb=2
i_snr=3
i_fwhmsnr=4
i_fwhm2snr=5
i_anglesnr=6
i_poserr=7
i_pos2err=8
i_zeros=9

c1=pl.genfromtxt('c40-1n.cfgnew-constFlux_2.parErrs.txt')
c2=pl.genfromtxt('c40-2n.cfgnew-constFlux_2.parErrs.txt')
c3=pl.genfromtxt('c40-3n.cfgnew-constFlux_2.parErrs.txt')
c4=pl.genfromtxt('c40-4n.cfgnew-constFlux_2.parErrs.txt')
c5=pl.genfromtxt('c40-5n.cfgnew-constFlux_2.parErrs.txt')
c6=pl.genfromtxt('c40-6n.cfgnew-constFlux_2.parErrs.txt')
c7=pl.genfromtxt('c40-7n2.cfgnew-constFlux_2.parErrs.txt')
c8=pl.genfromtxt('c40-8n2.cfgnew-constFlux_2.parErrs.txt')
c9=pl.genfromtxt('c40-9n2.cfgnew-constFlux_2.parErrs.txt')

lam=3e8/(100e9)
b1=tv.getbaselines('c40-1n.cfg',lam=lam)
b2=tv.getbaselines('c40-2n.cfg',lam=lam)
b3=tv.getbaselines('c40-3n.cfg',lam=lam)
b4=tv.getbaselines('c40-4n.cfg',lam=lam)
b5=tv.getbaselines('c40-5n.cfg',lam=lam)
b6=tv.getbaselines('c40-6n.cfg',lam=lam)
b7=tv.getbaselines('c40-7n2.cfg',lam=lam)
b8=tv.getbaselines('c40-8n2.cfg',lam=lam)
b9=tv.getbaselines('c40-9n2.cfg',lam=lam)

sbfac = sp.ones(c1[:,0].size)
# or
sbfac = sp.ones(c1[:,0].size) * (c1[:,i_theta]/c1[5,i_pb])**2 * 10

##########################
#
# Begin SNR plot (const flux)
#

i_plot = i_snr
ymax=20

i_plot=i_fwhmsnr
ymax=5

pl.close()
pl.ion()
pl.plot(c1[:,i_theta]/c1[:,i_pb],c1[:,i_plot],'b')
pl.plot(c2[:,i_theta]/c2[:,i_pb],c2[:,i_plot],'g')
pl.plot(c3[:,i_theta]/c3[:,i_pb],c3[:,i_plot],'r')
pl.plot(c4[:,i_theta]/c4[:,i_pb],c4[:,i_plot],'c')
pl.plot(c5[:,i_theta]/c5[:,i_pb],c5[:,i_plot],'m')
pl.plot(c6[:,i_theta]/c6[:,i_pb],c6[:,i_plot],'y')
pl.plot(c7[:,i_theta]/c7[:,i_pb],c7[:,i_plot],'k')
#pl.plot(c8[:,i_theta]/c8[:,i_pb],c8[:,i_plot],'b-')
pl.plot(c9[:,i_theta]/c9[:,i_pb],c9[:,i_plot],'g-')
## or -
pl.semilogy(c1[:,i_theta]/c1[:,i_pb],c1[:,i_plot],'b')
pl.semilogy(c2[:,i_theta]/c2[:,i_pb],c2[:,i_plot],'g')
pl.semilogy(c3[:,i_theta]/c3[:,i_pb],c3[:,i_plot],'r')
pl.semilogy(c4[:,i_theta]/c4[:,i_pb],c4[:,i_plot],'c')
pl.semilogy(c5[:,i_theta]/c5[:,i_pb],c5[:,i_plot],'m')
pl.semilogy(c6[:,i_theta]/c6[:,i_pb],c6[:,i_plot],'y')
pl.semilogy(c7[:,i_theta]/c7[:,i_pb],c7[:,i_plot],'k')
#pl.semilogy(c8[:,i_theta]/c8[:,i_pb],c8[:,i_plot],'b-')
pl.semilogy(c9[:,i_theta]/c9[:,i_pb],c9[:,i_plot],'g-')
#pl.loglog(c9[:,i_theta]/c9[:,i_pb],c9[:,i_plot],'g-')
# limits and labels-
pl.xlim([1e-4,0.6])
pl.ylim([0.1,ymax])
pl.xlabel('FWHMsrc/FWHM_PB')
pl.title('Constant Flux')
pl.legend(['C40-1n','C40-2n','C40-3n','C40-4n','C40-5n','C40-6n','C40-7n2','C40-9n2'])

# oplot the nominal values for comparison
#  LAS = 0.5 lambda/b_x  x={min,10%,15%)
#      = 0.5 / q_x  (q=b/lam)

xval=0.5/(b1['qspecial'])[0]/c1[5,i_pb]
pl.semilogy( [xval,xval],[0.1,ymax],'b-.')
xval=0.5/(b1['qspecial'])[1]/c1[5,i_pb]
pl.semilogy([xval,xval],[0.1,ymax],'b-')

xval=0.5/(b2['qspecial'])[0]/c1[5,i_pb]
pl.semilogy( [xval,xval],[0.1,ymax],'g-.')
xval=0.5/(b2['qspecial'])[1]/c1[5,i_pb]
pl.semilogy([xval,xval],[0.1,ymax],'g-')

xval=0.5/(b3['qspecial'])[0]/c1[5,i_pb]
pl.semilogy( [xval,xval],[0.1,ymax],'r-.')
xval=0.5/(b3['qspecial'])[1]/c1[5,i_pb]
pl.semilogy([xval,xval],[0.1,ymax],'r-')

xval=0.5/(b4['qspecial'])[0]/c1[5,i_pb]
pl.semilogy( [xval,xval],[0.1,ymax],'c-.')
xval=0.5/(b4['qspecial'])[1]/c1[5,i_pb]
pl.semilogy([xval,xval],[0.1,ymax],'c-')

xval=0.5/(b5['qspecial'])[0]/c1[5,i_pb]
pl.semilogy( [xval,xval],[0.1,ymax],'m-.')
xval=0.5/(b5['qspecial'])[1]/c1[5,i_pb]
pl.semilogy([xval,xval],[0.1,ymax],'m-')

xval=0.5/(b6['qspecial'])[0]/c1[5,i_pb]
pl.semilogy( [xval,xval],[0.1,ymax],'y-.')
xval=0.5/(b6['qspecial'])[1]/c1[5,i_pb]
pl.semilogy([xval,xval],[0.1,ymax],'y-')

xval=0.5/(b7['qspecial'])[0]/c1[5,i_pb]
pl.semilogy( [xval,xval],[0.1,ymax],'k-.')
xval=0.5/(b7['qspecial'])[1]/c1[5,i_pb]
pl.semilogy([xval,xval],[0.1,ymax],'k-')

xval=0.5/(b9['qspecial'])[0]/c1[5,i_pb]
pl.semilogy( [xval,xval],[0.1,ymax],'g-.')
xval=0.5/(b9['qspecial'])[1]/c1[5,i_pb]
pl.semilogy([xval,xval],[0.1,ymax],'g-')

# SNR plot -
pl.text(0.19,10,'$0.5 \lambda/b_{10}$')
pl.text(0.35,10,'$0.5 \lambda/b_{min}$')
pl.ylabel('Ssrc/Uncertainty(Ssrc)')

# SNRfwhm plot -
pl.text(0.19,4,'$0.5 \lambda/b_{10}$')
pl.text(0.35,4,'$0.5 \lambda/b_{min}$')
pl.ylabel('FWHM/Uncertainty(FWHM)')

#
# End SNR plot
#####################

#################
# begin const SB plot


i_plot = i_snr
ymax=1

i_plot=i_fwhmsnr
ymax=1

ymin=0.001

pl.close()
pl.ion()
## or -
sbfac = sp.ones(c1[:,0].size) * (c1[:,i_theta]/c1[5,i_pb])**2 * 10
pl.semilogy(c1[:,i_theta]/c1[:,i_pb],sbfac*c1[:,i_plot],'b')
sbfac = sp.ones(c1[:,0].size) * (c2[:,i_theta]/c1[5,i_pb])**2 * 10
pl.semilogy(c2[:,i_theta]/c2[:,i_pb],sbfac*c2[:,i_plot],'g')
sbfac = sp.ones(c1[:,0].size) * (c3[:,i_theta]/c1[5,i_pb])**2 * 10
pl.semilogy(c3[:,i_theta]/c3[:,i_pb],sbfac*c3[:,i_plot],'r')
sbfac = sp.ones(c1[:,0].size) * (c4[:,i_theta]/c1[5,i_pb])**2 * 10
pl.semilogy(c4[:,i_theta]/c4[:,i_pb],sbfac*c4[:,i_plot],'c')
sbfac = sp.ones(c1[:,0].size) * (c5[:,i_theta]/c1[5,i_pb])**2 * 10
pl.semilogy(c5[:,i_theta]/c5[:,i_pb],sbfac*c5[:,i_plot],'m')
sbfac = sp.ones(c1[:,0].size) * (c6[:,i_theta]/c1[5,i_pb])**2 * 10
pl.semilogy(c6[:,i_theta]/c6[:,i_pb],sbfac*c6[:,i_plot],'y')
sbfac = sp.ones(c1[:,0].size) * (c7[:,i_theta]/c1[5,i_pb])**2 * 10
pl.semilogy(c7[:,i_theta]/c7[:,i_pb],sbfac*c7[:,i_plot],'k')
sbfac = sp.ones(c1[:,0].size) * (c8[:,i_theta]/c1[5,i_pb])**2 * 10
#pl.semilogy(c8[:,i_theta]/c8[:,i_pb],sbfac*c8[:,i_plot],'b-')
sbfac = sp.ones(c1[:,0].size) * (c9[:,i_theta]/c1[5,i_pb])**2 * 10
pl.semilogy(c9[:,i_theta]/c9[:,i_pb],sbfac*c9[:,i_plot],'g-')
#pl.loglog(c9[:,i_theta]/c9[:,i_pb],sbfac*c9[:,i_plot],'g-')
# limits and labels-
#pl.xlim([1e-4,0.6])
pl.ylim([ymin,ymax])
pl.xlabel('FWHMsrc/FWHM_PB')
pl.title('Constant Surface Brightness')
pl.legend(['C40-1n','C40-2n','C40-3n','C40-4n','C40-5n','C40-6n','C40-7n2','C40-9n2'])

# oplot the nominal values for comparison
#  LAS = 0.5 lambda/b_x  x={min,10%,15%)
#      = 0.5 / q_x  (q=b/lam)

xval=0.5/(b1['qspecial'])[0]/c1[5,i_pb]
pl.semilogy( [xval,xval],[ymin,ymax],'b-.')
xval=0.5/(b1['qspecial'])[1]/c1[5,i_pb]
pl.semilogy([xval,xval],[ymin,ymax],'b-')

xval=0.5/(b2['qspecial'])[0]/c1[5,i_pb]
pl.semilogy( [xval,xval],[ymin,ymax],'g-.')
xval=0.5/(b2['qspecial'])[1]/c1[5,i_pb]
pl.semilogy([xval,xval],[ymin,ymax],'g-')

xval=0.5/(b3['qspecial'])[0]/c1[5,i_pb]
pl.semilogy( [xval,xval],[ymin,ymax],'r-.')
xval=0.5/(b3['qspecial'])[1]/c1[5,i_pb]
pl.semilogy([xval,xval],[ymin,ymax],'r-')

xval=0.5/(b4['qspecial'])[0]/c1[5,i_pb]
pl.semilogy( [xval,xval],[ymin,ymax],'c-.')
xval=0.5/(b4['qspecial'])[1]/c1[5,i_pb]
pl.semilogy([xval,xval],[ymin,ymax],'c-')

xval=0.5/(b5['qspecial'])[0]/c1[5,i_pb]
pl.semilogy( [xval,xval],[ymin,ymax],'m-.')
xval=0.5/(b5['qspecial'])[1]/c1[5,i_pb]
pl.semilogy([xval,xval],[ymin,ymax],'m-')

xval=0.5/(b6['qspecial'])[0]/c1[5,i_pb]
pl.semilogy( [xval,xval],[ymin,ymax],'y-.')
xval=0.5/(b6['qspecial'])[1]/c1[5,i_pb]
pl.semilogy([xval,xval],[ymin,ymax],'y-')

xval=0.5/(b7['qspecial'])[0]/c1[5,i_pb]
pl.semilogy( [xval,xval],[ymin,ymax],'k-.')
xval=0.5/(b7['qspecial'])[1]/c1[5,i_pb]
pl.semilogy([xval,xval],[ymin,ymax],'k-')

xval=0.5/(b9['qspecial'])[0]/c1[5,i_pb]
pl.semilogy( [xval,xval],[ymin,ymax],'g-.')
xval=0.5/(b9['qspecial'])[1]/c1[5,i_pb]
pl.semilogy([xval,xval],[ymin,ymax],'g-')

# SNR plot -
pl.text(0.25,0.7,'$0.5 \lambda/b_{5}$')
pl.text(0.35,0.7,'$0.5 \lambda/b_{min}$')
pl.ylabel('Ssrc/Uncertainty(Ssrc)')

# SNRfwhm plot -
pl.text(0.25,0.7,'$0.5 \lambda/b_{5}$')
pl.text(0.35,0.7,'$0.5 \lambda/b_{min}$')
pl.ylabel('FWHM/Uncertainty(FWHM)')

##################

#########################
# new ones
#  snr good except for when it gets to very low values - then jumps way up
#  fwhmsnr good except for minor glitches at small theta
#   -seems promising for imaging criterion
#  anglesnr seems just goofy


pl.close()
pl.plot(c1[:,i_theta]/c1[:,i_pb],c1[:,i_plot]*(c1[:,i_theta]/c1[:,i_pb]/0.1)**2,'b')
pl.plot(c2[:,i_theta]/c2[:,i_pb],c2[:,i_plot]*(c2[:,i_theta]/c2[:,i_pb]/0.1)**2,'g')
pl.plot(c3[:,i_theta]/c3[:,i_pb],c3[:,i_plot]*(c3[:,i_theta]/c3[:,i_pb]/0.1)**2,'r')
pl.plot(c4[:,i_theta]/c4[:,i_pb],c4[:,i_plot]*(c4[:,i_theta]/c4[:,i_pb]/0.1)**2,'c')
pl.plot(c5[:,i_theta]/c5[:,i_pb],c5[:,i_plot]*(c5[:,i_theta]/c5[:,i_pb]/0.1)**2,'m')
pl.plot(c6[:,i_theta]/c6[:,i_pb],c6[:,i_plot]*(c6[:,i_theta]/c6[:,i_pb]/0.1)**2,'y')
pl.plot(c7[:,i_theta]/c7[:,i_pb],c7[:,i_plot]*(c7[:,i_theta]/c7[:,i_pb]/0.1)**2,'k')

# SB SNR plots - these are very nearly congruent with each other! ie follow each other
#  closely when SB normalization plus peak renorm is used.
pl.close()
pl.plot(c1[:,i_theta]/c1[:,i_pb],c1[:,i_plot]*(c1[:,i_theta]/c1[:,i_pb]/0.1)**2/6.62,'b')
pl.plot(c2[:,i_theta]/c2[:,i_pb],c2[:,i_plot]*(c2[:,i_theta]/c2[:,i_pb]/0.1)**2/4.42,'g')
pl.plot(c3[:,i_theta]/c3[:,i_pb],c3[:,i_plot]*(c3[:,i_theta]/c3[:,i_pb]/0.1)**2/2.73,'r')
pl.plot(c4[:,i_theta]/c4[:,i_pb],c4[:,i_plot]*(c4[:,i_theta]/c4[:,i_pb]/0.1)**2/1.80,'c')
pl.plot(c5[:,i_theta]/c5[:,i_pb],c5[:,i_plot]*(c5[:,i_theta]/c5[:,i_pb]/0.1)**2/1.31,'m')
pl.plot(c6[:,i_theta]/c6[:,i_pb],c6[:,i_plot]*(c6[:,i_theta]/c6[:,i_pb]/0.1)**2,/0.97'y')
pl.plot(c7[:,i_theta]/c7[:,i_pb],c7[:,i_plot]*(c7[:,i_theta]/c7[:,i_pb]/0.1)**2,'k')
#pl.plot(c8[:,i_theta]/c8[:,i_pb],c8[:,i_plot]*(c8[:,i_theta]/c8[:,i_pb]/0.1)**2,'b-')
#pl.plot(c9[:,i_theta]/c9[:,i_pb],c9[:,i_plot]*(c9[:,i_theta]/c9[:,i_pb]/0.1)**2,'g-')
pl.xlim([0,0.8])

pl.close()
pl.semilogy(c1[:,i_theta]/c1[:,i_pb],c1[:,i_plot],'b')
pl.semilogy(c2[:,i_theta]/c2[:,i_pb],c2[:,i_plot],'g')
pl.semilogy(c3[:,i_theta]/c3[:,i_pb],c3[:,i_plot],'r')
pl.semilogy(c4[:,i_theta]/c4[:,i_pb],c4[:,i_plot],'c')
pl.semilogy(c5[:,i_theta]/c5[:,i_pb],c5[:,i_plot],'m')
pl.semilogy(c6[:,i_theta]/c6[:,i_pb],c6[:,i_plot],'y')
pl.semilogy(c7[:,i_theta]/c7[:,i_pb],c7[:,i_plot],'k')
pl.semilogy(c8[:,i_theta]/c8[:,i_pb],c8[:,i_plot],'b-')
pl.semilogy(c9[:,i_theta]/c9[:,i_pb],c9[:,i_plot],'g-')
pl.xlim([0,0.3])

