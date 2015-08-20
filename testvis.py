
import scipy as sp
import scipy.integrate as spi

def vis_igrnd(l,m,beamfwhm,u,v,norm,l0,m0,fwhm_1,fwhm_2,axis_angle):
    '''
    visibility of gaussian source. args
    0,1: ell,emm coords (inverse u/v units)
    2,3: u,v (conjugate to ell,emm)
    4: norm
    5,6: l0,m0 (center of source) [units of ell,emm]
    7,8: axis1,axis2 of source [same units as ell,emm]
    9: axis_angle [deg.]
    10: primary beam fwhm [same units as ell,emm]
    '''
    beamfactor  = sp.exp( -0.5* (l**2+m**2) / (beamfwhm/2.3548)**2)
    new_l = (l-l0) * sp.cos(0.0174533*axis_angle) + (m-m0) * \
            sp.sin(0.0174533*axis_angle)
    new_m = (m-m0) * sp.cos(0.0174533*axis_angle) - (l-l0) * \
            sp.sin(0.0174533*axis_angle)
    modelfactor = sp.exp( -0.5 * (new_l*2.3548/fwhm_1)**2 + \
                          (new_m*2.3548/fwhm_2)**2)
    return beamfactor * norm * modelfactor * \
        sp.exp( -2.0*sp.pi*sp.sqrt(-1.0)*(u*l + v*m))


# BSM - TBF - dblquad() only supports real values. need to integrate real and im. separately.
# use FFT and/or break up integral into re and im parts, then return via numpy COMPLEX type...
def visval(beamfwhm,u,v,norm,l0,m0,fwhm_1,fwhm_2,axis_angle):
    '''
    visibility for a single (u,v)
    '''
    lm_lim= 5.0 * (fwhm_1 * fwhm_2)**0.5
    ret_val = sp.zeros(1,dtype=complex)
    real_part =  ...
    return spi.dblquad(vis_igrnd,-1.0*lm_lim,lm_lim,lambda x: -1.0*lm_lim,lambda x:lm_lim, \
                       args=(beamfwhm,u,v,norm,l0,m0,fwhm_1,fwhm_2,axis_angle))

def getbaselines(infile_name,lam=1.0):
    '''
    read an antenna position file and turn it into a snapshot uv coverage.
    '''
    # read standard ALMA antenna position configuration files-
    ant_pos = sp.genfromtxt(infile_name,comments='#')
    n_ant=(ant_pos[:,0]).size
    ant_x=ant_pos[:,0]
    ant_y=ant_pos[:,1]
    # note - build this up including the reflected
    #  visibilities (just use even elements if you do not want)-
    nb=int(n_ant * (n_ant-1) * 0.5)
    u=sp.zeros(nb)
    v=sp.zeros(nb)
    ii=0
    print n_ant
    print nb
    for i in range(0,n_ant-1):
        for j in range(i+1,n_ant):
            u[ii]=(ant_x[i]-ant_x[j])/lam
            v[ii]=(ant_y[i]-ant_y[j])/lam
            #u[ii+1]= -1.0*u[ii]
            #v[ii+1]= -1.0*v[ii]
            ii += 1
    bl={'u':u,'v':v}
    return bl

# good defaults are 1e-4 for both vis_err [Jy?] and frac_stepsize 
def make_fisher_mx(bl,vis_err,deriv_frac_stepsize,beamfwhm,param_vec):
    # param_vec = [norm,l0,m0,fwhm_1,fwhm_2,axis_angle,beamfwhm]
    nb=(bl['u']).size
    nparam=param_vec.size
    # for single Gaussian = 6 free params...
    #nparam=6
    fish_mx=sp.zeros((nparam,nparam))
    for i in range(nparam):
        v1plus=param_vec
        v1minus=param_vec
        v1plus[i]=param_vec[i]*(1.0+deriv_frac_stepsize)
        v1minus[i]=param_vec[i]*(1.0-deriv_frac_stepsize)        
        delta_1 = param_vec[i]*2.0*deriv_frac_stepsize
        for j in range(i,nparam):
            v2plus=param_vec
            v2minus=param_vec
            v2plus[j]=param_vec[j]*(1.0+deriv_frac_stepsize)
            v2minus[j]=param_vec[j]*(1.0-deriv_frac_stepsize)
            delta_2= param_vec[j]*2.0*deriv_frac_stepsize
            this_sum = 0.0
            for k in range(nb):
                u=(bl['u'])[k]
                v=(bl['v'])[k]
                first_term = (visval(beamfwhm,u,v,v1plus[0],v1plus[1],v1plus[2],v1plus[3],v1plus[4],v1plus[5]) -
                              visval(beamfwhm,u,v,v1minus[0],v1minus[1],v1minus[2],v1minus[3],v1minus[4],v1minus[5]))/delta_1
                second_term = (visval(beamfwhm,u,v,v2plus[0],v2plus[1],v2plus[2],v2plus[3],v2plus[4],v2plus[5]) -
                              visval(beamfwhm,u,v,v2minus[0],v2minus[1],v2minus[2],v2minus[3],v2minus[4],v2minus[5]))/delta_2
                this_sum += first_term * second_term
            fish_mx[i,j]=this_sum/ vis_err**2
            fish_mx[j,i]=fish_mx[i,j]

    return fish_mx
