
import scipy as sp
import scipy.integrate as spi
import scipy.linalg as spla

def rot_mx(phi):
    '''
    2d rotation matrix, phi is in radians
    '''
    mx=sp.zeros((2,2))
    mx[0,0]=sp.cos(phi)
    mx[0,1]=sp.sin(phi)
    mx[1,0]= (-1.0)*sp.sin(phi)
    mx[1,1]=sp.cos(phi)
    return mx

def beam_mx(sig_1,sig_2,phi):
    mx=sp.zeros((2,2))
    mx[0,0]=sig_1**2
    mx[1,1]=sig_2**2
    if phi != 0.0:
        mx=mx.dot(rot_mx(phi))
    return mx

def vis_igrnd(l,m,beamfwhm,is_imaginary,u,v,norm,l0,m0,fwhm_1,fwhm_2,axis_angle):
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
    modelfactor = sp.exp( -0.5 * ((new_l*2.3548/fwhm_1)**2 + \
                          (new_m*2.3548/fwhm_2)**2))
    if is_imaginary:
        trigfactor = sp.sin(-2.0*sp.pi*(u*l + v*m))
    else:
        trigfactor = sp.cos(-2.0*sp.pi*(u*l + v*m))
    #print new_l,new_m,beamfactor,modelfactor,trigfactor
    return beamfactor * norm * modelfactor * trigfactor
        
# BSM - TBF - dblquad() only supports real values. need to integrate real and im. separately.
#  use FFT and/or break up integral into re and im parts, then return via numpy COMPLEX type...
# implicitly assumes the telescope primary beam is the phase center (ell,m)=(0,0).
def visval(beamfwhm,u,v,norm,l0,m0,fwhm_1,fwhm_2,axis_angle,brute_force=True):
    '''
    computation of complex visibility for a single (u,v).
    two methods of calculating it are available:
    brute_force=True - a not very clever, brute force calculation
    brute_force=False - use the analytic form
    '''
    lm_lim= 5.0 * beamfwhm
    ret_val = sp.zeros(1,dtype=complex)
    if brute_force:
        real_part = spi.dblquad(vis_igrnd,-1.0*lm_lim,lm_lim,lambda l: -1.0*lm_lim,lambda l:lm_lim, \
                                args=(beamfwhm,False,u,v,norm,l0,m0,fwhm_1,fwhm_2,axis_angle),epsrel=1e-4)
        imag_part = spi.dblquad(vis_igrnd,-1.0*lm_lim,lm_lim,lambda l: -1.0*lm_lim,lambda l:lm_lim, \
                                args=(beamfwhm,True,u,v,norm,l0,m0,fwhm_1,fwhm_2,axis_angle),epsrel=1e-4)
        ret_val = real_part[0] + sp.sqrt(-1.0)*imag_part[0]
    else:
        # uv in vector form-
        u_vec=sp.zeros((2,1))
        u_vec[0]=u
        u_vec[1]=v
        # signal centroid position in vector form-
        x_vec=sp.zeros((2,1))
        x_vec[0]=l0
        x_vec[1]=m0
        # matrix describing signal shape-
        sig_a = beam_mx(fwhm_1/2.3548,fwhm_2/2.3548,axis_angle*0.0174533)
        sig_a_inv=spla.inv(sig_a)
        # matrix describing beam shape-
        sig_b = beam_mx(beamfwhm/2.3548,beamfwhm/2.3548,0.0)
        sig_b_inv=spla.inv(sig_b)
        # combined quantities-
        sig_c_inv=sig_a_inv+sig_b_inv
        sig_c=spla.inv(sig_c_inv)
        sig_ab_inv=spla.inv(sig_a+sig_b)
        mc = sig_c.dot(sig_a_inv.dot(x_vec))
        pre_factor = 2.0*sp.pi* (spla.det(sig_c))**0.5 * \
                     sp.exp(-0.5*sp.transpose(x_vec).dot(sig_ab_inv.dot(x_vec))) * \
                     sp.exp(-4.0* sp.pi**2 * sp.transpose(u_vec).dot(sig_c.dot(u_vec)))
        real_part = norm*pre_factor * sp.cos( -2.0*sp.pi* sp.transpose(mc).dot(u_vec))
        imag_part = norm*pre_factor * sp.sin( -2.0*sp.pi* sp.transpose(mc).dot(u_vec))
        ret_val = real_part + sp.sqrt(-1.0)*imag_part
    return ret_val

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
def make_fisher_mx(bl,vis_err,deriv_frac_stepsize,beamfwhm,param_vec,brute_force=True):
    # param_vec = [norm,l0,m0,fwhm_1,fwhm_2,axis_angle,beamfwhm]
    nb=(bl['u']).size
    nparam=param_vec.size
    # for single Gaussian = 6 free params...
    #nparam=6
    fish_mx=sp.zeros((nparam,nparam),dtype=sp.float128)
    for i in range(nparam):
        v1plus=sp.copy(param_vec)
        v1minus=sp.copy(param_vec)
        v1plus[i]=param_vec[i]*(1.0+deriv_frac_stepsize)
        v1minus[i]=param_vec[i]*(1.0-deriv_frac_stepsize)        
        delta_1 = param_vec[i]*2.0*deriv_frac_stepsize
        #print v1plus-v1minus
        #print delta_1
        for j in range(i,nparam):
            v2plus=sp.copy(param_vec)
            v2minus=sp.copy(param_vec)
            v2plus[j]=param_vec[j]*(1.0+deriv_frac_stepsize)
            v2minus[j]=param_vec[j]*(1.0-deriv_frac_stepsize)
            delta_2= param_vec[j]*2.0*deriv_frac_stepsize
            #print v2plus-v2minus
            #print delta_2
            this_sum = 0.0
            print i,j
            for k in range(nb):
                u=(bl['u'])[k]
                v=(bl['v'])[k]
                first_term = (visval(beamfwhm,u,v,v1plus[0],v1plus[1],v1plus[2],v1plus[3],v1plus[4],v1plus[5],brute_force) -
                              visval(beamfwhm,u,v,v1minus[0],v1minus[1],v1minus[2],v1minus[3],v1minus[4],v1minus[5],brute_force))/delta_1
                second_term = (visval(beamfwhm,u,v,v2plus[0],v2plus[1],v2plus[2],v2plus[3],v2plus[4],v2plus[5],brute_force) -
                               visval(beamfwhm,u,v,v2minus[0],v2minus[1],v2minus[2],v2minus[3],v2minus[4],v2minus[5],brute_force))/delta_2
                this_sum += sp.real(first_term) * sp.real(second_term)
                this_sum += sp.imag(first_term) * sp.imag(second_term)
                #print "   ",k,u,v,first_term,second_term,delta_1,delta_2
            fish_mx[i,j]=this_sum/ vis_err**2
            fish_mx[j,i]=fish_mx[i,j]
            # looks like the second_term is giving some NaNs in the off-diagonals...
            # i think the trouble is my stepsize is fractional and some of the parameter values are zero.
    return fish_mx
