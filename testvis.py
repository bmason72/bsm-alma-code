
import scipy as sp

def vis_igrnd(l,m,u,v,norm,l0,m0,fwhm_1,fwhm_2,axis_angle,beamfwhm):
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


def visval(u,v,norm,l0,m0,fwhm_1,fwhm_2,axis_angle,beamfwhm):
    '''
    visibility for a single (u,v)
    '''

    lm_lim= 5.0 * (fwhm_1 * fwhm_2)**0.5

    return sp.dblquad(vis_igrnd,-1.0*lm_lim,lm_lim,lambda x: -1.0*lm_lim,lambda x:lm_lim)

def getbaselines(infile_name,lam=1.0):
    # read standard ALMA antenna position configuration files-
    ant_pos = sp.genfromtxt(infile_name,comments='#')
    n_ant=(ant_pos[:,0]).size
    ant_x=ant_pos[:,0]
    ant_y=ant_pos[:,1]
    # note - build this up including the reflected
    #  visibilities (just use even elements if you do not want)-
    nb=int(n_ant * (n_ant-1))
    u=sp.zeros(nb)
    v=sp.zeros(nb)
    ii=0
    print n_ant
    print nb
    for i in range(0,n_ant-1):
        for j in range(i+1,n_ant):
            u[ii]=(ant_x[i]-ant_x[j])/lam
            v[ii]=(ant_y[i]-ant_y[j])/lam
            u[ii+1]= -1.0*u[ii]
            v[ii+1]= -1.0*v[ii]
            ii += 2
    bl={'u':u,'v':v}
    return bl

