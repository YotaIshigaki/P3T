from math import *
import numpy as np
import sys
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PI = acos(-1.)
L = 149597870700
M = 1.9884e30


def main():
    dir = sys.argv[1]
    xrange = [0.97,1.03]
    yrange = [-0.002,0.025]
    mass_c = 1.e20/M
    mass_s = 1.e20/M
    mass_b = 1.e20/M
    size = 10
    ax = "a-ei"
    pngfile = dir + "/massDistribution.png"
    
    if len(sys.argv) > 5 :
        xrange = [float(sys.argv[2]), float(sys.argv[3])] 
        yrange = [float(sys.argv[4]), float(sys.argv[5])]

    if len(sys.argv) > 8 :
        mass_c = float(sys.argv[6])/M
        mass_s = float(sys.argv[7])/M
        mass_b = float(sys.argv[8])/M

    if len(sys.argv) > 9 :
        ax = sys.argv[9]
    
    number_min, number_max = makeFigureAll(dir, xrange, yrange, ax, mass_c, mass_s, mass_b, size)

    if len(sys.argv) > 11 :
        number_min = float(sys.argv[10])
        number_max = float(sys.argv[11])

    file_list = []
    number4 = int( (number_max-number_min)/4 + number_min )
    number2 = int( (number_max-number_min)/2 + number_min )
    file_list.append( dir + "/snap" + "{:0>6}".format(number4) + ".dat" )
    file_list.append( dir + "/snap" + "{:0>6}".format(number2) + ".dat" )
    file_list.append( dir + "/snap" + "{:0>6}".format(number_max) + ".dat" ) 

    makeCumulativeNumber(file_list, pngfile, mass_b=mass_b)
    
    return

class Particle:
    def __init__(self, mass, x, y, z, vx, vy, vz):
        self.mass = mass
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz
        
    def pos(self, target=None):
        if target is None:
            return [self.x, self.y, self.z]
        else:
            return [self.x-target.x, self.y-target.y, self.z-target.z]

    def pos2(self, target=None):
        x, y, z = self.pos(target)
        return x**2 + y**2 + z**2
        
    def vel(self, target=None):
        if target is None:
            return [self.vx, self.vy, self.vz]
        else:
            return [self.vx-target.vx, self.vy-target.vy, self.vz-target.vz]

    def vel2(self, target=None):
        vx, vy, vz = self.vel(target)
        return vx**2 + vy**2 + vz**2

    def getOrbitalElement(self, mu=1):
        r = sqrt(self.pos2())
        v2 = self.vel2()
        rv = self.x*self.vx + self.y*self.vy + self.z*self.vz
        rxv = [self.y*self.vz - self.z*self.vy,\
               self.z*self.vx - self.x*self.vz,\
               self.x*self.vy - self.y*self.vx\
        ]
        
        ax  = 1./(2./r - v2/mu)
        ecc = sqrt( (1.-r/ax)**2 + (rv)**2/(mu*ax) )
        inc = atan2(sqrt(rxv[0]**2+rxv[1]**2), rxv[2])

        return [ax, ecc, inc]

class Header:
    def __init__(self, time, n_body, \
                 etot0, ekin0, ephi_sun0, ephi_planet0, edisp0, \
                 etot1, ekin1, ephi_sun1, ephi_planet1, edisp1):
        self.time   = time
        self.n_body = n_body
        self.etot0        = etot0
        self.ekin0        = ekin0
        self.ephi_sun0    = ephi_sun0
        self.ephi_planet0 = ephi_planet0
        self.edisp0       = edisp0
        self.etot1        = etot1
        self.ekin1        = ekin1
        self.ephi_sun1    = ephi_sun1
        self.ephi_planet1 = ephi_planet1
        self.edisp1       = edisp1

        
#########################
###     FUNCTIONS     ###
#########################

def makeFigureAll(dir, xrange, yrange, \
                  axis="a-e", \
                  mass_c=1.e20/M, \
                  mass_s=1.e20/M, \
                  mass_b=1.e20/M, \
                  size=10):
    """
    Make snapshot figures and other figures in directory

    Arguments
    dir:     character string indicating a directory
    xrange:  x direction range of a snapshot figure
    yrange:  y direction range of a snapshot figure
    axis:    character string indicating elements of x,y in a snapshot figure (default "a-e")
             "a-e" -- semimajor axis - eccentricity
             "a-i" -- semimajor axis - inclination
             "a-ei" -- semimajor axis - root mean squared of eccentricity and inclination
             "e-i" -- eccentricity - inclination
    mass_c:  minimum particle mass to draw in a snapshot figure (default 1.e20 kg)
    mass_s:  mass as a reference for determining sizes of markers (default 1.e20 kg)
    mass_b:  mass_b:  width of mass bin (default 1.e20 kg) 
    size:    marker size of particles with mass mass_s (default 10)

    Return Values
    number_min: minimum number of snapshot file
    number_max: maximum number of snapshot file
    """

    figure_dir = dir + "/figure"
    
    if not os.path.exists(figure_dir):
        os.makedirs(figure_dir)

    files_all = os.listdir(dir)
    files = [f for f in files_all if os.path.isfile(os.path.join(dir, f))]
    snaps = [f for f in files if f[0:4]=="snap"]
    snaps.sort()

    datfile = dir + "/outcomes.dat" 
    fout = open(datfile, 'w')

    time_v = []
    mass_mx_v = []
    mass_av_v = []
    ecc_rms_v = []
    inc_rms_v = []
    n_body_v = []
    number_max = 0
    number_min = int(snaps[0][4:-4])
    for snap in snaps :
        number = int(snap[4:-4])
        if number > number_max : number_max = number
        if number < number_min : number_min = number
        
        snapfile = dir + "/" + snap
        if axis in ["a-e", "a-i", "a-ei", "e-i"] :
            tail = ""
            if axis == "a-e" : tail = "a_e"
            elif axis == "a-i" : tail = "a_i"
            elif axis == "a-ei" : tail = "a_ei"
            elif axis == "e-i" : tail = "e_i"
            figurefile = figure_dir + "/" + snap[0:-4] + tail + ".png"

        else : figurefile = None
        
        header, mass_mx, mass_av, ecc_rms, inc_rms \
            = makeSnap(snapfile, figurefile, xrange, yrange, axis, mass_c, mass_s, size)
        
        time_v.append(header.time)
        mass_mx_v.append(mass_mx)
        mass_av_v.append(mass_av)
        ecc_rms_v.append(ecc_rms)
        inc_rms_v.append(inc_rms)
        n_body_v.append(header.n_body)

        dammy0, dammy1, dammy2, a, p = calcCumulativeNumber(snapfile, mass_b=mass_b)

        fout.write(str(header.time) + "\t" \
                   + str(mass_mx) + "\t" \
                   + str(mass_av) + "\t" \
                   + str(ecc_rms) + "\t" \
                   + str(inc_rms) + "\t" \
                   + str(a) + "\t" \
                   + str(p) + "\n" )

    fout.close()

    plt.rcParams['font.family'] = "Times New Roman"
    plt.rcParams['font.size'] = 13
    
    plt.figure()
    plt.xscale("log")
    plt.yscale("log")
    plt.plot([time/(2.*PI) for time in time_v], [mass*M for mass in mass_mx_v], \
             linewidth=1, \
             label="$M_{\mathsf{max}}$", color='b', ls='-')
    plt.plot([time/(2.*PI) for time in time_v], [mass*M for mass in mass_av_v], \
             linewidth=1, \
             label="<$m$>", color='b', ls='--')
    plt.ylabel("$M_{\mathsf{max}}$, <$m$> (kg)", fontsize=15)
    plt.xlabel("$t$ (year)", fontsize=15)
    plt.legend(fontsize=13, loc='upper left')
    plt.xlim([min(time_v)/(2.*PI),max(time_v)/(2.*PI)])
    plt.ylim([mass_b*M/5,max(mass_mx_v)*M*5])
    
    pngfile = dir + "/massEvolution.png"
    plt.savefig(pngfile, dpi=100)
    plt.close()

    plt.figure()
    plt.xscale("log")
    plt.yscale("log")    
    plt.yticks([1.e-4, 2.e-4, 5.e-4, 1.e-3, 2.e-3, 5.e-3, \
                1.e-2, 2.e-2, 5.e-2, 1.e-1, 2.e-1, 5.e-1, 1.e0], \
                ["$1{\mathsf{x}}10^{-4}$","$2{\mathsf{x}}10^{-4}$","$5{\mathsf{x}}10^{-4}$",\
                 "$1{\mathsf{x}}10^{-3}$","$2{\mathsf{x}}10^{-3}$","$5{\mathsf{x}}10^{-3}$",\
                 "$1{\mathsf{x}}10^{-2}$","$2{\mathsf{x}}10^{-2}$","$5{\mathsf{x}}10^{-2}$",\
                 "$1{\mathsf{x}}10^{-1}$","$2{\mathsf{x}}10^{-1}$","$5{\mathsf{x}}10^{-1}$","$1{\mathsf{x}}10^{0}$"])
    plt.plot([time/(2.*PI) for time in time_v], ecc_rms_v, \
             linewidth=1, \
             label="<$e^2$>${}^{1/2}$", color='b', ls='-')
    plt.plot([time/(2.*PI) for time in time_v], inc_rms_v, \
             linewidth=1, \
             label="<$i^2$>${}^{1/2}$", color='r', ls='-.')
    plt.ylabel("<$e^2$>${}^{1/2}$, <$i^2$>${}^{1/2}$", fontsize=15)
    plt.xlabel("$t$ (year)", fontsize=15)
    plt.legend(fontsize=13, loc='upper left')
    plt.xlim([min(time_v)/(2.*PI),max(time_v)/(2.*PI)])
    plt.ylim([min(ecc_rms_v+inc_rms_v)/2,max(ecc_rms_v+inc_rms_v)*2])
    
    pngfile = dir + "/eccincEvolution.png"
    plt.savefig(pngfile, dpi=100)
    plt.close()

    return [number_min, number_max]


def calcCumulativeNumber(datfile,\
                         xrange=None,\
                         mass_b=1.e20/M):
    """
    Calculate parameter a,p of cumulative number distributuion as Nc~a*mass^p from snapshot files.
    
    Arguments
    datfile: character strings indicating a snapshot filename
    xrange:  mass range to fit (default None)
    mass_b:  width of mass bin (default 1.e20 kg) 

    Return Values
    header: Header object indicating snapshot file's header
    mass_d: list of particle mass to make figure
    nc_d:   list of cumulative number to make figure
    a:      parameter a of cumulative number distributuion
    p:      parameter p of cumulative number distributuion
    """
    
    header, pp = readSnap(datfile)

    mass_v = []
    i = 0
    for ptcl in pp.values():
        if xrange is None :
            mass_v.append(ptcl.mass)
        else :
            if xrange[0] <= ptcl.mass and ptcl.mass <= xrange[1] :
                mass_v.append(ptcl.mass)

    mass_v.sort()
    mass_v.reverse()

    mass_d = []
    nc_d = []
    n_body = len(mass_v)
    mass = ceil(mass_v[0]/mass_b) * mass_b
    nc = 0
    while nc < n_body and mass > 0:
        mass_d.append(mass)
        while nc < n_body and mass_v[nc] >= mass:
            nc = nc + 1
        nc_d.append(nc)
        mass = mass - mass_b

    sum_log_nc = 0.
    sum_log_mass = 0.
    sum_log_mass2 = 0.
    sum_log_nc_log_mass = 0.
    n = 0
    for i in range(len(mass_d)):
        if nc_d[i] <= 0. :
            nc_d[i] = 1.e-30
        else :
            sum_log_nc += log(nc_d[i])
            sum_log_mass += log(mass_d[i])
            sum_log_mass2 += log(mass_d[i])**2
            sum_log_nc_log_mass += log(nc_d[i])*log(mass_d[i])
            n = n + 1

    if n*sum_log_mass2 - sum_log_mass**2 == 0. :
        a = None
        p = None
    else :
        a = (n*sum_log_nc_log_mass - sum_log_mass*sum_log_nc) \
            / (n*sum_log_mass2 - sum_log_mass**2)
        b = (sum_log_mass2*sum_log_nc - sum_log_nc_log_mass*sum_log_mass) \
            / (n*sum_log_mass2 - sum_log_mass**2)
        p = exp(b)

    return [header, mass_d, nc_d, a, p]


def makeCumulativeNumber(datfiles, pngfile,
                         xrange=None, \
                         yrange=None, \
                         mass_b=1.e20/M, \
                         labels=None, \
                         colors=["royalblue", "mediumblue", "midnightblue"], \
                         styles=["-.","--","-"]):
    """
    Make a figure of cumulative number distributuion from snapshot files.
    
    Arguments
    datfiles: list of character strings indicating snapshot filenames
    pngfile:  character string indicating a snapshot figure filename
    xrange:   x direction range of a snapshot figure (default None)
    yrange:   y direction range of a snapshot figure (default None)
    mass_b:   width of mass bin (default 1.e20 kg)
    labels:   list of characters indicating labels (default None)
    colors:   list of characters indicating marker colors (default ["royalblue","mediumblue","midnightblue"])
    styles:   list of characters indicating line styles (default ["-.","--","-"])

    Return Values
    header:  Header object indicating snapshot file's header
    a_v:     list of parameter a of cumulative number distributuion
    p_v:     list of parameter p of cumulative number distributuion
    """
    
    plt.rcParams['font.family'] = "Times New Roman"
    plt.rcParams['font.size'] = 13

    if labels == None :
        flag = True
        labels = []

    plt.figure()
    plt.xscale("log")
    plt.yscale("log")

    i = 0
    n_body = 0
    p_v = []
    a_v = []
    mass_d_mx = []
    mass_d_mn = []
    for datfile in datfiles:
        header, mass_d, nc_d, a, p = calcCumulativeNumber(datfile, mass_b=mass_b)
        if flag : labels.append("%0.f yr"%(header.time/(2*PI)))
        if header.n_body > n_body : n_body = header.n_body
        p_v.append(p)
        a_v.append(a)

        mass_d_mx.append(mass_d[0])
        mass_d_mn.append(mass_d[-2])
        
        plt.plot([mass*M for mass in mass_d], nc_d, linewidth=1, label=labels[i], \
                 color=colors[i], ls=styles[i])

        i = i + 1

    plt.xlabel("$m$ (kg)")
    plt.ylabel("$N_c$")
    plt.legend(fontsize=13)
    if xrange is None:
        plt.xlim([min(mass_d_mn)*M/2,max(mass_d_mx)*M*2])
    else :
        plt.xlim(xrange)
    if yrange is None:
        plt.ylim([0.5,n_body])
    else :
        plt.ylim(ylim)

    plt.savefig(pngfile, dpi=100)
    plt.close()

    return [header, a_v, p_v]

    
def makeSnap(datfile, pngfile, xrange, yrange, \
             axis="a-e",\
             mass_c=1.e20/M,\
             mass_s=1.e20/M,\
             size=10,\
             color="blue"):
    """
    Make a snapshot figure from a snapshot file.

    Arguments
    datfile: character string indicating a snapshot filename
    pngfile: character string indicating a snapshot figure filename
    xrange:  x direction range of a snapshot figure
    yrange:  y direction range of a snapshot figure
    axis:    character string indicating elements of x,y in a snapshot figure (default "a-e")
             "a-e" -- semimajor axis - eccentricity
             "a-i" -- semimajor axis - inclination
             "a-ei" -- semimajor axis - root mean squared of eccentricity and inclination
             "e-i" -- eccentricity - inclination
    mass_c:  minimum particle mass to draw in a snapshot figure (default 1.e20 kg)
    mass_s:  mass as a reference for determining sizes of markers (default 1.e20 kg)
    size:    marker size of particles with mass mass_s (default 10)
    color:   character indicating marker color (default "blue")

    Return Values
    header:  Header object indicating a snapshot file's header
    mass_mx: maximum particle mass
    mass_av: mean particle mass
    ecc_rms: root mean squared of eccentricities
    inc_rms: root mean squared of inclinations
    """
    
    header, pp = readSnap(datfile)

    time_str = "%0.f yr"%(header.time/(2*PI))

    mass_v = []
    ax_v   = []
    ecc_v  = []
    inc_v  = []

    mass_mx = 0.
    mass_av = 0.
    ecc_rms = 0.
    inc_rms = 0.

    for ptcl in pp.values():
        ax, ecc, inc = ptcl.getOrbitalElement()
        if ptcl.mass >= mass_c :
            mass_v.append(ptcl.mass)
            ax_v.append(ax)
            ecc_v.append(ecc)
            inc_v.append(inc)

        if mass_mx < ptcl.mass :
            mass_mx = ptcl.mass
        mass_av += ptcl.mass
        ecc_rms += ecc*ecc
        inc_rms += inc*inc

    mass_av /= header.n_body
    ecc_rms = sqrt(ecc_rms/header.n_body)
    inc_rms = sqrt(inc_rms/header.n_body)

    # Make Figure
    if axis in ["a-e", "a-i", "a-ei", "e-i"] and pngfile is not None: 
        plt.rcParams['font.family'] = "Times New Roman"
        plt.rcParams['font.size'] = 13
        
        plt.figure()
        if axis == "a-e" :
            plt.scatter(ax_v, ecc_v,
                        s=[(mass/mass_s)**(2./3.)*size for mass in mass_v],\
                        c=color, alpha=0.2, linewidths=0.5 )
            plt.xlabel("$a$ (AU)", fontsize=15)
            plt.ylabel("$e$", fontsize=15)
            
        elif axis == "a-i" :
            plt.scatter(ax_v, inc_v,
                        s=[(mass/mass_s)**(2./3.)*size for mass in mass_v],\
                        c=color, alpha=0.2, linewidths=0.5 )
            plt.xlabel("$a$ (AU)", fontsize=15)
            plt.ylabel("$i$", fontsize=15)
            
        elif axis == "a-ei" :
            plt.scatter(ax_v, [sqrt(ecc**2 + inc**2) for (ecc, inc) in zip(ecc_v, inc_v)],
                        s=[(mass/mass_s)**(2./3.)*size for mass in mass_v],\
                        c=color, alpha=0.2, linewidths=0.5 )
            plt.xlabel("$a$ (AU)", fontsize=15)
            plt.ylabel("$\sqrt{e^2+i^2}$", fontsize=15)
        elif axis == "e-i" :
            plt.scatter(ecc_v, inc_v,
                        s=[(mass/mass_s)**(2./3.)*size for mass in mass_v],\
                        c=color, alpha=0.2, linewidths=0.5 )
            plt.xlabel("$e$", fontsize=15)
            plt.ylabel("$i$", fontsize=15)
            
        plt.title(time_str)
        plt.xlim(xrange)
        plt.ylim(yrange)
            
        plt.savefig(pngfile, dpi=130)
        plt.close()
    
    return [header, mass_mx, mass_av, ecc_rms, inc_rms]

def readSnapForParticles(filename, ids):
    """
    Read paticular particles' data from a snapshot file.

    Arguments
    filename: character string indicating a snapshot filename
    ids:      list of particle id

    Return Values
    header: Header object indicating snapshot file's header
    pp:     dectionary {id : Particle object} indicating particles 
    """

    i = 0
    pp = {}
    
    with open(filename) as f:
        for line in f:
            part = [p.strip() for p in line.split("\t")]

            if i==0 :
                header = Header(float(part[0]),\
                                int(part[1]),\
                                float(part[2]),\
                                float(part[3]),\
                                float(part[4]),\
                                float(part[5]),\
                                float(part[6]),\
                                float(part[7]),\
                                float(part[8]),\
                                float(part[9]),\
                                float(part[10]),\
                                float(part[11]) )
            else :
                idx = int(part[0])
                if idx in ids :
                    ptcl = Particle(float(part[1]),\
                                    float(part[2]),\
                                    float(part[3]),\
                                    float(part[4]),\
                                    float(part[5]),\
                                    float(part[6]),\
                                    float(part[7]) )

                    pp[idx] = ptcl
                    
            i += 1

    return [header, pp]

def readSnap(filename):
    """
    Read particles' data from a snapshot file.

    Arguments
    filename: character string indicating snapshot filename

    Return Values
    header: Header object indicating snapshot file's header
    pp:     dectionary {id : Particle object} indicating particles 
    """
    
    i = 0
    pp = {}

    with open(filename) as f:
        for line in f:
            part = [p.strip() for p in line.split("\t")]

            if i==0 :
                header = Header(float(part[0]),\
                                int(part[1]),\
                                float(part[2]),\
                                float(part[3]),\
                                float(part[4]),\
                                float(part[5]),\
                                float(part[6]),\
                                float(part[7]),\
                                float(part[8]),\
                                float(part[9]),\
                                float(part[10]),\
                                float(part[11]) )
                
            else :
                idx = int(part[0])
                ptcl = Particle(float(part[1]),\
                                float(part[2]),\
                                float(part[3]),\
                                float(part[4]),\
                                float(part[5]),\
                                float(part[6]),\
                                float(part[7]) )
                pp[idx] = ptcl;

            i += 1

    return [header, pp]

def calcEnergy(pp, m_sun=1., eps=0.):
    """
    Calculate particle system's energy.

    Arguments
    pp:    dectionary {id : Particle object} indicating particles 
    m_sun: mass of a central star (default 1.0)
    eps:   softening parameter (default 0.0)

    Return Values
    etot:        total energy
    ekin:        kinetic energy
    ephi_sun:    potential energy by central star
    ephi_planet: potential energy by interaction of particles
    """

    ekin = 0.
    ephi_sun = 0.
    ephi_planet = 0.
    
    for id, ptcl in pp.items():
        dr2 = ptcl.pos2() + eps**2
        dv2 = ptcl.vel2()
        
        ekin     += ptcl.mass * dv2
        ephi_sun += ptcl.mass / sqrt(dr2)

        for id2, ptcl2 in pp.items():
            if id < id2 :
                dr2 = ptcl.pos2(ptcl2) + eps**2               
                ephi_planet += ptcl.mass * ptcl2.mass / sqrt(dr2)
                
    ekin *= 0.5
    ephi_sun *= m_sun
    etot = ekin + ephi_sun + ephi_planet

    return [etot, ekin, ephi_sun, ephi_planet]

if __name__ == '__main__':
    main()


        
        
        

    

    

    
         
    
        

                
                        
                        
