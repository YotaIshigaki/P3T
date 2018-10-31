from analysis import *

PI = acos(-1.)
L = 149597870700
M = 1.9884e30

def readSnapPENTACLE(filename):
    """
    Read particles' data from a snapshot file in PENTACLE format.
    
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
            if i==0 :
                part = [p.strip() for p in line.split("\t")]
                header = Header(float(part[1]),\
                                int(part[0]),
                                float(part[5]),\
                                float(part[2]),\
                                float(part[3]),\
                                float(part[4]),\
                                float(part[6]),\
                                float(part[10]),\
                                float(part[7]),\
                                float(part[8]),\
                                float(part[9]),\
                                float(part[11]) )
                
            else :
                part = [p.strip() for p in line.split(" ")]
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


def changeFormatToPENTACLE(filename, filename_p) :
    """
    Make snapshot in PENTACLE format.
    
    Arguments
    filename:   character string indicating snapshot filename
    filename_p: character string indicating snapshot filename in PENTACLE format
    """
    
    header, pp = readSnap(filename)
    
    fout = open(filename_p, 'w')
    fout.write(str(header.n_body) + "\t" \
               + str(header.time) + "\t" \
               + str(header.ekin0) + "\t" \
               + str(header.ephi_sun0) + "\t" \
               + str(header.ephi_planet0) + "\t" \
               + str(header.etot0) + "\t" \
               + str(header.edisp0) + "\t" \
               + str(header.ekin1) + "\t" \
               + str(header.ephi_sun1) + "\t" \
               + str(header.ephi_planet1) + "\t" \
               + str(header.etot1) + "\t" \
               + str(header.edisp1) + "\n" )

    for id, ptcl in pp.items():
        fout.write(str(id) + " " \
                   + str(ptcl.mass) + " " \
                   + str(ptcl.x) + " " \
                   + str(ptcl.y) + " " \
                   + str(ptcl.z) + " " \
                   + str(ptcl.vx) + " " \
                   + str(ptcl.vy) + " " \
                   + str(ptcl.vz) + " " \
                   + str(0.) + " " \
                   + str(0.) + " " \
                   + str(0.) + " " \
                   + str(0.) + " " \
                   + str(0) + "\n" )

    fout.close()

    return

def changeFormatFromPENTACLE(filename_p, filename) :
    """
    Make snapshot from snapshot in  PENTACLE.
    
    Arguments
    filename_p: character string indicating snapshot filename in PENTACLE format
    filename:   character string indicating snapshot filename
    """

    header, pp = readSnapPENTACLE(filename_p)
    
    fout = open(filename, 'w')
    fout.write(str(header.time) + "\t" \
               + str(header.n_body) + "\t" \
               + str(header.etot0) + "\t" \
               + str(header.ekin0) + "\t" \
               + str(header.ephi_sun0) + "\t" \
               + str(header.ephi_planet0) + "\t" \
               + str(header.edisp0) + "\t" \
               + str(header.etot1) + "\t" \
               + str(header.ekin1) + "\t" \
               + str(header.ephi_sun1) + "\t" \
               + str(header.ephi_planet1) + "\t" \
               + str(header.edisp1) + "\n" )

    for id, ptcl in pp.items():
        fout.write(str(id) + " " \
                   + str(ptcl.mass) + " " \
                   + str(ptcl.x) + " " \
                   + str(ptcl.y) + " " \
                   + str(ptcl.z) + " " \
                   + str(ptcl.vx) + " " \
                   + str(ptcl.vy) + " " \
                   + str(ptcl.vz) + "\n")

    fout.close()

    return
