#!/usr/bin/env pypy
# -*- coding: utf-8 -*-

"""
© Univ. Paris Diderot & Inserm, 2015
guillaume.postic@univ-paris-diderot.fr
This software is a computer program whose purpose is to assign membrane
boundaries to a protein three-dimensional structure, by using the spatial
coordinates of the alpha carbons. It can also process coarse-grained
models of protein structures. The algorithm follows an approach that
treats the problem of membrane assignment as a binary classification.
The software is implemented in Python and requires the PyPy interpreter.
It also requires Naccess 2.1.1 (S.J. Hubbard, 1996) to calculate the
atomic accessible surface.
This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 
As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 
In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 
The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
"""


# Notes:
# - in the source code, the 'C value' is called Q
# - membrane thickness is sometimes called 'width'


def fsplit_path_name(filename):
    mylist = filename.split('/')
    mypath = ''

    for element in mylist:
        if element != mylist[-1]:
            mypath+=element
            mypath+='/'

    return mypath, mylist[-1]
# fsplit_path_name()


def freadpdb(infile, af=False):
    """ Read a PDB file and return a list of all lines (without newlines)
    Two types of behavior: 
    1. Default ANVIL
    2. [modified by A.G. 5.06.2022] 
    For AlphaFold2 files: containes ENDMDL and END as last lines -> 
    thus, need to add membrane definition inside the model, otherwise 
    it will not be displayed with convenctional PDB 3D viewers
    """
    file_name=infile
    try:
        with open( file_name, 'r') as data:
            lines=data.readlines()
        for i,item in enumerate(lines):
            item=item.replace('\n','')
            lines[i]=item
        # [A.G. 5.06.22]
        # input precondition: PDB file contains one model and is well formatted (END)
        # otherwise can not insert the membrane to the output file because
        # for displaying PDB mebrane (heteroatoms) need to be a part of the model.
        if af:
            assert ((lines[-1] == 'END                                                                             ') and (lines[-2] == 'ENDMDL                                                                          ')), "ENDMDL and END terminators at the end of an AlphaFold2 PDB file were expected."
    except:
        sys.exit('Error while reading {}'.format(file_name) )

    return lines
# freadpdb()


def fnaccess(infile, naccessdir, outputdir, asatype, atomtype): 
    """ Parse NACCESS data
    """
    path, infilemv = fsplit_path_name(infile)
    infilemv=os.path.splitext(infilemv) # To avoid problems such as *~
    infilemv=infilemv[0]
    datafile=outputdir+infilemv+'.rsa'

    infileCA = infile + atomtype

    try:
        os.system('grep " {} " {} > {}'.format(atomtype, infile, infileCA))
    except:
        sys.exit('Error while extracting (grep) atoms')

    try:
        with open(datafile,'r') as data:
            lines = data.readlines()
    except:
        try:
            # probe radius = 4 Å
            os.system('{}naccess {} -p 4 > /dev/null'.format(naccessdir, infileCA))
        except:
            sys.exit('Error while running NACCESS')
        try:
            os.system('mv {}.rsa {} 2> /dev/null'.format(infilemv, outputdir))
            os.system('mv {}.log {} 2> /dev/null'.format(infilemv, outputdir))
            os.system('mv {}.asa {} 2> /dev/null'.format(infilemv, outputdir))
            os.system('mv {} {}'.format(infileCA, outputdir))
            with open(datafile,'r') as data:
                lines = data.readlines()
        except:
            sys.exit('Error: Cannot move output files')

    naccess=[]
    # naccess[0] = chain; [1] = num res; [2] = relative SA value

    try:
        for i, item in enumerate(lines):
            if item[0:3] == 'RES':
                val=[  item[8].replace(' ',''), int(item[9:13].replace(' ','')) ]
                if asatype == 'relative':
                    val.append( float(item[35:42].replace(' ','')) )
                else:
                    val.append( float(item[28:36].replace(' ','')) )
                naccess.append(val)
    except:
        pass

    return naccess
# fnaccess()


def fcoords_atoms(lines, atomtype):
    """ Extract the atomic coordinates from the PDB file
        Works with fparsePDBline()
    """
    coord=[] #
    for i in lines:
        if i[0:4] == 'ATOM' and (i[12:16].replace(' ','')) == atomtype: # HERE
            coord.append(fparsePDBline(i))

    return coord
# fcoords_atoms()


def fparsePDBline(i):
    """ Parse one line of a PDB file (record: ATOM)
    """
    x=float(i[30:38])
    y=float(i[38:46])
    z=float(i[46:54])
    name=i[76:78].replace(' ','')
    num=int(i[22:26].replace(' ',''))
    chain=i[21].replace(' ','')
    amino_acid=i[17:20].replace(' ','')     
    if name != '':
        return [x,y,z, name, chain, num, amino_acid]
    else:
        name=i[12:16].replace(' ','')
        name=name[0]
        return [x,y,z, name, chain, num, amino_acid]
# fparsePDBline()


def fcenterofmass(lines):
    """ Calculate the center of mass of the protein molecule
    """
    mass = 0
    x = 0.0
    y = 0.0
    z = 0.0
    for i in lines:
        x += 12.0107*i[0] # atomic mass of carbon = 12.0107
        y += 12.0107*i[1]
        z += 12.0107*i[2]
        mass+=12.0107
    return [x/mass,y/mass,z/mass]
# fcenterofmass()


def fdistance(A,B):
    """ Distance between two points
    """
    return ( (A[0]-B[0])*(A[0]-B[0]) + (A[1]-B[1])*(A[1]-B[1]) + (A[2]-B[2])*(A[2]-B[2]) )**0.5
# fdistance()


def fradius(lines,center):
    """ Return the distance between the center of mass and the most peripheral atom of the protein
    """
    max_dist=0
    for i in lines:
        temp = fdistance(center,i[0:3])
        if temp > max_dist:
            max_dist = temp
    return max_dist
# fradius()


def fsphere_points(N,r,center):
    """ Calculate the coordinates of N points distributed on a sphere
    """
    points = []
    k=0
    old_phi=0
    while k<N:
        k+=1
        h = -1 + 2*(k-1)/float(N-1)
        theta = acos(h)
        if k==1 or k==N:
            phi = 0
        else:
            phi = (old_phi + 3.6/(N*(1-h*h))**0.5) % (2*pi)
        points.append((r*sin(phi)*sin(theta)+center[0], 
        r*cos(theta)+center[1], 
        r*cos(phi)*sin(theta)+center[2]))
        old_phi = phi
    return points
# fsphere_points()


def fjoin_Coords_Naccess(coords, naccess):
    """ Join NACCESS and structural data into a file containing: 
        x, y and z coordinates,
        atom name, (NOT USED)
        ASA value, 
        name res,
        num res,
        chain
    """
    data=[]
    error=0
    jmax= len(naccess)-1
    if jmax > 0:
        for i, item in enumerate(coords):
            j=0
            while (naccess[j][0] != item[4] or naccess[j][1] != item[5] ) and j < jmax:
                j+=1
            if naccess[j][0] == item[4] and naccess[j][1] == item[5]:
                data.append([ coords[i][0], coords[i][1], coords[i][2],
                coords[i][3], naccess[j][2], coords[i][6], naccess[j][1], 
                naccess[j][0]  ])
                nb=0
            else:
                error+=1
                if error >= len(naccess)*10:
                    sys.exit('Critical failure! Too many errors with the NACCESS file')

    else:
        sys.exit('Error: NACCESS *.rsa file is empty')

    #Debug
    #if error >0:
    #    sys.stderr.write("Warning: NACCESS file: {} residues not taken into account\n".format(error))

    return data
# fjoin_Coords_Naccess()


def fvect(A,B):
    """ Calculate the AB vector (not BA)
    """
    AB=[0,0,0]
    AB[0]=B[0]-A[0]
    AB[1]=B[1]-A[1]
    AB[2]=B[2]-A[2]
    return AB
# fvect()


def fnorm(vect):
    """ Calculate the norm of a vector
    """
    n=(vect[0]*vect[0] + vect[1]*vect[1] + vect[2]*vect[2])**0.5
    return n
# fnorm()


def fthales(vector, d, A):
    """ Returns coordinates of point C, located on 'vector' at a distance d from point A (also on 'vector')
    """ 
    coeff = d / fnorm(vector)
    C=[]
    C.append( vector[0]*coeff + A[0])
    C.append( vector[1]*coeff + A[1])
    C.append( vector[2]*coeff + A[2])
    return C
# fthales()


def fis_in_space(diam_vect, C1, C2, point):
    """ Use spatial coordinates to determine if a point belongs to the space
        between the two planes which
        - have an orthogonal vector 'diam_vect'
        - and pass through points C1 ant C2
    """
    d1 = -(diam_vect[0]*C1[0] + diam_vect[1]*C1[1] + diam_vect[2]*C1[2])
    d2 = -(diam_vect[0]*C2[0] + diam_vect[1]*C2[1] + diam_vect[2]*C2[2])
    
    d = -(diam_vect[0]*point[0] + diam_vect[1]*point[1] + diam_vect[2]*point[2])
    
    var1=min(d1,d2)
    var2=max(d1,d2)

    if d > var1 and d < var2:
        return True
    else:
        return False
# fis_in_space()


def fhphob_hphil(data, asapound, afilter, dict_scale):
    """ Return the total number of M and S residues
    """
    hphob=0
    hphil=0

    for item in data:
        if item[4] > afilter:
            if asapound:
                accessibility=item[4]
            else:
                accessibility=1

            val=dict_scale[item[5]] #item[5] = residue type (3 letters)

            if val > 0:
                hphob+=val*accessibility

            else:
                val=fabs(val)
                hphil+=val*accessibility

    return hphil, hphob
# fhphob_hphil()


def fQ_value(hphil, hphob, hphiltot, hphobtot):
    """ Calculate the C value
    """
 
    if hphobtot <1:
        hphobtot=0.1
    if hphiltot <1:
        hphiltot+=1

    Qval = 0

    tot=hphobtot+hphiltot

    part_tot=hphob+hphil # numbers of exposed residues in the slice

    try:
        Qval = (hphob*(hphiltot-hphil) - hphil*(hphobtot-hphob))/ \
        ( part_tot*hphobtot*hphiltot*(tot-part_tot) )**0.5
    except:
        Qval = 0
    # A.G. 27.06
    if isinstance(Qval, complex):
        Qval = 0
    return Qval
# fQ_value()


def fscalar_product(u, v):
    """ Return the scalar product of two vectors u and v
    """
    scalar = u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
    return scalar
# fscalar_product()


def fproxim_axes(sphere_pts, Qbest, N, r):
    """ This function is aimed at finding proximate axes
    """
    imax=len(sphere_pts)
    j=4
    sphere_pts2=[]
    while len(sphere_pts2) < N:
        d=2*r/N+j # arbitrary
        sphere_pts2=[] 
        i=0   
        while i < imax:
            if fdistance(sphere_pts[i], Qbest[0]) < d:
                sphere_pts2.append(sphere_pts[i])
            i+=1
        j+=0.2
    return sphere_pts2
# fproxim_axes()


def frenaming(filename):
    """ Return the name of the structure file created by ANVIL, which
        contains the coordinates of the lipid bilayer
    """
    mylist = filename.split('.')
    newname = ''

    for element in mylist:
        if element == mylist[-2]:
            newname+=element
            newname+='_ANVIL.'
        elif element == mylist[-1]:
            newname+=element
        else:
            newname+=element
            newname+='.'

    return newname
# frenaming()


def fdraw_membrane(Qbest, infile, center, outputdir, sphere_pts, radius, density, af=False):
    """ Draw the lipid bilayer
    """
    
    newname = frenaming(infile)
    # [A.G.] TODO: reimplement os calls with native Python ...
    os.system('cp {} {}'.format(infile, newname))

    # [A.G. 5.06.22] #############################
    # drop END and ENDMDL at the end
    if af:
        #os.system("sed -i '$d' {} | sed '$d'".format(newname))
        # drop last 2 lines
        os.system("sed -i -e :a -e '$d;N;2,{}ba' -e 'P;D' {}".format(2, newname))
    ##############################################
    
    outfile=open(newname, 'a')
    
    # Remember:
    # Qbest = [point, C1, C2, hphob, hphil, tot, Qmax, center]

    point=Qbest[0]
    diam_vect=fvect(point, center)
    diam_norm=fnorm(diam_vect)
 
    array=[Qbest[1], Qbest[2]]
    for layer in array:
        d = -(diam_vect[0]*layer[0] + diam_vect[1]*layer[1] + diam_vect[2]*layer[2])
        i = -1000.0
        while i < 1000:
            j=-1000.0
            while j <1000:
                try:
                    # 'memb_atom': spatial coordinates of one membrane atom
                    memb_atom=[i,j,-(d +i*diam_vect[0] + j*diam_vect[1])/diam_vect[2]]
                except:
                    memb_atom=[i,j,0]

                if (memb_atom[0]-layer[0])**2 + (memb_atom[1]-layer[1])**2 + (memb_atom[2]-layer[2])**2 <= radius**2:
                    x=str(round(memb_atom[0],3))
                    y=str(round(memb_atom[1],3))
                    z=str(round(memb_atom[2],3))

                    # s1, s2, s3 = spacings
                    s1 = ' '*(8-len(x))
                    s2 = ' '*(8-len(y))
                    s3 = ' '*(32-len(z))

                    if layer == Qbest[1]:
                        outfile.write("HETATM90001  S   DUM  9001    {}{}{}{}{}{}\n".format(x,s1,y,s2,z,s3))
                    else:
                        outfile.write("HETATM90002  S   DUM  9002    {}{}{}{}{}{}\n".format(x,s1,y,s2,z,s3))

                j+=density
            i+=density
    
    # [A.G. 5.06.22] ##########
    if af:
        outfile.write('ENDMDL                                                                          ')
        outfile.write('\nEND                                                                             ')
    ###########################

    outfile.close()

    os.system('mv {} {}'.format(newname, outputdir))
# fdraw_membrane()


def ffind_all_axis(var):
    """ Explore the protein structure to find the membrane that maximizes the C value
    """

    center=var[0]
    sphere_pts=var[1]
    data=var[2]
    dict_scale=var[3]
    minthick=var[4]
    maxthick=var[5]
    step=var[6]
    hphiltot=var[7]
    hphobtot=var[8]
    asapound=var[9]
    afilter=var[10]

    membraneKept = []
    Qmax=0

    for point in sphere_pts:
        diam_vect=fvect(point, center)

        for i,rr in enumerate(diam_vect):
            diam_vect[i]=2*rr

        i=0
        diam_norm=fnorm(diam_vect)

        Qvartemp=[]

        while (i+step)  <  diam_norm:
            d_pointC1=i
            d_pointC2=i+step
            
            C1=fthales(diam_vect, d_pointC1, point )
            C2=fthales(diam_vect, d_pointC2, point )
            hphil=0
            hphob=0
            tot=0
            for item in data:
                if fis_in_space(diam_vect, C1, C2, item[:3]):
                    if item[4] >= afilter:
                        tot+=1

                        hphobtemp=0
                        hphiltemp=0

                        if dict_scale[item[5]] > 0: # item[5] = residue
                            hphobtemp=dict_scale[item[5]]
                        else:
                            hphiltemp=fabs(dict_scale[item[5]])

                        if asapound:
                            hphil+=hphiltemp*item[4] # item[4] = accessibility
                            hphob+=hphobtemp*item[4]
                        else:
                            hphil+=hphiltemp
                            hphob+=hphobtemp
 
            Qvartemp.append([C1, C2, hphob, hphil, tot])
            i+=step

        jmax = int((minthick/step)-1)

        width2 = 0

        while width2 <= maxthick:
            imax=(len(Qvartemp)-1)-jmax
            i=0
            while i < imax:
                C1=Qvartemp[i][0]
                C2=Qvartemp[(i+jmax)][1]
                j=0
                hphob=0
                hphil=0
                tot=0

                while j < jmax:
                    if j == 0:
                        hphil+=Qvartemp[i+j][3]/2
                        hphob+=Qvartemp[i+j][2]/2
                    elif j == jmax-1:
                        hphil+=Qvartemp[i+j][3]/2
                        hphob+=Qvartemp[i+j][2]/2
                    else:
                        hphil+=Qvartemp[i+j][3]
                        hphob+=Qvartemp[i+j][2]

                    tot+=Qvartemp[i+j][4]

                    j+=1

                if hphob != 0:
                    Qvaltest=fQ_value(hphil, hphob, hphiltot, hphobtot)
                    if Qvaltest >= Qmax:
                        Qmax=Qvaltest
                        membraneKept = [point, C1, C2, hphob, hphil, tot, Qmax, center]
                i+=1

            jmax+=1
            width2 = (jmax+1)*step

    return membraneKept
# ffind_all_axis()


def fwmembrane(Qbest, data, center, warn):
    """ Return the intra-membrane residues
    """
    # Between two intramembrane residues, if an amino acid is missing in the PDB file, it will be
    # considered as intramembrane

    result={}
    nonresult={}
    chains=fname_chains(data)
    diam_vect=fvect(Qbest[0], center)

    for name in chains:
        result[name]=[]
        nonresult[name]=[]
    for item in data:
        if fis_in_space(diam_vect, Qbest[1], Qbest[2], item[:3]):
            result[item[7]].append(item[6]) # item[7] = chain; item[6] = num res
        else:
            nonresult[item[7]].append(item[6])


    result2={}
    noname=[]

    for name in chains:
        try:
            result[name].sort()
            result2[name]=[[str(result[name][0])]]
        except:
            if result[name] == []:
                noname.append(name)
            else:
                sys.stderr.write("Error with fwmembrane()\n")
                return []
    if noname != []:
        if warn != 'mute':
             sys.stderr.write("No TM segment in chain(s) {}\n".format(noname))

    for name in chains:
        try:
            j=result[name][0]-1
            new=0
            for i, item in enumerate(result[name]):
                if item != j+1:
                    ok=False
                    for item2 in nonresult[name]:
                        if item2 == j+1:
                            ok = True
                    if ok:
                        j=item
                        result2[name][new].append(str(result[name][i-1]))
                        result2[name].append([str(j)])
                        new+=1
                    else:
                        j=item
                else:
                    j+=1
            if len(result2[name][new]) ==1:
                result2[name][new].append(str(result[name][i]))
        except:
            pass
            #if warn != 'mute':
                 #sys.stderr.write("Warning: \n")

    return result2
# fwmembrane()


def fname_chains(data):
    """ Get the names of all chains
    """
    chains=[]
    for item in data:
        name=item[7]
        new=True
        for ch in chains:
            if ch == name:
                new=False
        if new:
            chains.append(name)
    return chains
# fname_chains()


def ftransmembrane(Qbest, rawresults, data, center, adjust):
    """ Return the number and sequence positions of transmembrane segments
    """
    transmembrane={}
    vect=fvect(center, Qbest[0])
    d1= -( vect[0]*Qbest[1][0]+vect[1]*Qbest[1][1]+vect[2]*Qbest[1][2] )
    d2= -( vect[0]*Qbest[2][0]+vect[1]*Qbest[2][1]+vect[2]*Qbest[2][2] )
    dmax=max(d1, d2)
    dmin=min(d1, d2)
    
    for chain in rawresults:
        for segment in rawresults[chain]:
            # Residues outside the membrane
            segment[0]=int(segment[0])-1
            segment[1]=int(segment[1])+1
    
    for chain in rawresults:
        for segment in rawresults[chain]:
            coords1=coords2=0

            for item in data:
                if chain == item[7]:
                    if segment[0] == item[6]: # item[6]: num res
                        coords1=item[0:3]
                    elif segment[1] == item[6]:
                        coords2=item[0:3]

            if coords1 != 0 and coords2 !=0:
                d3= -( vect[0]*coords1[0]+vect[1]*coords1[1]+vect[2]*coords1[2])
                d4= -( vect[0]*coords2[0]+vect[1]*coords2[1]+vect[2]*coords2[2])
                if min(d3,d4) < dmin and max(d3,d4) > dmax:
                    segment[0]+=1
                    segment[1]-=1
                    try:
                        transmembrane[chain].append(segment)
                    except:
                        transmembrane[chain]=[segment]
    n=0
    for chain in transmembrane:
        n+=len(transmembrane[chain])

    reject = 0

    segments = [] # NB: different than segment[] without 's'
    # segments[i][0] = chain of segment i
    # segments[i][1] = extremity of segment i (res num)
    # segments[i][2] = the other extremity
    
    # Below: modification of 'transmembrane' dictionary
    for chain in transmembrane:
        for i,item in enumerate(transmembrane[chain]):
            if (abs(transmembrane[chain][i][0]-transmembrane[chain][i][1])+1) < adjust:
                reject = 1

            segments.append([chain, transmembrane[chain][i][0], transmembrane[chain][i][1]])

            transmembrane[chain][i][0]=str(transmembrane[chain][i][0])
            transmembrane[chain][i][1]=str(transmembrane[chain][i][1])
            transmembrane[chain][i]='-'.join(transmembrane[chain][i])

        transmembrane[chain]=']['.join(transmembrane[chain])
        transmembrane[chain]='['+transmembrane[chain]+']'
 
    return n, transmembrane, reject, segments
# ftransmembrane()


def fcompare(lines):
    """ Process the results from OPM or MEMEMBED:
        the thickness of the membrane and normal vector of the membrane plane
    """
    dataOPM=[]
    try:
        for item in lines:
            if item[0:6] == 'HETATM' and item[17:20] == 'DUM':
                dataOPM.append(item)    
        n=len(dataOPM)//2
        x1=float(dataOPM[n][30:38])
        y1=float(dataOPM[n][38:46])
        z1=float(dataOPM[n][46:54])
        if dataOPM[n][7:12] ==dataOPM[n-1][7:12]:
            n2=n-1
        elif dataOPM[n][7:12] ==dataOPM[n+1][7:12]:
            n2=n+1
        else:
            print('Cannot use data from OPM/MEMEMBED')
            n2=False
        if n2 != False:
            x2=float(dataOPM[n2][30:38])
            y2=float(dataOPM[n2][38:46])
            z2=float(dataOPM[n2][46:54])
            widthOPM = fdistance([x1,y1,z1], [x2,y2,z2])
            vectOPM  = fvect([x1,y1,z1], [x2,y2,z2])
            return widthOPM, vectOPM, n, n2, [x1,y1,z1], [x2,y2,z2]
        else:
            return 0,[0,0,0],0,0,[0,0,0],[0,0,0]
    except:
        return 0,[0,0,0],0,0,[0,0,0],[0,0,0]
# fcompare()


def fcalc_Qvar(nb_process, sphere_pts, center, data, dict_scale, minthick, maxthick, \
    step, hphil, hphob, asapound, afilter):
    """ To run the ffind_all_axis() function in multiprocessing (experimental!)
    """

    if nb_process == 1:
        Qvar = ffind_all_axis([center, sphere_pts, data, dict_scale, minthick, \
        maxthick, step, hphil, hphob, asapound, afilter])

    else:
        if not nb_process:
            nb_process = (cpu_count()/2)

        nsph=len(sphere_pts)//nb_process

        parall=[]
        i=0
        while i < nb_process-1:
            parall.append([center, sphere_pts[i*nsph:(i+1)*nsph],\
            data, dict_scale, minthick, maxthick, step, hphil, \
            hphob, asapound, afilter])
            i+=1
        parall.append([center, sphere_pts[i*nsph:], data, dict_scale, minthick, \
        maxthick, step, hphil, hphob, asapound, afilter])

        pool = Pool(processes=nb_process)
        result = pool.map_async(ffind_all_axis, parall)
        pool.close()
        pool.join()
        profile = result.get()
        pool.terminate()

        Qvar=[]
        Qmax=0

        for item in profile:
            if item[6] > Qmax:
                Qmax = item[6]
                Qvar = item

    return Qvar
# fcalc_Qvar()


def fdict_coords(lines, atomtype):
    coords={}

    for i in lines:
        if i[0:4] == 'ATOM' and (i[12:16].replace(' ','')) == atomtype:
            key=str(int(i[22:26].replace(' ','')))+'-'+(i[21].replace(' ','')) # num+'-'+chain
            value=[float(i[30:38]), float(i[38:46]), float(i[46:54])] # x, y, z
            coords[key]=value

    return coords
# fdict_coords()


def fquery_tilt(selseg, dict_coords, Qbest):

    print('')

    for query in selseg.split('_'):

        chain = query[0]
        myrange = query[1:]
        key1, key2 = myrange.split('-')
        key1+=('-'+chain)
        key2+=('-'+chain)

        if key1 in dict_coords and key2 in dict_coords:
            value1 = dict_coords[key1] # corresponding x, y, z coordinates (a list)
            value2 = dict_coords[key2]
            seg_vec = fvect(value1, value2) # vector representing the TM segment
            seg_norm = fdistance(value1, value2) # same as fnorm()
        
            # Remember: Qbest[0] and Qbest[7] define the best axis
            # i.e., orthogonal to the membrane plane
            scal = fscalar_product(fvect(Qbest[0], Qbest[7]), seg_vec)
            tiltangle = degrees(acos(scal/(fdistance(Qbest[0], Qbest[7])*seg_norm)))
            if tiltangle > 90:
                tiltangle = 180-tiltangle

            print('Query: {} - Segment: {} - Tilt: {}'.format(chain, myrange, round(tiltangle, 1)))

        else:
            print('Error: Segment {} does not exist'.format(myrange))
    
# fquery_tilt()


def fcalc_tilt(segments, dict_coords, Qbest):

    print('')

    # for each TM segment:
    for i, item in enumerate(segments):
        key1 = str(segments[i][1])+'-'+segments[i][0] # key = num+'-'+chain (e.g., 79-A or 42-B...)
        key2 = str(segments[i][2])+'-'+segments[i][0]
        value1 = dict_coords[key1] # corresponding x, y, z coordinates (a list)
        value2 = dict_coords[key2]
        seg_vec = fvect(value1, value2) # vector representing the TM segment
        seg_norm = fdistance(value1, value2) # same as fnorm()

        # Remember: Qbest[0] and Qbest[7] define the best axis
        # i.e., orthogonal to the membrane plane
        scal = fscalar_product(fvect(Qbest[0], Qbest[7]), seg_vec)
        tiltangle = degrees(acos(scal/(fdistance(Qbest[0], Qbest[7])*seg_norm)))
        if tiltangle > 90:
            tiltangle = 180-tiltangle

        myrange = str(segments[i][1])+'-'+str(segments[i][2])
        print('{} - Segment: {} - Tilt: {}'.format(segments[i][0], myrange, round(tiltangle, 1)))

    # This function returns nothing!

# fcalc_tilt()


def get_args():

    parser = argparse.ArgumentParser(description = 'ANVIL: Assignment aNd VIsualization of the Lipid bilayer')
    parser.add_argument('-i', dest='infile', type=str, help="input file (PDB file)", required=True)
    parser.add_argument('-n', dest='naccessdir', type=str, help="directory of the NACCESS program (default=./naccess2.1.1/)")
    parser.add_argument('-o', dest='outputdir', type=str, help="directory of output files (default=./)")
    parser.add_argument('--beta', dest="beta", action='store_true', help="for beta-strand protein structures")
    parser.add_argument('--martini', dest="martini", action='store_true', help="for MARTINI coarse-grained models")
    parser.add_argument('--inclin', dest="inclination", action='store_true', help="angle with OPM/MEMEMBED membrane")
    parser.add_argument('--min', dest="minthick", type=float, help="minimum thickness (default=20.0 Å)")
    parser.add_argument('--max', dest="maxthick", type=float, help="maximum thickness (default=40.0 Å)")
    parser.add_argument('--step', dest="step", type=float, help="step for the membrane search (default=1.0 Å)")
    parser.add_argument('--tilt', dest="tilt", action='store_true', help="calculate tilt for TM segments found")
    parser.add_argument('--selseg', dest="selseg", type=str, help="selected segments for tilt calculation (e.g., A532-544_B15-35)")
    parser.add_argument('-r', dest="radius", type=float, help="manually set the radius (Å) of the membrane disks")
    parser.add_argument('-d', dest="density", type=float, help="atom spacing in the membrane disks (default=2.0 Å)")
    # [A.G. 5.06.22]
    parser.add_argument('-af', dest="af", action='store_true', help="flag for alpha fold generated pdbs")

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    if args.outputdir:
        if args.outputdir[-1] != '/':
            args.outputdir+='/'

        if os.path.exists(args.outputdir) == False:
            sys.exit('Error: Outpout directory does not exist')

    if args.naccessdir:
        if os.path.exists(args.naccessdir) == False:
            sys.exit('Error: NACCESS directory does not exist')

    return args.naccessdir, args.outputdir, \
    args.infile, \
    args.inclination, args.minthick, args.maxthick, args.step, args.beta, args.selseg, args.tilt, args.radius, args.density, args.martini, \
    args.af
# get_args()


def main():

    asapound = False
    inclination = False
    adjust = 14
    nb_process = 1 # multiprocess is experimental
    afilter = 40
    asatype = 'relative'
    N = 350
    atomtype = 'CA'

    naccessdir, outputdir, infile, inclination, minthick, maxthick, step, beta, selseg, tilt, radius, density, martini, af = get_args()

    if not naccessdir:
        naccessdir = './naccess2.1.1/'

    if not outputdir:
        outputdir = './'

    if not minthick:
        minthick = 20.0

    if not maxthick:
        maxthick = 40.0

    if not step:
        step = 1.0

    if not density:
        density = 2.0

    if beta:
        adjust = 5

    if martini:
        atomtype = 'BB'

    print("\nProcessing: {}".format(infile))

    lines = freadpdb(infile, af=af) 
    naccess = fnaccess(infile, naccessdir, outputdir, asatype, atomtype)
    coordsATM = fcoords_atoms(lines, atomtype)

    print("{}s = {}".format(atomtype, len(coordsATM)))

    if len(coordsATM) == 0:
        sys.exit('Stop')

    dict_coords = fdict_coords(lines, atomtype) # key: num+'-'+chain, value:[x, y, z]

    center = fcenterofmass(coordsATM)
    r = fradius(coordsATM, center)

    if not radius:
        radius = r

    r*=1.2 # sphere radius will be 1.2x longer than protein radius
    data=fjoin_Coords_Naccess(coordsATM, naccess)

    dict_scale = {'ALA':+1,'ARG':-1,'ASN':-1,'ASP':-1,'CYS':+1,'GLN':-1,'GLU':-1, \
                  'GLY':+1,'HIS':+1,'ILE':+1,'LEU':+1,'LYS':-1,'MET':+1,'PHE':+1, \
                  'PRO':-1,'SER':+1,'THR':-1,'TRP':+1,'TYR':-1,'VAL':+1}

    hphil, hphob = fhphob_hphil(data, asapound, afilter, dict_scale)
    sphere_pts = fsphere_points(N,r,center) # returns N points

    # Try to find the best Q value
    Qbest = fcalc_Qvar(nb_process, sphere_pts, center, data, dict_scale, minthick, \
    maxthick, step, hphil, hphob, asapound, afilter)

    if Qbest == []:
        sys.exit('Error: Possible corruption of the structure file {}\n'.format(infile))
    else:
        print('Membrane found\nRefining inclination...')

    # Try to improve the results by finding a better axis
    sphere_pts = fsphere_points(30000,r,center)
    sphere_pts = fproxim_axes(sphere_pts, Qbest, N, r) # find N proximate axes among 30000
    Qbest2 = fcalc_Qvar(nb_process, sphere_pts, center, data, dict_scale, minthick, \
    maxthick, step, hphil, hphob, asapound, afilter)

    if Qbest2[6] >= Qbest[6]:
        Qbest=Qbest2
        print('done')
    else:
        print('better inclination not found')

    if adjust:
        nb_process = 1
        step = 0.3
        maxthick = fdistance(Qbest[1], Qbest[2])

        print('Adjusting thickness...')

        b = fwmembrane(Qbest, data, center, 'mute')
        cmax, d, e, f = ftransmembrane(Qbest, b, data, center, adjust)

        Qbestfinal = []

        while maxthick > minthick:
            Qtemp = fcalc_Qvar(nb_process, [Qbest[0]], center, data, dict_scale, minthick, \
            maxthick, step, hphil, hphob, asapound, afilter)
            if (len(Qtemp)) > 0:
                b = fwmembrane(Qtemp, data, center, 'mute')
                c, d, e, f = ftransmembrane(Qtemp, b, data, center, adjust)
                if c > cmax and e == 0:
                    cmax = c
                    Qbestfinal = Qtemp
            maxthick-=step

        if len(Qbestfinal) > 0:
            Qbest = Qbestfinal
            print('done')
        else:
            print('better thickness not found')
    

    rawresults = fwmembrane(Qbest, data, center, 'warn')
    a = 0 # number of TM segments
    a, transmembr, reject, segments = ftransmembrane(Qbest, rawresults, data, center, adjust)

    prettydict = str(transmembr).replace("{","").replace("}","").replace("'","").replace(",","\n")
    print("\nNumber of TM segments = {}".format(a))
    print("\nSequence positions:\n {}".format(prettydict))
    print("\nMembrane width = {} Å".format(fdistance(Qbest[1], Qbest[2])))

    # Write a PDB file that includes the membrane coordinates
    fdraw_membrane(Qbest, infile, center, outputdir, sphere_pts, radius, density, af=af)

    if tilt:
        if a > 0:
            fcalc_tilt(segments, dict_coords, Qbest)
        else:
            print('Tilt = n/a')

    if selseg:
        try:
            fquery_tilt(selseg, dict_coords, Qbest)
        except:
            sys.stderr.write("Error with --selseg: Check the input format")

    # Compare with OPM/MEMEMBED
    if inclination:
        widthOPM, vectOPM, n, n2, C1, C2 = fcompare(lines)
        Qopm = [C1,C1,C2,0,0,0,0,C2] # the zeros won't be used anyway
        if n != n2:
            scal=fscalar_product(fvect(Qbest[0], Qbest[7]), vectOPM)
            angleOPM = round(degrees( acos( scal/(fdistance(Qbest[0], Qbest[7])*widthOPM ))), 1)

            while angleOPM > 90:
                angleOPM = 180.0-angleOPM

            print("\nAngle with OPM/MEMEMBED = {}°".format(angleOPM))

        else:
            print("\nNo data from OPM/MEMEMBED in {}".format(infile))
    
    print("\nDone")
# main()


if __name__ == '__main__':

    from math import pi, acos, sin, cos, fabs, degrees
    #from multiprocessing import Pool, cpu_count # experimental
    import sys, os, argparse

    main()
