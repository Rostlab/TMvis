       `shhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh+       
         :dd:.```./mh/-M+.`.yMmNs-``+MM/.`-dm: .odo-`:MMh.`-MmNy:.`.sNs.        
          .Mo. . `+M/ -M+`  `+moM:  +MNd. .dd    dd..sNNh. -Mo`Mo. .Ms          
          dm. .o. `sN--M+` .``:mM:  +MoMo` /M+  :M+`.Nymh. -Mo Mo. -Ms          
         +M: `sMy. .ym/M+` -+` -h:  +M:yN- .yN. hm.`+M-mh. -Mo Mo. -Ms          
 `......-Ny-./MhMs-.:mNMo.`:Mo` .-  +M:.Ny` -Ny:M+`.dh mh. -Mo Mo. -Ms       `-.
`-+hdNNdhhhhhhhyhhyhhhNMms.:MMy. `  +M: oM: `oMmm. /M: mh. -Mo Mo. -Ms    `-ohM+
    `:+mdyo:-``  ``-ohdN+` :Mhmh-   +M: `mh. .dMo`.hm  mh. -Ms`Mo. -My-/oyhy+/M+
      -Ns:odmyo/-/ymd:`//` .oddMd:-:oM:  /M+` /m- :M+  mh-`.ymhMo. -dhyo/:`` :M+
     `dd. -mh.:syh/sN/:/oossyhmmdhyysy.   hm- .:``sN`  ydhhyyymMs-````       :M+
     oN: `/No`      shso/:--...``         .Ns` . -Ns    ``.--/+oyyhy+:.`     :M+
    -No`  `/mh.                            sN:  `oM.             `.-+syho:.` :M+
   `dd. `-+hds-                            `mh. .my                   `./yds:-M+
   oN::ohho-`                               /M/`/M-                       .+hdM+
  -Mmdho-`                                   dm-hd                          `/h+
 `dh+.`                                      -MyM/                            ``
 ``                                           yMm                    Version 1.0
                                              `No                               
                                               +`                               


            ANVIL: Assignment aNd VIsualization of the Lipid bilayer           



============
Description:
============
Briefly, ANVIL is aimed at assigning membrane boundaries to a protein 3D structure,
by using the spatial coordinates of the alpha carbons. It can also process coarse-
grained models of protein structures. The algorithm follows an approach that treats
the problem of membrane assignment as a binary classification. The software is
implemented in Python and intended to be run with the PyPy interpreter. In output,
ANVIL writes a PDB file that contains the coordinates of the membrane, which can
be visualized with softwares, such as RasMol or PyMOL. The ANVIL source code is
freely available under a CeCILL License.


====================
System requirements:
====================
Operating system: Linux, Mac OS X

The PyPy Python interpreter
e.g., under a Debian-like Linux distribution:
$ sudo apt-get install pypy

Naccess 2.1.1 (S.J. Hubbard, 1996)
http://www.bioinf.manchester.ac.uk/naccess/
This program requires to install Fortran 77 and the C shell:
$ sudo apt-get install fort77
$ sudo apt-get install csh


=================
How to run ANVIL?
=================
In a command-line interface:

ANVIL.py [-h] -i INFILE [-n NACCESSDIR] [-o OUTPUTDIR] [--beta]
         [--martini] [--inclin] [--min MINTHICK] [--max MAXTHICK]
         [--step STEP] [--tilt] [--selseg SELSEG] [-r RADIUS]
         [-d DENSITY]

optional arguments:
  -h, --help       show this help message and exit
  -i INFILE        input file (PDB file)
  -n NACCESSDIR    directory of the NACCESS program (default=./naccess2.1.1/)
  -o OUTPUTDIR     directory of output files (default=./)
  --beta           for beta-strand protein structures
  --martini        for MARTINI coarse-grained models
  --inclin         angle with OPM/MEMEMBED membrane
  --min MINTHICK   minimum thickness (default=20.0 Å)
  --max MAXTHICK   maximum thickness (default=40.0 Å)
  --step STEP      step for the membrane search (default=1.0 Å)
  --tilt           calculate tilt for TM segments found
  --selseg SELSEG  selected segments for tilt calculation (e.g., A532-544_B15-35)
  -r RADIUS        manually set the radius (Å) of the membrane disks
  -d DENSITY       atom spacing in the membrane disks (default=2.0 Å)


=========
Features:
=========
If you are not satisfied with the aspect of the lipid bilayer in output,
you can run ANVIL using the -d option with a value >2.0 (e.g., 1.0 or 0.5)


=================
Versions history:
=================
v1.0 (G. Postic)
- martini
- TM segments maximization
- tilt-related functions
- membrane aspect functions
- PyPy implementation

v0.1 (Y. Ghouzam & V. Guiraud)
- exploration procedure
- membrane inclination
- parallelization



Guillaume Postic, Ph.D.
guillaume.postic@univ-paris-diderot.fr
