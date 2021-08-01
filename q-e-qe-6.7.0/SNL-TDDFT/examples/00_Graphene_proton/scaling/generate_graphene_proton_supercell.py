# 2021-08-01: first version of script created by Andrew Baczewski @ SNL (adbacze@sandia.gov)
import numpy as np
import sys

# this script will generate the input file for a SCF calculation
# consisting of a graphene sheet tesselated a specific number of times with a proton above it
#
# usage is: python3 generate_graphene_proton_supercell.py <# of primitive cells on direction A> <# of primitive cells on direction B>
# this will create a file, graphene_proton_nCellsA_nCellsB-scf.in (nCellsA = # of primitive cells on direction A, etc.)
# once the wave function and rho cutoffs are specified, you can execute the scf calculation specified by this file
# note that it is easy to edit the cutoffs without opening the file in something like vim using the sed utility (sed = 'stream editor'),
# e.g., sed -i "s/ECUTWFC = 100/ECUTWFC = 200/" graphene_proton_nCellsA_nCellsB-scf.in sets the wave function cutoff to 200 in the specified input file

# note that the pseudo directory in the input file is in the directory ABOVE the working directory
# i.e., '../pseudo/'
# feel empowered to change this with hard coding (or anything else for that matter)

# the two command-line inputs are the number of 2-atom primitive cells along the A and B lattice vectors
nCellsA = np.int(sys.argv[1])
nCellsB = np.int(sys.argv[2])

# the lattice constant for graphene (in this 2-atom/cell picture)
rGrapheneLattice = 4.92
# the a and b lattice constants for the supercell will be proportionally larger
rA = nCellsA*rGrapheneLattice
rB = nCellsB*rGrapheneLattice
# note: the c lattice constant is the amount of vacuum between images...I've hard coded it to be 20 Angstroms, but you can change that
rC = 20.0

# the total number of atoms is determined by the number of primitive cells in the tesselation (remember that there are 2 atoms per primitive cell)
# plus 1 for the proton!
nAtoms = 2*nCellsA*nCellsB + 1

# the A and B coordinates of the 1st atom in any given cell
# (2nd atom is 2x this)
# ...easiest to just look at the part of the code where we write the C atom positions
rAtomAEps = 1.0/(3.0*nCellsA)
rAtomBEps = 1.0/(3.0*nCellsB)

# this is the filename of the file that we're creating
sFilename = 'graphene_proton_'+str(nCellsA)+'x'+str(nCellsB)+'-scf.in'

# open the file...
fSCFInput = open(sFilename, 'w')

# the control part of the input namelist
fSCFInput.write("&control\n")
fSCFInput.write("    prefix = 'graphene_proton_"+str(nCellsA)+"x"+str(nCellsB)+"'\n")
fSCFInput.write("    calculation = 'scf'\n")
fSCFInput.write("    restart_mode = 'from_scratch'\n")
fSCFInput.write("    pseudo_dir = '../pseudo/'\n")
fSCFInput.write("    outdir = './scratch/'\n") 
fSCFInput.write("    tstress = .false.\n")
fSCFInput.write("    tprnfor = .true.\n")
fSCFInput.write("/\n")

# the system part of the input namelist
fSCFInput.write("&system\n")
fSCFInput.write("    ibrav = 12\n")
fSCFInput.write("    a = "+str(rA)+"\n")
fSCFInput.write("    b = "+str(rB)+"\n")
fSCFInput.write("    c = "+str(rC)+"\n")
fSCFInput.write("    cosab = -0.5\n")
fSCFInput.write("    nat = "+str(nAtoms)+"\n")
fSCFInput.write("    ntyp = 2\n")
fSCFInput.write("    ecutwfc = 100.0\n")
fSCFInput.write("    ecutrho = 400.0\n")
fSCFInput.write("    occupations = 'smearing'\n")
fSCFInput.write("    smearing = 'gaussian'\n")
fSCFInput.write("    degauss = 0.01\n")
fSCFInput.write("/\n")

# the electrons part of the input namelist
fSCFInput.write("&electrons\n")
fSCFInput.write("    diagonalization = 'davidson'\n")
fSCFInput.write("    conv_thr = 1.d-6\n")
fSCFInput.write("/\n")

# the atomic species part of the input namelist
fSCFInput.write("ATOMIC_SPECIES\n")
fSCFInput.write("H    1   H.upf\n")
fSCFInput.write("C   12   C.upf\n")
fSCFInput.write("\n\n")

# the k point part of the input namelist
fSCFInput.write("K_POINTS automatic\n")
fSCFInput.write("1 1 1    0 0 0\n")
fSCFInput.write("\n")

# the atomic positions part of the input namelist
fSCFInput.write("ATOMIC_POSITIONS crystal\n")
fSCFInput.write("H   0.500   0.500   0.000\n")
# loop over the individual cells (A = slow, B = fast)
for iAIdx in np.arange(nCellsA):
        # slow loop is over cells in A direction
        # compute the A positions of the two atoms in each cell
	rAtom1A = (iAIdx/nCellsA) + rAtomAEps
	rAtom2A = (iAIdx/nCellsA) + 2.0*rAtomAEps

	for iBIdx in np.arange(nCellsB):
                # fast loop is over cells in B direction
		# compute the B positions of the two atoms in each cell
		rAtom1B = (iBIdx/nCellsB) + rAtomBEps
		rAtom2B = (iBIdx/nCellsB) + 2.0*rAtomBEps

		# then write the two atomic positions per individual cell
		fSCFInput.write("C   "+str(rAtom1A)+"   "+str(rAtom1B)+"   0.500\n")
		fSCFInput.write("C   "+str(rAtom2A)+"   "+str(rAtom2B)+"   0.500\n")

# close the file...
fSCFInput.close()
