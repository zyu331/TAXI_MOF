# NOTE: This script can be modified for different atomic structures, 
# units, etc. See in.elastic for more info.
#

# Define the finite deformation size. Try several values of this
# variable to verify that results do not depend on it.
variable up equal 1.0e-2
 
# Define the amount of random jiggle for atoms
# This prevents atoms from staying on saddle points
variable atomjiggle equal 1.0e-5

# Uncomment one of these blocks, depending on what units
# you are using in LAMMPS and for output


# real units, elastic constants in GPa
units		real
atom_style      full
variable cfac equal 1.01325e-4
variable cunits string GPa

# Define minimization parameters
variable etol equal 1.0e-11 
variable ftol equal 1.0e-10
variable maxiter equal 10000
variable maxeval equal 100000
variable dmax equal 1.0e-2

# generate the box and atom positions using a diamond lattice

boundary	p p p
box tilt        large
read_data       data.XXXX
change_box all triclinic

