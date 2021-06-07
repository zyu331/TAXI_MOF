
# real units, elastic constants in GPa
variable up equal 1.0e-2
variable atomjiggle equal 1.0e-5
units		real
variable cfac equal 1.01325e-4
variable cunits string GPa
atom_style      full

# Define minimization parameters
variable etol equal 1.0e-10 
variable ftol equal 1.0e-10
variable maxiter equal 100000
variable maxeval equal 10000000
variable dmax equal 1.0e-3
bond_style      harmonic
angle_style     hybrid cosine/periodic fourier
dihedral_style  harmonic
improper_style  fourier
box tilt large
 read_data       data.MOF
change_box all triclinic