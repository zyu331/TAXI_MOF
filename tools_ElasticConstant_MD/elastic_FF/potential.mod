pair_style      lj/cut 12.5
bond_style      harmonic
angle_style     hybrid cosine/periodic fourier
dihedral_style  harmonic
improper_style  fourier


special_bonds   lj 0.0 0.0 1.0
pair_modify     tail yes mix arithmetic
include data.pair

# Setup neighbor style
neighbor 1.0 nsq
neigh_modify once no every 1 delay 0 check yes

# Setup minimization style
min_style	     cg
min_modify	     dmax ${dmax} line quadratic

# Setup output
thermo		1
thermo_style custom step temp pe press pxx pyy pzz pxy pxz pyz lx ly lz vol
thermo_modify norm no
