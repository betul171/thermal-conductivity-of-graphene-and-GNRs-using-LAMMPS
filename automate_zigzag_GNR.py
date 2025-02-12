# The LAMMPS script runs for 10 different random seed numbers, and generates a HCACF file for each seed.

from lammps import lammps

# List of 10 random seed numbers
seeds = [302184, 174592, 537801, 248304, 766334, 845243, 91215300, 4036723, 127366487, 754621503]


# LAMMPS initialization
lmp = lammps()

# LAMMPS command template
command_template = """
units       metal
variable    T equal 300
variable    dt equal 0.001
variable    p equal 200        # correlation length
variable    s equal 10          # sample interval
variable    d equal $p*$s       # dump interval
variable    r equal 10000     # Run
variable    ir equal 35000     # Equilibration run 

variable    kB equal 1.3806504e-23    # [J/K] Boltzmann
variable    eV2J equal 1.60217663e-19
variable    A2m equal 1.0e-10
variable    ps2s equal 1.0e-12
variable    convert equal ${eV2J}*${eV2J}/${ps2s}/${A2m}

dimension 3
boundary p p p 
atom_style atomic
region box block -15 30 0 101 -20 20
create_box 1 box
mass 1 12.011
lattice custom 1.42 a1 3 0 0 a2 0 1.732 0 a3 0 0 20 &
basis 0 0 0 &
basis 0.333 0 0 &
basis 0.5 0.5 0 &
basis 0.833 0.5 0
variable x_0 equal -1
variable x_f equal 15
variable y_0 equal ylo
variable y_f equal yhi
variable x equal ${x_f}-${x_0}
variable y equal ${y_f}-${y_0}
variable fixed1_yf equal ${y_0}+1 # fixed atoms' final y value
variable fixed2_yi equal ${y_f}-1 # fixed atoms' final y value
variable thickness equal 3.35
variable V equal $x*$y*${thickness}
region graphene block ${x_0} ${x_f} ${y_0} ${y_f} -0.1 0.1 units box
create_atoms 1 region graphene

region fixed_atoms1 block ${x_0} ${x_f} ${y_0} ${fixed1_yf} -0.1 0.1 units box 
region fixed_atoms2 block ${x_0} ${x_f} ${fixed2_yi} ${y_f} -0.1 0.1 units box 
group fixed_atoms1 region fixed_atoms1
group fixed_atoms2 region fixed_atoms2
group fixed_atoms union fixed_atoms1 fixed_atoms2

pair_style airebo 3.0
pair_coeff * * CH.airebo C
minimize 1.0e-25 1.0e-25 5000 10000
timestep ${dt}
thermo $d
thermo_style custom step temp pe ke etotal
dump 1 all custom 10000 nrr.dump id type x y z vx vy vz

velocity all create 10.0 ${seed}    # Random seed
fix tmp1 all nvt temp 10.0 10.0 0.001
run ${ir}
unfix tmp1
fix tmp2 all nvt temp 10.0 300 0.001
run ${ir}
unfix tmp2
fix tmp3 all nvt temp 300 300 0.001
run ${ir}
undump 1

reset_timestep 0
dump myDump all atom 10000 dump_reset_nr_fixed.lammpstrj

compute ke all ke/atom
compute pe all pe/atom
compute stress all stress/atom NULL virial
compute flux all heat/flux ke pe stress
variable     Jx equal c_flux[1]/$V
variable     Jy equal c_flux[2]/$V
variable     Jz equal c_flux[3]/$V

fix          JJ all ave/correlate $s $p $d &
             c_flux[1] c_flux[2] c_flux[3] type auto file HCACF${index}.dat ave running

variable     scale equal ${convert}/${kB}/$T/$T/$V*$s*${dt}
variable     k11 equal trap(f_JJ[3])*${scale}
variable     k22 equal trap(f_JJ[4])*${scale}
variable     k33 equal trap(f_JJ[5])*${scale}

thermo_style custom step temp v_k11 v_k22 v_k33
run          $r
   
"""

# Loop through all seed values and create unique file names
for i, seed in enumerate(seeds, start=1):
    lmp = lammps()
    command_with_seed = command_template.replace("${seed}", str(seed)).replace("${index}", str(i))
    lmp.commands_string(command_with_seed)
    lmp.close()

lmp.close()
