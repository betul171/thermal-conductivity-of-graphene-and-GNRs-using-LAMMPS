# Calculate the thermal conductivity of silicon bulk.
# This script is adapted from an example input script from the official LAMMPS documentation.

from lammps import lammps

lmp = lammps()

command = """ 

units       metal
variable    T equal 450
variable    V equal vol
variable    dt equal 0.001
variable    p equal 200         # correlation length
variable    s equal 10          # sample interval
variable    d equal $p*$s       # dump interval
# variable    r equal 100000     # Run
# variable    ir equal 8000     # Equilibration run 
variable    r equal 10000     # Run
variable    ir equal 8000     # Equilibration run 


variable    kB equal 1.3806504e-23    # [J/K] Boltzmann
variable    eV2J equal 1.60217663e-19
variable    A2m equal 1.0e-10
variable    ps2s equal 1.0e-12
variable    convert equal ${eV2J}*${eV2J}/${ps2s}/${A2m}

dimension    3
boundary     p p p
lattice      diamond 5.43 orient x 1 0 0 orient y 0 1 0 orient z 0 0 1
region       box block 0 5 0 5 0 5
create_box   1 box
create_atoms 1 box
mass         1 28.0855
pair_style   sw
pair_coeff   * * Si.sw Si
minimize 1.0e-6 1.0e-8 1000 10000
timestep ${dt}
thermo $d
thermo_style custom step temp pe ke etotal
dump 1 all custom 100 si.dump id type x y z vx vy vz


velocity all create 10.0 48279
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
dump myDump all atom 100 dump_reset_silicon.lammpstrj


compute ke all ke/atom
compute pe all pe/atom
compute stress all stress/atom NULL virial
compute flux all heat/flux ke pe stress
#fix temp_avg all ave/time 10 10 100 temp file temp_avg.txt
variable     Jx equal c_flux[1]/vol
variable     Jy equal c_flux[2]/vol
variable     Jz equal c_flux[3]/vol

fix          JJ all ave/correlate $s $p $d &
             c_flux[1] c_flux[2] c_flux[3] type auto file J0Jt.dat ave running

variable     scale equal ${convert}/${kB}/$T/$T/$V*$s*${dt}
variable     k11 equal trap(f_JJ[3])*${scale}
variable     k22 equal trap(f_JJ[4])*${scale}
variable     k33 equal trap(f_JJ[5])*${scale}

thermo_style custom step temp v_k11 v_k22 v_k33
run          $r

variable     k equal (v_k11+v_k22+v_k33)/3.0
print        "average thermal conductivity: $k[W/mK] @ $T K"
print        "kxx: ${k11}[W/mK]"
print        "kyy: ${k22}[W/mK]"
print        "kzz: ${k33}[W/mK]"     
"""

lmp.commands_string(command)

lmp.close()
