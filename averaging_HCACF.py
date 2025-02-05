# Averege over 10 HCACF data files with different initial conditions (different random seed numbers).
# The aim is to create "avereged_HCACF.dat" file.
### Finally, compute thermal conductivity.

fID1 = open("HCACF.dat")

# Read all lines from the file
lines1 = fID1.readlines()

time = []
Jy_1 = []
# Process only the last 200 lines
for line in lines1[-200:]:
    s = line.split()
    time.append(float(s[1]))
    Jy_1.append(float(s[4]))

fID1.close()

############################################################################################

fID2 = open("HCACF2.dat")

# Read all lines from the file
lines2 = fID2.readlines()

time = []
Jy_2 = []
# Process only the last 200 lines
for line in lines2[-200:]:
    s = line.split()
    time.append(float(s[1]))
    Jy_2.append(float(s[4]))

fID2.close()

############################################################################################

fID3 = open("HCACF3.dat")

# Read all lines from the file
lines3 = fID3.readlines()

time = []
Jy_3 = []
# Process only the last 200 lines
for line in lines3[-200:]:
    s = line.split()
    time.append(float(s[1]))
    Jy_3.append(float(s[4]))

fID3.close()

############################################################################################

fID4 = open("HCACF4.dat")

# Read all lines from the file
lines4 = fID4.readlines()

time = []
Jy_4 = []
# Process only the last 200 lines
for line in lines4[-200:]:
    s = line.split()
    time.append(float(s[1]))
    Jy_4.append(float(s[4]))

fID4.close()

############################################################################################

fID5 = open("HCACF5.dat")

# Read all lines from the file
lines5 = fID5.readlines()

time = []
Jy_5 = []
# Process only the last 200 lines
for line in lines5[-200:]:
    s = line.split()
    time.append(float(s[1]))
    Jy_5.append(float(s[4]))

fID5.close()

############################################################################################

fID6 = open("HCACF6.dat")

# Read all lines from the file
lines6 = fID6.readlines()

time = []
Jy_6 = []
# Process only the last 200 lines
for line in lines6[-200:]:
    s = line.split()
    time.append(float(s[1]))
    Jy_6.append(float(s[4]))

fID6.close()

############################################################################################

fID7 = open("HCACF7.dat")

# Read all lines from the file
lines7 = fID7.readlines()

time = []
Jy_7 = []
# Process only the last 200 lines
for line in lines7[-200:]:
    s = line.split()
    time.append(float(s[1]))
    Jy_7.append(float(s[4]))

fID7.close()

############################################################################################

fID8 = open("HCACF8.dat")

# Read all lines from the file
lines8 = fID8.readlines()

time = []
Jy_8 = []
# Process only the last 200 lines
for line in lines8[-200:]:
    s = line.split()
    time.append(float(s[1]))
    Jy_8.append(float(s[4]))

fID8.close()

############################################################################################

fID9 = open("HCACF9.dat")

# Read all lines from the file
lines9 = fID9.readlines()

time = []
Jy_9 = []
# Process only the last 200 lines
for line in lines9[-200:]:
    s = line.split()
    time.append(float(s[1]))
    Jy_9.append(float(s[4]))

fID9.close()

############################################################################################

fID10 = open("HCACF10.dat")

# Read all lines from the file
lines10 = fID10.readlines()

time = []
Jy_10 = []
# Process only the last 200 lines
for line in lines10[-200:]:
    s = line.split()
    time.append(float(s[1]))
    Jy_10.append(float(s[4]))

fID10.close()

##################################################################################
##################################################################################
##################################################################################

Jy = []
for i in range(200):
    Jy_i = Jy_1[i] + Jy_2[i] + Jy_3[i] + Jy_4[i] + Jy_5[i] + Jy_6[i] + Jy_7[i] + Jy_8[i] + Jy_9[i] + Jy_10[i]
    Jy_k = Jy_i / 10
    Jy.append(Jy_k)



fID = open("avereged_HCACF.dat","w")

time_in_ps = [t * 0.001 for t in time]  ## in ps right???

for i in range(200):
    fID.write("{} {}".format(time_in_ps[i], Jy[i]))

fID.close()


T = 300
kB = 1.3806504e-23  # Boltzmann constant [J/K]
eV2J = 1.60217663e-19  # eV to J conversion factor
A2m = 1.0e-10  # Angstrom to meter conversion factor
ps2s = 1.0e-12  # picoseconds to seconds conversion factor


dt = 0.001  # time step
p = 200  # correlation length
s = 10  # sample interval
d = p * s  # dump interval
r = 10000  # Run
ir = 35000  # Equilibration run

x_0 = -1
x_f = 15
y_0 = 0
y_f = 25
x = x_f - x_0
y = y_f - y_0
thickness = 3.35
V = x * y * thickness

convert = eV2J * eV2J / (ps2s * A2m)
scale = convert / (kB * T * T * V * s * dt)

# Numerically integrate Jy using the trapezoidal rule
integrated_Jy = np.trapz(Jy, time_in_ps)

thermal_conductivity = integrated_Jy * scale # along y



