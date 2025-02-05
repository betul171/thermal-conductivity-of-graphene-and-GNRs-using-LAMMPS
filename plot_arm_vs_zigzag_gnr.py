# PLOT
#Thermal conductivity of rectangular GNRs with armchair (solid line) and zigzag (dashed line) long edges versus their length at 300K.

import matplotlib.pyplot as plt

data = {
    "zigzag": {
        2.5: 32,
        5: 33,
        7.5: 22,
        10.1: 26,
        12.6: 74,
        14.9: 248,
        17.5: 17,
        20: 223,
        22.5: 55,
        25.1: 18,
        27.4: 115,
        30: 73
    },
    "armchair": {
        2.5: 78,
        5: 41,
        7.6: 92,
        9.8: 126,
        12.3: 3,
        14.9: 36,
        17.4: 249,
        20: 208,
        22.5: 160,
        25.1: 74,
        27.2: 233,
        29.8: 41
    }
}

length_zigzag = list(data["zigzag"].keys())
tc_zigzag = list(data["zigzag"].values())

length_armchair = list(data["armchair"].keys())
tc_armchair = list(data["armchair"].values())

plt.xlabel('Length (nm)', fontsize=11, fontweight='bold')
plt.ylabel('Thermal conductivity (W / m·K)', fontsize=11, fontweight='bold' )
plt.plot(length_zigzag , tc_zigzag, label='Zigzag', linestyle='--', marker='o', color='blue')
plt.plot(length_armchair , tc_armchair, label='Armchair', linestyle='-', marker='s', color='red')
plt.minorticks_on()
plt.legend()
plt.xticks(fontsize=10, fontweight='bold')
plt.yticks(fontsize=10, fontweight='bold')
plt.savefig('/Users/betulalkan/Desktop/thermal_conductivity_plot.png', dpi=300)
plt.show()

