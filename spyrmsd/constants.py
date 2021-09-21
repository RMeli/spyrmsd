from typing import Dict

connectivity_tolerance: float = 0.4

# Values from QCElemental
anum_to_mass: Dict = {
    1: 1.00782503223,
    2: 4.00260325413,
    3: 7.0160034366,
    4: 9.012183065,
    5: 11.00930536,
    6: 12.0,
    7: 14.00307400443,
    8: 15.99491461957,
    9: 18.99840316273,
    10: 19.9924401762,
    11: 22.989769282,
    12: 23.985041697,
    13: 26.98153853,
    14: 27.97692653465,
    15: 30.97376199842,
    16: 31.9720711744,
    17: 34.968852682,
    18: 39.9623831237,
    19: 38.9637064864,
    20: 39.962590863,
    21: 44.95590828,
    22: 47.94794198,
    23: 50.94395704,
    24: 51.94050623,
    25: 54.93804391,
    26: 55.93493633,
    27: 58.93319429,
    28: 57.93534241,
    29: 62.92959772,
    30: 63.92914201,
    31: 68.9255735,
    32: 73.921177761,
    33: 74.92159457,
    34: 79.9165218,
    35: 78.9183376,
    36: 83.9114977282,
    37: 84.9117897379,
    38: 87.9056125,
    39: 88.9058403,
    40: 89.9046977,
    41: 92.906373,
    42: 97.90540482,
    43: 97.9072124,
    44: 101.9043441,
    45: 102.905498,
    46: 105.9034804,
    47: 106.9050916,
    48: 113.90336509,
    49: 114.903878776,
    50: 119.90220163,
    51: 120.903812,
    52: 129.906222748,
    53: 126.9044719,
    54: 131.9041550856,
    55: 132.905451961,
    56: 137.905247,
    57: 138.9063563,
    58: 139.9054431,
    59: 140.9076576,
    60: 141.907729,
    61: 144.9127559,
    62: 151.9197397,
    63: 152.921238,
    64: 157.9241123,
    65: 158.9253547,
    66: 163.9291819,
    67: 164.9303288,
    68: 165.9302995,
    69: 168.9342179,
    70: 173.9388664,
    71: 174.9407752,
    72: 179.946557,
    73: 180.9479958,
    74: 183.95093092,
    75: 186.9557501,
    76: 191.961477,
    77: 192.9629216,
    78: 194.9647917,
    79: 196.96656879,
    80: 201.9706434,
    81: 204.9744278,
    82: 207.9766525,
    83: 208.9803991,
    84: 208.9824308,
    85: 209.9871479,
    86: 222.0175782,
    87: 223.019736,
    88: 226.0254103,
    89: 227.0277523,
    90: 232.0380558,
    91: 231.0358842,
    92: 238.0507884,
    93: 237.0481736,
    94: 244.0642053,
    95: 243.0613813,
    96: 247.0703541,
    97: 247.0703073,
    98: 251.0795886,
    99: 252.08298,
    100: 257.0951061,
    101: 258.0984315,
    102: 259.10103,
    103: 266.11983,
    104: 267.12179,
    105: 268.12567,
    106: 271.13393,
    107: 270.13336,
    108: 269.13375,
    109: 278.15631,
    110: 281.16451,
    111: 282.16912,
    112: 285.17712,
    113: 286.18221,
    114: 289.19042,
    115: 289.19363,
    116: 293.20449,
    117: 294.21046,
}

# Data from QCElemental (in angstroms)
anum_to_covalentradius: Dict = {
    1: 0.31,
    2: 0.28,
    3: 1.28,
    4: 0.96,
    5: 0.84,
    6: 0.76,
    7: 0.71,
    8: 0.66,
    9: 0.57,
    10: 0.58,
    11: 1.66,
    12: 1.41,
    13: 1.21,
    14: 1.11,
    15: 1.07,
    16: 1.05,
    17: 1.02,
    18: 1.06,
    19: 2.03,
    20: 1.76,
    21: 1.7,
    22: 1.6,
    23: 1.53,
    24: 1.39,
    25: 1.61,
    26: 1.52,
    27: 1.5,
    28: 1.24,
    29: 1.32,
    30: 1.22,
    31: 1.22,
    32: 1.2,
    33: 1.19,
    34: 1.2,
    35: 1.2,
    36: 1.16,
    37: 2.2,
    38: 1.95,
    39: 1.9,
    40: 1.75,
    41: 1.64,
    42: 1.54,
    43: 1.47,
    44: 1.46,
    45: 1.42,
    46: 1.39,
    47: 1.45,
    48: 1.44,
    49: 1.42,
    50: 1.39,
    51: 1.39,
    52: 1.38,
    53: 1.39,
    54: 1.4,
    55: 2.44,
    56: 2.15,
    57: 2.07,
    58: 2.04,
    59: 2.03,
    60: 2.01,
    61: 1.99,
    62: 1.98,
    63: 1.98,
    64: 1.96,
    65: 1.94,
    66: 1.92,
    67: 1.92,
    68: 1.89,
    69: 1.9,
    70: 1.87,
    71: 1.87,
    72: 1.75,
    73: 1.7,
    74: 1.62,
    75: 1.51,
    76: 1.44,
    77: 1.41,
    78: 1.36,
    79: 1.36,
    80: 1.32,
    81: 1.45,
    82: 1.46,
    83: 1.48,
    84: 1.4,
    85: 1.5,
    86: 1.5,
    87: 2.6,
    88: 2.21,
    89: 2.15,
    90: 2.06,
    91: 2.0,
    92: 1.96,
    93: 1.9,
    94: 1.87,
    95: 1.8,
    96: 1.69,
}
