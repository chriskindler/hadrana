# This dictionary collects all configurations per replica and is based on the collection of Wolfgang Soeldner
# /glurch/scratch/analysis/clsmeasurement2pt
# extended by Marcel Rodekamp and Christian Kindler

_CONFIGURATION_NUMBERS = {
    "A650": {
        "r000": {
            "end_cfg": 630,
            "broken_cfg": []
        },
        "r001": {
            "end_cfg": 4432,
            "broken_cfg": []
        }
    },
    "A653": {
        "r000": {
            "start_cfg": 1,
            "step_cfg": 2, # only odd configs
            "end_cfg": 5050,
            "broken_cfg": []
        }
    },
    "A654": {
        "r000": {
            "start_cfg": 1,
            "step_cfg": 2, # only odd configs
            "end_cfg": 5068,
            "broken_cfg": []
        }
    },
    "B450": {
        "r000": {
            "end_cfg": 1612,
            "broken_cfg": []
        }
    },
    "B452": {
        "r000": {
            "end_cfg": 1944,
            "broken_cfg": []
        }
    },
    "C101": {
        "r014": {
            "end_cfg": 2000,
            "broken_cfg": []
        }
    },
    "C102": {
        "r003": {
            "end_cfg": 500,
            "broken_cfg": []
        },
        "r004": {
            "end_cfg": 500,
            "broken_cfg": []
        },
        "r005": {
            "end_cfg": 500,
            "broken_cfg": []
        }
    },
    "D101": {
        "r005": {
            "end_cfg": 286,
            "broken_cfg": []
        },
        "r006": {
            "end_cfg": 323,
            "broken_cfg": []
        }
    },
    "D150": {
        "r000": {
            "end_cfg": 404,
            "broken_cfg": [147]
        },
        "r001": {
            "end_cfg": 199,
            "broken_cfg": [10, 47, 170, 174]
        }
    },
    "D200": {
        "r000": {
            "end_cfg": 2001,
            "broken_cfg": []
        }
    },
    "D201": {
        "r001": {
            "end_cfg": 1078,
            "broken_cfg": [257, 302, 310, 318, 370, 50, 89]
        }
    },
    "D250": {
        "r000": {
            "end_cfg": 1000,
            "broken_cfg": []
        }
    },
    "D251": {
        "r000": {
            "end_cfg": 403,
            "broken_cfg": []
        },
        "r001": {
            "end_cfg": 1610,
            "broken_cfg": []
        }
    },
    "D450": {
        "r010": {
            "end_cfg": 500,
            "broken_cfg": []
        },
        "r011": {
            "end_cfg": 1000,
            "broken_cfg": []
        }
    },
    "D451": {
        "r000": {
            "end_cfg": 1028,
            "broken_cfg": []
        }
    },
    "D452": {
        "r002": {
            "end_cfg": 1000,
            "broken_cfg": []
        }
    },
    "E250": {
        "r001": {
            "end_cfg": 503,
            "broken_cfg": []
        }
    },
    "E300": {
        "r001": {
            "end_cfg": 1139,
            "broken_cfg": []
        },
        # "r002": { # has only run9
        #     "end_cfg": 537,
        #     "broken_cfg": []
        # },
        # "r003": { # has only run9
        #     "end_cfg": 483,
        #     "broken_cfg": []
        # }
    },
    "H101": {
        "r000": {
            "end_cfg": 1007,
            "broken_cfg": []
        },
        "r001": {
            "end_cfg": 1009,
            "broken_cfg": []
        }
    },
    "H102": {
        # We combine both replicas here in one even though these should be analysed seperatly
        "r001": {
            "end_cfg": 1029,
            "broken_cfg": []
        },
        "r002": {
            "end_cfg": 1008,
            "broken_cfg": []
        }
    },
    "H105": {
        "r001": {
            "end_cfg": 1027,
            "broken_cfg": [ 20, 40, 60, 80 ]
        },
        "r002": {
            "end_cfg": 1042,
            "broken_cfg": []
        }
    },
    "H106": {
        "r001": {
            "end_cfg": 661,
            "broken_cfg": []
        },
        "r002": {
            "end_cfg": 589,
            "broken_cfg": []
        },
        "r003": {
            "end_cfg": 256,
            "broken_cfg": []
        },
        "r004": {
            "end_cfg": 47,
            "broken_cfg": []
        }
    },
    "H107": {
        "r000": {
            "end_cfg": 150,
            "broken_cfg": []
        },
        "r001": {
            "end_cfg": 150,
            "broken_cfg": []
        },
        "r002": {
            "end_cfg": 248,
            "broken_cfg": []
        },
        "r003": {
            "end_cfg": 252,
            "broken_cfg": []
        },
        "r004": {
            "end_cfg": 388,
            "broken_cfg": []
        },
        "r005": {
            "end_cfg": 376,
            "broken_cfg": []
        }
    },
    "J250": {
        "r000": {
            "end_cfg": 903,
            "broken_cfg": []
        }
    },
    "J303": {
        "r003": {
            "end_cfg": 998,
            "broken_cfg": []
        }
    },
    "J304": {
        "r000": {
            "end_cfg": 830,
            "broken_cfg": []
        },
        "r001": {
            "end_cfg": 804,
            "broken_cfg": []
        }
    },
    "J305": {
        "r000": {
            "end_cfg": 817,
            "broken_cfg": []
        },
        "r001": {
            "end_cfg": 404,
            "broken_cfg": []
        }
    },
    "J500": {
        "r004": {
            "end_cfg": 789,
            "broken_cfg": []
        },
        "r005": {
            "end_cfg": 655,
            "broken_cfg": []
        },
        "r006": {
            "end_cfg": 431,
            "broken_cfg": []
        }
    },
    "J501": {
        "r001": {
            "end_cfg": 1635,
            "broken_cfg": []
        },
        "r002": {
            "end_cfg": 1142,
            "broken_cfg": []
        },
        "r003": {
            "end_cfg": 1150,
            "broken_cfg": []
        }
    },
    "N101": {
        "r001": {
            "end_cfg": 280,
            "broken_cfg": list(range(211,281))
        },
        "r003": {
            "end_cfg": 404,
            "broken_cfg": list(range(1,46))
        },
        "r004": {
            "end_cfg": 240,
            "broken_cfg": [ 24 ]
        },
        "r005": {
            "end_cfg": 352,
            "broken_cfg": []
        },
        "r006": {
            "end_cfg": 320,
            "broken_cfg": []
        }
    },
    "N200": {
        "r000": {
            "end_cfg": 856,
            "broken_cfg": []
        },
        "r001": {
            "end_cfg": 856,
            "broken_cfg": []
        }
    },
    "N201": {
        "r001": {
            "end_cfg": 1522,
            "broken_cfg": []
        }
    },
    "N202": {
        "r001": {
            "end_cfg": 899,
            "broken_cfg": []
        },
        "r002": {
            "end_cfg": 1003,
            "broken_cfg": []
        }
    },
    "N203": {
        "r000": {
            "end_cfg": 756,
            "broken_cfg": []
        },
        "r001": {
            "end_cfg": 787,
            "broken_cfg": []
        }
    },
    "N204": {
        "r000": {
            "end_cfg": 1500,
            "broken_cfg": []
        }
    },
    "N300": {
        "r002": {
            "end_cfg": 1540,
            "broken_cfg": []
        }
    },
    "N302": {
        "r001": {
            "end_cfg": 2201,
            "broken_cfg": []
        }
    },
    "N304": {
        "r000": {
            "end_cfg": 974,
            "broken_cfg": []
        },
        "r001": {
            "end_cfg": 752,
            "broken_cfg": []
        }
    },
    "N306": {
        "r000": {
            "end_cfg": 1507,
            "broken_cfg": []
        },
    },
    "N401": {
        "r000": {
            "end_cfg": 1100,
            "broken_cfg": [381, 390, 464, 468, 580, 629, 690, 946, 947, 948, 1007, 1014, 1044, 1089]
        }
    },
    "N450": {
        "r000": {
            "end_cfg": 808,
            "broken_cfg": []
        },
        "r001": {
            "end_cfg": 387,
            "broken_cfg": []
        }
    },
    "N451": {
        "r000": {
            "end_cfg": 1011,
            "broken_cfg": []
        }
    },
    "S100": {
        "r003": {
            "end_cfg": 499,
            "broken_cfg": []
        },
        "r004": {
            "end_cfg": 485,
            "broken_cfg": []
        }
    },
    "S201": {
        "r002": {
            "end_cfg": 1112,
            "broken_cfg": [ 920, 921, 922, 923, 924, 925, 926, 927, 928, 929, 930, 931, 932, 933, 934, 935, 936, 937, 938, 939, 940, 941, 942, 943, 944, 945, 946, 947, 948, 949, 950, 951, 952, 953, 954, 955, 956, 957, 958, 959, 960, 961, 962, 963, 964, 965, 966, 967, 968, 969, 970, 971, 972, 973, 974 ]
        },
        "r003": {
            "end_cfg": 1036,
            "broken_cfg": []
        }
    },
    "S400": {
        "r000": {
            "end_cfg": 872,
            "broken_cfg": []
        },
        "r001": {
            "end_cfg": 2001,
            "broken_cfg": []
        }
    },
    "U102": {
        "r001": {
            "end_cfg": 1299,
            "broken_cfg": []
        },
        "r002": {
            "end_cfg": 3562,
            "broken_cfg": []
        }
    },
    "U103": {
        "r001": {
            "start_cfg": 1,
            "step_cfg": 2, # only odd configs
            "end_cfg": 2121,
            "broken_cfg": []
        },
        "r002": {
            "start_cfg": 1,
            "step_cfg": 2, # only odd configs
            "end_cfg": 1781,
            "broken_cfg": []
        },
        "r003": {
            "start_cfg": 1,
            "step_cfg": 2, # only odd configs
            "end_cfg": 1819,
            "broken_cfg": [1781, 1783, 1793, 1795, 1799, 1801]
        }
    },
    "X151": {
        "r000": {
            "end_cfg": 1000,
            "broken_cfg": []
        }
    },
    "X250": {
        "r000": {
            "end_cfg": 345,
            "broken_cfg": []
        },
        "r001": {
            "end_cfg": 1155,
            "broken_cfg": []
        }
    },
    "X251": {
        "r000": {
            "end_cfg": 436,
            "broken_cfg": []
        },
        "r001": {
            "end_cfg": 1064,
            "broken_cfg": []
        }
    },
    "X252": {
        "r000": {
            "end_cfg": 3508,
            "broken_cfg": []
        }
    },
    "X253": {
        "r000": {
            "end_cfg": 3005,
            "broken_cfg": []
        }
    },
    "X450": {
        "r001": {
            "end_cfg": 401,
            "broken_cfg": [ 401 ]
        },
        "r002": {
            "end_cfg": 659,
            "broken_cfg": []
        }
    },
    "rqcd021": {
        "r021": {
            "start_cfg": 500,
            "end_cfg": 2057,
            "broken_cfg": [ 1037, 1038, 1043, 1044, 1045, 1047, 1052, 1053, 1056, 1066, 1067, 1069, 1071, 1219, 1220, 1222, 1269 ]
        }
    },
    "rqcd030": {
        "r030": {
            "start_cfg": 300,
            "end_cfg": 1525,
            "broken_cfg": [ 1369, 1371 ]
        }
    },

}

