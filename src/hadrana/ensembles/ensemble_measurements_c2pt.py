# This dictionary was originally created by Simon Weishaeupl and Wolfgang Soeldner. Later extended by Marcel Rodekamp and Christian Kindler
# It contains currently available data of nucleon two-point functions generated within the CLS effort
# Updated: 24/02/2026

_MEASUREMENTS_C2PT = {
    'A650': {
        'c2pt': {
            'source_set0': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run13_nog2
            'source_set1': { 'src_ids': [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23], 'src_ids is t_src': False }, # run9
            'source_set2': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], 'src_ids is t_src': False }, # run13b_nog2

        }
    },
    'A653': {
        'c2pt': {
            'source_set0': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run13
            'source_set1': { 'src_ids': [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23], 'src_ids is t_src': False }, # run9
            'source_set2': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], 'src_ids is t_src': False }, # run13b_nog2

        }
    },
    'A654': {
        'c2pt': {
            # 'source_set0': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run11
            'source_set2': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run13_nog2
            'source_set3': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], 'src_ids is t_src': False }, # run13b_nog2
            'source_set4': { 'src_ids': [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23], 'src_ids is t_src': False }, # run9

        }
    },
    'B450': {
        'c2pt': {
            # 'source_set0': { 'src_ids': [], 'src_ids is t_src': False }, # run6_BAD_NAMING_THREEPT
            # 'source_set1': { 'src_ids': [], 'src_ids is t_src': False }, # run6
            'source_set2': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run13_nog2
            'source_set3': { 'src_ids': [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30], 'src_ids is t_src': False }, # run9
            'source_set4': { 'src_ids': [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30], 'src_ids is t_src': False }, # run9_nocharm
            'source_set5': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], 'src_ids is t_src': False }, # run13b_nog2

        }
    },
    'B452': {
        'c2pt': {
            # 'source_set0': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run11
            'source_set2': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run13_nog2
            'source_set3': { 'src_ids': [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31], 'src_ids is t_src': False }, # run9
            'source_set4': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], 'src_ids is t_src': False }, # run13b_nog2

        }
    },
    'C101': {
        'c2pt': {
            # 'source_set0': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run11
            'source_set2': { 'src_ids': [30, 36, 44, 51, 53, 55], 'src_ids is t_src': True }, # run67
            'source_set3': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run13b_nog2
            'source_set4': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7], 'src_ids is t_src': False }, # run13c_nog2
            'source_set5': { 'src_ids': [32, 34, 38, 40, 42, 46, 47, 48, 49, 57, 59, 61, 63, 65], 'src_ids is t_src': True }, # run9+7

        }
    },
    'C102': {
        'c2pt': {
            # 'source_set0': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run11
            # 'source_set1': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run11b
            'source_set2': { 'src_ids': [0, 1, 2, 3, 4], 'src_ids is t_src': False }, # run13_nog2
            'source_set3': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run13b_nog2
            'source_set4': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7], 'src_ids is t_src': False }, # run13c_nog2
            'source_set6': { 'src_ids': [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35], 'src_ids is t_src': False }, # run9

        }
    },
    'D101': {
        'c2pt': {
            'source_set2': { 'src_ids': [30, 42, 54, 64, 75, 87], 'src_ids is t_src': True }, # run6

        }
    },
    'D150': {
        'c2pt': {
            'source_set2': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run6

        }
    },
    'D200': {
        'c2pt': {
            # 'source_set0': { 'src_ids': [], 'src_ids is t_src': False }, # run11
            'source_set2': { 'src_ids': [40, 50, 60, 68, 71, 73], 'src_ids is t_src': True }, # run6
            'source_set3': { 'src_ids': [0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38], 'src_ids is t_src': False }, # run9

        }
    },
    'D201': {
        'c2pt': {
            # 'source_set0': { 'src_ids': [], 'src_ids is t_src': False }, # run11
            # 'source_set1': { 'src_ids': [], 'src_ids is t_src': False }, # run11b
            'source_set2': { 'src_ids': [40, 50, 60, 68, 71, 73], 'src_ids is t_src': True }, # run6
        }
    },
    'D250': {
        'c2pt': {
            'source_set0': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], 'src_ids is t_src': False }, # run13_nog2

        }
    },
    'D251': {
        'c2pt': {
            'source_set0': { 'src_ids': [0, 2, 4, 6, 8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62], 'src_ids is t_src': False}, # run9
            'source_set1': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7], 'src_ids is t_src': False }, # run13_nog2

        }
    },
    # D450 has older replicas r00{4,5,6,9} which have little statistic.
    # For simplicity we ignore the, here and only consider replicas 10,11
    #'D450': {
    #    'c2pt': {
    #        # 'source_set0': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run11
    #        'source_set2': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run13
    #        'source_set3': { 'src_ids': [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62], 'src_ids is t_src': False }, # run9
    #        # source_set4 (run13b_nog2) not done, but for replicas r010,r011 (see below)
    #    }
    #},
    'D450': {
        'c2pt': {
            'source_set2': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7], 'src_ids is t_src': False }, # run13
            'source_set3': { 'src_ids': [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62], 'src_ids is t_src': False }, # run9
            'source_set4': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7], 'src_ids is t_src': False }, # run13b_nog2
        }
    },
    'D451': {
        'c2pt': {
            # 'source_set0': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run11
            'source_set2': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run13
            'source_set3': { 'src_ids': [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62], 'src_ids is t_src': False }, # run9
            'source_set4': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], 'src_ids is t_src': False }, # run13b_nog2

        }
    },
    'D452': {
        'c2pt': {
            # 'source_set0': { 'src_ids': [], 'src_ids is t_src': False }, # run11
            # 'source_set1': { 'src_ids': [], 'src_ids is t_src': False }, # run11b
            'source_set2': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], 'src_ids is t_src': False }, # run13_nog2
            'source_set3': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31], 'src_ids is t_src': False }, # run13b_nog2

        }
    },
    'E250': {
        'c2pt': {
            # 'source_set0': { 'src_ids': [0, 2, 4, 6, 8, 10, 12, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62], 'src_ids is t_src': False }, # run9
            'source_set2': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run13

        }
    },
    'E300': {
        'c2pt': {
            # 'source_set0': { 'src_ids': [], 'src_ids is t_src': False }, # run11
            'source_set2': { 'src_ids': [0, 1, 2, 3, 4, 5], 'src_ids is t_src': False }, # run13_nog2

        }
    },
    'H101': {
        'c2pt': {
            'source_set0': { 'src_ids': [30, 51, 53, 55, 57], 'src_ids is t_src': True }, # run6
            # 'source_set1': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], 'src_ids is t_src': False }, # run9_nosmearedlocal
            'source_set2': { 'src_ids': [30, 32, 34, 36, 38, 40, 42, 44, 46, 47, 48, 49, 51, 53, 55, 57, 59, 61, 63, 65], 'src_ids is t_src': True }, # run9
            'source_set3': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21], 'src_ids is t_src': False }, # run13_nog2

        }
    },
    'H102': {
        'c2pt': {
            # 'source_set0': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run11
            # 'source_set1': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run11b
            'source_set2': { 'src_ids': [30, 36, 44, 51, 53, 55], 'src_ids is t_src': True }, # run7
            'source_set3': { 'src_ids': [32, 34, 38, 40, 42, 46, 47, 48, 49, 57, 59, 61, 63, 65], 'src_ids is t_src': True }, # run7
            # 'source_set4': { 'src_ids': [0, 1, 2, 3, 4], 'src_ids is t_src': False }, # run12
            'source_set5': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], 'src_ids is t_src': False }, # run13_nog2

        }
    },
    # 'H102r002': {
    #     'c2pt': {
    #         # 'source_set0': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run11
    #         # 'source_set1': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run11b
    #         'source_set2': { 'src_ids': [30, 36, 44, 51, 53, 55], 'src_ids is t_src': True }, # run7
    #         'source_set3': { 'src_ids': [32, 34, 38, 40, 42, 46, 47, 48, 49, 57, 59, 61, 63, 65], 'src_ids is t_src': True }, # run7
    #         'source_set4': { 'src_ids': [0, 1, 2, 3, 4], 'src_ids is t_src': False }, # run12
    #         'source_set5': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], 'src_ids is t_src': False }, # run13_nog2

    #     }
    # },
    'H105': {
        'c2pt': {
            # 'source_set0': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run11
            'source_set2': { 'src_ids': [30, 36, 44, 51, 53, 55], 'src_ids is t_src': True }, # run7
            'source_set3': { 'src_ids': [32, 34, 38, 40, 42, 46, 47, 48, 49, 57, 59, 61, 63, 65], 'src_ids is t_src': True }, # run7
            'source_set4': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], 'src_ids is t_src': False }, # run13_nog2

        }
    },
    'H106': {
        'c2pt': {
            # 'source_set0': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run11
            # 'source_set1': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run11b
            'source_set2': { 'src_ids': [0, 1, 2, 3, 4], 'src_ids is t_src': False }, # run13
            'source_set3': { 'src_ids': [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35], 'src_ids is t_src': False }, # run9
            'source_set4': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], 'src_ids is t_src': False }, # run13b_nog2

        }
    },
    'H107': {
        'c2pt': {
            # 'source_set0': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run11
            'source_set2': { 'src_ids': [0, 1, 2, 3, 4], 'src_ids is t_src': False }, # run13
            'source_set3': { 'src_ids': [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35], 'src_ids is t_src': False }, # run9
            'source_set4': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], 'src_ids is t_src': False }, # run13b_nog2

        }
    },
    'J250': {
        'c2pt': {
            'source_set2': { 'src_ids': [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23], 'src_ids is t_src': False }, # run13_nog2
        }
    },
    'J303': {
        'c2pt': {
            # 'source_set0': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run11
            'source_set2': { 'src_ids': [40, 64, 80, 88, 110, 127], 'src_ids is t_src': True }, # run6
            'source_set3': { 'src_ids': [50, 64, 80, 98, 100, 120, 130], 'src_ids is t_src': True }, # run6b
            'source_set4': { 'src_ids': [40, 46, 52, 58, 64, 70, 76, 82, 88, 94, 96, 100, 106, 112, 118, 124, 130, 136, 142, 148], 'src_ids is t_src': True }, # run8

        }
    },
    'J304': {
        'c2pt': {
            # 'source_set0': { 'src_ids': [], 'src_ids is t_src': False }, # run11
            'source_set2': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run13

        }
    },
    'J500': {
        'c2pt': {
            'source_set0': { 'src_ids': [50, 69, 80, 89, 109, 114, 119], 'src_ids is t_src': True }, # run6
            'source_set1': { 'src_ids': [48, 53, 58, 63, 68, 73, 78, 83, 88, 93, 98, 103, 108, 113, 118, 123, 128, 133, 138, 143], 'src_ids is t_src': True }, # run9

        }
    },
    'J501': {
        'c2pt': {
            'source_set1': { 'src_ids': [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38], 'src_ids is t_src': False }, # run9
            'source_set2': { 'src_ids': [50, 69, 80, 89, 109, 114, 119], 'src_ids is t_src': True }, # run6
        }
    },
    'N101': {
        'c2pt': {
            'source_set2': { 'src_ids': [30, 42, 54, 64, 75, 87], 'src_ids is t_src': True }, # run7
        }
    },
    'N200': {
        'c2pt': {
            # 'source_set0': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run11/twopt
            # 'source_set1': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run11b
            'source_set2': { 'src_ids': [40, 49, 61, 69, 72, 75], 'src_ids is t_src': True }, # run7
            'source_set3': { 'src_ids': [37, 43, 46, 52, 55, 58, 63, 64, 66, 78, 81, 84, 87, 90], 'src_ids is t_src': True }, # run7
            'source_set4': { 'src_ids': [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19], 'src_ids is t_src': False }, # run13_nog2
        }
    },
    'N201': {
        'c2pt': {
            # 'source_set0': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run11
            'source_set2': { 'src_ids': [0, 1, 2, 3, 4], 'src_ids is t_src': False }, # run13
            'source_set3': { 'src_ids': [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47], 'src_ids is t_src': False }, # run9

        }
    },
    'N202': {
        'c2pt': {
            'source_set0': { 'src_ids': [40, 50, 60, 68, 71, 73], 'src_ids is t_src': True }, # run6
            'source_set1': { 'src_ids': [37, 40, 43, 46, 49, 52, 55, 58, 61, 63, 64, 66, 69, 72, 75, 78, 81, 84, 87, 90], 'src_ids is t_src': True }, # run9

        }
    },
    'N203': {
        'c2pt': {
            # 'source_set0': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run11
            # 'source_set1': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run11b
            'source_set2': { 'src_ids': [40, 49, 61, 69, 72, 75], 'src_ids is t_src': True }, # run7
            'source_set3': { 'src_ids': [37, 43, 46, 52, 55, 58, 63, 64, 66, 78, 81, 84, 87, 90], 'src_ids is t_src': True }, # run7
        }
    },
    'N204': {
        'c2pt': {
            # 'source_set0': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run11
            'source_set2': { 'src_ids': [0, 1, 2, 3, 4], 'src_ids is t_src': False }, # run13_nog2
            # run12 has only low momentum data, and is thus exlcuded
            # 'source_set3': { 'src_ids': [0, 1, 2, 3, 4], 'src_ids is t_src': False }, # run12
            'source_set4': { 'src_ids': [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47], 'src_ids is t_src': False }, # run14

        }
    },
    'N300': {
        'c2pt': {
            'source_set0': { 'src_ids': [40, 45, 55, 63, 66, 70], 'src_ids is t_src': True }, # run6
            'source_set1': { 'src_ids': [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34], 'src_ids is t_src': False }, # run9
            'source_set3': { 'src_ids': [0, 1, 2, 3, 4, 5], 'src_ids is t_src': False }, # run6

        }
    },
    'N302': {
        'c2pt': {
            # 'source_set0': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run11
            'source_set2': { 'src_ids': [40, 46, 55, 63, 66, 69], 'src_ids is t_src': True }, # run7
            'source_set3': { 'src_ids': [37, 43, 49, 52, 58, 61, 64, 72, 75, 78, 81, 84, 87, 90], 'src_ids is t_src': True }, # run7

        }
    },
    'N304': {
        'c2pt': {
            # 'source_set0': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run11
            'source_set2': { 'src_ids': [0, 1, 2, 3, 4], 'src_ids is t_src': False }, # run13
            'source_set3': { 'src_ids': [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47], 'src_ids is t_src': False }, # run14

        }
    },
    'N306': {
        'c2pt': {
            'source_set1': { 'src_ids': [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38], 'src_ids is t_src': False }, # run9
            'source_set2': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29], 'src_ids is t_src': False }, # run13_nog2

        }
    },
    'N401': {
        'c2pt': {
            # 'source_set0': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run11
            'source_set2': { 'src_ids': [40, 49, 61, 69, 72, 75], 'src_ids is t_src': True }, # run7
            'source_set3': { 'src_ids': [37, 43, 46, 52, 55, 58, 63, 64, 66, 78, 81, 84, 87, 90], 'src_ids is t_src': True }, # run7

        }
    },
    'N450': {
        'c2pt': {
            # 'source_set0': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run11
            'source_set2': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run13
            'source_set3': { 'src_ids': [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63], 'src_ids is t_src': False }, # run9

        }
    },
    'N451': {
        'c2pt': {
            # 'source_set0': { 'src_ids': [], 'src_ids is t_src': False }, # run11
            # 'source_set1': { 'src_ids': [], 'src_ids is t_src': False }, # run11b
            'source_set2': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run13_nog2
            'source_set3': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27], 'src_ids is t_src': False }, # run13b_nog2

        }
    },
    'S100': {
        'c2pt': {
            'source_set2': { 'src_ids': [30, 42, 54, 64, 75, 87], 'src_ids is t_src': True }, # run7
            'source_set3': { 'src_ids': [32, 34, 36, 38, 40, 44, 46, 48, 50, 52, 56, 58, 60, 62, 63, 65, 67, 69, 71, 73, 77, 79, 81, 83, 85, 89, 91, 93, 95, 97], 'src_ids is t_src': True }, # run7
            # 'source_set4': { 'src_ids': list(range(1,136,2)), 'src_ids is t_src':True }, # run15
            # 'source_set5': { 'src_ids': list(range(1,136,2)), 'src_ids is t_src':True }, # run15b
        }
    },
    'S201': {
        'c2pt': {
            'source_set2': { 'src_ids': [40, 49, 61, 69, 72, 75], 'src_ids is t_src': True }, # run6
            # have no smeared smeared correlators
            #'source_set3': { 'src_ids': list(range(1,95,2)), 'src_ids is t_src': True }, # run15
            #'source_set4': { 'src_ids': list(range(1,95,2)), 'src_ids is t_src': True }, # run15b
            'source_set6': { 'src_ids': [0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38], 'src_ids is t_src':False }, # run9
        }
    },
    'S400': {
        'c2pt': {
            # 'source_set0': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run11
            'source_set2': { 'src_ids': [40, 49, 61, 69, 72, 75], 'src_ids is t_src': True }, # run7
            'source_set3': { 'src_ids': [37, 43, 46, 52, 55, 58, 63, 64, 66, 78, 81, 84, 87, 90], 'src_ids is t_src': True }, # run7
            # run 12 has only low momenta and thus is ignored here
            # 'source_set4': { 'src_ids': [0, 1, 2, 3, 4], 'src_ids is t_src': False }, # run12
            'source_set5': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], 'src_ids is t_src': False }, # run13_nog2

        }
    },
    'U102': {
        'c2pt': {
            # run7 has only 1597/3562 configs
            # 'source_set2': { 'src_ids': [30, 42, 54, 64, 75, 87], 'src_ids is t_src': True }, # run7
            # run15{,b,c} have only measurements on odd configs
            # 'source_set3': { 'src_ids': list(range(1,136,2)), 'src_ids is t_src': True }, # run15
            # 'source_set4': { 'src_ids': list(range(1,136,2)), 'src_ids is t_src': True }, # run15b
            # 'source_set5': { 'src_ids': list(range(1,136,2)), 'src_ids is t_src': True }, # run15c
            'source_set6': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], 'src_ids is t_src': True }, # run13_nog2
        }
    },
    'U103': {
        'c2pt': {
            'source_set0': { 'src_ids': [30, 46, 63, 67, 81, 83, 85, 87], 'src_ids is t_src': True }, # run7
            'source_set1': { 'src_ids': [32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96], 'src_ids is t_src': True }, # run10
            'source_set5': { 'src_ids': [30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97], 'src_ids is t_src': True }, # run15c
            'source_set6': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], 'src_ids is t_src': False }, # run13_nog2

        }
    },
    'X151': {
        'c2pt': {
            'source_set0': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], 'src_ids is t_src': False }, # run13_nog2
            'source_set1': { 'src_ids': [0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62], 'src_ids is t_src': False }, # run9
        }
    },
    'X250': {
        'c2pt': {
            'source_set0': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run13_nog2
            'source_set1': { 'src_ids': [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63], 'src_ids is t_src': False }, # run14
            'source_set2': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], 'src_ids is t_src': False }, # run13b_nog2

        }
    },
    'X251': {
        'c2pt': {
            'source_set0': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run13_nog2
            'source_set1': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run13b_nog2
            'source_set2': { 'src_ids': [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63], 'src_ids is t_src': False }, # run14
            'source_set3': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], 'src_ids is t_src': False }, # run13c_nog2

        }
    },
    'X252': {
        'c2pt': {
            'source_set0': { 'src_ids': [0, 2, 4, 6, 8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62], 'src_ids is t_src': False }, # run9
            'source_set1': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23], 'src_ids is t_src': False }, # run13_nog2

        }
    },
    'X253': {
        'c2pt': {
            'source_set0': { 'src_ids': [0, 2, 4, 6, 8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62], 'src_ids is t_src': False }, # run9
            'source_set1': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23], 'src_ids is t_src': False }, # run13_nog2

        }
    },
    'X450': {
        'c2pt': {
            'source_set2': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run13_nog2
            'source_set3': { 'src_ids': [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63], 'src_ids is t_src': False }, # run14
            'source_set4': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], 'src_ids is t_src': False }, # run13b_nog2

        }
    },
    'rqcd021': {
        'c2pt': {
            'source_set0': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run13_nog2
            'source_set1': { 'src_ids': [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31], 'src_ids is t_src': False }, # run14
            'source_set2': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], 'src_ids is t_src': False }, # run13b_nog2

        }
    },
    'rqcd030': {
        'c2pt': {
            'source_set0': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run13_nog2
            'source_set1': { 'src_ids': [0, 1, 2, 3], 'src_ids is t_src': False }, # run13b_nog2
            'source_set2': { 'src_ids': [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31], 'src_ids is t_src': False }, # run9
            'source_set3': { 'src_ids': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23], 'src_ids is t_src': False }, # run13c_nog2

        }
    },

}
