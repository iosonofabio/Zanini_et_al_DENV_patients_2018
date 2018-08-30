# vim: fdm=indent
'''
author:     Fabio Zanini
date:       21/06/18
content:    Support module with criteria for cell subtypes
'''
cell_subtypes = {
    # NOTE: CD3+ CD56+ T/NK cells are probably gd-Tcells rather than NKT,
    # but I have to double check
    'T cell': {
        'helper': '(CD4 >= 100) & (CD8A < 100) & (CD8B < 100)',
        'killer': '(CD4 < 100) & ((CD8A >= 100) | (CD8B >= 100))',
        'cytolytic': '(PRF1 >= 100) | (GZMB >= 100) | (GZMA >= 100)',
        'all': '',
        },
    'B cell': {
        'naive': '(TCL1A >= 100) & ((IGHM >= 100) | (IGHD >= 100))',
        'isoswitched': 'isotype not in ("M", "D")',
        'isonaive': 'isotype in ("M", "D")',
        'plasma': '(MS4A1 < 100) & (PRDM1 >= 100)',
        'all': '',
        },
    'NK cell': {
        # The first 4 are relatively well described subsets
        'CD56+': '(NCAM1 >= 100)',
        'CD16+': '(FCGR3A >= 100)',
        'CD62L+': '(SELL >= 100)',
        'CD57+': '(B3GAT1 >= 100)',
        'KIR2DL3+': '(KIR2DL3 >= 100)',
        'KLRB1+': '(KLRB1 >= 100)',
        'all': '',
        },
    'monocyte': {
        'classical': '(CD14 >= 100) & (FCGR3A < 100)',
        'nonclassical': '(CD14 < 100) & (FCGR3A >= 100)',
        'double_positive': '(CD14 >= 100) & (FCGR3A >= 100)',
        'all': '',
        },
    'pDC': {
        'all': '',
        },
    'all': {
        'all': '',
        }
    }
