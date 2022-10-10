#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Custom color palettes for matplotlib"""

import matplotlib.colors


_cdict = {'red':   [(0.0, 0.2, 0.2),
                   (0.1825395, 1.0, 1.0),
                   (0.5, 1.0, 1.0),
                   (4.0/6, 0.0, 0.0),
                   (1.0,   0.0, 0.0)],
         'green': [(0.0, 0.0, 0.0),
                   (0.1825395, 0.0, 0.0),
                   (0.373016, 0.9, 0.9),
                   (0.5, 1.0, 1.0),
                   (5.0/6, 0.0, 0.0),
                   (1.0,   0.0, 0.0)],
         'blue':  [(0.0,   0.0, 0.0),
                   (0.373016, 0.0, 0.0),
                   (0.5, 1.0, 1.0),
                   (5.0/6, 1.0, 1.0),
                   (1.0,   0.0, 0.0)]}
rwb=matplotlib.colors.LinearSegmentedColormap('rwb', _cdict, N=256)
rwb.set_under((0.2,0,0))
rwb.set_over((0,0,0)) #noir


_cdict = {'red':   [(0.0, 1.0, 1.0),
                   (1.0/3, 0.0, 0.0),
                   (1.0,   0.0, 0.0)],
         'green': [(0.0, 1.0, 1.0),
                   (2.0/3, 0.0, 0.0),
                   (1.0,   0.0, 0.0)],
         'blue':  [(0.0, 1.0, 1.0),
                   (2.0/3, 1.0, 1.0),
                   (1.0,   0.0, 0.0)]}
wb=matplotlib.colors.LinearSegmentedColormap('wb', _cdict, N=256)
wb.set_under((1,1,1)) #blanc
wb.set_over((0,0,0)) #noir
