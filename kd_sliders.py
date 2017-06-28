#!/usr/bin/env python

import numpy as np
# import matplotlib.pyplot as plt


def get_distribution(Kd, conc):
    """
    Calculates distribution of monomers and dimers.

    :param float Kd: dissociation constant
    :param float conc: initial concentration expressed as amount of monomers

    :return mono, dimer: concentrations of monomers and dimers
    :rtype tuple[float]
    """
    mono = 0.25 * (np.sqrt(Kd**2 + 8 * Kd * conc) - Kd)
    dimer = mono ** 2 / Kd
    return mono, dimer


def test_get_distribution():
    """
    Testing method for get_distribution().
    """
    for Kd, conc in [(0.001, 0.001), (0.001, 1000), (1, 1), (1000, 0.001), (1000, 1), (1000, 1000)]:
        mono, dimer = get_distribution(Kd, conc)
        print mono**2 / dimer - Kd
        assert abs(mono**2 / dimer - Kd) < 10**-10
    print 'OK - all tests passed.'

