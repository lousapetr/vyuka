#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


def get_distribution(Kd, conc):
    """
    Calculates distribution of monomers and dimers.

    :param float Kd: dissociation constant
    :param conc: initial concentration expressed as amount of monomers

    :return mono, dimer: concentrations of monomers and dimers
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


def get_percentage(Kd, conc):
    """
    Calculate, what part of all monomers is free ('mono') or bound into dimer ('dimer').
    :param float Kd: dissociation constant
    :param conc: initial concentration of all monomers
    :return: part of monomers dict['mono'], part of dimers dict['dimer']
    :rtype: dict
    """
    mono_abs, dimer_abs = get_distribution(Kd, conc)
    mono_ratio = mono_abs/conc
    dimer_ratio = 2*dimer_abs/conc
    return {'mono': mono_ratio, 'dimer': dimer_ratio}



def sliders():
    x = np.logspace(-3, 3)
    kd0 = 1
    conc0 = 1

    r_at_conc = get_percentage(kd0, conc0)
    percentage = get_percentage(kd0, x)

    fig, ax = plt.subplots()
    plt.xscale('log')

    mono,  = plt.plot(x, percentage['mono'], 'r-', label='monomer')
    dimer, = plt.plot(x, percentage['dimer'], 'g-', label='dimer')
    # plt.legend(loc='best')

    left, right = 0.25, 0.85
    plt.subplots_adjust(left=left, right=right, bottom=0.3)
    axkd = plt.axes([left, 0.15, right - left, 0.03])
    axconc = plt.axes([left, 0.05, right - left, 0.03])

    kd_log = float(np.log10(kd0))
    conc_log = float(np.log10(conc0))
    skd = Slider(axkd, 'Kd', kd_log - 3, kd_log + 3, valinit=kd_log)
    sconc = Slider(axconc, 'Conc.', conc_log - 3, conc_log + 3, valinit=conc_log)

    def set_number_format(number):
        exp = 10**number
        if exp >= 10:
            return '%.1f nM' % exp
        else:
            return '%.3f nM' % exp

    skd.valtext.set_text(set_number_format(skd.val))
    sconc.valtext.set_text(set_number_format(sconc.val))

    def update(_):
        skd.valtext.set_text(set_number_format(skd.val))
        sconc.valtext.set_text(set_number_format(sconc.val))
        kd = 10**skd.val
        percentage = get_percentage(kd, x)
        mono.set_ydata(percentage['mono'])
        dimer.set_ydata(percentage['dimer'])


    skd.on_changed(update)
    sconc.on_changed(update)

    plt.show()


# noinspection PyUnresolvedReferences
def dilution_sliders(data, L_tot, fit_params):
    """
    Plot experimental data and curve based on values
    that can be interactively changed.

    :param data: experimental data
    :param L_tot: initial concentration of labeled
    :param tuple fit_params: fitted Kd, mono_signal, dimer_signal
    """
    from matplotlib.widgets import Slider

    conc_unlab, signal, errors = data
    # r_exp = conc_unlab / L_tot
    sig_max, sig_min = max(signal[0], signal[-1]), min(signal[0], signal[-1])

    kd0, mono_signal0, dimer_signal0 = fit_params
    prediction0 = get_signal(L_tot, conc_unlab, *fit_params)

    fig, ax = plt.subplots()
    ax.set_ylim([sig_min - 5, sig_max + 5])
    # r2_text = ax.text(0.1, sig_min - 3, 'R^2 = %.4f' % r_squared(mean_signal, prediction0))
    # ax.text(0.04, sig_min - 4, 'best R^2 = %.4f' % r_squared(mean_signal, prediction0))
    r2_text = ax.text(0.01, sig_min - 3, 'chi^2 = %.4f' % get_reduced_chi(signal, errors, prediction0))
    ax.text(0.01, sig_min - 4, 'best chi^2 = %.4f' % get_reduced_chi(signal, errors, prediction0))
    plt.subplots_adjust(bottom=0.3)
    plt.xscale('log')

    plt.errorbar(conc_unlab, signal, errors, fmt='ro')
    # plot curve
    conc_unlab_pred = np.logspace(-4, 4, 50)
    l, = plt.plot(conc_unlab_pred, get_signal(L_tot, conc_unlab_pred, kd0, mono_signal0, dimer_signal0))

    # Sizes of axes
    axkd = plt.axes([0.25, 0.15, 0.65, 0.03])
    axmono = plt.axes([0.25, 0.1, 0.65, 0.03])
    axdimer = plt.axes([0.25, 0.05, 0.65, 0.03])

    kd_log = float(np.log10(kd0))
    skd = Slider(axkd, 'Kd', kd_log - 2, kd_log + 2, valinit=kd_log)
    skd.valtext.set_text('%.5f' % kd0)
    smono = Slider(axmono, 'Mono', mono_signal0 - 15, mono_signal0 + 15, valinit=mono_signal0)
    sdimer = Slider(axdimer, 'Dimer', dimer_signal0 - 15, dimer_signal0 + 15, valinit=dimer_signal0)

    # noinspection PyUnresolvedReferences
    def update(_):
        kd = 10**skd.val
        skd.valtext.set_text('%.3f' % kd)
        mono_signal = smono.val
        dimer_signal = sdimer.val
        l.set_ydata(get_signal(L_tot, conc_unlab_pred, kd, mono_signal, dimer_signal))
        prediction = get_signal(L_tot, conc_unlab, kd, mono_signal, dimer_signal)
        # r2_text.set_text('R^2 = %.4f' % r_squared(mean_signal, prediction))
        r2_text.set_text('chi^2 = %.4f' % get_reduced_chi(signal, errors, prediction))
        fig.canvas.draw_idle()
    skd.on_changed(update)
    smono.on_changed(update)
    sdimer.on_changed(update)
    plt.show()


def main():
    sliders()

if __name__ == '__main__':
    main()
