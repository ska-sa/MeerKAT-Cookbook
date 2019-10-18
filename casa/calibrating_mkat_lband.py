#!/usr/bin/python3
# Bandpass, gain and flux calibration
# https://casa.nrao.edu/docs/TaskRef/fixvis-task.html
# https://casa.nrao.edu/docs/TaskRef/clearcal-task.html
# https://casa.nrao.edu/docs/TaskRef/setjy-task.html
# https://casa.nrao.edu/docs/TaskRef/bandpass-task.html
# https://casa.nrao.edu/docs/TaskRef/gaincal-task.html
# https://casa.nrao.edu/docs/TaskRef/fluxscale-task.html
# https://casa.nrao.edu/docs/TaskRef/applycal-task.html

from os import F_OK
from taskinit import *
from tasks import *

import argparse
import casac
import os

flux_calibrators = {'J0408-6545': [17.1, 0, 0, 0],
                    'J1939-6342': 'Stevens-Reynolds 2016',
                    'J1331+3030': 'Perley-Butler 2013',
                    }

def print_msg(msg):
    """
    Intended for status messages within CASA during calibration.
    Extra padding to make it stand out amongst all the other CASA output
    """
    border = '#' * ((len(msg) + 8))
    msg_text = '\n{}'.format(border)
    msg_text += '### {} ###'.format(msg)
    msg_text += '{}\n'.format(border)
    print(msg_text)


def cli():
    usage = "%%prog [options]"
    description = 'antenna based calibration corrections'

    parser = argparse.ArgumentParser(
            usage=usage,
            description=description,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    group = parser.add_argument_group(
            title="required arguments",
            description="arguments required by the script to run")
    group.add_argument(
            '--msfile',
            type=str,
            required=True,
            help='filename of measurement set',
            )
    group.add_argument(
            '-r', '--ref-ant',
            type=str,
            required=True,
            help='reference antenna full name, m0xx',
            )
    group.add_argument(
            '-t', '--target',
            type=str,
            required=True,
            help='comma separated list of observation target(s)',
            )
    group.add_argument(
            '-g', '--gain',
            type=str,
            help='comma separated list of gain calibrator(s)',
            )
    flux_cals = ['J0408-6545', 'J1939-6342', 'J1331+3030']
    group.add_argument(
            '-f', '--flux',
            type=str,
            required=True,
            help='comma separated list of flux calibrator(s) from list: {}'.format(flux_cals),
            )
    group = parser.add_argument_group(
            title="additional optional arguments",
            description="arguments only specified when necessary, logical defaults will be supplied")
    group.add_argument(
            '--prefix',
            type=str,
            help='the prefix for all output files',
            )
    group.add_argument(
            '--fixvis',
            action='store_true',
            help='recalculates (u, v, w) and or phase centre',
            )
    group.add_argument(
            '--setjy',
            type=str,
            choices=flux_cals,
            metavar='flux_cal',
            help='specific flux calibrator when more than one is used in the observation, picked from list: {}'.format(flux_cals),
            )
    group.add_argument(
            '-b', '--bandpass',
            type=str,
            help='comma separated list of bandpass calibrator(s)',
            )
    group.add_argument(
            '-c', '--ref-chans',
            type=str,
            default='',
            help='channels to use for preliminary gain calibration, \'spw:start_chan,end_chan\'',
            )
    group.add_argument(
            '--applycal',
            action='store_true',
            help='applycal on interim solutions',
            )
    parser.add_argument(
            '-v', '--verbose',
            action='store_true',
            help='display calibration results as output during calibration',
            )
    args = parser.parse_args()

    # set defaults for optional parameters
    if args.prefix is None:
        args.prefix = os.path.splitext(os.path.basename(args.msfile))[0]
    if args.bandpass is None:
        args.bandpass = args.flux
    if args.setjy is None:
        args.setjy = args.flux

    return args


def primary_calibrators(msfile,
                        ref_ant,
                        f_cal,
                        b_cal,
                        prelim_gcal=False,
                        ref_chans='',
                        verbose=False,
                        interimcal=False,
                        ):

    gaintable_list=['']
    if prelim_gcal:
        print_msg('Doing preliminary gain calibration')
        gtable0 = msfile + '.G0'
        if os.access(gtable0, F_OK):
            print('Deleting preliminary gain table before gaincal\n')
            rmtables(gtable0)
        gaincal(vis=msfile,
                caltable=gtable0,
                field=f_cal,
                gaintype='G',
                solint='inf',
                refant=ref_ant,
                spw=ref_chans,
                calmode='p',
                minsnr=3.0,
                solnorm=True,
                gaintable=[''])
        gaintable_list=[gtable0]

    print_msg('Doing delay calibration')
    ktable = msfile + '.K'
    if os.access(ktable, F_OK):
        print('Deleting delay table before gaincal\n')
        rmtables(ktable)
    gaincal(vis=msfile,
            caltable=ktable,
            field=f_cal,
            gaintype='K',
            solint='inf',
            refant=ref_ant,
            combine='scan',
            solnorm=False,
            minsnr=3.0,
            gaintable=gaintable_list)
    gaintable_list.append(ktable)

    print_msg('Doing bandpass calibration')
    btable = msfile + '.B'
    if os.access(btable, F_OK):
        print('Deleting bandpass table before gaincal\n')
        rmtables(btable)
    bandpass(vis=msfile,
             caltable=btable,
             field=b_cal,
             bandtype='B',
             solint='inf',
             refant=ref_ant,
             combine='scan',
             solnorm=True,
             minsnr=3.0,
             gaintable=gaintable_list)

    if interimcal:
        applycal(vis=msfile,
                 gaintable=[btable, ktable],
                 calwt=False,
                 applymode='calflag')

    print_msg('Gain calibration on fluxscale calibrator')
    gtable = msfile + '.G'
    if os.access(gtable, F_OK):
        print('Deleting gain table before gaincal\n')
        rmtables(gtable)
    gaincal(vis=msfile,
            caltable=gtable,
            field=f_cal,
            gaintype='G',
            calmode='ap',
            solint='int',
            refant=ref_ant,
            combine='scan',
            solnorm=False,
            minsnr=1.0,
            gaintable=[btable, ktable])

    if verbose:
        # Plotcal does not release table lock, so can't append to caltable.
        # Make plots as the last thing.
        plotcal(caltable=gtable0,
                xaxis='time',
                yaxis='phase',
                showgui=False,
                plotrange=[-1, -1, -180, 180],
                figfile=gtable0 + '_phase.png'
                )
        plotcal(caltable=ktable,
                xaxis='antenna',
                yaxis='delay')
        plotcal(caltable=btable,
                xaxis='chan',
                yaxis='phase',
                showgui=False,
                plotrange=[-1, -1, -180, 180],
                figfile=btable + '_phase.png'
                )
        plotcal(caltable=btable,
                xaxis='chan',
                yaxis='amp',
                showgui=False,
                figfile=btable + '_amp.png'
                )
        plotcal(caltable=gtable,
                field=f_cal,
                xaxis='time',
                yaxis='phase',
                showgui=False,
                plotrange=[-1, -1, -180, 180],
                figfile=f_cal + '_phase.png'
                )
        plotcal(caltable=gtable,
                field=f_cal,
                xaxis='time',
                yaxis='amp',
                showgui=False,
                plotrange=[-1, -1, -180, 180],
                figfile=f_cal + '_amp.png'
                )

    return [ktable, btable, gtable]

def secondary_calibrators(msfile,
                          ref_ant,
                          f_cal,
                          b_cal,
                          g_cal,
                          ktable,
                          btable,
                          gtable,
                          ref_chans='',
                          avg_time='int',
                          verbose=False,
                          ):

    cals = b_cal.split(',') + g_cal.split(',')
    cals = list(set(cals))
    for cal in f_cal.split(','):
        cals.remove(cal)
    print_msg('Gain calibration for secondary cal')
    for cal in cals:
        print('calibrating field {}'.format(cal))
        gaincal(vis=msfile,
                caltable=gtable,
                field=cal,
                gaintype='G',
                calmode='ap',
                solint=avg_time,
                refant=ref_ant,
                combine='scan',
                solnorm=False,
                minsnr=1.0,
                append=True,
                gaintable=[btable, ktable])

    print_msg('Fluxscale calibration for secondary cal')
    ftable = msfile + '.flux'
    if os.access(ftable, F_OK):
        print('Deleting flux table before gaincal\n')
        rmtables(ftable)
    fluxscale(vis=msfile,
              caltable=gtable,
              fluxtable=ftable,
              reference=f_cal,
              transfer=','.join(cals))

    if verbose:
        # Plotcal does not release table lock, so can't append to caltable.
        # Make plots as the last thing.
        for cal in cals:
            plotcal(caltable=gtable,
                    field=cal,
                    xaxis='time',
                    yaxis='phase',
                    showgui=False,
                    plotrange=[-1, -1, -180, 180],
                    figfile=cal + '_phase.png'
                    )
            plotcal(caltable=gtable,
                    field=cal,
                    xaxis='time',
                    yaxis='amp',
                    showgui=False,
                    plotrange=[-1, -1, -180, 180],
                    figfile=cal + '_amp.png'
                    )

    return ftable

def calibrate(msfile,
              ref_ant,
              f_cal,
              b_cal,
              g_cal,
              setjy_cal=None,
              standard=None,
              prelim_gcal=False,
              ref_chans='',
              avg_time='int',
              verbose=False,
              interimcal=False,
              ):

    clearcal(msfile)

    # To convert correlation coefficients to absolute flux densities
    print_msg('Setting fluxes on flux calibrator')
    if setjy_cal is None:
        setjy_cal=f_cal.split(',')[0]
    if standard is None:
        standard=flux_calibrators[setjy_cal]
    if type(standard) is list:
        setjy(vis=msfile,
              field=setjy_cal,
              scalebychan=True,
              standard='manual',
              fluxdensity=standard)
    else:
        setjy(vis=msfile,
              field=setjy_cal,
              scalebychan=True,
              standard=standard,
              fluxdensity=-1)

    [ktable, btable, gtable] = primary_calibrators(msfile,
                                                   ref_ant,
                                                   f_cal,
                                                   b_cal,
                                                   prelim_gcal=prelim_gcal,
                                                   ref_chans=ref_chans,
                                                   verbose=verbose,
                                                   interimcal=interimcal,
                                                   )

    if interimcal:
        clearstat()
        applycal(vis=msfile,
                 gaintable=[gtable, btable, ktable],
                 gainfield=[f_cal, b_cal, f_cal],
                 interp=['', 'nearest', ''],
                 calwt=False,
                 applymode='calflag')

    if g_cal is not None:
        ftable = secondary_calibrators(msfile,
                                       ref_ant,
                                       f_cal,
                                       b_cal,
                                       g_cal,
                                       ktable,
                                       btable,
                                       gtable,
                                       ref_chans='',
                                       avg_time='int',
                                       verbose=False,
                                       )
    clearstat()
    applycal(vis=msfile,
             gaintable=[ftable, btable, ktable],
             interp=['', 'nearest', ''],
             applymode='calflag')


if __name__ == '__main__':
    args = cli()

    msfile = args.msfile
    prefix = args.prefix
    target = args.target
    f_cal = args.flux
    b_cal = args.bandpass
    g_cal = args.gain
    ref_ant = args.ref_ant
    ref_chans = args.ref_chans
    setjy_cal = args.setjy

    if args.verbose:
        print('Processing with CASA')
        print('  msfile={}'.format(msfile))
        print('  prefix={}'.format(prefix))
        print('  ref_ant={}'.format(ref_ant))
        print('  target={}'.format(target))
        print('  g_cal={}'.format(g_cal))
        print('  b_cal={}'.format(b_cal))
        print('  f_cal={}'.format(f_cal))

    if args.fixvis:
        fixvis(vis=msfile, outputvis=prefix+'_fixvis.ms')
        msfile = prefix+'_fixvis.ms'

    calibrate(msfile,
              ref_ant,
              f_cal,
              b_cal,
              g_cal,
              setjy_cal=setjy_cal,
              prelim_gcal=True,
              ref_chans=ref_chans,
#               avg_time='50sec',
              verbose=args.verbose,
              interimcal=args.applycal,
              )

# -fin-
