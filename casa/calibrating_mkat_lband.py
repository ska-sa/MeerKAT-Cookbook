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

## -- global parameters for script --
flux_calibrators = {'J0408-6545': [17.1, 0, 0, 0],
                    'J1939-6342': 'Stevens-Reynolds 2016',
                    'J1331+3030': 'Perley-Butler 2013',
                    }
## -- global parameters for script --


def _get_prefix(msfile):
    return os.path.splitext(os.path.basename(msfile))[0]


def _str2list(string_):
    list_ = [str_.strip() for str_ in string_.split(',')]
    return list(set(list_))


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
    usage = "%%prog [options] --msfile <filename.ms> -f <J...>"
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
    flux_cals = flux_calibrators.keys()
    group.add_argument(
            '-f', '--fluxcal',
            type=str,
            required=True,
            help='flux calibrator: '
                 'suggested flux cals for MeerKAT {}'.format(flux_cals),
            )
    group = parser.add_argument_group(
            title="additional optional arguments",
            description="arguments only specified when necessary, logical defaults will be supplied")
    group.add_argument(
            '--fixvis',
            action='store_true',
            help='recalculates (u, v, w) and or phase centre',
            )
    group.add_argument(
            '--standard',
            type=str,
            help='standard flux model to apply, '
                 'leave blank to apply MeerKAT standard if available',
            )
    group.add_argument(
            '-b', '--bpcal',
            type=str,
            help='comma separated list of bandpass calibrator(s)',
            )
    group.add_argument(
            '-d', '--delaycal',
            type=str,
            help='comma separated list of delay calibrator(s)',
            )
    group.add_argument(
            '-g', '--gaincal',
            type=str,
            help='comma separated list of gain calibrator(s)',
            )
    group.add_argument(
            '-t', '--target',
            type=str,
            help='comma separated list of observation target(s)',
            )
    group.add_argument(
            '-r', '--ref-ant',
            type=str,
            default='',
            help='reference antenna full name, m0xx',
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
            help='apply calibration solution to measurement set',
            )
    args = parser.parse_args()

    # set defaults for optional parameters
    if args.standard is None:
        if args.fluxcal in flux_cals:
            args.standard = flux_calibrators[args.fluxcal]
        else:
            msg = ("Not a standard MeerKAT calibrator, "
                   "please specify flux model: "
                   "--standard <flux_model> or --standard 'I, Q, U, V'")
            raise RuntimeError(msg)
    else:
        _eval_stokes = args.standard.split(',')
        if len(_eval_stokes) == 4:
            args.standard = [float(stokes) for stokes in _eval_stokes]

    if args.bpcal is None:
        args.bpcal = args.fluxcal
    if args.delaycal is None:
        args.delaycal = args.bpcal
    if args.debug:
        global DEBUG
        DEBUG = True

    return args


def primary_calibrators(msfile,
                        f_cal,
                        bp_cal,
                        delay_cal,
                        ref_ant='',
                        prelim_gcal=False,
                        ref_chans='',
                        ):

    prefix = _get_prefix(msfile)

    print_msg('Doing delay calibration using flux calibrator')
    ktable = prefix + '.K'
    if os.access(ktable, F_OK):
        print('Deleting delay table before gaincal\n')
        rmtables(ktable)
    gaincal(vis=msfile,
            caltable=ktable,
            field=delay_cal,
            gaintype='K',
            solint='inf',
            refant=ref_ant,
            combine='scan,field',
            solnorm=False,
            minsnr=3.0,
            gaintable=[])

    gaintable_list=[ktable]
    if prelim_gcal:
        print_msg('Doing preliminary gain calibration')
        gtable0 = prefix + '.G0'
        if os.access(gtable0, F_OK):
            print('Deleting preliminary gain table before gaincal\n')
            rmtables(gtable0)
        gaincal(vis=msfile,
                caltable=gtable0,
                field=bp_cal,
                gaintype='G',
                solint='inf',
                refant=ref_ant,
                spw=ref_chans,
                calmode='p',
                minsnr=3.0,
                solnorm=True,
                gaintable=[ktable])
        gaintable_list=[gtable0,ktable]

    print_msg('Doing bandpass calibration')
    btable = prefix + '.B'
    if os.access(btable, F_OK):
        print('Deleting bandpass table before gaincal\n')
        rmtables(btable)
    bandpass(vis=msfile,
             caltable=btable,
             field=f_cal,
             bandtype='B',
             solint='inf',
             refant=ref_ant,
             combine='scan',
             solnorm=True,
             minsnr=3.0,
             gaintable=gaintable_list)

    print_msg('Gain calibration for primary calibrators')
    gtable = prefix + '.G'
    if os.access(gtable, F_OK):
        print('Deleting gain table before gaincal\n')
        rmtables(gtable)
    gaincal(vis=msfile,
            caltable=gtable,
            field=bp_cal,
            gaintype='G',
            calmode='ap',
            solint='inf',
            refant=ref_ant,
            combine='spw',
            solnorm=False,
            minsnr=1.0,
            gaintable=[btable, ktable])

    return [ktable, btable, gtable]


def secondary_calibrators(msfile,
                          ktable,
                          btable,
                          gtable,
                          g_cal,
                          ref_ant='',
                          ):

    print_msg('Appending gain calibration for secondary calibrators')
    gaincal(vis=msfile,
            caltable=gtable,
            field=g_cal,
            gaintype='G',
            calmode='ap',
            solint='inf',
            refant=ref_ant,
            combine='spw',
            solnorm=False,
            minsnr=1.0,
            append=True,
            gaintable=[btable, ktable])

    return gtable


def flux_calibration(msfile,
                     gtable,
                     flux_cals,
                     ):
    prefix = _get_prefix(msfile)
    print_msg('Fluxscale calibration for secondary cal')
    ftable = prefix + '.flux'
    if os.access(ftable, F_OK):
        print('Deleting flux table before gaincal\n')
        rmtables(ftable)
    fluxscale(vis=msfile,
              caltable=gtable,
              fluxtable=ftable,
              reference=f_cal,
              transfer=flux_cals)

    return ftable


def apply_calibration(msfile,
                      cal_tables,
                      delay_cal,
                      f_cal,
                      g_cals=None,
                      targets=None,
                      ):

    print_msg('Apply calibration results')
    clearstat()

    # apply calibration to flux calibrators
    cmd = _apply_cmd(msfile, f_cal, cal_tables, )
    applycal(vis=msfile,
             field=f_cal,
             gaintable=cal_tables,
             gainfield=[f_cal, f_cal, delay_cal],
             interp=['', 'nearest', ''],
             calwt=False,
             applymode='calflag')

    # apply calibration to other calibrators
    if g_cals is not None:
        for cal in g_cals:
            applycal(vis=msfile,
                     field=cal,
                     gaintable=cal_tables,
                     gainfield=[cal, f_cal, delay_cal],
                     interp=['', 'nearest', ''],
                     calwt=False,
                     applymode='calflag')

    # target calibration
    if target is not None:
        for tgt in targets:
            applycal(vis=msfile,
                     field=tgt,
                     gaintable=cal_tables,
                     gainfield=[g_cal, f_cal, delay_cal],
                     interp=['', 'nearest', ''],
                     calwt=False,
                     applymode='calflag')


def calibrate(msfile,
              f_cal,
              bp_cal,
              delay_cal,
              g_cal=None,
              ref_ant='',
              standard=None,
              prelim_gcal=False,
              ref_chans='',
              applycal=False,
              ):

    print_msg('Clear existing calibration results')
    clearcal(msfile)

    # To convert correlation coefficients to absolute flux densities
    print_msg('Apply flux model to flux calibrator')
    if type(standard) is list:
        setjy(vis=msfile,
              field=f_cal,
              scalebychan=True,
              standard='manual',
              fluxdensity=standard)
    else:
        setjy(vis=msfile,
              field=f_cal,
              scalebychan=True,
              standard=standard,
              fluxdensity=-1)

    [ktable, btable, gtable] = primary_calibrators(msfile,
                                                   f_cal,
                                                   bp_cal,
                                                   delay_cal,
                                                   ref_ant=ref_ant,
                                                   prelim_gcal=prelim_gcal,
                                                   ref_chans=ref_chans,
                                                   )

    if g_cal is not None:
        g_cals = _str2list(g_cal)
        primcals = _str2list('{},{}'.format(bp_cal, f_cal))
        for cal in primcals:
            if cal in g_cals:
                g_cals.remove(cal)
        gtable = secondary_calibrators(msfile,
                                       ktable,
                                       btable,
                                       gtable,
                                       g_cal=','.join(g_cals),
                                       ref_ant=ref_ant,
                                       )
    cal_tables = [gtable, btable, ktable]

    secondcals = _str2list('{},{}'.format(bp_cal, g_cal))
    for cal in f_cal.split(','):
        if cal in secondcals:
            secondcals.remove(cal)
    if len(secondcals) > 0:
        ftable = flux_calibration(msfile,
                                  gtable,
                                  flux_cals=','.join(secondcals))
        if not os.path.isdir(ftable) and not DEBUG:
            raise RuntimeError(msg('flux'))
        cal_tables = [ftable, btable, ktable]

    if applycal:
        bp_cals = _str2list(bp_cal)
        for cal in f_cal.split(','):
            if cal in bp_cals:
                bp_cals.remove(cal)
        targets = _str2list(target)
        apply_calibration(msfile,
                          cal_tables,
                          delay_cal,
                          f_cal,
                          g_cals=secondcals,
                          targets=targets)


if __name__ == '__main__':
    args = cli()

    msfile = args.msfile
    bp_cal = args.bpcal
    delay_cal = args.delaycal
    f_cal = args.fluxcal
    g_cal = args.gaincal
    target = args.target
    ref_ant = args.ref_ant

    if args.fixvis:
        prefix = _get_prefix(msfile)
        print_msg('Correct phase center with fixvis')
        fixvis(vis=msfile, outputvis=prefix+'_fixvis.ms')
        msfile = prefix+'_fixvis.ms'

    calibrate(msfile,
              f_cal,
              bp_cal,
              delay_cal,
              g_cal=g_cal,
              ref_ant=ref_ant,
              standard=args.standard,
              prelim_gcal=True,
              ref_chans=args.ref_chans,
              applycal=args.applycal,
              )

# -fin-
