import os, sys, configparser, ast
from optparse import OptionParser
import namib

from namib.utils import Posteriors, Plots
from namib.utils import create_directory, save_posteriors_to_txt, save_output_medians
from namib.options import usage


def main():

    parser = OptionParser(usage)
    parser.add_option('--config-file', type='string', metavar = 'config_file', default = None)
    (opts, args) = parser.parse_args()

    config_file = opts.config_file
    if not config_file: parser.error('Please specify a config file.\n')
    if not os.path.exists(config_file): parser.error('Config file {} not found.\n'.format(config_file))

    Config = configparser.ConfigParser()
    Config.read(config_file)

    input_pars = {

        # [input]

        'samp-dir'           : '',
        'output'             : '',
        'screen-output'      : 0,

        'stack-mode'         : 'event',
        'compare'            : '',

        'parameters'         : ['m1', 'm2', 'chi1', 'chi2'],
        'bounds'             : [],
        'ordering'           : [],
        'compare-ordering'   : [],
        
        'compare-hard'       : 0,
        'evidence'           : 0,
        'include-prior'      : 0,
        'include-IMR'        : 0,
        'truths'             : [],

        'modes'              : [(2,2)],
        'ds-scaling'         : 0,
        'qnms-pyRing'        : 1,
        'IMR-fits'           : 'JimenezForteza_TEOBPM',
        'IMR-fits-IMR'       : 'NRSur7dq4Remnant',

        'save-post'          : 0,
        'save-medians'       : 0,
        'downsample'         : 1,

        # [plots]

        'corner'             : 0,
        'violin'             : 0,
        'ridgeline'          : 0,
        'TGR-plot'           : 0,
        'corner-sns'         : 1,

        'corner-settings'    : {'figsize':  8,       'figname': 'corner',    'alpha': 0.5, 'smooth': 0,    'linewidth': 1},
        'violin-settings'    : {'figsize': (15, 25), 'figname': 'violin',    'alpha': 0.5, 'rotation': 0,  'pad': -0.5},
        'ridgeline-settings' : {'figsize': (20, 10), 'figname': 'ridgeline', 'alpha': 0.5, 'overlap': 0.5, 'fade': 0, 'borderaxespad': 0.5},
        'label-sizes'        : {'xtick': 15, 'ytick': 15, 'legend': 17, 'axes': 17},
        'palette'            : 'crest',

        'extra-row'          : '',
        'BF-comparison'      : 0,
        'evidence-top'       : 0,
        'time-percentiles'   : [],
        'automatic-bounds'   : 0,
        'IMR-posteriors'     : 0,

        'horizontal-legend'  : 0,
        'event-name'         : '',
        'remove-xticks'      : 0,
        'remove-legend'      : 0,
        'fix-dimensions'     : 0,
        
        'single-prior'       : '',
        'prior-color'        : '#828F61',
        'truth-color'        : 'k',

    }
    
    for key in input_pars.keys():

        if ('samp-dir' in key) or ('output' in key) or ('stack-mode' in key) or ('compare' in key):
            try: input_pars[key] = Config.get('input', key)
            except: pass
        if ('screen-output' in key) or ('compare-hard' in key) or ('evidence' in key) or ('save-post' in key) or ('include-prior' in key) or ('include-IMR' in key) or ('ds-scaling' in key) or ('screen-medians' in key) or ('save-medians' in key) or ('qnms-pyRing' in key):
            try: input_pars[key] = Config.getboolean('input', key)
            except: pass
        if ('downsample' in key):
            try: input_pars[key] = Config.getfloat('input', key)
            except: pass
        if ('parameters' in key) or ('modes' in key) or ('ordering' in key) or ('bounds' in key) or ('compare-ordering' in key ) or ('truths' in key) or ('IMR-fits' in key) or ('IMR-fits-IMR' in key):
            try: input_pars[key] = ast.literal_eval(Config.get('input', key))
            except: pass
        if ('corner' in key) or ('violin' in key) or ('ridgeline' in key) or ('TGR-plot' in key) or ('BF-comparison' in key) or ('evidence-top' in key) or ('remove-xticks' in key) or ('remove-legend' in key) or ('horizontal-legend' in key) or ('fix-dimensions' in key) or ('corner-sns' in key) or ('automatic-bounds' in key) or ('IMR-posteriors' in key):
            try: input_pars[key] = Config.getboolean('plots', key)
            except: pass
        if ('extra-row' in key) or ('single-prior' in key) or ('prior-color' in key) or ('event-name' in key) or ('truth-color' in key):
            try: input_pars[key] = Config.get('plots', key)
            except: pass
        if ('palette' in key) or ('time-percentiles' in key) or ('corner-settings' in key) or ('violin-settings' in key) or ('ridgeline-settings' in key) or ('label-sizes' in key):
            try: input_pars[key] = ast.literal_eval(Config.get('plots', key))
            except: pass

    # Deviate stdout and stderr to file
    if input_pars['screen-output']:
        sys.stdout = open(os.path.join(input_pars['output'], 'stdout_namib.txt'), 'w')
        sys.stderr = open(os.path.join(input_pars['output'], 'stderr_namib.txt'), 'w')
    else: pass
    print('\nn a m i b\n')
    print(('Reading config file:\n{}'.format(config_file)))

    if not os.path.exists(input_pars['samp-dir']):
        raise ValueError('\nSamples directory {} not found.\n'.format(input_pars['samp-dir']))
    else: print('\nPosteriors are read from:\n{}'.format(input_pars['samp-dir']))

    # Read the posteriors and create the .txt files with the reduced posteriors
    PostOutput = Posteriors(input_pars)
    SampDataFrame, PriorDataFrame, EvidenceDataFrame, IMRDataFrame = PostOutput.return_samples_dict()

    if input_pars['save-post']:
        red_post_dir = create_directory(input_pars['output'], 'reduced_posteriors')
        save_posteriors_to_txt(input_pars, red_post_dir, SampDataFrame)
    if input_pars['evidence'] and input_pars['save-medians']:
        output_medians_dir = create_directory(input_pars['output'], 'output_medians')
        save_output_medians(input_pars, SampDataFrame, EvidenceDataFrame, output_medians_dir)

    if  input_pars['corner'] or input_pars['violin'] or input_pars['ridgeline'] or input_pars['TGR-plot']:
        input_pars['plots-dir'] = create_directory(input_pars['output'], '')
        Plots(input_pars, SampDataFrame, PriorDataFrame, EvidenceDataFrame, IMRDataFrame)
        print('\nPlots are saved in:\n{}'.format(input_pars['plots-dir']))
    
    print('\nFinished.\n')


if __name__=='__main__':
    main()
