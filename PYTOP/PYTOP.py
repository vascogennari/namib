import os, configparser, ast
from optparse import OptionParser

from utils import Posteriors, Plots, save_posteriors_to_txt
from utils import create_directory


if __name__=='__main__':

    parser = OptionParser()
    parser.add_option('--config-file', type='string', metavar = 'config_file', default = None)
    (opts, args) = parser.parse_args()

    cwd = os.getcwd()
    pardir_path = os.path.abspath(os.path.join(cwd, os.pardir))
    config_path = os.path.join(os.path.join(pardir_path, 'config_files'), opts.config_file)

    Config = configparser.ConfigParser()
    Config.read(config_path)

    input_pars = {
        'samp-dir'           : 'samples',
        'output'             : 'posteriors',
        'file-path'          : '', 
        'parameters'         : ['m1', 'm2', 'chi1', 'chi2'],
        'bounds'             : [],
        'modes'              : [(2,2)],
        'stack-mode'         : 'event',
        'compare'            : '',
        'compare-hard'       : 0,
        'evidence'           : 0,
        'save-post'          : 0,
        'ordering'           : [],
        'compare-ordering'   : [],
        'include-prior'      : 0,
        'ds-scaling'         : 0,

        'corner'             : 0,
        'violin'             : 0,
        'ridgeline'          : 0,
        'TGR-plot'           : 0,
        'plot-cpnest'        : '',
        'BF-comparison'      : 0,
        'evidence-top'       : 0,
        'time-percentiles'   : [],
        'palette'            : 'crest',
        'single-prior'       : '',
        'prior-color'        : '#828F61',
        'corner-settings'    : {'figsize': (15, 15), 'smooth': 0},
        'violin-settings'    : {'figsize': (15, 25), 'alpha': 0.5, 'rotation': 0},
        'ridgeline-settings' : {'figsize': (20, 10), 'alpha': 0.5, 'overlap': 0.5, 'fade': 0},
    }
    for key in input_pars.keys():

        if ('samp-dir' in key) or ('output' in key) or ('stack-mode' in key) or ('compare' in key):
            try: input_pars[key] = Config.get('input', key)
            except: pass
        if ('compare-hard' in key) or ('evidence' in key) or ('save-post' in key) or ('include-prior' in key) or ('ds-scaling' in key):
            try: input_pars[key] = Config.getboolean('input', key)
            except: pass
        if ('parameters' in key) or ('bounds' in key) or ('modes' in key) or ('ordering' in key) or ('compare-ordering' in key):
            try: input_pars[key] = ast.literal_eval(Config.get('input', key))
            except: pass
        if ('corner' in key) or ('violin' in key) or ('ridgeline' in key) or ('TGR-plot' in key) or ('BF-comparison' in key) or ('evidence-top' in key):
            try: input_pars[key] = Config.getboolean('plots', key)
            except: pass
        if ('plot-cpnest' in key) or ('single-prior' in key) or ('prior-color' in key):
            try: input_pars[key] = Config.getboolean('plots', key)
            except: pass
        if ('palette' in key) or ('time-percentiles' in key) or ('corner-settings' in key) or ('violin-settings' in key) or ('ridgeline-settings' in key):
            try: input_pars[key] = ast.literal_eval(Config.get('plots', key))
            except: pass

    input_pars['parent-dir'] = pardir_path
    samp_dir                 = create_directory(input_pars['parent-dir'], 'samples')
    input_pars['samp-dir']   = os.path.join(samp_dir, input_pars['samp-dir'])
    print('\nPosteriors are read from:\n{}'.format(input_pars['samp-dir']))

    res_dir = create_directory(input_pars['parent-dir'], 'results')
    out_dir = create_directory(res_dir                 , input_pars['output'])

    # Read the posteriors and create the .txt files with the reduced posteriors
    PostOutput = Posteriors(input_pars)
    SampDataFrame, PriorDataFrame, EvidenceDataFrame = PostOutput.return_samples_dict()

    if input_pars['save-post']:
        red_post_dir = create_directory(out_dir, 'reduced_posteriors')
        save_posteriors_to_txt(input_pars, red_post_dir, SampDataFrame)
        if input_pars['evidence'] == True:
            evidence_path = os.path.join(out_dir, 'evidence.txt')
            EvidenceDataFrame.to_csv(evidence_path, sep='\t', index=False)
            print('\nEvidences are saved in:\n{}\n'.format(evidence_path))

    if  input_pars['corner'] or input_pars['violin'] or input_pars['ridgeline'] or input_pars['TGR-plot']:
        input_pars['plots-dir'] = create_directory(out_dir, '')
        Plots(input_pars, SampDataFrame, PriorDataFrame, EvidenceDataFrame)
        print('\nPlots are saved in:\n{}'.format(input_pars['plots-dir']))
    
    print('\nFinished.\n')