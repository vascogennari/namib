from cgi import parse_header
import os
import h5py, pandas as pd
import numpy as np
from tqdm import tqdm

import pyRing.waveform as wf
import plots as plots


def from_IMR_to_RD_samples(df, pars):

    # granite
    if (set(['m1_detect', 'm2_detect']) <= set(df.keys())):
        df.rename(columns = {'m1_detect' : 'm1', 'm2_detect' : 'm2'}, inplace = True)
    if (set(['s1z', 's2z']) <= set(df.keys())):
        df.rename(columns = {'s1z' : 'chi1', 's2z' : 'chi2'}, inplace = True)
    if (set(['spin1', 'spin2']) <= set(df.keys())):
        df.rename(columns = {'spin1' : 'chi1', 'spin2' : 'chi2'}, inplace = True)
    if (set(['Mc', 'q']) <= set(df.keys())) or (set(['mc', 'q']) <= set(df.keys())):
        df = compute_progenitors_from_IMR(df)
    if (set(['Mc', 'q']) <= set(pars['parameters'])) or (set(['mc', 'q']) <= set(pars['parameters'])) or (set(['q']) <= set(pars['parameters'])) and (set(['m1', 'm2']) <= set(df.keys())):
        df = compute_progenitors_from_IMR(df, inverse = True)

    # Compute remnant pars
    if (set(['Mf', 'af']) <= set(pars['parameters'])) and (set(['m1', 'm2', 'chi1', 'chi2']) <= set(df.keys())):
        df = compute_Mf_af_from_IMR(df)
    if (set(['f_22', 'tau_22']) <= set(pars['parameters'])) and (set(['Mf', 'af']) <= set(df.keys())):
        df = compute_qnms_from_Mf_af(df,  pars['modes'])
    if (set(['f_22', 'tau_22']) <= set(pars['parameters'])) and (set(['m1', 'm2', 'chi1', 'chi2']) <= set(df.keys())) and not (set(['Mf', 'af']) <= set(pars['parameters'])):
        df = compute_Mf_af_from_IMR(df)
        df = compute_qnms_from_Mf_af(df, pars['modes'])

    # Damped-sinusoids
    if (set(['f_22', 'tau_22']) <= set(pars['parameters'])) and (set(['f_t_0', 'tau_t_0']) <= set(df.keys())):
        if pars['ds-scaling'] and (set(['f_t_0', 'tau_t_0']) <= set(df.keys())): df.tau_t_0 *= 1000  # Set time in [ms]
        df.rename(columns = {'f_t_0' : 'f_22', 'tau_t_0' : 'tau_22'}, inplace = True)
    if 'A2220' in set(pars['parameters']) and 'logA_t_0' in set(df.keys()):
        df['logA_t_0'] = df['logA_t_0'].apply(lambda x: np.exp(x))
        if pars['ds-scaling'] and 'logA_t_0' in set(df.keys()): df.logA_t_0 *= 1e10  # Scale amplitude as [1e-21]
        df.rename(columns = {'logA_t_0' : 'A2220'}, inplace = True)    
    if 'A2330' in set(pars['parameters']) and 'logA_t_1' in set(df.keys()):
        df['logA_t_1'] = df['logA_t_1'].apply(lambda x: np.exp(x))
        if pars['ds-scaling'] and 'logA_t_1' in set(df.keys()): df.logA_t_1 *= 1e10  # Scale amplitude as [1e-21]
        df.rename(columns = {'logA_t_1' : 'A2330'}, inplace = True)

    # Extrinsic parameters
    if (set(['distance']) <= set(pars['parameters'])) and (set(['logdistance']) <= set(df.keys())):
        df.insert(0, 'distance', np.exp(df.logdistance))
    if (set(['iota']) <= set(pars['parameters'])) and (set(['costheta_jn']) <= set(df.keys())):
        df.insert(0, 'iota', np.arccos(df.costheta_jn))
    if (set(['iota']) <= set(pars['parameters'])) and (set(['cosiota']) <= set(df.keys())):
        df.insert(0, 'iota', np.arccos(df.cosiota))

    # Set positive spins
    if 'chi1' in set(pars['parameters']):
        df['chi1'] = df['chi1'].apply(lambda x: np.abs(x))
    if (set(['chi2']) in set(pars['parameters'])):
        df['chi2'] = df['chi2'].apply(lambda x: np.abs(x))

    # TGR analysis
    if pars['TGR-plot']:
        if (set(['f_22'])   <= set(pars['parameters'])) and (set(['domega_220']) <= set(df.keys())):
            try:    df.f_22   *= 1. + df.domega_220
            except: pass
        if (set(['f_33'])   <= set(pars['parameters'])) and (set(['domega_330']) <= set(df.keys())):
            try:    df.f_33   *= 1. + df.domega_330
            except: pass
        if (set(['tau_22']) <= set(pars['parameters'])) and (set(['dtau_220'])   <= set(df.keys())):
            try:    df.tau_22 *= 1. + df.dtau_220
            except: pass

    if not (set(pars['parameters']).difference(df.keys()) == set()):
        additional_pars = set(pars['parameters']).difference(df.keys())
        for additional_par in additional_pars:
            df.insert(0, additional_par, np.nan)
            df[additional_par][0] = 0

    return df

def read_posteriors_event(file_path, pars):
    '''
    Read the posteriors distribution of a single file.
    The posterior distributions for the passed parameters are returned in a Pandas DF.
    '''
    filename = os.path.basename(os.path.normpath(file_path))
    if ('.txt' in filename) or ('.dat' in filename):
        load = np.genfromtxt(file_path, names = True)
    if '.h5' in filename:
        with h5py.File(file_path, 'r') as f:
            try:
                tmp = f['EXP1']['posterior_samples']
            except:
                tmp = f['combined']['posterior_samples']   # IMPROVE ME: Check the CPNest version
                if pars['include-prior']:
                    try:    tmpp = f['combined']['prior_samples']
                    except: raise ValueError('Invalid option for prior reading: cannot find prior samples. Exiting...')
            load = np.array(tmp)
            if pars['include-prior']: loadp = np.array(tmpp)

    df = pd.DataFrame(load)
    df = from_IMR_to_RD_samples(df, pars)
    df = df.filter(items = pars['parameters'])

    if pars['include-prior']:
        dfp = pd.DataFrame(loadp)
        dfp = from_IMR_to_RD_samples(dfp, pars)
        dfp = dfp.filter(items = pars['parameters'])
    else:
        dfp = pd.DataFrame()

    nsamp = len(df)

    return df, dfp, nsamp

def read_evidence_event(pars, dir_path, file_path):

    evt_evidence = {}

    # Read the noise evidence
    filename = os.path.basename(os.path.normpath(file_path))
    if   '.txt' in filename: root_file = filename.replace('.txt', '')
    elif '.dat' in filename: root_file = filename.replace('.dat', '')
    elif '.h5'  in filename: root_file = filename.replace('.h5' , '')

    noise_file = root_file + '_noise.txt'
    noise_evt_path = os.path.join(dir_path, 'noise_evidences', noise_file)

    tmp  = np.genfromtxt(noise_evt_path, names = True)
    evt_evidence['lnZ_noise'] = np.array(tmp['lnZ_noise'])
    if '.dat' in filename:
        evt_evidence['lnZ'] = np.array(tmp['lnZ_signal'])

    # Read the signal evidence
    filename = os.path.basename(os.path.normpath(file_path))
    if ('.txt' in filename) or ('.dat' in filename):
        load = np.genfromtxt(file_path, names = True)
    if '.h5' in filename:
        with h5py.File(file_path, 'r') as f:
            try:
                evt_evidence['lnZ']       = np.array(f['EXP1']['logZ'])
                evt_evidence['lnZ_error'] = np.array(f['EXP1']['logZ_error'])
                evt_evidence['H']         = np.array(f['EXP1']['information'])
            except:
                evt_evidence['lnZ']       = np.array(f['combined']['logZ'])
                evt_evidence['lnZ_error'] = np.array(f['combined']['logZ_error'])
                evt_evidence['H']         = np.array(f['combined']['information'])
    evt_evidence['lnB'] = evt_evidence['lnZ'] - evt_evidence['lnZ_noise']

    if pars['save-medians']:
        try:
            # Read the SNR samples
            SNR_file = root_file + '_SNR.dat'
            SNR_evt_path = os.path.join(dir_path, 'SNR_samples', SNR_file)
            tmp  = np.genfromtxt(SNR_evt_path, names = True)
            evt_evidence['SNR_net_opt'] = np.array(tmp['Network'])
        except:
            raise ValueError('SNR samples have not been found. Exiting.')

    df = pd.DataFrame([evt_evidence])

    return df

def stack_posteriors(stack_list, post_list):
    '''
    Stack the posteriors DF into a single dictionary.
    '''
    samp_dict = {}
    for model,post in zip(stack_list, post_list):
        samp_dict[model] = post
    
    return samp_dict

def add_parameters_to_event(evt_df, new_params_list, new_params_samps):
    '''
    Add new parameters to the single event DF.
    '''
    df = evt_df
    for par,samp in zip(new_params_list, new_params_samps):
        df[par] = samp
    
    return df

def compute_Mf_af_from_IMR(df):
    '''
    Compute Mf and af of the remnant BH from IMR parameters, using the Jimenez-Forteza fits implemented in pyRing.
    '''
    nsamp = len(df)
    Mf = np.zeros(nsamp)
    af = np.zeros(nsamp)
    for i in range(nsamp):
        m1, m2, chi1, chi2 = df.m1[i], df.m2[i], df.chi1[i], df.chi2[i]
        tmp   = wf.TEOBPM(0, m1, m2, chi1, chi2, {}, 100, 0, 0, [], {})
        Mf[i] = tmp.JimenezFortezaRemnantMass()   # [M_\odot]
        af[i] = tmp.JimenezFortezaRemnantSpin()   # []
    df.insert(0, 'Mf', Mf)
    df.insert(0, 'af', af)
    
    return df

def compute_qnms_from_Mf_af(df, modes):
    '''
    Compute QNMs frequency and damping time from Mf and af for one mode (l,m), using the Berti's fits implemented in pyRing.
    '''
    nsamp = len(df)
    for mode in modes:
        l, m = mode[0], mode[1]
        omg = np.zeros(nsamp)
        tau = np.zeros(nsamp)
        for i in range(nsamp):
            Mf, af = df.Mf[i], df.af[i]
            omg[i] = wf.QNM_fit(l, m, 0).f(Mf, af)            # [Hz]
            tau[i] = wf.QNM_fit(l, m, 0).tau(Mf, af) * 1000   # [ms]
        df.insert(0, 'f_{}{}'.format(l,m),   omg)
        df.insert(0, 'tau_{}{}'.format(l,m), tau)

    return df

def compute_progenitors_from_IMR(df, inverse = False):

    def McQ2Masses(mc, q):
        factor = mc * np.power(1. + q, 1.0/5.0)
        m1     = factor * np.power(q, -3.0/5.0)
        m2     = factor * np.power(q, +2.0/5.0)
        return m1, m2
    
    def Masses2McQ(m1, m2):
        q  = m2/m1
        mc = (m1*m2)**(3./5.)/(m1+m2)**(1./5.)
        return mc, q
    
    if not inverse: df['m1'], df['m2'] = McQ2Masses(df['mc'], df['q'])
    else:           df['mc'], df['q']  = Masses2McQ(df['m1'], df['m2'])

    return df

def clean_empty_keys_violin(samp_L, samp_R):

    tmp = np.array([1000,1001])
    for samp in [samp_L, samp_R]:
        for i,elem in enumerate(samp):
            if   not list(elem):       samp[i] = tmp
            elif np.isnan(elem).any(): samp[i] = tmp

    return samp_L, samp_R

def clean_empty_keys_corner(samp):

    tmp = np.array([0,1])
    for i,elem in enumerate(samp.T):
        if np.isnan(elem).any(): samp.T[i] = np.zeros(len(elem))
    
    return samp

def create_directory(parent_path, name):

    dir_path = os.path.join(parent_path, name)
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

    return dir_path

def compute_bayes_factor(pars, df):

    keys      = plots.sort_times_list(pd.unique(df[pars['stack-mode']]))
    comp_pars = pd.unique(df[pars['compare']])
    if not pars['compare-ordering'] == []:
        if ((set(pars['compare-ordering']) <= set(comp_pars))) and (len(pars['compare-ordering']) == len(comp_pars)): comp_pars = pars['compare-ordering']
        else: raise ValueError('Invalid option for {compare} ordering.'.format(compare = pars['compare-ordering']))
    df.insert(0, 'Bayes_factor', 0)
    df.insert(0, 'Bayes_factor_error', 0)

    for key in keys:
        evidence_L = (df[pars['stack-mode']] == key) & (df[pars['compare']] == comp_pars[0])
        evidence_R = (df[pars['stack-mode']] == key) & (df[pars['compare']] == comp_pars[1])
        df.loc[df.index[df[pars['stack-mode']] == key].tolist(), 'Bayes_factor'] =               float(df[evidence_R]['lnZ'].iloc[0])       - float(df[evidence_L]['lnZ'].iloc[0])
        try:    df.loc[df.index[df[pars['stack-mode']] == key].tolist(), 'Bayes_factor_error'] = float(df[evidence_R]['lnZ_error'].iloc[0]) + float(df[evidence_L]['lnZ_error'].iloc[0])
        except: df.loc[df.index[df[pars['stack-mode']] == key].tolist(), 'Bayes_factor_error'] = 0.1

    return df

def save_posteriors_to_txt(pars, path, df):

    keys = plots.sort_times_list(pd.unique(df[pars['stack-mode']]))
    if not pars['compare'] == '': comp_pars = pd.unique(df[pars['compare']])
    else:                         comp_pars = 'rand'
    for comp in comp_pars:
        if not pars['compare'] == '': df_comp = df[df[pars['compare']] == comp]
        else:                         df_comp = df
        for key in keys:
            if not pars['compare'] == '': filename = '{path}/posteriors_{mode}_{key}_{comp}.txt'.format(path = path, mode = pars['stack-mode'], key = key, comp = comp)
            else:                         filename = '{path}/posteriors_{mode}_{key}.txt'.format(       path = path, mode = pars['stack-mode'], key = key)
            df_filt = df_comp[df_comp[pars['stack-mode']] == key]
            df_filt = df_filt[pars['parameters']]
            df_filt.to_csv(filename, sep='\t', index=False)

    print('\nProcessed posteriors and are saved in:\n{}'.format(path))

    return 0

def save_output_medians(pars, df, df_evidence, out_dir):

    # Median values of the selected parameters
    output_path_pars = os.path.join(out_dir, 'parameters.txt')
    if os.path.isfile(output_path_pars): os.remove(output_path_pars)
    keys = plots.sort_times_list(pd.unique(df[pars['stack-mode']]))
    if not pars['compare'] == '': comp_pars = pd.unique(df[pars['compare']])
    else:                         comp_pars = 'rand'
    pars_header = '{}\t'.format(pars['stack-mode'])
    for par in pars['parameters']: pars_header += '{}\t\t\t'.format(par)

    with open(output_path_pars, 'a') as f:
        for comp in comp_pars:
            if not pars['compare'] == '': df_comp = df[df[pars['compare']] == comp]
            else:                         df_comp = df
            if comp == 'rand': tmp = ''
            else:              tmp = comp
            f.write('{}\n'.format(tmp))
            f.write('{}\n'.format(pars_header))

            for key in keys:
                df_filt = df_comp[df_comp[pars['stack-mode']] == key]
                pars_line = '{}\t'.format(key)
                for par in pars['parameters']:
                    median, low_perc, upp_perc = np.percentile(df_filt[par], np.array([0.5, 0.1, 0.9]))
                    low_err = median - low_perc
                    upp_err = upp_perc - median
                    pars_line += '{1:.{0}f}\t+{2:.{0}f}\t-{3:.{0}f}\t'.format(2, median, upp_err, low_err)
                f.write('{}\n'.format(pars_line))
            f.write('\n')
    print('\nMedian values of the parameters are saved in:\n{}\n'.format(output_path_pars))

    # Bayes factors, informations, median values of the network optimal SNR
    output_path = os.path.join(out_dir, 'SNR-H-BF.txt')
    if os.path.isfile(output_path): os.remove(output_path)
    pars_header  = '{}\t'.format(pars['stack-mode'])
    pars_header += 'SNR\t\t\tH\tBF\t'.format(pars['stack-mode'])

    with open(output_path, 'a') as f:
        for comp in comp_pars:
            if not pars['compare'] == '': df_comp = df_evidence[df_evidence[pars['compare']] == comp]
            else:                         df_comp = df_evidence
            if comp == 'rand': tmp = ''
            else:              tmp = comp
            f.write('{}\n'.format(tmp))
            f.write('{}\n'.format(pars_header))

            for key in keys:
                df_filt = df_comp[df_comp[pars['stack-mode']] == key]
                median, low_perc, upp_perc = np.percentile(df_filt['SNR_net_opt'].values[0], np.array([0.5, 0.1, 0.9]))
                low_err = median - low_perc
                upp_err = upp_perc - median
                if pars['compare-hard']:
                    BF     = df_filt['Bayes_factor'].values[0]
                    BF_err = df_filt['Bayes_factor_error'].values[0]
                else:
                    BF     = df_filt['lnB'].values[0]
                    BF_err = df_filt['lnZ_error'].values[0]

                pars_line =  '{}\t'.format(key)
                pars_line += '{1:.{0}f}\t+{2:.{0}f}\t-{3:.{0}f}\t'.format(2, median, upp_err, low_err)
                pars_line += '{1:.{0}f}\t'.format(2, df_filt['H'].values[0])
                if not pars['compare'] == '' and comp == comp_pars[0]: pass
                else:                                                  pars_line += '{1:.{0}f}\t+-{2:.{0}f}'.format(2, BF, BF_err)
                f.write('{}\n'.format(pars_line))
            f.write('\n')
    print('Median values of SNR, information and BF are saved in:\n{}\n'.format(output_path))

    return 0


class Posteriors:
    '''
    Returns a data frame with all the posterior samples, for a given list of events and models.
    Also, it computes the remnant properties Mf and af, and the QNMs frequencies and damping times when possible.
    '''
    def __init__(self, pars):

        dir_path = pars['samp-dir']
        if not (pars['file-path'] == ''): dir_path = pars['file-path']

        self.SampDataFrame     = pd.DataFrame(columns = pars['parameters'])
        self.PriorDataFrame    = pd.DataFrame(columns = pars['parameters'])
        self.EvidenceDataFrame = pd.DataFrame()
        single_evt_keys = {'event': str(), 'pipeline': str(), 'model': str(), 'submodel': str(), 'time': str(), 'GR_tag': str()}

        for file in tqdm(os.listdir(dir_path), desc = 'Reading Posteriors'):
            if not ((file == '.DS_Store') or (file == 'noise_evidences') or (file == 'ignore') or (file == 'SNR_samples')):

                file_path = os.path.join(dir_path, file)
                keys = file.split('_')
                keys[-1] = keys[-1].split('.')[0]
                for i,key in enumerate(single_evt_keys.keys()):
                    single_evt_keys[key] = keys[i]

                EventDataFrame, EventPriorDataFrame, _ = read_posteriors_event(file_path, pars)
                if not pars['stack-mode'] == '':
                    EventDataFrame = EventDataFrame.assign(par = single_evt_keys[pars['stack-mode']])
                    EventDataFrame.rename(columns={'par': pars['stack-mode']}, inplace = True)
                    if pars['include-prior']:
                        EventPriorDataFrame = EventPriorDataFrame.assign(par = single_evt_keys[pars['stack-mode']])
                        EventPriorDataFrame.rename(columns={'par': pars['stack-mode']}, inplace = True)
                if not pars['compare'] == '':
                    EventDataFrame = EventDataFrame.assign(par = single_evt_keys[pars['compare']])
                    EventDataFrame.rename(columns={'par': pars['compare']}, inplace = True)                

                self.SampDataFrame      = pd.concat([self.SampDataFrame,  EventDataFrame],      ignore_index=True)
                if pars['include-prior']:
                    self.PriorDataFrame = pd.concat([self.PriorDataFrame, EventPriorDataFrame], ignore_index=True)
                if pars['evidence']:
                    EventEvidenceDataFrame = read_evidence_event(pars, dir_path, file_path)
                    EventEvidenceDataFrame.insert(0, pars['stack-mode'], single_evt_keys[pars['stack-mode']])
                    if not pars['compare'] == '': EventEvidenceDataFrame.insert(0, pars['compare'], single_evt_keys[pars['compare']])
                    self.EvidenceDataFrame = pd.concat([self.EvidenceDataFrame, EventEvidenceDataFrame], ignore_index=True)

        if (not pars['compare'] == '') and pars['BF-comparison']: self.EvidenceDataFrame = compute_bayes_factor(pars, self.EvidenceDataFrame)

    def return_samples_dict(self):
        return self.SampDataFrame, self.PriorDataFrame, self.EvidenceDataFrame


class Plots:
    '''
    Compute corner and violin plots.
    '''
    def __init__(self, pars, SampDataFrame, PriorDataFrame, EvidenceDataFrame):
        
        if pars['corner']:    plots.corner_plots(   pars, SampDataFrame, PriorDataFrame)
        if pars['violin']:    plots.violin_plots(   pars, SampDataFrame, PriorDataFrame, EvidenceDataFrame)
        if pars['ridgeline']: plots.ridgeline_plots(pars, SampDataFrame, PriorDataFrame)
        if pars['TGR-plot']:  plots.TGR_plots(      pars, SampDataFrame)
