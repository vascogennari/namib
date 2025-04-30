from multiprocessing import Value
import os, h5py, pandas as pd
import fnmatch
import numpy as np
from tqdm import tqdm
import qnm
from astropy import constants as const

import namib.plots as plots

import warnings
warnings.filterwarnings('ignore')


def Adapt_Samples(df, pars, IMR_flag = False):
    '''
        Adapt the input samples to namib internal conventions.
        Compute remnant paramters and QNMs if needed.
        Remove unnecessary paramters from the data frame.
    '''
    if IMR_flag: IMR_fits = pars['IMR-fits-IMR']
    else:        IMR_fits = pars['IMR-fits']

    def LVK_conventions(df, pars):
        if (set(['mass_1', 'mass_2']) <= set(df.keys())):
            df.rename(columns = {'mass_1' : 'm1', 'mass_2' : 'm2'}, inplace = True)
        if (set(['distance']) <= set(pars['parameters'])) and (set(['luminosity_distance']) <= set(df.keys())):
            df.insert(0, 'distance', df.luminosity_distance / 1000)     # [Gpc]
        if (set(['cosiota']) <= set(pars['parameters'])) and (set(['cos_theta_jn']) <= set(df.keys())):
            df.rename(columns = {'cos_theta_jn' : 'cosiota'}, inplace = True)
        if (set(['iota']) <= set(pars['parameters'])) and (set(['cos_theta_jn']) <= set(df.keys())):
            df.insert(0, 'iota', np.arccos(df.cos_theta_jn))
        if (set(['a_1', 'a_2']) <= set(df.keys())):
            df.rename(columns = {'a_1' : 'chi1', 'a_2' : 'chi2'}, inplace = True)

    def granite_conventions(df, pars):
        if (set(['m1_detect', 'm2_detect']) <= set(df.keys())):
            df.rename(columns = {'m1_detect' : 'm1', 'm2_detect' : 'm2'}, inplace = True)
        if (set(['s1z', 's2z']) <= set(df.keys())):
            df.rename(columns = {'s1z' : 'chi1', 's2z' : 'chi2'}, inplace = True)
        if (set(['spin1', 'spin2']) <= set(df.keys())):
            df.rename(columns = {'spin1' : 'chi1', 'spin2' : 'chi2'}, inplace = True)

    def compute_remnant_from_IMR(df, pars):
        if not (set(['Mf', 'af']) <= set(df.keys())) and not (set(['f_t_0', 'tau_t_0']) <= set(df.keys())):
            if IMR_fits == 'IMRPhenomXPrecessing':
                if not (set(['eta', 'chi_p']) <= set(df.keys())):
                    df = compute_progenitors_from_IMR(df, func = 'SymmetricMassRatio')
                    df = compute_progenitors_from_IMR(df, func = 'ChiSymmetric')
            if IMR_fits == 'NRSur7dq4Remnant':
                df = compute_progenitors_from_IMR(df, func = 'MassRatio')
                df = remove_mass_ratio_over_threshold(df)
            df = compute_Mf_af_from_IMR(df, pars, IMR_fits)
        return df

    def compute_phase_amplitude_from_IMR(df, pars):
        if not ((set(['A2220']) <= set(df.keys())) or (set(['phi2220']) <= set(df.keys())) or (set(['A2220_1']) <= set(df.keys())) or (set(['phi2220_1']) <= set(df.keys()))) and not (set(['f_t_0', 'tau_t_0']) <= set(df.keys())):
            if not (set(['eta', 'chi_p', 'chi_a']) <= set(df.keys())):
                df = compute_progenitors_from_IMR(df, func = 'SymmetricMassRatio')
                df = compute_progenitors_from_IMR(df, func = 'ChiSymmetric')
                df = compute_progenitors_from_IMR(df, func = 'ChiAntiymmetric')
            df = compute_phase_amplitude_from_progenitors(df, pars['modes'])
        return df
    
    def compute_phase_amplitude_test(df,pars):
        if (set(['A2220']) <= set(df.keys()) or set(['phi2220']) <= set(df.keys())):
            for mode in pars['modes']:
                l, m, n = mode[0], mode[1], mode[2]
                if not (mode == (2,2,0)):

                    df[f'AR{l}{m}{n}'] = df[f'A2{l}{m}{n}']/df['A2220']
                    df[f'deltaphi{l}{m}{n}'] = (m/2.)*df['phi2220'] - df[f'phi2{l}{m}{n}']

                    nsamp = len(df)
                    deltaphi = df[f'deltaphi{l}{m}{n}']
                    if (m % 2 == 0): 
                        cc=2.
                    else:
                        cc = 1.
                    for i in range(nsamp):
                        while (deltaphi[i] > cc*np.pi) or (deltaphi[i] < 0):
                            if deltaphi[i] > cc*np.pi: deltaphi[i] = deltaphi[i] - cc*np.pi
                            if deltaphi[i] < 0:        deltaphi[i] = deltaphi[i] + cc*np.pi
                    df[f'deltaphi{l}{m}{n}'] = deltaphi
        return df

    def compute_qnms_from_remnant(df, pars):
        if (set(['f_22', 'tau_22']) <= set(pars['parameters'])) and not (set(['f_22', 'tau_22']) <= set(df.keys())) and not (set(['f_t_0', 'tau_t_0']) <= set(df.keys())):
            if not (set(['Mf', 'af']) <= set(df.keys())):
                df = compute_remnant_from_IMR(df, pars)
            df = compute_qnms_from_Mf_af(df, pars['modes'], pars, scaling = 1)
        if (set(['f_t_0', 'tau_t_0']) <= set(pars['parameters'])) and not (set(['f_t_0', 'tau_t_0']) <= set(df.keys())):
            if not (set(['Mf', 'af']) <= set(df.keys())):
                df = compute_remnant_from_IMR(df, pars)
            df = compute_qnms_from_Mf_af(df, [(2,2,0)], pars, scaling = 0)
            df.rename(columns = {'f_22' : 'f_t_0', 'tau_22' : 'tau_t_0'}, inplace = True)
            # if (set(['f_t_0', 'tau_t_0']) <= set(df.keys())):
            # df.rename(columns = {'f_t_0' : 'f_22', 'tau_t_0' : 'tau_22'}, inplace = True)
        return df

    def pyring_damped_sinusoids_conventions(df, pars):
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

    def extrinsic_parameters_conventions(df, pars):
        if (set(['distance']) <= set(pars['parameters'])) and (set(['logdistance']) <= set(df.keys())):
            df.insert(0, 'distance', np.exp(df.logdistance) / 1000)     # [Gpc]
        if (set(['iota']) <= set(pars['parameters'])) and (set(['costheta_jn']) <= set(df.keys())):
            df.insert(0, 'iota', np.arccos(df.costheta_jn))
        if (set(['iota']) <= set(pars['parameters'])) and (set(['cosiota']) <= set(df.keys())):
            df.insert(0, 'iota', np.arccos(df.cosiota))

    def set_positive_spins(df, pars):
        if (set(['chi1']) <= set(pars['parameters'])) and (set(['chi1']) <= set(df.keys())):
            df['chi1'] = df['chi1'].apply(lambda x: np.abs(x))
        if (set(['chi2']) <= set(pars['parameters'])) and (set(['chi2']) <= set(df.keys())):
            df['chi2'] = df['chi2'].apply(lambda x: np.abs(x))
    
    def compute_dependent_parameters(df, pars):
        if (set(['q'])     <= set(pars['parameters'])) and (set(['m1', 'm2']) <= set(df.keys())):
            df = compute_progenitors_from_IMR(df, func = 'MassRatio')
        if (set(['mc'])    <= set(pars['parameters'])) and (set(['m1', 'm2']) <= set(df.keys())):
            df = compute_progenitors_from_IMR(df, func = 'Masses2McQ')
        if (set(['eta'])   <= set(pars['parameters'])) and (set(['m1', 'm2']) <= set(df.keys())):
            df = compute_progenitors_from_IMR(df, func = 'SymmetricMassRatio')
        if (set(['chi_s']) <= set(pars['parameters'])) and (set(['m1', 'm2', 'chi1', 'chi2']) <= set(df.keys())):
            df = compute_progenitors_from_IMR(df, func = 'ChiSymmetric')
        if (set(['chi_a']) <= set(pars['parameters'])) and (set(['m1', 'm2', 'chi1', 'chi2']) <= set(df.keys())):
            df = compute_progenitors_from_IMR(df, func = 'ChiAntiymmetric')

    # FIXME: Implement as a separate option non related to the TGR plot, for all the possible modes.
    def TGR_plot(df, pars):
        if (set(['f_22'])   <= set(pars['parameters'])) and (set(['domega_220']) <= set(df.keys())):
            try:    df.f_22   *= 1. + df.domega_220
            except: pass
        if (set(['f_33'])   <= set(pars['parameters'])) and (set(['domega_330']) <= set(df.keys())):
            try:    df.f_33   *= 1. + df.domega_330
            except: pass
        if (set(['tau_22']) <= set(pars['parameters'])) and (set(['dtau_220'])   <= set(df.keys())):
            try:    df.tau_22 *= 1. + df.dtau_220
            except: pass

    LVK_conventions(                      df, pars)
    granite_conventions(                  df, pars)
    if (set(['Mf', 'af']) <= set(pars['parameters'])):
        df = compute_remnant_from_IMR(    df, pars)
    df = compute_phase_amplitude_from_IMR(df, pars)
    df = compute_phase_amplitude_test(    df, pars)
    df = compute_qnms_from_remnant(       df, pars)
    pyring_damped_sinusoids_conventions(  df, pars)
    extrinsic_parameters_conventions(     df, pars)
    set_positive_spins(                   df, pars)
    compute_dependent_parameters(         df, pars)
    if pars['TGR-plot']: TGR_plot(        df, pars)

    if not (set(pars['parameters']).difference(df.keys()) == set()):
        additional_pars = set(pars['parameters']).difference(df.keys())
        for additional_par in additional_pars:
            df.insert(0, additional_par, np.nan)
            df[additional_par][0] = 0

    return df

def read_posteriors_event(file_path, pars, IMR_flag = False):
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
                tmp = f['C01:IMRPhenomPv2']['posterior_samples']
            except:
                tmp = f['combined']['posterior_samples']   # IMPROVE ME: Check the CPNest version
                if pars['include-prior']:
                    try:    tmpp = f['combined']['prior_samples']
                    except: raise ValueError('Invalid option for prior reading: cannot find prior samples. Exiting...')
            load = np.array(tmp)
            if pars['include-prior']: loadp = np.array(tmpp)
    if '.hdf5' in filename:
        with h5py.File(file_path, 'r') as f:
            try:
                tmp = f['posterior']
            except: raise ValueError('Invalid structure for samples reading: cannot find posterior samples. Exiting...')
            load = {key: np.array(tmp[key]) for key in tmp.keys()}

    df = pd.DataFrame(load)
    df = downsampling(df, pars)    # Downsample the df if required
    df = Adapt_Samples(df, pars, IMR_flag = IMR_flag)
    df = df.filter(items = pars['parameters'])

    if pars['include-prior']:
        dfp = pd.DataFrame(loadp)
        dfp = Adapt_Samples(dfp, pars)
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
    if ('.dat' in filename) or ('.txt' in filename):
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
                evt_evidence['max_logL']  = np.max(np.array(f['EXP1']['posterior_samples']['logL']))
            except:
                evt_evidence['lnZ']       = np.array(f['combined']['logZ'])
                evt_evidence['lnZ_error'] = np.array(f['combined']['logZ_error'])
                evt_evidence['H']         = np.array(f['combined']['information'])
                evt_evidence['max_logL']  = np.max(np.array(f['combined']['posterior_samples']['logL']))
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

def downsampling(df, pars):
    '''
    Return the data frame downsampled according to the required probability.
    downsample = 1 takes the 100% of the data, i.e. no downsampling
    '''
    if not pars['downsample'] == 1:
        if (pars['downsample'] <= 0) or (pars['downsample'] > 1):
            raise ValueError('Invalid option for the downsampling. Its value needs to be in the interval [0, 1].')
        
        new_nsamp = int(len(df.index) * pars['downsample'])
        df = df.sample(new_nsamp)
        df = df.reset_index()

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

def compute_Mf_af_NRSur(df):
    try:    from pesummary.gw.conversions.nrutils import NRSur_fit
    except: raise ValueError('Unable to find the NRSur remnant fits. Please either install pesummary and  make sure that the LAL_DATA_PATH is properly set, or use a different option for the remnant fits.')

    if not (set(['m1', 'm2', 'chi1', 'chi2', 'tilt_1', 'tilt_2', 'phi_12', 'phi_jl', 'theta_jn', 'phase']) <= set(df.columns)):
        raise ValueError('The IMR samples are not compatible with the selected remnant fits. Please make sure they are consistent.')

    fits = NRSur_fit(df.m1, df.m2, df.chi1, df.chi2, df.tilt_1, df.tilt_2, df.phi_12, df.phi_jl, df.theta_jn, df.phase,
                        20.0, np.full_like(df.m1, 20.0),
                        model       = 'NRSur7dq4Remnant',
                        approximant = 'IMRPhenomXPHM')

    Mf = fits['final_mass']
    af = fits['final_spin']
    
    df.insert(0, 'Mf', Mf)
    df.insert(0, 'af', af)

    return df

def compute_Mf_af_from_IMR(df, pars, IMR_fits):
    '''
    Compute Mf and af of the remnant BH from IMR parameters. Both aligned-spin and precessing fits are implemented.
    Current options are: JimenezForteza_TEOBPM, UIB2016, NRSur7dq4Remnant, IMRPhenomXPrecessing.
    '''
    # Aligned-spin fits.
    if   IMR_fits == 'JimenezForteza_TEOBPM':
        try:    import pyRing.waveform as wf
        except: raise ValueError('Unable to find the pyRing remnant fits. Please either install pyRing or use a different option for the remnant fits.')

        if not (set(['m1', 'm2', 'chi1', 'chi2']) <= set(df.columns)):
            raise ValueError('The IMR samples are not compatible with the selected remnant fits. Please make sure they are consistent.')

        nsamp = len(df)
        Mf = np.zeros(nsamp)
        af = np.zeros(nsamp)

        for i in range(nsamp):
            tmp   = wf.TEOBPM(0, df.m1[i], df.m2[i], df.chi1[i], df.chi2[i], {}, 100, 0, 0, [], {})
            Mf[i] = tmp.JimenezFortezaRemnantMass()
            af[i] = tmp.JimenezFortezaRemnantSpin()

    # Aligned-spin fits.
    elif IMR_fits == 'UIB2016':
        try:    from lalinference.imrtgr.nrutils import bbh_final_mass_projected_spins, bbh_final_spin_projected_spins, bbh_Kerr_trunc_opts
        except: raise ValueError('Unable to find the UIB2016 remnant fits. Please either install lalinference or use a different option for the remnant fits.')

        if not (set(['m1', 'm2', 'chi1', 'chi2']) <= set(df.columns)):
            raise ValueError('The IMR samples are not compatible with the selected remnant fits. Please make sure they are consistent.')

        nsamp = len(df)
        Mf = np.zeros(nsamp)
        af = np.zeros(nsamp)

        for i in range(nsamp):
            # Adapt to final state fits conventions
            if df.chi1[i] < 0: tilt1 = np.pi
            else:              tilt1 = 0.0
            if df.chi2[i] < 0: tilt2 = np.pi
            else:              tilt2 = 0.0
            chi1_abs = np.abs(df.chi1[i])
            chi2_abs = np.abs(df.chi2[i])

            Mf[i] = bbh_final_mass_projected_spins(df.m1[i], df.m2[i], chi1_abs, chi2_abs, tilt1, tilt2, 'UIB2016')
            af[i] = bbh_final_spin_projected_spins(df.m1[i], df.m2[i], chi1_abs, chi2_abs, tilt1, tilt2, 'UIB2016', truncate = bbh_Kerr_trunc_opts.trunc)

    # Precessing fits.
    elif IMR_fits == 'NRSur7dq4Remnant':

        # This option requires the LAL_DATA_PATH to be set correctly, cloning https://git.ligo.org/lscsoft/lalsuite-extra and pointing to lalsuite-extra/data/lalsimulation:
        # export LAL_DATA_PATH=~/lalsuite-extra/data/lalsimulation
        try:    from pesummary.gw.conversions.nrutils import NRSur_fit
        except: raise ValueError('Unable to find the NRSur remnant fits. Please either install pesummary and  make sure that the LAL_DATA_PATH is properly set, or use a different option for the remnant fits.')

        if not (set(['m1', 'm2', 'chi1', 'chi2', 'tilt_1', 'tilt_2', 'phi_12', 'phi_jl', 'theta_jn', 'phase']) <= set(df.columns)):
            raise ValueError('The IMR samples are not compatible with the selected remnant fits. Please make sure they are consistent.')

        fits = NRSur_fit(df.m1, df.m2, df.chi1, df.chi2, df.tilt_1, df.tilt_2, df.phi_12, df.phi_jl, df.theta_jn, df.phase,
                            20.0, np.full_like(df.m1, 20.0),
                            model       = 'NRSur7dq4Remnant',
                            approximant = 'IMRPhenomXPHM')

        Mf = fits['final_mass']
        af = fits['final_spin']

    # Precessing fits.
    elif IMR_fits == 'IMRPhenomXPrecessing':
        try:    from lalsimulation import SimIMRPhenomXPrecessingFinalSpin2017, SimIMRPhenomXFinalMass2017
        except: raise ValueError('Unable to find the IMRPhenomXPrecessing remnant fits. Please either install lalsimulation or use a different option for the remnant fits.')

        if not (set(['eta', 'chi1', 'chi2', 'chi_p']) <= set(df.columns)):
            raise ValueError('The IMR samples are not compatible with the selected remnant fits. Please make sure they are consistent.')
        
        Mf = SimIMRPhenomXFinalMass2017(          df.eta, df.chi1, df.chi2)
        af = SimIMRPhenomXPrecessingFinalSpin2017(df.eta, df.chi1, df.chi2, df.chi_p)

    # Precessing fits.
    elif IMR_fits == 'HBR2016':
        try:    from pesummary.gw.conversions.nrutils import bbh_final_mass_non_precessing_UIB2016, bbh_final_spin_precessing_HBR2016
        except: raise ValueError('Unable to find the HBR2016 remnant fits. Please either install pesummary and  make sure that the LAL_DATA_PATH is properly set, or use a different option for the remnant fits.')

        if not (set(['m1', 'm2', 'chi1', 'chi2', 'a_1', 'a_2', 'tilt_1', 'tilt_2', 'phi_12']) <= set(df.columns)):
            raise ValueError('The IMR samples are not compatible with the selected remnant fits. Please make sure they are consistent.')

        fits_m = bbh_final_mass_non_precessing_UIB2016(df.m1, df.m2, df.chi1, df.chi2)
        fits_a = bbh_final_spin_precessing_HBR2016(    df.m1, df.m2, df.a_1,  df.a_2, df.tilt_1, df.tilt_2, df.phi_12)

        Mf = fits_m['final_mass']
        af = fits_a['final_spin']

    else:
        raise ValueError('Option not available for the NR fits to compute remnant parameters from IMR samples. Exiting.')

    df.insert(0, 'Mf', Mf)
    df.insert(0, 'af', af)

    return df

def compute_qnms_from_Mf_af(df, modes, pars, scaling = 1):
    '''
    Compute QNMs frequency and damping time from Mf and af for one mode (l,m,n)
    using the qnm python package [https://github.com/duetosymmetry/qnm]
    '''
    nsamp = len(df)

    for mode in modes:
        l, m, n = mode[0], mode[1], mode[2]
        omg = np.zeros(nsamp)
        tau = np.zeros(nsamp)
        for i in range(nsamp):
            Mf, af = df.Mf[i], df.af[i]
            if not pars['qnms-pyRing']:
                Warning('Using qnm fits to compute the remnant samples [Mf, af]. This option is still experimental and it is currently very slow: we suggest to use the option "qnms-pyRing".')
                omg[i], tau[i] = get_qnms(Mf, af, l, m, n)
            else:
                try:
                    import pyRing.waveform as wf
                except:
                    raise ValueError('Unable to find the pyRing installation for the QNMs fits. Please either install pyRing or disactivate the option "qnms-pyRing".')
                omg[i] = wf.QNM_fit(l, m, n).f(Mf, af)                # [Hz]
                if scaling == 1:
                    tau[i] = wf.QNM_fit(l, m, n).tau(Mf, af) * 1000   # [ms]
                else:
                    tau[i] = wf.QNM_fit(l, m, n).tau(Mf, af)          # [ms]

        df.insert(0, 'f_{}{}'.format(l,m),   omg)
        df.insert(0, 'tau_{}{}'.format(l,m), tau)

    return df

def get_qnms(Mf, af, l, m, n = 0):
    '''
        Return frequency f [Hz] and damping time tau [ms] of the QNM
    '''
    T_MSUN = const.M_sun.value * const.G.value / const.c.value**3

    qnms = qnm.modes_cache(-2, l, m, n)
    tmp, _, _ = qnms(a = af)
    f     =     np.real(tmp) /( 2*np.pi) / (Mf * T_MSUN)
    gamma = abs(np.imag(tmp)) / (Mf * T_MSUN)
    tau   = 1./gamma * 1000

    return f, tau

def compute_progenitors_from_IMR(df, func = None):

    def McQ2Masses(mc, q):
        factor = mc * np.power(1. + q, 1.0/5.0)
        m1     = factor * np.power(q, -3.0/5.0)
        m2     = factor * np.power(q, +2.0/5.0)
        return m1, m2
    
    def Masses2McQ(m1, m2):
        q  = m2/m1
        mc = (m1*m2) ** (3./5.) / (m1+m2) ** (1./5.)
        return mc, q

    def MassRatio(m1, m2):
        q  = m2/m1
        return q

    def SymmetricMassRatio(m1, m2):
        # Taken from Bilby
        # https://git.ligo.org/lscsoft/bilby/-/blob/c6bcb81649b7ebf97ae6e1fd689e8712fe028eb0/bilby/gw/conversion.py#L1030
        eta = np.minimum((m1*m2) / (m1+m2) ** 2., 1./4.)
        return eta

    # FIXME: Check which spins needs to be passed.
    def ChiPrecessing(m1, m2, chi1, chi2):
        # Taken from Bilby
        # https://git.ligo.org/lscsoft/bilby/-/blob/c6bcb81649b7ebf97ae6e1fd689e8712fe028eb0/bilby/gw/conversion.py#L2082
        q = m2/m1
        chi_p = np.maximum(chi1, (4 * q + 3) / (3 * q + 4) * q * chi2)
        return chi_p

    def ChiSymmetric(m1, m2, chi1, chi2):
        # See App.C.2 of [https://arxiv.org/pdf/2001.09082]
        chi_s = (m1 * chi1 + m2 * chi2) / (m1 + m2)
        return chi_s

    def ChiAntiymmetric(m1, m2, chi1, chi2):
        # See App.C.2 of [https://arxiv.org/pdf/2001.09082]
        chi_a = (m1 * chi1 - m2 * chi2) / (m1 + m2)
        return chi_a

    if   func == 'McQ2Masses'        : df['m1'], df['m2'] = McQ2Masses(        df['mc'], df['q'])
    elif func == 'Masses2McQ'        : df['mc'], df['q']  = Masses2McQ(        df['m1'], df['m2'])
    elif func == 'MassRatio'         : df['q']            = MassRatio(         df['m1'], df['m2'])
    elif func == 'SymmetricMassRatio': df['eta']          = SymmetricMassRatio(df['m1'], df['m2'])
    elif func == 'ChiPrecessing'     : raise ValueError('Chi Precessing conversion in not yet implemented.')
    elif func == 'ChiSymmetric'      : df['chi_s']        = ChiSymmetric(      df['m1'], df['m2'], df['chi1'], df['chi2'])
    elif func == 'ChiAntiymmetric'   : df['chi_a']        = ChiAntiymmetric(   df['m1'], df['m2'], df['chi1'], df['chi2'])
    else:
        raise ValueError('The selected sample conversion does not exist. Please check the available options.')

    return df

def phase_amplitude_fits(eta,chi_p,chi_m,mode):
    '''
    Compute amplitude and phase fits from eta, chi_p and chi_m
    uses the fits of Cheung et al [https://arxiv.org/abs/2310.04489]
    '''
    delta = np.sqrt(1-4*eta)

    Amps = {(2,2,0): 4.004 + 1.349*chi_p + 0.333*chi_m - 1.325*eta**2 - 1.369*eta*chi_m + 2.622*chi_p*chi_m - 32.74*eta**2*chi_p + 4.313*eta*chi_p**2 - 25.18*eta*chi_p*chi_m + 83.37*eta**3*chi_p - 13.39*eta**2*chi_p**2 + 58.01*eta**2*chi_p*chi_m - 0.3837 *eta*chi_p**3 - 0.2075* chi_p**4,
            
            (2,2,1): 15.46 - 407*eta**2 + 55.43*eta*chi_p - 413.5*eta*chi_m + 14.82*chi_p**2 - 65.08*chi_p*chi_m + 17.99*chi_m**2 + 1731*eta**3 + 4245*eta**2*chi_m + 876.8*eta*chi_p*chi_m - 72.06*eta*chi_m**2 + 11.46*chi_p**3 + 101.2*chi_p*chi_m**2 -2.499*chi_m**3 - 10310*eta**3*chi_m - 2485*eta**2*chi_p*chi_m - 400*eta*chi_p*chi_m**2,
            
            (2,1,0): 0.9376*abs(chi_m) + delta*(6.697 - 148.3*eta - 1.035*chi_m + 1603*eta**2 - 0.96*eta*chi_p + 3.022*chi_p*chi_m - 4.27*chi_m**2 - 7388*eta**3 - 37.87*eta**2*chi_m - 15.85*eta*chi_p*chi_m + 12060* eta**4 - 13.17*eta*chi_p*chi_m**2 + 11.61*eta*chi_m**3 - 2.666*chi_p**2*chi_m**2 + 4.661*chi_m**4),
            
            (3,3,0): 0.2115*abs(chi_m) + delta*(1.82 + 0.6007*chi_p + 0.4653*chi_m + 16.49*eta**2 + 0.9369*chi_p*chi_m - 0.2701*chi_m**2 - 53.16*eta**3 - 4.201*eta**2*chi_m + 2.18*eta*chi_p**2 - 6.289*eta*chi_p*chi_m ),
            
            (3,2,0): 0.7695 - 3.308*eta - 1.446*eta*chi_p - 61.87*eta**3 + 72.14*eta**2*chi_p - 127.1*eta**2*chi_m - 2.769*eta*chi_p*chi_m + 0.3681*eta*chi_m**2 - 0.5065*chi_p*chi_m**2 + 0.5483*chi_m**3 + 293.4*eta**4 - 527.6*eta**3*chi_p + 1110*eta**3*chi_m + 11.14*eta**2*chi_p*chi_m + 2.18*eta*chi_p*chi_m**2 - 2.023*eta*chi_m**3 + 1014*eta**4*chi_p - 2407*eta**4*chi_m,
            
            (4,4,0): 0.6505 + 2.978*eta*chi_m + 0.4262*chi_p*chi_m + 106.1*eta**3 + 67.45*eta**2*chi_p - 12.08*eta**2*chi_m - 1.738*eta*chi_p*chi_m - 2041*eta**4 - 614.2*eta**3*chi_p + 5974*eta**5 + 1387*eta**4*chi_p
            
            }

    phases =   {(2,2,0): 0.,
    
                (2,2,1): 3.918 + 30.68*eta + 1.65*chi_p + 2.251*chi_m - 196.8*eta**2 - 15.94*eta*chi_p - 35.86*eta*chi_m - 0.2809*chi_p**2 - 2.797*chi_p*chi_m + 324.6*eta**3 + 32.04*eta**2*chi_p + 107*eta**2*chi_m + 11.19*eta*chi_p*chi_m - 0.2427*chi_p**3,
                
                (2,1,0): 4.282 + 2.075*eta - 0.8584*chi_p - 5.04*eta*chi_m - 1.626*chi_p*chi_m - 4.319*eta**2*chi_p + 21.01*eta**2*chi_m - 2.27*eta*chi_p**2 + 5.414*eta*chi_p*chi_m,
                
                (3,3,0): 0.08988 + 1.049*eta*chi_p + 40.79*eta**3,
                 
                (3,2,0): -32.08 + 889.7*eta - 81.88*chi_p + 93.05*chi_m - 9292*eta**2 + 1584*eta*chi_p - 1817*eta*chi_m - 0.3888*chi_m**2 + 40350*eta**3 - 9588*eta**2*chi_p + 10930*eta**2*chi_m - 6.121*eta*chi_p**2 - 60250*eta**4 + 18190*eta**3*chi_p - 20600*eta**3*chi_m,
                 
                (4,4,0): 153.6 - 6463*eta + 114700*eta**2 - 1053000*eta**3 + 5278000*eta**4 + 478.4*eta**3*chi_p - 13680000*eta**5 - 1960*eta**4*chi_p + 65.4*eta**4*chi_m + 14320000*eta**6
                 
                }

    Amp = eta*Amps[mode]
    phi = phases[mode]

    return Amp, phi

def compute_phase_amplitude_from_progenitors(df,modes):

    nsamp = len(df)
    A   = np.zeros(nsamp)
    phi = np.zeros(nsamp)

    for mode in modes:
        l, m, n = mode[0], mode[1], mode[2]

        for i in range(nsamp):

            eta,chi_p,chi_m = df.eta[i], df.chi_s[i], df.chi_a[i]
            A[i], phi[i] = phase_amplitude_fits(eta,chi_p,chi_m,mode)

        df.insert(0, 'A2{}{}{}'.format(l,m,n), A)
        df.insert(0, 'phi2{}{}{}'.format(l,m,n), phi)

    return df
        

def remove_mass_ratio_over_threshold(df):
    
    if not (set(['q']) <= set(df.columns)):
        raise ValueError('The IMR samples do not contain the mass ratio. Please make sure they are consistent.')
    
    nsamp0 = len(df)
    df = df[df['q'] > 1./6.]
    df = df.reset_index()

    if len(df) < nsamp0:    print(f'{nsamp0-len(df)} incompatible samples with NRSur fit were removed')

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

def set_keys_and_comp_pars(pars, df):
    
    keys = pd.unique(df[pars['stack-mode']])
    if pars['stack-mode'] == 'time': keys = plots.sort_times_list(keys)
    if not pars['ordering'] == []:
        if ((set(pars['ordering']) <= set(keys))) and (len(pars['ordering']) == len(keys)): keys = pars['ordering']
        else: raise ValueError('Invalid option for {stack_mode} ordering.'.format(stack_mode = pars['stack-mode']))
    if not (pars['compare'] == ''):
        comp_pars = pd.unique(df[pars['compare']])
        if not pars['compare-ordering'] == []:
            if ((set(pars['compare-ordering']) <= set(comp_pars))) and (len(pars['compare-ordering']) == len(comp_pars)): comp_pars = pars['compare-ordering']
            else: raise ValueError('Invalid option for {compare} ordering.'.format(compare = pars['compare-ordering']))
    else: comp_pars = 'a'

    return keys, comp_pars

def compute_bayes_factor(pars, df):

    keys = pd.unique(df[pars['stack-mode']])
    if pars['stack-mode'] == 'time': keys = plots.sort_times_list(keys)
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

    keys = pd.unique(df[pars['stack-mode']])
    if pars['stack-mode'] == 'time': keys = plots.sort_times_list(keys)
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
    keys, comp_pars = set_keys_and_comp_pars(pars, df)

    pars_header = '{}\t\t'.format(pars['stack-mode'])
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
                    median, low_perc, upp_perc = np.percentile(df_filt[par], np.array([50, 10, 90]))
                    low_err = median - low_perc
                    upp_err = upp_perc - median
                    pars_line += '{1:.{0}f}\t+{2:.{0}f}\t-{3:.{0}f}\t'.format(1, median, upp_err, low_err)
                f.write('{}\n'.format(pars_line))
            f.write('\n')
    print('\nMedian values of the parameters are saved in:\n{}\n'.format(output_path_pars))

    # Bayes factors, informations, median values of the network optimal SNR
    output_path = os.path.join(out_dir, 'SNR-H-BF.txt')
    if os.path.isfile(output_path): os.remove(output_path)
    pars_header  = '{}\t'.format(pars['stack-mode'])
    pars_header += '\tSNR\t\t\tH\tmaxL\tBF\t'.format(pars['stack-mode'])

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
                median, low_perc, upp_perc = np.percentile(df_filt['SNR_net_opt'].values[0], np.array([50, 10, 90]))
                low_err = median - low_perc
                upp_err = upp_perc - median
                if pars['BF-comparison']:
                    BF     = df_filt['Bayes_factor'].values[0]
                    BF_err = df_filt['Bayes_factor_error'].values[0]
                else:
                    BF     = df_filt['lnB'].values[0]
                    BF_err = df_filt['lnZ_error'].values[0]

                pars_line =  '{}\t'.format(key)
                pars_line += '{1:.{0}f}\t+{2:.{0}f}\t-{3:.{0}f}\t'.format(1, median, upp_err, low_err)
                pars_line += '{1:.{0}f}\t'.format(1, df_filt['H'].values[0])
                pars_line += '{1:.{0}f}\t'.format(0, df_filt['max_logL'].values[0])
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
        self.SampDataFrame     = pd.DataFrame(columns = pars['parameters'])
        self.PriorDataFrame    = pd.DataFrame(columns = pars['parameters'])
        self.IMRDataFrame      = pd.DataFrame(columns = pars['parameters'])
        self.EvidenceDataFrame = pd.DataFrame()
        single_evt_keys = {'event': str(), 'pipeline': str(), 'model': str(), 'submodel': str(), 'time': str(), 'GR_tag': str()}
        IMR_keys        = {'event': str(), 'pipeline': str(), 'model': str()}
        stack_mode_keys = list()

        if not (pars['include-IMR'] == 0 and pars['stack-mode'] in IMR_keys.keys()) and pars['ridgeline'] == 1:
            for file in os.listdir(dir_path):
                if not ((file == '.DS_Store') or (file == 'noise_evidences') or (file == 'ignore') or (file == 'SNR_samples') or (fnmatch.fnmatch(file, '*IMR*'))):
                    keys = file.split('_')
                    keys[-1] = keys[-1].split('.')[0]
                    for i,key in enumerate(single_evt_keys.keys()):
                        if (key == pars['stack-mode'] and (keys[i] not in stack_mode_keys)):
                            stack_mode_keys.append(keys[i])
                    
        for file in tqdm(os.listdir(dir_path), desc = 'Reading Posteriors'):
            if not ((file == '.DS_Store') or (file == 'noise_evidences') or (file == 'ignore') or (file == 'SNR_samples') or (fnmatch.fnmatch(file, '*_IMR*'))):

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

            # Case when there is one IMR analysis to compare separately.
            if  (fnmatch.fnmatch(file, '*_IMR*') and not pars['include-IMR'] == 0):

                file_path = os.path.join(dir_path, file)
                keys = file.split('_')
                keys[-1] = keys[-1].split('.')[0]
                for i,key in enumerate(IMR_keys.keys()):
                    IMR_keys[key] = keys[i]
                EventDataFrame0, _, _ = read_posteriors_event(file_path, pars, IMR_flag = True)
                if not pars['stack-mode'] == '':
                    if pars['stack-mode'] in IMR_keys.keys():
                        EventDataFrame = EventDataFrame0.assign(par = IMR_keys[pars['stack-mode']])
                        EventDataFrame.rename(columns={'par': pars['stack-mode']}, inplace = True)
                        if not pars['compare'] == '':
                            if pars['compare'] in IMR_keys.keys():
                                EventDataFrame = EventDataFrame.assign(par = IMR_keys[pars['compare']])
                            else:
                                EventDataFrame = EventDataFrame.assign(par = 'IMR')
                            EventDataFrame.rename(columns={'par': pars['compare']}, inplace = True)         
                        if not pars['IMR-posteriors']:
                            self.IMRDataFrame  = pd.concat([self.IMRDataFrame,  EventDataFrame], ignore_index=True)
                        else:
                            self.SampDataFrame = pd.concat([self.SampDataFrame, EventDataFrame], ignore_index=True)
                    else:
                        # This block is to avoid re-computing the fits multiple times.
                        if pars['ridgeline'] == 1:
                            for key in stack_mode_keys:
                                EventDataFrame = EventDataFrame0.assign(par = key)
                                EventDataFrame.rename(columns={'par': pars['stack-mode']}, inplace = True)
                                if not pars['compare'] == '':
                                    if pars['compare'] in IMR_keys.keys():
                                        EventDataFrame = EventDataFrame.assign(par = IMR_keys[pars['compare']])
                                    else:
                                        EventDataFrame = EventDataFrame.assign(par = 'IMR')
                                    EventDataFrame.rename(columns={'par': pars['compare']}, inplace = True) 
                                if not pars['IMR-posteriors']:
                                    self.IMRDataFrame  = pd.concat([self.IMRDataFrame,  EventDataFrame], ignore_index=True)
                                else:
                                    self.SampDataFrame = pd.concat([self.SampDataFrame, EventDataFrame], ignore_index=True)
                        else:
                            EventDataFrame = EventDataFrame0.assign(par = 'IMR')
                            EventDataFrame.rename(columns={'par': pars['stack-mode']}, inplace = True)
                            if not pars['compare'] == '':
                                if pars['compare'] in IMR_keys.keys():
                                    EventDataFrame = EventDataFrame.assign(par = IMR_keys[pars['compare']])
                                else:
                                    EventDataFrame = EventDataFrame.assign(par = 'IMR')
                                EventDataFrame.rename(columns={'par': pars['compare']}, inplace = True)         
                            if not pars['IMR-posteriors']:
                                self.IMRDataFrame  = pd.concat([self.IMRDataFrame,  EventDataFrame], ignore_index=True)
                            else:
                                self.SampDataFrame = pd.concat([self.SampDataFrame, EventDataFrame], ignore_index=True)
                else:
                    if not pars['compare'] == '':
                        if pars['compare'] in IMR_keys.keys():
                            EventDataFrame = EventDataFrame0.assign(par = IMR_keys[pars['compare']])
                        else:
                            EventDataFrame = EventDataFrame0.assign(par = 'IMR')
                        EventDataFrame.rename(columns={'par': pars['compare']}, inplace = True)         
                    if not pars['IMR-posteriors']:
                        self.IMRDataFrame  = pd.concat([self.IMRDataFrame,  EventDataFrame], ignore_index=True)
                    else:
                        self.SampDataFrame = pd.concat([self.SampDataFrame, EventDataFrame], ignore_index=True)
            if (not pars['compare'] == '') and pars['BF-comparison']:
                if not pars['evidence']: raise ValueError('Please activate the evidence option to compute the Bayes factor.')
                self.EvidenceDataFrame = compute_bayes_factor(pars, self.EvidenceDataFrame)

    def return_samples_dict(self):

        return self.SampDataFrame, self.PriorDataFrame, self.EvidenceDataFrame, self.IMRDataFrame


class Plots:
    '''
    Compute corner and violin plots.
    '''
    def __init__(self, pars, SampDataFrame, PriorDataFrame, EvidenceDataFrame, IMRDataFrame):
        
        if pars['corner']:
            if pars['corner-sns']: plots.corner_plots_sns(pars, SampDataFrame, PriorDataFrame, IMRDataFrame)
            else:                  plots.corner_plots(    pars, SampDataFrame, PriorDataFrame, IMRDataFrame)
        if pars['violin']:         plots.violin_plots(    pars, SampDataFrame, PriorDataFrame, IMRDataFrame, EvidenceDataFrame)
        if pars['ridgeline']:      plots.ridgeline_plots( pars, SampDataFrame, PriorDataFrame, IMRDataFrame)
        if pars['TGR-plot']:       plots.TGR_plots(       pars, SampDataFrame)
