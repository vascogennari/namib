from matplotlib import rcParams
from distutils.spawn import find_executable

if find_executable('latex'): rcParams["text.usetex"] = True
rcParams["xtick.labelsize"] = 15
rcParams["ytick.labelsize"] = 15
rcParams["xtick.direction"] = "in"
rcParams["ytick.direction"] = "in"
rcParams["legend.fontsize"] = 15
rcParams["legend.frameon"]  = False
rcParams["legend.loc"]      = "best"
rcParams["axes.labelsize"]  = 15
rcParams["axes.grid"]       = True
rcParams["grid.alpha"]      = 0.6
rcParams["grid.linestyle"]  = "dotted"
rcParams["lines.linewidth"] = 0.7

def rc_labelsizes(pars):
    
    rcParams["xtick.labelsize"] = pars['label-sizes']['xtick']
    rcParams["ytick.labelsize"] = pars['label-sizes']['ytick']
    rcParams["legend.fontsize"] = pars['label-sizes']['legend']
    rcParams["axes.labelsize"]  = pars['label-sizes']['axes']

    Warning('\nLabelsize of axes is not currently read from the label-sizes dictionary. Please manually change the value at the beginning of labels_palettes.py\n')

    return 0

def labels_parameters(pars_list):

    labels      = []
    labels_dict = {}
    for par in pars_list:

        if   par == 'f_t_0':        string = '$f_{1}\ [Hz]$'
        elif par == 'tau_t_0':      string = '$\\tau_{1}\ [ms]$'
        elif par == 'logA_t_0':     string = '$lnA_{1}$'
        elif par == 'f_t_1':        string = '$f_{2}\ [Hz]$'
        elif par == 'tau_t_1':      string = '$\\tau_{2}\ [ms]$'
        elif par == 'logA_t_1':     string = '$lnA_{2}$'
        elif par == 'f_t_2':        string = '$f_{3}\ [Hz]$'
        elif par == 'tau_t_2':      string = '$\\tau_{3}\ [ms]$'

        elif par == 'Mf':           string = '$M_f\ [M_{\odot}]$'
        elif par == 'af':           string = '$a_f$'
        elif par == 'A2220':        string = '$A_{220}$'
        elif par == 'A2330':        string = '$A_{330}$'
        elif par == 'A2210':        string = '$A_{210}$'
        elif par == 'f_22':         string = '$f_{22}\ [Hz]$'
        elif par == 'tau_22':       string = '$\\tau_{22}\ [ms]$'
        elif par == 'f_33':         string = '$f_{33}\ [Hz]$'
        elif par == 'tau_33':       string = '$\\tau_{33}\ [ms]$'
        elif par == 'f_44':         string = '$f_{44}\ [Hz]$'
        elif par == 'tau_44':       string = '$\\tau_{4}\ [ms]$'

        elif par == 'm1':           string = '$m_1\ [M_{\odot}]$'
        elif par == 'm2':           string = '$m_2\ [M_{\odot}]$'
        elif par == 'chi1':         string = '$\\chi_1$'
        elif par == 'chi2':         string = '$\\chi_2$'
        elif par == 'cosiota':      string = '$cos\\iota$'
        elif par == 'iota':         string = '$\\iota\ [rad]$'

        elif par == 'logdistance':  string = '$ln d_L\ [Mpc]$'
        elif par == 'distance':     string = '$d_L\ [Gpc]$'
        elif par == 'psi':          string = '$\\psi$'
        elif par == 'phase_22':     string = '$\\phi_{22}$'
        elif par == 'phase_33':     string = '$\\phi_{33}$'

        elif par == 'mc':           string = '$M_c\ [M_{\odot}]$'
        elif par == 'q':            string = '$q$'

        elif par == 'domega_220':   string = '$\\delta\\omega_{22}$'
        elif par == 'domega_330':   string = '$\\delta\\omega_{33}$'
        elif par == 'domega_221':   string = '$\\delta\\omega_{221}$'
        elif par == 'dtau_220':     string = '$\\delta\\tau_{22}$'
        elif par == 'dtau_330':     string = '$\\delta\\tau_{33}$'
        elif par == 'domega_220':   string = '$\\delta\\omega_{22}$'
        elif par == 'ell':          string = '$l\ [km]$'

        # Cosmology
        elif par == 'H0':           string = '$H_{0}\ [km/s/Mpc]$'
        elif par == 'Om0':          string = '$\\Omega_{m_0}$'
        elif par == 'alpha':        string = '$\\alpha$'
        elif par == 'alpha_a':      string = '$\\alpha_a$'
        elif par == 'alpha_b':      string = '$\\alpha_b$'
        elif par == 'alpha_c':      string = '$\\alpha_c$'
        elif par == 'alpha_z0':     string = '$\\alpha_{z_0}$'
        elif par == 'alpha_z1':     string = '$\\alpha_{z_1}$'
        elif par == 'alpha_a_z0':   string = '$\\alpha_{a_{z_0}}$'
        elif par == 'alpha_a_z1':   string = '$\\alpha_{a_{z_1}}$'
        elif par == 'alpha_b_z0':   string = '$\\alpha_{b_{z_0}}$'
        elif par == 'alpha_b_z1':   string = '$\\alpha_{b_{z_1}}$'
        elif par == 'alpha_c_z0':   string = '$\\alpha_{c_{z_0}}$'
        elif par == 'alpha_c_z1':   string = '$\\alpha_{c_{z_1}}$'
        elif par == 'mmin':         string = '$m_{min}\ [M_{\odot}]$'
        elif par == 'mmin_z0':      string = '$m_{min}\ [M_{\odot}]$'#'$m_{min_{z_0}}\ [M_{\odot}]$'
        elif par == 'mmin_z1':      string = '$m_{min_{z_1}}\ [M_{\odot}]$'
        elif par == 'mmin_a_z0':    string = '$m_{min_{a_{z_0}}}\ [M_{\odot}]$'
        elif par == 'mmin_a_z1':    string = '$m_{min_{a_{z_1}}}\ [M_{\odot}]$'
        elif par == 'mmin_b_z0':    string = '$m_{min_{b_{z_0}}}\ [M_{\odot}]$'
        elif par == 'mmin_b_z1':    string = '$m_{min_{b_{z_1}}}\ [M_{\odot}]$'
        elif par == 'mmin_c_z0':    string = '$m_{min_{c_{z_0}}}\ [M_{\odot}]$'
        elif par == 'mmin_c_z1':    string = '$m_{min_{c_{z_1}}}\ [M_{\odot}]$'
        elif par == 'mmin_a':       string = '$m_{min_{a}}\ [M_{\odot}]$'
        elif par == 'mmin_b':       string = '$m_{min_{b}}\ [M_{\odot}]$'
        elif par == 'mmin_c':       string = '$m_{min_{c}}\ [M_{\odot}]$'
        elif par == 'mmax':         string = '$m_{max}\ [M_{\odot}]$'
        elif par == 'mmax_a_z0':    string = '$m_{max_{a_{z_0}}}\ [M_{\odot}]$'
        elif par == 'mmax_b_z0':    string = '$m_{max_{b_{z_0}}}\ [M_{\odot}]$'
        elif par == 'mmax_c_z0':    string = '$m_{max_{c_{z_0}}}\ [M_{\odot}]$'
        elif par == 'mmax_c_z1':    string = '$m_{max_{c_{z_1}}}\ [M_{\odot}]$'
        elif par == 'mu_g':         string = '$\\mu_{g}\ [M_{\odot}]$'#'$\\mu\ [M_{\odot}]$'
        elif par == 'sigma_g':      string = '$\\sigma_{g}\ [M_{\odot}]$'#'$\\sigma\ [M_{\odot}]$'
        elif par == 'mu_z0':        string = '$\\mu\ [M_{\odot}]$'#'$\\mu_{z_0}\ [M_{\odot}]$'
        elif par == 'mu_z1':        string = '$\\mu_{z_1}\ [M_{\odot}]$'
        elif par == 'sigma_z0':     string = '$\\sigma\ [M_{\odot}]$'#'$\\sigma_{z_0}\ [M_{\odot}]$'
        elif par == 'sigma_z1':     string = '$\\sigma_{z_1}\ [M_{\odot}]$'
        elif par == 'mu_a_z0':      string = '$\\mu_{a}\ [M_{\odot}]$'
        elif par == 'mu_b_z0':      string = '$\\mu_{b}\ [M_{\odot}]$'
        elif par == 'sigma_a_z0':   string = '$\\sigma_{a}\ [M_{\odot}]$'
        elif par == 'sigma_b_z0':   string = '$\\sigma_{b}\ [M_{\odot}]$'
        elif par == 'lambda_peak':  string = '$\\lambda_{p}$'
        elif par == 'mix_z0':       string = '$mix_{z_0}$'#'$mix$'
        elif par == 'mix_z1':       string = '$mix_{z_1}$'
        elif par == 'mix_alpha_z0': string = '$mix_{\\alpha_{z_0}}$'
        elif par == 'mix_alpha_z1': string = '$mix_{\\alpha_{z_1}}$'
        elif par == 'mix_beta_z0':  string = '$mix_{\\beta_{z_0}}$'
        elif par == 'mix_beta_z1':  string = '$mix_{\\beta_{z_1}}$'
        elif par == 'delta_m':      string = '$\\delta_{m}$'
        elif par == 'delta_m_a':    string = '$\\delta_{m_a}$'
        elif par == 'delta_m_b':    string = '$\\delta_{m_b}$'
        elif par == 'delta_m_c':    string = '$\\delta_{m_c}$'
        elif par == 'mu_q':         string = '$\\mu_{q}$'
        elif par == 'sigma_q':      string = '$\\sigma_{q}$'
        elif par == 'alpha_q':      string = '$\\alpha_{q}$'
        elif par == 'gamma':        string = '$\\gamma$'
        elif par == 'kappa':        string = '$\\kappa$'
        elif par == 'zp':           string = '$z_{p}$'
        elif par == 'R0':           string = '$R_0\ [Gpc^{-3}yr^{-1}]$'
        elif par == 'beta':         string = '$\\beta$'
        elif par == 'm_b':          string = '$m_{b}\ [M_{\odot}]$'
        elif par == 'delta':        string = '$\\delta$'
        elif par == 'gamma':        string = '$\\gamma$'
        elif par == 'a_gamma':      string = '$a_{q}$'
        elif par == 'theta':        string = '$\\theta_{q}$'
        elif par == 'a_johnson':    string = '$a_{\\mathrm{j}}$'
        elif par == 'b_johnson':    string = '$b_{\\mathrm{j}}$'
        elif par == 'loc_johnson':  string = '$loc_{\\mathrm{j}}$'
        elif par == 'scale_johnson':string = '$scale_{\\mathrm{j}}$'
        elif par == 'a':            string = '$a_{z}$'
        elif par == 'b':            string = '$b_{z}$'
        elif par == 'loc':          string = '$loc_{z}$'
        elif par == 'scale':        string = '$scale_{z}$'

        else:
            raise ValueError('At least one of the selected parameters does not have its corresponding label. Please, fix it in labels_palettes.py')
        
        labels.append(string)
        labels_dict[par] = string

    return labels, labels_dict

def labels_parameters_evidence(par):

    label = ''
    if   par == 'BF_comparison': label = '$lnB$'
    elif par == 'bayes-factor':  label = '$lnB_{s/n}$'
    elif par == 'information':   label = '$H$'
    elif par == 'likelihood':    label = '$H$'
    else:
        raise ValueError('The selected option for the evidence does not have its corresponding label. Please, fix it in labels_palettes.py')
    
    return label

def labels_legend(par):

    label = ''
    try:
        if   par == '1DS':       label = '$1 \mathrm{DS}$'
        elif par == '2DS':       label = '$2 \mathrm{DS}$'
        elif par == '3DS':       label = '$3 \mathrm{DS}$' 
        elif par == '1mode':     label = '$1 \mathrm{DS}$'
        elif par == '2mode':     label = '$2 \mathrm{DS}$'
        elif par == '3mode':     label = '$3 \mathrm{DS}$'

        elif par == '22':        label = '$(2,2)$'
        elif par == '22-33':     label = '$(2,2),(3,3)$'
        elif par == '22-21':     label = '$(2,2),(2,1)$'
        elif par == '22-21-33':  label = '$(2,2),(2,1),(3,3)$'
        elif par == '22-21-33-44-55':  label = '$(2,2),(2,1),(3,3),(4,4),(5,5)$'
        elif par == 'HM':        label = '$\mathrm{HM}$'

        elif par == '220':       label = '$(2,2,0)$'
        elif par == '220-330':   label = '$(2,2,0),(3,3,0)$'
        elif par == '220-221':   label = '$(2,2,0),(2,2,1)$'
        elif par == '220-210':   label = '$(2,2,0),(2,1,0)$'
        elif par == '220-200':   label = '$(2,2,0),(2,0,0)$'
        elif par == '220-320':   label = '$(2,2,0),(3,2,0)$'
        elif par == '220-440':   label = '$(2,2,0),(4,4,0)$'

        elif par == 'GR':        label = '$\mathrm{GR}$'
        elif par == 'nGR':       label = '$\mathrm{nGR}$'
        elif par == 'do22':      label = '$\\delta\\omega_{22}$'
        elif par == 'do33':      label = '$\\delta\\omega_{33}$'
        elif par == 'dt22':      label = '$\\delta\\tau_{22}$'
        elif par == 'do22-dt22': label = '$\\delta\\omega_{22},\ \\delta\\tau_{22}$'

        elif par == '22-do22':   label = '$\\delta\\omega_{22}$'
        elif par == '22-dt22':   label = '$\\delta\\tau_{22}$'
        elif par == '22-do22-dt22': label = '$(2,2),\ \{\\delta\\omega_{22},\ \\delta\\tau_{22}\}$'
        elif par == 'HM-do22-dt22':  label = '$\mathrm{HM},\ \{\\delta\\omega_{22},\ \\delta\\tau_{22}\}$'

        elif par == 'EdGB':      label = '$\mathrm{EdGB}$'

        elif par == 'DS':        label = '$\mathrm{DS}$'
        elif par == 'TEOBPM':    label = '$\mathrm{TEOBPM}$'
        elif par == 'Kerr':      label = '$\mathrm{Kerr}$'
        elif par == 'Kerr-220':  label = '$\mathrm{Kerr\ 220}$'
        elif par == 'Kerr-221':  label = '$\mathrm{Kerr\ 221}$'
        elif par == 'MMRDNP':    label = '$\mathrm{MMRDNP}$'
        elif par == 'LVK-RD':    label = '$\mathrm{LVK\ RD}$'
        elif par == 'LVK-IMR':   label = '$\mathrm{LVK\ IMR}$'
        elif par == 'IMR':       label = '$\mathrm{IMR}$'
        elif par == 'KerrBinary':     label = '$\mathrm{KerrBinary}$'
        elif par == 'KerrPostmerger': label = '$\mathrm{KerrPostmerger}$'

        elif par == 'NRSur7dq4':  label = '$\mathrm{NRSur7dq4}$'
        elif par == 'SEOBv5PHM':  label = '$\mathrm{SEOBv5PHM}$'

        # Cosmology
        elif par == 'beta':               label = '$\\beta\ \mathrm{function}$'
        elif par == 'MD':                 label = '$\mathrm{Madau-Dickinson}$'
        elif par == 'Mass2-PowerLaw':     label = '$m_2\sim\mathrm{PL}$'
        elif par == 'MassRatio-Gaussian': label = '$q\sim\mathrm{G}$'
        elif par == 'MassRatio-PowerLaw': label = '$q\sim\mathrm{PL}$'
        elif par == 'PowerLaw':           label = '$\mathrm{PL\:\:rate}$'
        elif par == 'MadauDickinson':     label = '$\mathrm{MD\:\:rate}$'

        elif par == 'PL-PL':              label = '$\mathrm{PL\ PL}$'
        elif par == 'PL-G':               label = '$\mathrm{PL\ G}$'
        elif par == 'PL-Gz':              label = '$\mathrm{PL\ G(z)}$'
        elif par == 'PLz-G':              label = '$\mathrm{PL(z)\ G}$'
        elif par == 'PLz-Gz':             label = '$\mathrm{PL(z)\ G(z)}$'
        elif par == 'PL-PL-PL':           label = '$\mathrm{PL\ PL\ PL}$'
        elif par == 'PL-PL-G':            label = '$\mathrm{PL\ PL\ G}$'
        elif par == 'PL-G-G':             label = '$\mathrm{PL\ G\ G}$'
        elif par == 'PLz-PL-PL':          label = '$\mathrm{PL(z)\ PL\ PL}$'
        elif par == 'PL-PLz-PL':          label = '$\mathrm{PL\ PL(z)\ PL}$'
        elif par == 'PL-PL-PLz':          label = '$\mathrm{PL\ PL\ PL(z)}$'
        elif par == 'PLz-PLz-PLz':        label = '$\mathrm{PL(z)\ PL(z)\ PL(z)}$'
        elif par == 'PLz-PL-PL-mixture':  label = '$\mathrm{PL(z)\ PL\ PL,\ \mathrm{mix(z)}}$'
        elif par == 'PL-PLz-PL-mixture':  label = '$\mathrm{PL\ PL(z)\ PL,\ \mathrm{mix(z)}}$'
        elif par == 'PL-PL-PLz-mmax':     label = '$\mathrm{PL\ PL\ PL(z),\ m_{max}(z)}$'
        elif par == 'PL-PL-GW190521':     label = '$\mathrm{PL\ PL,\ \mathrm{GW}190521}$'
        elif par == 'PL-PL-PL-GW190521':  label = '$\mathrm{PL\ PL\ PL,\ \mathrm{GW}190521}$'

        elif par == 'alpha8-sigma6':      label = '$\mathrm{PL\ G,\ narrow\ priors}$'
        elif par == 'alpha12-sigma6':     label = '$\mathrm{PL\ G,\ narrow\ priors}$'
        elif par == 'alpha12-sigma10':    label = '$\mathrm{PL\ G,\ medium\ priors}$'
        elif par == 'alpha12-sigma30':    label = '$\mathrm{PL\ G,\ medium\ priors}$'
        elif par == 'alpha200-sigma30':   label = '$\mathrm{PL\ G,\ wide\ priors}$'
        elif par == 'alpha100-sigma30':   label = '$\mathrm{PL\ G,\ wide\ priors}$'

        elif par == 'stationary':         label = '$\mathrm{Stationary}$'
        elif par == 'evolving':           label = '$\mathrm{Evolving}$'
        elif par == 'sharp':              label = '$\mathrm{PL\ PL\ PL}$'
        elif par == 'smoothed':           label = '$\mathrm{PL\ PL\ PL,\ Smoothing}$'

        elif par == '90yrs-Johnson':        label = '$\mathrm{JD\ (90yrs)}$'
        elif par == '90yrs-DoublePowerlaw': label = '$\mathrm{DPL\ (90yrs)}$'

        # FIXME: to remove
        elif par == 'injected_pop':       label = '$\mathrm{Injected\ population}$'
        elif par == 'data_histogram':     label = '$\mathrm{MBH\ Q3d\ catalog}$'

        else:
            raise ValueError('Unknown legend parameter.')
    except: label = par

    return label

def labels_events(par):

    label = ''
    if   par == 'GW150914':        label = '$\mathrm{GW}150914$'
    elif par == 'GW170729':        label = '$\mathrm{GW}170729$'
    elif par == 'GW190521':        label = '$\mathrm{GW}190521$'
    elif par == 'GW190521_074359': label = '$\mathrm{GW}190521\_074359$'
    elif par == 'GW191109_010717': label = '$\mathrm{GW}191109\_010717$'
    elif par == 'GW200129_065458': label = '$\mathrm{GW}200129\_065458$'
    elif par == 'GW200224_222234': label = '$\mathrm{GW}200224\_222234$'
    else:
        raise ValueError('The selected option for the event name does not have its corresponding label. Please, fix it in labels_palettes.py')
    
    return label

def labels_curves(pars_list, log10_PDF = False):

    labels_dict = {}
    for par in pars_list:

        if   par == 'PrimaryMassDistribution':
            if not log10_PDF:                                        dict = {'x': '$m_1\ [M_{\odot}]$',       'y': '$p(m_1)$'}
            else:                                                    dict = {'x': '$log_{10}(m_1)$',          'y': '$p(log_{10}(m_1))$'}
        elif par == 'PrimaryMassDistribution_NoSelectionEffects':    dict = {'x': '$m_1\ [M_{\odot}]$',       'y': '$p(m_1)$'}
        elif par == 'PrimaryMassDistribution_DetectorFrame':         dict = {'x': '$m_{1,det}\ [M_{\odot}]$', 'y': '$p(m_{1,det})$'}
        elif par == 'SecondaryMassDistribution':
            if not log10_PDF:                                        dict = {'x': '$q$',                      'y': '$p(q)$'}  # dict = {'x': '$m_2\ [M_{\odot}]$',       'y': '$p(m_2)$'}
            else:                                                    dict = {'x': '$log_{10}(q)$',            'y': '$p(log_{10}(q))$'}  # dict = {'x': '$m_2\ [M_{\odot}]$',       'y': '$p(m_2)$'}
        elif par == 'SecondaryMassDistribution_NoSelectionEffects':  dict = {'x': '$q_{det}$',                'y': '$p(q_{det})$'}  # dict = {'x': '$m_2\ [M_{\odot}]$',       'y': '$p(m_2)$'}
        elif par == 'SecondaryMassDistribution_DetectorFrame':       dict = {'x': '$m_{2,det}\ [M_{\odot}]$', 'y': '$p(m_{2,det})$'}
        elif par == 'RateEvolutionFunction':                         dict = {'x': '$z$',                      'y': '$R(z)/R_0$'}
        elif par == 'RateEvolutionDistribution_Probability':         dict = {'x': '$z$',                      'y': '$p(z)$'}
        elif par == 'RedshiftDistribution_NoSelectionEffects':       dict = {'x': '$z$',                      'y': '$p(z)$'}
        elif par == 'LuminosityDistranceDistribution_DetectorFrame': dict = {'x': '$d_L\ [Mpc]$',             'y': '$p(d_L)$'}

        else:
            raise ValueError('At least one of the selected parameters does not have its corresponding label. Please, fix it in labels_palettes.py')

        labels_dict[par] = dict

    return labels_dict

def palettes(pars, colormap, number_colors, corner_plot = False):

    # To take inspiration for choosing the palette:
    # https://r02b.github.io/seaborn_palettes/
    # https://matplotlib.org/stable/tutorials/colors/colormaps.html

    import seaborn as sns
    from matplotlib import cm
    from mycolorpy import colorlist as mcp

    if type(pars['palette']) == str:
        try:
            if colormap: colors = sns.color_palette(pars['palette'], as_cmap = True)
            else:        colors = sns.color_palette(pars['palette'], number_colors)
        except:
            if colormap: colors = cm.get_cmap(pars['palette'])
            else:        colors = mcp.gen_color(cmap = pars['palette'], n = number_colors)
    else:
        colors = pars['palette']
    
    # if corner_plot:
    #     try:    colors = sns.color_palette(pars['palette'], number_colors)
    #     except: colors = mcp.gen_color(cmap = pars['palette'], n = number_colors)

    return colors