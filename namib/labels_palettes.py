from matplotlib import rcParams
from distutils.spawn import find_executable

if find_executable('latex'): rcParams["text.usetex"] = True
rcParams["xtick.labelsize"] = 15
rcParams["ytick.labelsize"] = 15
rcParams["xtick.direction"] = "in"
rcParams["ytick.direction"] = "in"
rcParams["legend.fontsize"] = 17
rcParams["legend.frameon"]  = False
rcParams["legend.loc"]      = "best"
rcParams["axes.labelsize"]  = 20
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

        if   par == 'f_t_0':       string = '$f_{1}\ [Hz]$'
        elif par == 'tau_t_0':     string = '$\\tau_{1}\ [ms]$'
        elif par == 'logA_t_0':    string = '$lnA_{1}$'
        elif par == 'f_t_1':       string = '$f_{2}\ [Hz]$'
        elif par == 'tau_t_1':     string = '$\\tau_{2}\ [ms]$'
        elif par == 'logA_t_1':    string = '$lnA_{2}$'

        elif par == 'Mf':          string = '$M_f\ [M_{\odot}]$'
        elif par == 'af':          string = '$a_f$'
        elif par == 'A2220':       string = '$A_{220}$'
        elif par == 'A2330':       string = '$A_{330}$'
        elif par == 'f_22':        string = '$f_{22}\ [Hz]$'
        elif par == 'tau_22':      string = '$\\tau_{22}\ [ms]$'
        elif par == 'f_33':        string = '$f_{33}\ [Hz]$'
        elif par == 'tau_33':      string = '$\\tau_{33}\ [ms]$'

        elif par == 'm1':          string = '$m_1\ [M_{\odot}]$'
        elif par == 'm2':          string = '$m_2\ [M_{\odot}]$'
        elif par == 'chi1':        string = '$\\chi_1$'
        elif par == 'chi2':        string = '$\\chi_2$'
        elif par == 'cosiota':     string = '$cos\\iota$'
        elif par == 'iota':        string = '$\\iota\ [rad]$'

        elif par == 'logdistance': string = '$ln d_L\ [Mpc]$'
        elif par == 'distance':    string = '$d_L\ [Gpc]$'
        elif par == 'psi':         string = '$\\psi$'
        elif par == 'phase_22':    string = '$\\phi_{22}$'
        elif par == 'phase_33':    string = '$\\phi_{33}$'

        elif par == 'mc':          string = '$M_c\ [M_{\odot}]$'
        elif par == 'q':           string = '$q$'

        elif par == 'domega_220':  string = '$\\delta\\omega_{22}$'
        elif par == 'domega_330':  string = '$\\delta\\omega_{33}$'
        elif par == 'dtau_220':    string = '$\\delta\\tau_{22}$'
        elif par == 'ell':         string = '$l\ [km]$'

        # Cosmology
        elif par == 'H0':          string = '$H_{0}\ [km/s/Mpc]$'
        elif par == 'Om0':         string = '$\\Omega_{m_0}$'
        elif par == 'alpha':       string = '$\\alpha$'
        elif par == 'mmin':        string = '$m_{min}\ [M_{\odot}]$'
        elif par == 'mmax':        string = '$m_{max}\ [M_{\odot}]$'
        elif par == 'mu_g':        string = '$\\mu_{g}\ [M_{\odot}]$'
        elif par == 'sigma_g':     string = '$\\sigma_{g}\ [M_{\odot}]$'
        elif par == 'mu_z0':       string = '$\\mu_{z_0}\ [M_{\odot}]$'
        elif par == 'mu_z1':       string = '$\\mu_{z_1}$'
        elif par == 'sigma_z0':    string = '$\\sigma_{z_0}\ [M_{\odot}]$'
        elif par == 'sigma_z1':    string = '$\\sigma_{z_1}$'
        elif par == 'lambda_peak': string = '$\\lambda_{p}$'
        elif par == 'mix_z0':      string = '$mix_{z_0}$'
        elif par == 'mix_z1':      string = '$mix_{z_1}$'
        elif par == 'delta_m':     string = '$\\delta_{m}$'
        elif par == 'mu_q':        string = '$\\mu_{q}$'
        elif par == 'sigma_q':     string = '$\\sigma_{q}$'
        elif par == 'alpha_q':     string = '$\\alpha_{q}$'
        elif par == 'gamma':       string = '$\\gamma$'
        elif par == 'kappa':       string = '$\\kappa$'
        elif par == 'zp':          string = '$z_{p}$'
        elif par == 'R0':          string = '$R_0\ [Gpc^{-3}yr^{-1}]$'

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
        if   par == '22':        label = '$(2,2)$'
        elif par == '22-33':     label = '$(2,2),(3,3)$'
        elif par == '22-21':     label = '$(2,2),(2,1)$'
        elif par == '22-21-33':  label = '$(2,2),(2,1),(3,3)$'

        elif par == '220':       label = '$(2,2,0)$'
        elif par == '220-330':   label = '$(2,2,0),(3,3,0)$'
        elif par == '220-221':   label = '$(2,2,0),(2,2,1)$'

        elif par == 'GR':        label = '$\mathrm{GR}$'
        elif par == 'nGR':       label = '$\mathrm{nGR}$'
        elif par == 'do22':      label = '$\\delta\\omega_{22}$'
        elif par == 'do33':      label = '$\\delta\\omega_{33}$'
        elif par == 'dt22':      label = '$\\delta\\tau_{22}$'
        elif par == 'do22-dt22': label = '$\\delta\\omega_{33},\ \\delta\\tau_{22}$'

        elif par == 'EdGB':       label = '$\mathrm{EdGB}$'

        elif par == 'DS':        label = '$\mathrm{DS}$'
        elif par == 'TEOBPM':    label = '$\mathrm{TEOBPM}$'
        elif par == 'Kerr-220':  label = '$\mathrm{Kerr\ 220}$'
        elif par == 'Kerr-221':  label = '$\mathrm{Kerr\ 221}$'
        elif par == 'MMRDNP':    label = '$\mathrm{MMRDNP}$'
        elif par == 'LVK-RD':    label = '$\mathrm{LVK\ RD}$'
        elif par == 'LVK-IMR':   label = '$\mathrm{LVK\ IMR}$'

        # Cosmology
        elif par == 'beta':      label = '$\\beta\ \mathrm{function}$'
        elif par == 'MD':        label = '$\mathrm{Madau-Dickinson}$'

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