from matplotlib import rcParams
from distutils.spawn import find_executable

tex_flag = False
if find_executable('latex'):
    rcParams["text.usetex"] = True
    tex_flag = True
rcParams["xtick.labelsize"] = 14
rcParams["ytick.labelsize"] = 14
rcParams["xtick.direction"] = "in"
rcParams["ytick.direction"] = "in"
rcParams["legend.fontsize"] = 12
rcParams["legend.frameon"]  = False
rcParams["legend.loc"]      = "best"
rcParams["axes.labelsize"]  = 16
rcParams["axes.grid"]       = True
rcParams["grid.alpha"]      = 0.6
rcParams["grid.linestyle"]  = "dotted"
rcParams["lines.linewidth"] = 0.7


def labels_parameters(pars_list):

    labels      = []
    labels_dict = {}
    for par in pars_list:
        if   par == 'f_t_0':
            labels.append('$f_{1}\ [Hz]$')
            labels_dict[par] = '$f_{1}\ [Hz]$'
        elif par == 'tau_t_0':
            labels.append('$\\tau_{1}\ [ms]$')
            labels_dict[par] = '$\\tau_{1}\ [ms]$'
        elif par == 'logA_t_0':
            labels.append('$lnA_{1}$')
            labels_dict[par] = '$lnA_{1}$'
        elif par == 'f_t_1':
            labels.append('$f_{2}\ [Hz]$')
            labels_dict[par] = '$f_{2}\ [Hz]$'
        elif par == 'tau_t_1':
            labels.append('$\\tau_{2}\ [ms]$')
            labels_dict[par] = '$\\tau_{2}\ [ms]$'
        elif par == 'logA_t_1':
            labels.append('$lnA_{2}$')
            labels_dict[par] = '$lnA_{2}$'

        elif par == 'Mf':
            labels.append('$M_f\ [M_{\odot}]$')
            labels_dict[par] = '$M_f\ [M_{\odot}]$'
        elif par == 'af':
            labels.append('$a_f$')
            labels_dict[par] = '$a_f$'
        elif par == 'A2220':
            labels.append('$A_{220}$')
            labels_dict[par] = '$A_{220}$'
        elif par == 'A2330':
            labels.append('$A_{330}$')
            labels_dict[par] = '$A_{330}$'
        elif par == 'f_22':
            labels.append('$f_{22}\ [Hz]$')
            labels_dict[par] = '$f_{22}\ [Hz]$'
        elif par == 'tau_22':
            labels.append('$\\tau_{22}\ [ms]$')
            labels_dict[par] = '$\\tau_{22}\ [ms]$'
        elif par == 'f_33':
            labels.append('$f_{33}\ [Hz]$')
            labels_dict[par] = '$f_{33}\ [Hz]$'
        elif par == 'tau_33':
            labels.append('$\\tau_{33}\ [ms]$')
            labels_dict[par] = '$\\tau_{33}\ [ms]$'

        elif par == 'm1':
            labels.append('$m_1\ [M_{\odot}]$')
            labels_dict[par] = '$m_1\ [M_{\odot}]$'
        elif par == 'm2':
            labels.append('$m_2\ [M_{\odot}]$')
            labels_dict[par] = '$m_2\ [M_{\odot}]$'
        elif par == 'chi1':
            labels.append('$\\chi_1$')
            labels_dict[par] = '$\\chi_1$'
        elif par == 'chi2':
            labels.append('$\\chi_2$')
            labels_dict[par] = '$\\chi_2$'
        elif par == 'cosiota':
            labels.append('$cos\\iota$')
            labels_dict[par] = '$cos\\iota$'
        elif par == 'iota':
            labels.append('$\\iota\ [rad]$')
            labels_dict[par] = '$\\iota\ [rad]$'
        elif par == 'logdistance':
            labels.append('$ln d_L\ [Mpc]$')
            labels_dict[par] = '$ln d_L\ [Mpc]$'
        elif par == 'distance':
            labels.append('$d_L\ [Mpc]$')
            labels_dict[par] = '$d_L\ [Mpc]$'
        elif par == 'psi':
            labels.append('$\\psi$')
            labels_dict[par] = '$\\psi$'
        elif par == 'phase_22':
            labels.append('$\\phi_{22}$')
            labels_dict[par] = '$\\phi_{22}$'

        elif par == 'mc':
            labels.append('$M_c\ [M_{\odot}]$')
            labels_dict[par] = '$M_c\ [M_{\odot}]$'
        elif par == 'q':
            labels.append('$q$')
            labels_dict[par] = '$q$'

        elif par == 'domega_220':
            labels.append('$\\delta\\omega_{22}$')
            labels_dict[par] = '$\\delta\\omega_{22}$'
        elif par == 'domega_330':
            labels.append('$\\delta\\omega_{33}$')
            labels_dict[par] = '$\\delta\\omega_{33}$'
        elif par == 'ell':
            labels.append('$l\ [km]$')
            labels_dict[par] = '$l\ [km]$'

    return labels, labels_dict

def labels_parameters_evidence(par):

    label = ''
    if   par == 'BF_comparison':
        label = '$lnB$'
    elif par == 'bayes-factor':
        label = '$lnB_{s/n}$'
    elif par == 'information':
        label = '$H$'
    elif par == 'likelihood':
        label = '$H$'
    
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
            if colormap: colors = sns.color_palette(pars['palette'], as_cmap=True)
            else:        colors = sns.color_palette(pars['palette'], number_colors)
        except:
            if colormap: colors = cm.get_cmap(pars['palette'])
            else:        colors = mcp.gen_color(cmap = pars['palette'], n = number_colors)
    else:
        colors = pars['palette']
    
    if corner_plot:
        try:    colors = sns.color_palette('crest', number_colors)
        except: colors = mcp.gen_color(cmap = 'crest', n = number_colors)

    return colors