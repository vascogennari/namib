from matplotlib import rcParams
from distutils.spawn import find_executable

if find_executable('latex'): rcParams["text.usetex"] = True
rcParams["xtick.direction"] = "inout"
rcParams["ytick.direction"] = "inout"
rcParams["legend.frameon"]  = False
rcParams["legend.loc"]      = "best"
rcParams["axes.grid"]       = True
rcParams["grid.alpha"]      = 0.6
rcParams["grid.linestyle"]  = "dotted"
rcParams["lines.linewidth"] = 0.7
#rcParams["axes.labelsize"]  = 20
rcParams["axes.titlepad"]   = 30.
rcParams["font.family"]     = "Computer Modern Roman"

def rc_labelsizes(pars):
    
    rcParams["xtick.labelsize"] = pars['label-sizes']['xtick']
    rcParams["ytick.labelsize"] = pars['label-sizes']['ytick']
    rcParams["legend.fontsize"] = pars['label-sizes']['legend']
    rcParams["axes.labelsize"]  = pars['label-sizes']['axes']
    rcParams["font.size"]       = pars['label-sizes']['font']

    return 0

def labels_parameters(pars):

    pars_list = pars['parameters']
    
    labels      = []
    labels_dict = {}
    for par in pars_list:

        if   par == 'f_t_0':
            if not pars['freq-log-scaling']: string = '$f_{1}\ [Hz]$'
            else:                            string = '$\\log f_{1}\ [Hz]$'
        elif par == 'tau_t_0':     string = '$\\tau_{1}\ [ms]$'
        elif par == 'logA_t_0':    string = '$\\log\ A_{1}$'
        elif par == 'f_t_1':       string = '$f_{2}\ [Hz]$'
        elif par == 'tau_t_1':     string = '$\\tau_{2}\ [ms]$'
        elif par == 'logA_t_1':    string = '$\\log\ A_{2}$'
        elif par == 'f_t_2':       string = '$f_{3}\ [Hz]$'
        elif par == 'tau_t_2':     string = '$\\tau_{3}\ [ms]$'
        elif par == 'logA_t_2':    string = '$\\log\ A_{3}$'

        elif par == 'Mf':          string = '$M_f^\\mathrm{det}\ [M_{\odot}]$'
        elif par == 'af':          string = '$a_f$'

        elif par == 'A2220':       string = '$A_{220}$'
        elif par == 'A2220_1':     string = '$A^{(1)}_{220}$'
        elif par == 'A2220_2':     string = '$A^{(2)}_{220}$'
        elif par == 'A2210':       string = '$A_{210}$'
        elif par == 'A2210_1':     string = '$A^{(1)}_{210}$'
        elif par == 'A2210_2':     string = '$A^{(2)}_{210}$'
        elif par == 'A2200':       string = '$A_{200}$'
        elif par == 'A2200_1':     string = '$A^{(1)}_{200}$'
        elif par == 'A2200_2':     string = '$A^{(2)}_{200}$'
        elif par == 'A2330':       string = '$A_{330}$'
        elif par == 'A2330_1':     string = '$A^{(1)}_{330}$'
        elif par == 'A2330_2':     string = '$A^{(2)}_{330}$'
        elif par == 'A2320':       string = '$A_{320}$'
        elif par == 'A2320_1':     string = '$A^{(1)}_{320}$'
        elif par == 'A2320_2':     string = '$A^{(2)}_{320}$'
        elif par == 'A2440':       string = '$A_{440}$'
        elif par == 'A2440_1':     string = '$A^{(1)}_{440}$'
        elif par == 'A2440_2':     string = '$A^{(2)}_{440}$'

        elif par == 'AR221':       
            if not pars['AR-log-scaling']: string = '$A^R_{221}$'
            else:                        string = '$\\log A^R_{221}$'
        elif par == 'AR210':       
            if not pars['AR-log-scaling']: string = '$A^R_{210}$'
            else:                        string = '$\\log A^R_{210}$'
        elif par == 'AR200':       
            if not pars['AR-log-scaling']: string = '$A^R_{200}$'
            else:                        string = '$\\log A^R_{200}$'
        elif par == 'AR330':       
            if not pars['AR-log-scaling']: string = '$A^R_{330}$'
            else:                        string = '$\\log A^R_{330}$'
        elif par == 'AR320':       
            if not pars['AR-log-scaling']: string = '$A^R_{320}$'
            else:                        string = '$\\log A^R_{320}$'
        elif par == 'AR440':       
            if not pars['AR-log-scaling']: string = '$A^R_{440}$'
            else:                        string = '$\\log A^R_{440}$'
        
        elif par == 'phi2220':     string = '$\\phi_{220}$'
        elif par == 'phi2220_1':   string = '$\\phi^{(1)}_{220}$'
        elif par == 'phi2220_2':   string = '$\\phi^{(2)}_{220}$'
        elif par == 'phi2210':     string = '$\\phi_{210}$'
        elif par == 'phi2210_1':   string = '$\\phi^{(1)}_{210}$'
        elif par == 'phi2210_2':   string = '$\\phi^{(2)}_{210}$'
        elif par == 'phi2200':     string = '$\\phi_{200}$'
        elif par == 'phi2200_1':   string = '$\\phi^{(1)}_{200}$'
        elif par == 'phi2200_2':   string = '$\\phi^{(2)}_{200}$'
        elif par == 'phi2330':     string = '$\\phi_{330}$'
        elif par == 'phi2330_1':   string = '$\\phi^{(1)}_{330}$'
        elif par == 'phi2330_2':   string = '$\\phi^{(2)}_{330}$'
        elif par == 'phi2320':     string = '$\\phi_{320}$'
        elif par == 'phi2320_1':   string = '$\\phi^{(1)}_{320}$'
        elif par == 'phi2320_2':   string = '$\\phi^{(2)}_{320}$'
        elif par == 'phi2440':     string = '$\\phi_{440}$'
        elif par == 'phi2440_1':   string = '$\\phi^{(1)}_{440}$'
        elif par == 'phi2440_2':   string = '$\\phi^{(2)}_{440}$'
        
        elif par == 'deltaphi221':     string = '$\\delta\\phi_{221}$'
        elif par == 'deltaphi210':     string = '$\\delta\\phi_{210}$'
        elif par == 'deltaphi200':     string = '$\\delta\\phi_{200}$'
        elif par == 'deltaphi330':     string = '$\\delta\\phi_{330}$'
        elif par == 'deltaphi320':     string = '$\\delta\\phi_{320}$'
        elif par == 'deltaphi440':     string = '$\\delta\\phi_{440}$'
        
        elif par == 'f_22':        string = '$f_{22}\ [Hz]$'
        elif par == 'tau_22':      string = '$\\tau_{22}\ [ms]$'
        elif par == 'f_220':       string = '$f_{220}\ [Hz]$'
        elif par == 'tau_220':     string = '$\\tau_{220}\ [ms]$'
        elif par == 'f_33':        string = '$f_{33}\ [Hz]$'
        elif par == 'tau_33':      string = '$\\tau_{33}\ [ms]$'
        elif par == 'f_44':        string = '$f_{44}\ [Hz]$'
        elif par == 'tau_44':      string = '$\\tau_{4}\ [ms]$'

        elif par == 'm1':          string = '$m_1\ [M_{\odot}]$'
        elif par == 'm2':          string = '$m_2\ [M_{\odot}]$'
        elif par == 'eta':         string = '$\\eta$'
        elif par == 'chi1':        string = '$\\chi_1$'
        elif par == 'chi2':        string = '$\\chi_2$'
        elif par == 'chi_s':       string = '$\\chi_s$'
        elif par == 'chi_a':       string = '$\\chi_a$'
        elif par == 'cosiota':     string = '$cos\\iota$'
        elif par == 'iota':        string = '$\\iota\ [rad]$'

        elif par == 'logdistance': string = '$ln d_L\ [Mpc]$'
        elif par == 'distance':    string = '$d_L\ [Gpc]$'
        elif par == 'psi':         string = '$\\psi$'
        elif par == 'phase_22':    string = '$\\phi_{22}$'
        elif par == 'phase_33':    string = '$\\phi_{33}$'

        elif par == 'mc':          string = '$M_c\ [M_{\odot}]$'
        elif par == 'q':           string = '$q$'
        elif par == 'eta':         string = '$\eta$'
        elif par == 'chi_s':       string = '$\chi_s$'
        elif par == 'chi_a':       string = '$\chi_a$'
        elif par == 'chi_p':       string = '$\chi_p$'

        elif par == 'domega_220':  string = '$\\delta\\omega_{220}$'
        elif par == 'domega_330':  string = '$\\delta\\omega_{330}$'
        elif par == 'domega_320':  string = '$\\delta\\omega_{320}$'
        elif par == 'domega_221':  string = '$\\delta\\omega_{221}$'
        elif par == 'dtau_220':    string = '$\\delta\\tau_{220}$'
        elif par == 'dtau_330':    string = '$\\delta\\tau_{330}$'
        elif par == 'dtau_320':    string = '$\\delta\\tau_{320}$'
        elif par == 'dtau_221':    string = '$\\delta\\tau_{221}$'
        elif par == 'ell':         string = '$l\ [km]$'

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
        if   par == '1DS':                           label = '$1 \mathrm{DS}$'
        elif par == '2DS':                           label = '$2 \mathrm{DS}$'
        elif par == '3DS':                           label = '$3 \mathrm{DS}$' 
        elif par == '1mode':                         label = '$1 \mathrm{DS}$'
        elif par == '2mode':                         label = '$2 \mathrm{DS}$'
        elif par == '3mode':                         label = '$3 \mathrm{DS}$'

        elif par == '22':                            label = '$(2,2)$'
        elif par == '22-33':                         label = '$(2,2),(3,3)$'
        elif par == '22-21':                         label = '$(2,2),(2,1)$'
        elif par == '22-44':                         label = '$(2,2),(4,4)$'
        elif par == '22-21-33':                      label = '$(2,2),(2,1),(3,3)$'
        elif par == '22-21-33-44-55':                label = '$(2,2),(2,1),(3,3),(4,4),(5,5)$'
        elif par == '22-21-33-44-55-domega220':      label = '$(2,2),(2,1),(3,3),(4,4),(5,5)+\\delta\\omega_{220}$'
        elif par == '22-21-33-44-55-dtau220':        label = '$(2,2),(2,1),(3,3),(4,4),(5,5)+\\delta\\tau_{220}$'
        elif par == '22-21-33-44-55-domega-dtau220': label = '$(2,2),(2,1),(3,3),(4,4),(5,5)+\\delta\\omega_{220}+\\delta\\tau_{220}$'
        elif par == '220':                           label = '$(2,2,0)$'
        elif par == '220-np':                        label = '$(2,2,0)-\\mathrm{n.p.}$'
        elif par == '220-330':                       label = '$(2,2,0),(3,3,0)$'
        elif par == '220-330-np':                    label = '$(2,2,0),(3,3,0)-\\mathrm{n.p.}$'
        elif par == '220-330-np-amp':                label = '$(2,2,0),(3,3,0)-\\mathrm{n.p.-amp}$'
        elif par == '220-221':                       label = '$(2,2,0),(2,2,1)$'
        elif par == '220-221-np':                    label = '$(2,2,0),(2,2,1)-\\mathrm{n.p.}$'
        elif par == '220-221-domega221':             label = '$(2,2,0),(2,2,1)+\\delta\\omega_{221}$'
        elif par == '220-221-dtau221':               label = '$(2,2,0),(2,2,1)+\\delta\\tau_{221}$'
        elif par == '220-221-domega-dtau221':        label = '$(2,2,0),(2,2,1)+\\delta\\omega_{221}+\\delta\\tau_{221}$'
        elif par == '221':                           label = '$(2,2,0),(2,2,1)$'
        elif par == '220-210':                       label = '$(2,2,0),(2,1,0)$'
        elif par == '220-210-np':                    label = '$(2,2,0),(2,1,0)-\\mathrm{n.p.}$'
        elif par == '220-210-np-amp':                label = '$(2,2,0),(2,1,0)-\\mathrm{n.p.-amp}$'
        elif par == '220-200':                       label = '$(2,2,0),(2,0,0)$'
        elif par == '220-200-np':                    label = '$(2,2,0),(2,0,0)-\\mathrm{n.p.}$'
        elif par == '220-200-np-amp':                label = '$(2,2,0),(2,0,0)-\\mathrm{n.p.-amp}$'
        elif par == '220-320':                       label = '$(2,2,0),(3,2,0)$'
        elif par == '220-320-domega320':             label = '$(2,2,0),(3,2,0)+\\delta\\omega_{320}$'
        elif par == '220-320-dtau320':               label = '$(2,2,0),(3,2,0)+\\delta\\tau_{320}$'
        elif par == '220-320-domega-dtau320':        label = '$(2,2,0),(3,2,0)+\\delta\\omega_{320}+\\delta\\tau_{320}$'
        elif par == '220-320-np':                    label = '$(2,2,0),(3,2,0)-\\mathrm{n.p.}$'
        elif par == '220-320-np-domega320':          label = '$(2,2,0),(3,2,0)-\\mathrm{n.p.}+\\delta\\omega_{320}$'
        elif par == '220-320-np-dtau320':            label = '$(2,2,0),(3,2,0)-\\mathrm{n.p.}+\\delta\\tau_{320}$'
        elif par == '220-320-np-domega-dtau320':     label = '$(2,2,0),(3,2,0)-\\mathrm{n.p.}+\\delta\\omega_{320}+\\delta\\tau_{320}$'
        elif par == '220-320-np-amp':                label = '$(2,2,0),(3,2,0)-\\mathrm{n.p.-amp}$'
        elif par == '220-440':                       label = '$(2,2,0),(4,4,0)$'
        elif par == '220-440-np':                    label = '$(2,2,0),(4,4,0)-\\mathrm{n.p.}$'
        elif par == '220-440-np-amp':                label = '$(2,2,0),(4,4,0)-\\mathrm{n.p.-amp}$'

        elif par == '220':                           label = '$(2,2,0)$'
        elif par == '220-330':                       label = '$(2,2,0),(3,3,0)$'
        elif par == '220-221':                       label = '$(2,2,0),(2,2,1)$'
        elif par == '221':                           label = '$(2,2,0),(2,2,1)$'
        elif par == '220-210':                       label = '$(2,2,0),(2,1,0)$'
        elif par == '220-200':                       label = '$(2,2,0),(2,0,0)$'
        elif par == '220-320':                       label = '$(2,2,0),(3,2,0)$'
        elif par == '220-440':                       label = '$(2,2,0),(4,4,0)$'

        elif par == 'GR':                            label = '$\mathrm{GR}$'
        elif par == 'nGR':                           label = '$\mathrm{nGR}$'
        elif par == 'do22':                          label = '$\\delta\\omega_{22}$'
        elif par == 'do33':                          label = '$\\delta\\omega_{33}$'
        elif par == 'dt22':                          label = '$\\delta\\tau_{22}$'
        elif par == 'do22-dt22':                     label = '$\\delta\\omega_{22},\ \\delta\\tau_{22}$'

        elif par == '22-do22':                       label = '$\\delta\\omega_{22}$'
        elif par == '22-dt22':                       label = '$\\delta\\tau_{22}$'
        elif par == '22-do22-dt22':                  label = '$\\delta\\omega_{22},\ \\delta\\tau_{22}$'

        elif par == 'EdGB':                          label = '$\mathrm{EdGB}$'

        elif par == 'DS':                            label = '$\mathrm{DS}$'
        elif par == 'KerrBinary':                    label = '$\mathrm{KerrBinary}$'
        elif par == 'KerrPostmerger':                label = '$\mathrm{KerrPostmerger}$'
        elif par == 'TEOB':                          label = '$\mathrm{TEOB}$'
        elif par == 'TEOBPM':                        label = '$\mathrm{TEOBPM}$'
        elif par == 'Kerr':                          label = '$\mathrm{Kerr}$'
        elif par == 'Kerr-np':                       label = '$\mathrm{Kerr-n.p.}$'
        elif par == 'Kerr-np-amp':                   label = '$\mathrm{Kerr-n.p.-amp}$'
        elif par == 'Kerr-220':                      label = '$\mathrm{Kerr\ 220}$'
        elif par == 'Kerr-221':                      label = '$\mathrm{Kerr\ 221}$'
        elif par == 'MMRDNP':                        label = '$\mathrm{MMRDNP}$'
        elif par == 'LVK-RD':                        label = '$\mathrm{LVK\ RD}$'
        elif par == 'LVK-IMR':                       label = '$\mathrm{LVK\ IMR}$'
        elif par == 'IMR':                           label = '$\mathrm{IMR}$'

        elif par == 'IMR-NRSur7dq4':                 label = '$\mathrm{IMR\ NRSur7dq4}$'
        elif par == 'IMR-combined':                  label = '$\mathrm{IMR\ combined}$'
        elif par == 'IMR-220':                       label = '$\mathrm{IMR\ }(\\ell,m,n) = (2,2,0)$'
        elif par == 'IMR-200':                       label = '$\mathrm{IMR\ }(\\ell,m,n) = (2,0,0)$'
        elif par == 'NRSur7dq4':                     label = '$\mathrm{NRSur7dq4}$'
        elif par == 'SEOBv5PHM':                     label = '$\mathrm{SEOBv5PHM}$'

        else:
            raise ValueError('Unknown legend parameter.')
    except: label = f'${par}$'

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