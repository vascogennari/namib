import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rcParams, colors
from corner import corner
import seaborn as sns
import numpy as np, os, pandas as pd
from statsmodels.graphics.boxplots import violinplot
from joypy import joyplot
import namib.utils as utils, namib.labels_palettes as lp
import matplotlib.colors as mcolors


def sort_times_list(input_keys, labels = False):

    num_keys = len(input_keys)
    # Convert the list into a numpy array and sort it
    tmp = np.empty(num_keys)
    for i, key in enumerate(input_keys):
        tmp[i] = float(key.strip('M'))
    sorted_array = np.sort(tmp)

    # Clean the array and re-add the M
    switch = False
    for key in input_keys:
        if '.' in key: switch = True
    keys = [0] * num_keys
    for i, key in enumerate(sorted_array):
        tmp = str(key)
        if not switch:
            if tmp.endswith('.0'): tmp = tmp[:-2]
        keys[i] = tmp
        if not labels: keys[i] += 'M'

    return keys

def bounds_dictionary(pars):

    bounds_dict = {}
    for i,par in enumerate(pars['parameters']):
        bounds_dict[par] = pars['bounds'][i]
    
    return bounds_dict

def convert_time_percentiles(val, locs, label_x):

    range_1 = max(locs)
    range_2 = max(label_x) - min(label_x)
    val = (val - min(label_x)) * range_1 / range_2

    return val

def hex_to_RGB(hex, alpha):

    RGB = colors.to_rgb(hex)
    RGB += (alpha,)

    return RGB

def get_sigma_bounds(df, pars, keys, comp_pars, par):

    # Set the parameter bounds as 3 times the 90% CI.
    dev = 3
    bounds_max = []
    bounds_min = []

    if pars['compare'] == '':
        for key in keys:
            samps = np.array(df[df[pars['stack-mode']] == key][par])
            samps = samps[~np.isnan(samps)]

            median = np.percentile(samps, 50)
            sigma  = np.percentile(samps, 90) - median
            bounds_max.append(median + dev * sigma)
            bounds_min.append(median - dev * sigma)
    else:
        for key in keys:
            for comp in comp_pars:
                tmp = df[df[pars['stack-mode']] == key]
                samps = np.array(tmp[tmp[pars['compare']] == comp][par])
                samps = samps[~np.isnan(samps)]

                median = np.percentile(samps, 50)
                sigma  = np.percentile(samps, 90) - median
                bounds_max.append(median + dev * sigma)
                bounds_min.append(median - dev * sigma)

    return [min(bounds_min), max(bounds_max)]
    
# --------------------------------------------------------- #

# ------------------------- #
# General plots:            #
# corner, violin, ridgeline #
# ------------------------- #

def corner_plots(pars, SampDataFrame, PriorDataFrame):

    if not pars['compare'] == '': comp_pars = pd.unique(SampDataFrame[pars['compare']])
    else:                         comp_pars = 'a'

    keys = pd.unique(SampDataFrame[pars['stack-mode']])
    if pars['stack-mode'] == 'time': keys = sort_times_list(keys)
    if not pars['ordering'] == []:
        if ((set(pars['ordering']) <= set(keys))) and (len(pars['ordering']) == len(keys)): keys = pars['ordering']
        else: raise ValueError('Invalid option for {stack_mode} ordering.'.format(stack_mode = pars['stack-mode']))

    if pars['truths'] == []: truths = None
    else:                    truths = pars['truths']

    levels = [0.5, 0.68, 0.90]
    colors = lp.palettes(pars, colormap = False, number_colors = len(keys), corner_plot = True)
    
    for comp in comp_pars:
        if not pars['compare'] == '': SampDataFrameComp = SampDataFrame[SampDataFrame[pars['compare']] == comp]
        else:                         SampDataFrameComp = SampDataFrame

        range = None
        if not pars['bounds'] == []: range = pars['bounds']
        labels, _ = lp.labels_parameters(pars['parameters'])

        flag = [(SampDataFrameComp[par] == 0).all() for par in pars['parameters']]
        if any(flag):
            params = [par for par,f in zip(pars['parameters'], flag) if f == False]
            labels = [lab for lab,f in zip(labels,             flag) if f == False]
            if not pars['bounds'] == []:
                range = [ran for ran,f in zip(range,           flag) if f == False]
        else: params = pars['parameters']

        fig = plt.figure(figsize = (pars['corner-settings']['figsize'], pars['corner-settings']['figsize']))
        for i,key in enumerate(keys):
            SampDataFrameFilt = SampDataFrameComp[SampDataFrameComp[pars['stack-mode']] == key]
            samp = np.column_stack(SampDataFrameFilt[par] for par in params)
            samp = utils.clean_empty_keys_corner(samp)  # FIXME: This is a hardfix that should be changed

            lab = key.replace('_', '\_')
            fig = corner(
                samp,
                fig              = fig,
                range            = range,
                levels           = levels,
                hist_kwargs      = {'density':True, 'label':'$\mathrm{'+lab+'}$'},
                labels           = labels,
                color            = colors[i],
                show_titles      = True,
                title_kwargs     = {"fontsize": 22},
                use_math_text    = True,
                no_fill_contours = True,
                smooth           = pars['corner-settings']['smooth'],
                truths           = truths,
                truth_color      = pars['truth-color'],
            )
        fig.axes[np.shape(samp)[-1]-1].legend(*fig.axes[0].get_legend_handles_labels(), loc = 'center', frameon = False)
        for axx in fig.axes: axx.grid(visible = True)

        # Plot prior samples if required
        if pars['include-prior']:
            Warning('The prior option is not fully implemented. Please be careful in using it.')
            if pars['single-prior'] == '':
                for i,key in enumerate(keys):
                    samp = np.column_stack(PriorDataFrame[par] for par in params)
                    fig = corner(
                        samp,
                        fig              = fig,
                        range            = range,
                        levels           = levels,
                        hist_kwargs      = {'density':True, 'label':'$\mathrm{'+lab+'}$'},
                        color            = pars['prior-color'],
                        show_titles      = False,
                        title_kwargs     = {"fontsize": 12},
                        use_math_text    = True,
                        no_fill_contours = True,
                        smooth           = pars['corner-settings']['smooth'],
                        truths           = truths,
                        truth_color      = pars['truth-color'],
                    )
            else:
                samp = np.column_stack(PriorDataFrame[PriorDataFrame[pars['stack-mode']] == pars['single-prior']][par] for par in params)
                fig = corner(
                    samp,
                    fig              = fig,
                    range            = range,
                    levels           = levels,
                    hist_kwargs      = {'density':True, 'label':'$\mathrm{'+lab+'}$'},
                    color            = pars['prior-color'],
                    show_titles      = False,
                    title_kwargs     = {"fontsize": 12},
                    use_math_text    = True,
                    no_fill_contours = True,
                    smooth           = pars['corner-settings']['smooth'],
                    truths           = truths,
                    truth_color      = pars['truth-color'],
                )
        
        utils.create_directory(pars['plots-dir'], 'PNG')
        for extension in ['pdf', 'png']:
            if extension == 'pdf': filename = os.path.join(pars['plots-dir'],        '{name}.{ext}'.format(name = pars['corner-settings']['figname'], ext = extension))
            if extension == 'png': filename = os.path.join(pars['plots-dir'], 'PNG', '{name}.{ext}'.format(name = pars['corner-settings']['figname'], ext = extension))
            fig.savefig(filename, bbox_inches = 'tight', transparent = True)



def corner_plots_sns(pars, SampDataFrame, PriorDataFrame):

    Warning('The corner_plots_sns option is not fully implemented. Please be careful in using it.')

    keys = pd.unique(SampDataFrame[pars['stack-mode']])
    if pars['stack-mode'] == 'time': keys = sort_times_list(keys)
    if not pars['ordering'] == []:
        if ((set(pars['ordering']) <= set(keys))) and (len(pars['ordering']) == len(keys)): keys = pars['ordering']
        else: raise ValueError('Invalid option for {stack_mode} ordering.'.format(stack_mode = pars['stack-mode']))

    _, labels_dict = lp.labels_parameters(pars['parameters'])

    comp_pars = keys
    colors = lp.palettes(pars, colormap = False, number_colors = len(comp_pars))
    SampDataFrame['ordering'] = pd.Categorical(SampDataFrame[pars['stack-mode']], categories = comp_pars, ordered = True)

    lp.rc_labelsizes(pars)    # Set label sizes of matplotlib RC parameters
    height = pars['corner-settings']['figsize'] / len(pars['parameters'])    # The figure size in seaborn is controlled by that of individual plots

    # Plot diagonal elements
    fig = sns.pairplot(SampDataFrame.sort_values('ordering'),
            corner    = True,
            hue       = 'ordering',
            diag_kind = 'kde',
            vars      = labels_dict,
            palette   = colors,
            height    = height,
            dropna    = 1,
            plot_kws  = dict(s = 0),
            diag_kws  = dict(alpha = pars['corner-settings']['alpha'], linewidth = pars['corner-settings']['linewidth'], common_norm = False, gridsize= 3000),
        )

    for i, var_x in enumerate(labels_dict):
        for j, var_y in enumerate(labels_dict):
            if j >= i: continue
            ax = fig.axes[i, j]
            for k, comp_par in enumerate(list(reversed(comp_pars))):
                colors_r = list(reversed(colors)) # Colors are already correctly assigned to samples both in pairplot and kdeplot, but the order is which they are plotted is different
                # Plot 2D level curves
                sns.kdeplot(
                    data        = SampDataFrame[SampDataFrame[pars['stack-mode']]==comp_par],
                    x           = var_y,
                    y           = var_x,
                    ax          = ax,
                    levels      = [0.1, 1],
                    fill        = False,
                    color       = colors_r[k],
                    common_norm = False
                )
                # Fill 2D level curves
                rgba = mcolors.to_rgba(colors_r[k], pars['corner-settings']['alpha'])
                cmap = mcolors.ListedColormap([rgba]) # Use same color for the 2D fill
                sns.kdeplot(
                    data        = SampDataFrame[SampDataFrame[pars['stack-mode']]==comp_par],
                    x           = var_y,
                    y           = var_x,
                    ax          = ax,
                    levels      = [0.1, 1],
                    fill        = True,
                    linewidth   = pars['corner-settings']['linewidth'],
                    cmap        = cmap,
                    common_norm = False
                )

    # Add legend
    fig._legend.remove()    # Remove default legend
    patch = [mpatches.Patch(facecolor = hex_to_RGB(colors[ci], pars['corner-settings']['alpha']), edgecolor = colors[ci], label = lp.labels_legend(comp_pars[ci])) for ci,c in enumerate(comp_pars)]
    fig.axes[0, 0].legend(handles = patch, loc = 'center', frameon = False, bbox_to_anchor = (len(pars['parameters'])-0.5, 0.5))

    # Add truths if required
    if not pars['truths'] == []:
        for pi,_ in enumerate(pars['parameters']):
            for qi,_ in enumerate(pars['parameters']):
                if pi == qi: # Diagonal elements
                    fig.axes[pi, qi].axvline(pars['truths'][qi], ls = '--', lw = 1, alpha = 0.8, color = pars['truth-color'])
                elif pi > qi:
                    fig.axes[pi, qi].axvline(pars['truths'][qi], ls = '--', lw = 1, alpha = 0.8, color = pars['truth-color'])
                    fig.axes[pi, qi].axhline(pars['truths'][pi], ls = '--', lw = 1, alpha = 0.8, color = pars['truth-color'])

    for pi,par in enumerate(pars['parameters']):
        # Set the bounds
        if not pars['bounds'] == []:
            fig.axes[pi, pi].set_xlim(pars['bounds'][pi])
            fig.axes[pi, pi].set_ylim(pars['bounds'][pi])

        # Remove the grid
        if pars['remove-grid']:
            for qi,_ in enumerate(pars['parameters']):
                if pi == qi or pi > qi: fig.axes[pi, qi].grid(visible = False)

        fig.axes[len(pars['parameters'])-1, pi].set_xlabel(labels_dict[par])
        if not pi==0: fig.axes[pi, 0].set_ylabel(labels_dict[par])

    utils.create_directory(pars['plots-dir'], 'PNG')
    for extension in ['pdf', 'png']:
        if extension == 'pdf': filename = os.path.join(pars['plots-dir'],        '{name}.{ext}'.format(name = pars['corner-settings']['figname'], ext = extension))
        if extension == 'png': filename = os.path.join(pars['plots-dir'], 'PNG', '{name}.{ext}'.format(name = pars['corner-settings']['figname'], ext = extension))
        fig.savefig(filename, bbox_inches = 'tight', transparent = True)



def violin_plots(pars, SampDataFrame, PriorDataFrame, EvidenceDataFrame):

    keys, comp_pars = utils.set_keys_and_comp_pars(pars, SampDataFrame)
    if pars['stack-mode'] == 'time': label_x = np.array(sort_times_list(keys, labels = True), dtype = float)
    else:                            label_x = keys

    positions = np.arange(len(keys))
    _, labels = lp.labels_parameters(pars['parameters'])
    if not pars['bounds'] == []: bounds_dict = bounds_dictionary(pars)

    if pars['compare'] == '':
        colors = lp.palettes(pars, colormap = False, number_colors = 1)
        plot_opts_L = {
            'violin_fc'      : colors[0],
            'label_fontsize' : 'small',
            'label_rotation' : pars['violin-settings']['rotation'],
        }
        plot_opts_R = plot_opts_L
        SampDataFrameComp_L = SampDataFrame
        SampDataFrameComp_R = SampDataFrame
    else:
        if not pars['compare-hard']:
            colors = lp.palettes(pars, colormap = False, number_colors = len(comp_pars))
        else:
            colors = lp.palettes(pars, colormap = False, number_colors = 2)
            plot_opts_L = {
                'violin_fc'      : colors[0],
                'label_fontsize' : 'small',
                'label_rotation' : pars['violin-settings']['rotation'],
            }
            plot_opts_R = {
                'violin_fc'      : colors[1],
                'label_fontsize' : 'small',
                'label_rotation' : pars['violin-settings']['rotation'],
            }

    params = pars['parameters']
    if not pars['compare'] == '':
        if pars['BF-comparison']:
            if pars['evidence-top']: params = ['BF_comparison'] + params
            else:                    params = params + ['BF_comparison']
        else:
            if not pars['extra-row'] == '':
                if pars['evidence-top']: params = [pars['extra-row']] + params
                else:                    params = params + [pars['extra-row']]

    if (not pars['compare'] == '') and pars['compare-hard']: comp_pars_loop = 'A'
    else:                                                    comp_pars_loop = comp_pars

    fig, ax = plt.subplots(len(params), figsize = pars['violin-settings']['figsize'], sharex = True)
    fig.subplots_adjust(wspace=0, hspace=0, top=0.95, bottom=0.18)

    for ci,comp_par in enumerate(comp_pars_loop):
        if   (not pars['compare'] == '') and (not pars['compare-hard']):
            plot_opts_L = {
                'violin_fc'      : colors[ci],
                'label_fontsize' : 'small',
                'label_rotation' : pars['violin-settings']['rotation'],
            }
            plot_opts_R = plot_opts_L
            SampDataFrameComp_L = SampDataFrame[SampDataFrame[pars['compare']] == comp_par]
            SampDataFrameComp_R = SampDataFrameComp_L
        elif (not pars['compare'] == '') and pars['compare-hard']:
            SampDataFrameComp_L = SampDataFrame[SampDataFrame[pars['compare']] == comp_pars[0]]
            SampDataFrameComp_R = SampDataFrame[SampDataFrame[pars['compare']] == comp_pars[1]]

        for pi,par in enumerate(params):
            if (not par == 'BF_comparison') and (not par == pars['extra-row']):
                SampDataFrameFilt_L = SampDataFrameComp_L[par]
                SampDataFrameFilt_R = SampDataFrameComp_R[par]
                if (not (SampDataFrameFilt_L == 0).all()) or (not (SampDataFrameFilt_R == 0).all()):
                    samp_L = [np.float_(SampDataFrameFilt_L[SampDataFrameComp_L[pars['stack-mode']] == key]) for key in keys]
                    samp_R = [np.float_(SampDataFrameFilt_R[SampDataFrameComp_R[pars['stack-mode']] == key]) for key in keys]
                    samp_L, samp_R = utils.clean_empty_keys_violin(samp_L, samp_R)  # FIXME: This is a hardfix which should be changed

                    violinplot(samp_L,
                            positions    = positions,
                            labels       = label_x,
                            show_boxplot = False,
                            side         = 'left',
                            ax           = ax[pi],
                            plot_opts    = plot_opts_L)
                    violinplot(samp_R,
                            positions    = positions,
                            labels       = label_x,
                            show_boxplot = False,
                            side         = 'right',
                            ax           = ax[pi],
                            plot_opts    = plot_opts_R)
                    ax[pi].set_ylabel(labels[par])
                    if not pars['bounds'] == []:   ax[pi].set_ylim(bounds_dict[par])
                    if not (pi == len(params )-1): ax[pi].xaxis.set_visible(False)
            else:
                if par == 'BF_comparison':
                    EvidenceDataFrame['ordering'] = pd.Categorical(EvidenceDataFrame[pars['stack-mode']], categories = keys, ordered = True)
                    value     = EvidenceDataFrame[EvidenceDataFrame[pars['compare']] == comp_pars[0]].sort_values('ordering').Bayes_factor
                    value_err = EvidenceDataFrame[EvidenceDataFrame[pars['compare']] == comp_pars[0]].sort_values('ordering').Bayes_factor_error
                    label_evidence = lp.labels_parameters_evidence(par)
                    ax[pi].errorbar(
                        keys                  ,
                        value                 ,
                        value_err             ,
                        fmt        = 'o'      ,
                        linewidth  = 1        ,
                        elinewidth = 1        ,
                        capsize    = 4        ,
                        ecolor     = 'k'      ,
                        mfc        = 'None'   ,
                        ms         = 8        ,
                        mew        = 1        ,
                        mec        = 'k'      ,
                        alpha      = pars['violin-settings']['alpha'],
                    )
                    cs = [colors[1] if b > 0 else colors[0] for b in value]
                    ax[pi].scatter(keys, value, s = 50, c = cs, alpha = pars['violin-settings']['alpha'])
                    ax[pi].set_ylabel(label_evidence)
                    if pars['remove-xticks']:  ax[pi].set_xticklabels([])
                elif par == pars['extra-row']:
                    EvidenceDataFrame['ordering'] = pd.Categorical(EvidenceDataFrame[pars['stack-mode']], categories = keys, ordered = True)
                    for c,comp in enumerate(comp_pars):
                        if   pars['extra-row'] == 'bayes-factor':
                            value     = EvidenceDataFrame[EvidenceDataFrame[pars['compare']] == comp].sort_values('ordering').lnB
                            value_err = EvidenceDataFrame[EvidenceDataFrame[pars['compare']] == comp].sort_values('ordering').lnZ_error
                        elif pars['extra-row'] == 'information':
                            value     = EvidenceDataFrame[EvidenceDataFrame[pars['compare']] == comp].sort_values('ordering').H
                            value_err = 0
                        elif pars['extra-row'] == 'likelihood':
                            raise ValueError('Maximum likelihood is not currently implemented.')
                            value     = EvidenceDataFrame[EvidenceDataFrame[pars['compare']] == comp].sort_values('ordering').maxL
                            value_err = 0
                        label_evidence = lp.labels_parameters_evidence(par)
                        ax[pi].errorbar(
                            keys                  ,
                            value                 ,
                            value_err             ,
                            fmt        = 'o'      ,
                            linewidth  = 1        ,
                            elinewidth = 1        ,
                            capsize    = 4        ,
                            ecolor     = 'k'      ,
                            mfc        = colors[ci],
                            ms         = 8        ,
                            mew        = 1        ,
                            mec        = 'k'      ,
                            alpha      = pars['violin-settings']['alpha'],
                        )
                        ax[pi].scatter(keys, value, s = 50, c = colors[ci], alpha = pars['violin-settings']['alpha'])
                    ax[pi].set_ylabel(label_evidence)
                if not pars['time-percentiles'] == []:
                    a = convert_time_percentiles(pars['time-percentiles'][0], ax[pi].get_xticks(), label_x)
                    b = convert_time_percentiles(pars['time-percentiles'][1], ax[pi].get_xticks(), label_x)
                    ax[pi].axvspan(a, b, alpha = 0.1, color = '#BA9934')

    [l.set_rotation(pars['violin-settings']['rotation']) for l in ax[len(params)-1].get_xticklabels()]
    if not pars['remove-xticks'] and pars['stack-mode'] == 'time': plt.xlabel('$Time\ [M_{f}]$')
    if not pars['compare'] == '':
        patch = [mpatches.Patch(facecolor = colors[ci], edgecolor = 'k', alpha = pars['violin-settings']['alpha'], label = lp.labels_legend(comp_pars[ci])) for ci,c in enumerate(comp_pars)]
        fig.axes[-1].legend(handles = patch, loc = 'best', frameon = False)

        # Add or remove grid
        if pars['remove-grid']:
            for axx in fig.axes: axx.grid(visible = False)
        else:
            for axx in fig.axes: axx.grid(visible = True)

    if not pars['event-name'] == '':
        if not '_' in pars['event-name']: a, b = 0.955, 0.88
        else:                             a, b = 0.930, 0.88
        fig.text(a, b, lp.labels_events(pars['event-name']), size = 22, horizontalalignment = 'center', verticalalignment = 'center', transform = fig.axes[0].transAxes)

    # Plot the truth values
    if not pars['truths'] == []:
        if not pars['evidence-top']:
            for xi, axx in enumerate(fig.axes):
                if not xi == len(params)-1: axx.axhline(pars['truths'][xi], ls = '--', lw = 1.5, alpha = 0.5, color = pars['truth-color'])
        else:
            for xi, axx in enumerate(fig.axes):
                if not xi == 0:             axx.axhline(pars['truths'][xi], ls = '--', lw = 1.5, alpha = 0.5, color = pars['truth-color'])

    utils.create_directory(pars['plots-dir'], 'PNG')
    for extension in ['pdf', 'png']:
        if extension == 'pdf': filename = os.path.join(pars['plots-dir'],        '{name}.{ext}'.format(name = pars['violin-settings']['figname'], ext = extension))
        if extension == 'png': filename = os.path.join(pars['plots-dir'], 'PNG', '{name}.{ext}'.format(name = pars['violin-settings']['figname'], ext = extension))
        if not pars['fix-dimensions']:
            fig.savefig(filename, bbox_inches = 'tight', transparent = True)
        else:
            plt.tight_layout(h_pad = pars['violin-settings']['pad'])
            fig.savefig(filename, transparent = True)



def ridgeline_plots(pars, SampDataFrame, PriorDataFrame):

    keys = pd.unique(SampDataFrame[pars['stack-mode']])
    if pars['stack-mode'] == 'time': keys = sort_times_list(keys)
    if not pars['ordering'] == []:
        if ((set(pars['ordering']) <= set(keys))) and (len(pars['ordering']) == len(keys)): keys = pars['ordering']
        else: raise ValueError('Invalid option for {stack_mode} ordering.'.format(stack_mode = pars['stack-mode']))

    _, labels_dict = lp.labels_parameters(pars['parameters'])
    if pars['stack-mode'] == 'time': label_y = np.array(sort_times_list(keys, labels = True), dtype = float)
    else:                            label_y = keys

    fig, ax = plt.subplots(len(keys), len(pars['parameters']), figsize = pars['ridgeline-settings']['figsize'])
    lp.rc_labelsizes(pars)  # Set label sizes of matplotlib RC parameters

    for pi,par in enumerate(pars['parameters']):

        if pars['bounds'] == []: bounds = None
        else:                    bounds = pars['bounds'][pi]
        if pi == 0: flag = True
        else:       flag = False

        if pars['compare'] == '':
            comp_pars = keys
            colors = lp.palettes(pars, colormap = True, number_colors = 0)
            SampDataFrame['ordering'] = pd.Categorical(SampDataFrame[pars['stack-mode']], categories = keys, ordered = True)
            if pars['stack-mode'] == 'time':
                SampDataFrame['ordering'] = SampDataFrame['ordering'].map(lambda x: x.replace('M', '$'))
                SampDataFrame['ordering'] = SampDataFrame['ordering'].map(lambda x: '$'+x)

            if ax.ndim == 1: subset = ax[pi]
            else:            subset = ax[:,pi]
            if pars['automatic-bounds']: bounds = get_sigma_bounds(SampDataFrame, pars, keys, comp_pars, par)
            joyplot(SampDataFrame.sort_values('ordering'),
                by        = 'ordering',
                column    = par,
                ylim      = 'own',
                legend    = flag,
                ylabels   = flag,
                colormap  = colors,
                alpha     = pars['ridgeline-settings']['alpha'],
                fade      = pars['ridgeline-settings']['fade'],
                overlap   = pars['ridgeline-settings']['overlap'],
                linewidth = 0.5,
                linecolor = 'k',
                x_range   = bounds,
                ax        = subset,
                xlabels   = labels_dict[par]
            )
        else:
            comp_pars = pd.unique(SampDataFrame[pars['compare']])
            if not pars['compare-ordering'] == []:
                if ((set(pars['compare-ordering']) <= set(comp_pars))) and (len(pars['compare-ordering']) == len(comp_pars)): comp_pars = pars['compare-ordering']
                else: raise ValueError('Invalid option for {compare} ordering.'.format(compare = pars['compare-ordering']))
            if not pars['compare-hard']:
                colors = lp.palettes(pars, colormap = False, number_colors = len(comp_pars))
            else:
                colors = lp.palettes(pars, colormap = False, number_colors = 2)

            for comp in comp_pars:
                try:    SampDataFrame = SampDataFrame.drop(columns={comp})
                except: pass
                SampDataFrame.insert(0, comp, np.float_(SampDataFrame[par]))
            for comp in comp_pars:
                tmp = comp_pars[:]
                tmp.remove(comp)
                for elems in tmp:
                    SampDataFrame.loc[SampDataFrame[pars['compare']] == comp, [elems]] = np.nan
            SampDataFrame['ordering'] = pd.Categorical(SampDataFrame[pars['stack-mode']], categories = keys, ordered = True)
            if pars['stack-mode'] == 'time':
                SampDataFrame['ordering'] = SampDataFrame['ordering'].map(lambda x: x.replace('M', '$'))
                SampDataFrame['ordering'] = SampDataFrame['ordering'].map(lambda x: '$'+x)

            if ax.ndim == 1: subset = ax[pi]
            else:            subset = ax[:,pi]
            if pi == 0: flag = True
            else:       flag = False
            
            if pars['automatic-bounds']: bounds = get_sigma_bounds(SampDataFrame, pars, keys, comp_pars, par)
            joyplot(SampDataFrame.sort_values('ordering'),
                by        = 'ordering',
                column    = comp_pars,
                ylim      = 'own',
                legend    = flag,
                ylabels   = flag,
                color     = colors,
                alpha     = pars['ridgeline-settings']['alpha'],
                fade      = pars['ridgeline-settings']['fade'],
                overlap   = pars['ridgeline-settings']['overlap'],
                linewidth = 0.5,
                linecolor = 'k',
                x_range   = bounds,
                ax        = subset,
                xlabels   = labels_dict[par]
            )

        if ax.ndim == 1:
            ax[pi].xaxis.set_visible(True)
            ax[pi].grid(visible = False)
            ax[pi].set_xlabel(labels_dict[par])
            ax[pi].tick_params(axis = 'both', which = 'major', labelsize = pars['label-sizes']['xtick'])
            if not pars['truths'] == []: ax[pi].axvline(pars['truths'][pi], ls = '--', lw = 1.5, alpha = 0.5, color = pars['truth-color'])  # Plot truth values
        else:
            ax[len(keys)-1][pi].xaxis.set_visible(True)
            ax[len(keys)-1][pi].grid(visible = False)
            ax[len(keys)-1][pi].set_xlabel(labels_dict[par])
            for ni in range(len(keys)):
                ax[ni][pi].tick_params(axis = 'both', which = 'major', labelsize = pars['label-sizes']['xtick'])
                if not pars['truths'] == []: ax[ni][pi].axvline(pars['truths'][pi], ls = '--', lw = 1.5, alpha = 0.5, color = pars['truth-color'])  # Plot truth values

    if pars['compare'] == '': colors = lp.palettes(pars, colormap = False, number_colors = len(comp_pars))
    patch = [mpatches.Patch(facecolor = colors[ci], edgecolor = 'k', alpha = pars['ridgeline-settings']['alpha'], label = lp.labels_legend(comp_pars[ci])) for ci,c in enumerate(comp_pars)]

    # Set multiple labels in columns if required
    if not pars['horizontal-legend']: ncol = 1
    else:                             ncol = len(comp_pars)
    fig.axes[0].legend(handles = patch, loc = 2, frameon = False, ncol = ncol, borderaxespad = pars['ridgeline-settings']['borderaxespad'])

    if pars['remove-legend']: fig.axes[0].get_legend().remove()
    if pars['stack-mode'] == 'time' and ax.ndim != 1:
        if ax.ndim == 1: ax[0].set_ylabel('$Time\ [M_{f}]$')
        else:            ax[round(len(keys)/2)][0].set_ylabel('$Time\ [M_{f}]$')

    utils.create_directory(pars['plots-dir'], 'PNG')
    for extension in ['pdf', 'png']:
        if extension == 'pdf': filename = os.path.join(pars['plots-dir'],        '{name}.{ext}'.format(name = pars['ridgeline-settings']['figname'], ext = extension))
        if extension == 'png': filename = os.path.join(pars['plots-dir'], 'PNG', '{name}.{ext}'.format(name = pars['ridgeline-settings']['figname'], ext = extension))
        plt.savefig(filename, bbox_inches = 'tight', transparent = True)



# --------------------- #
# Ringdown no-hair test #
# --------------------- #

def TGR_plots(pars, SampDataFrame):
    '''
        Similar plots in the literature can be found in:
        0309007, 1107.0854, 1111.5819, 2301.02267
    '''

    from matplotlib.lines import Line2D
    from scipy.optimize import newton
    from corner import hist2d
    from tqdm import tqdm
    import pyRing.waveform as wf

    def compute_QNMs(Mf, af, l, m, par):
        if af > 1. or af < 0.: return -np.inf
        if   par=='omg': res = wf.QNM_fit(l, m, 0).f(Mf, af)          #[Hz]
        elif par=='tau': res = wf.QNM_fit(l, m, 0).tau(Mf, af)*1000   #[ms]
        return res

    def find_level_curve_mass(level, l, m, par, spins):
        curve  = np.zeros(len(spins))
        for i, af in enumerate(spins):
            curve[i] = find_mass(level, af, l, m, par)
        return curve

    def find_mass(level, af, l, m, par):
        if   par=='omg': estimate = 10
        elif par=='tau': estimate = 10
        def objective(Mf, level, af, l, m, par):
            return level - compute_QNMs(Mf, af, l, m, par)
        return newton(objective, estimate, args=(level, af, l, m, par))

    af_range = [0, 1]
    num = 100

    omg_22 = SampDataFrame['f_22'][SampDataFrame[pars['stack-mode']]   == 'nGR'][1::num]
    omg_33 = SampDataFrame['f_33'][SampDataFrame[pars['stack-mode']]   == 'nGR'][1::num]
    tau_22 = SampDataFrame['tau_22'][SampDataFrame[pars['stack-mode']] == 'nGR'][1::num]

    spins  = np.linspace(af_range[0], af_range[1], 200)[:-1]

    print('\nComputing level curves for TGR plot:')
    masses_omg22 = np.array([find_level_curve_mass(freq, 2, 2, 'omg', spins) for freq in tqdm(omg_22, desc = 'freq 22')])
    masses_omg33 = np.array([find_level_curve_mass(freq, 3, 3, 'omg', spins) for freq in tqdm(omg_33, desc = 'freq 33')])
    masses_tau22 = np.array([find_level_curve_mass(freq, 2, 2, 'tau', spins) for freq in tqdm(tau_22, desc = 'tau  22')])

    percentiles = [pars['percentiles']['m'], pars['percentiles']['ll'], pars['percentiles']['l'], pars['percentiles']['h'], pars['percentiles']['hh']]
    p_f22 = {}
    p_f33 = {}
    p_t22 = {}
    for perc in percentiles:
        p_f22[perc] = np.percentile(masses_omg22, perc, axis = 0)
        p_f33[perc] = np.percentile(masses_omg33, perc, axis = 0)
        p_t22[perc] = np.percentile(masses_tau22, perc, axis = 0)
    
    samp_Mf_nGR = np.array(SampDataFrame['Mf'][SampDataFrame[pars['stack-mode']] == 'nGR'])
    samp_af_nGR = np.array(SampDataFrame['af'][SampDataFrame[pars['stack-mode']] == 'nGR'])
    samp_Mf     = np.array(SampDataFrame['Mf'][SampDataFrame[pars['stack-mode']] == 'GR'])
    samp_af     = np.array(SampDataFrame['af'][SampDataFrame[pars['stack-mode']] == 'GR'])

    fig, ax = plt.subplots()

    hist2d(samp_Mf,
           samp_af, 
           ax               = ax, 
           levels           = [0.5,0.9], 
           plot_datapoints  = False, 
           plot_density     = False, 
           no_fill_contours = False, 
           contour_kwargs   = {'linewidths': 0.6, 'colors': 'k'})
    hist2d(samp_Mf_nGR, 
           samp_af_nGR, 
           ax               = ax, 
           levels           = [0.5,0.9], 
           plot_datapoints  = False, 
           plot_density     = False, 
           no_fill_contours = False, 
           contour_kwargs   = {'linewidths': 0.6, 'colors': 'k', 'linestyles': 'dashed', 'alpha': 0.5})
    # CR
    # do (3,3)
    ax.fill_betweenx(spins, p_f33[pars['percentiles']['ll']], p_f33[pars['percentiles']['hh']], color = '#9C3F5C', alpha = 0.25)
    ax.fill_betweenx(spins, p_f33[pars['percentiles']['l' ]], p_f33[pars['percentiles']['h' ]], color = '#9C3F5C', alpha = 0.5)
    ax.plot(p_f33[50], spins, lw = 0.7, color = '#9C3F5C', label = '$\\omega_{33}$')

    # dt (2,2)
    ax.fill_betweenx(spins, p_t22[pars['percentiles']['ll']], p_t22[pars['percentiles']['hh']], color = '#E0AA07', alpha = 0.25)
    ax.fill_betweenx(spins, p_t22[pars['percentiles']['l' ]], p_t22[pars['percentiles']['h' ]], color = '#E0AA07', alpha = 0.5)
    ax.plot(p_t22[50], spins, lw = 0.7, color = '#E0AA07', label = '$\\tau_{22}$')

    # do (2,2)
    ax.fill_betweenx(spins, p_f22[pars['percentiles']['ll']], p_f22[pars['percentiles']['hh']], color = '#608F3A', alpha = 0.25)
    ax.fill_betweenx(spins, p_f22[pars['percentiles']['l' ]], p_f22[pars['percentiles']['h' ]], color = '#608F3A', alpha = 0.5)
    ax.plot(p_f22[50], spins, lw = 0.7, color = '#608F3A', label = '$\\omega_{22}$')

    ax.set_xlabel('$M_f\ [M_\\odot]$')
    ax.set_ylabel('$a_f$')
    ax.set_ylim(af_range)
    ax.grid(True, dashes=(1,3))

    handles, _     = ax.get_legend_handles_labels()
    new_object     = Line2D([0],[0], color = 'k', lw = 0.6, label = '$\mathrm{GR}$')
    new_object_nGR = Line2D([0],[0], color = 'k', lw = 0.6, label = '$\mathrm{nGR}$', ls = 'dashed')
    handles.append(new_object)
    handles.append(new_object_nGR)
    ax.legend(handles = handles, loc = 'lower right', frameon = False)
    if not pars['event-name'] == '':
        if not '_' in pars['event-name']: a, b = 0.885, 0.945
        else:                             a, b = 0.825, 0.945
        fig.text(a, b, lp.labels_events(pars['event-name']), size = 14, horizontalalignment = 'center', verticalalignment = 'center', transform = fig.axes[0].transAxes)

    filename = os.path.join(pars['plots-dir'], 'TGR_plot.pdf')
    fig.savefig(filename, bbox_inches = 'tight', transparent = True)



# ------------------------ #
# Population distributions #
# ------------------------ #

def reconstructed_distributions(pars, CurvesDictionary):
    '''
        Plot reconstructed population distributions that are not redshift dependent.
        The file structure should be compatible with the ICAROGW pipeline output.
    '''
    keys = CurvesDictionary.keys()
    if not pars['ordering'] == []:
            if ((set(pars['ordering']) <= set(keys))) and (len(pars['ordering']) == len(keys)): keys = pars['ordering']
            else: raise ValueError('Invalid option for {stack_mode} ordering.'.format(stack_mode = pars['stack-mode']))

    if not pars['bounds'] == []: bounds_dict = bounds_dictionary(pars)
    params = pars['parameters']
    labels = lp.labels_curves(params)
    colors = lp.palettes(pars, colormap = True, number_colors = len(keys))
    fig, ax = plt.subplots(len(pars['parameters']), figsize = (pars['corner-settings']['figsize'], pars['corner-settings']['figsize']))

    # Load the injected population.
    if not pars['injected-pop'] == '':
        try:
            injected_pop = np.loadtxt(pars['injected-pop'], unpack = True)
            m_array, pdf_array = injected_pop[0], injected_pop[1]
        except:
            raise ValueError('Injected population file not found: {}.'.format(pars['injected-pop']))

    # Loop on the models.
    for ki,key in enumerate(keys):
        CurvesDataFrame = CurvesDictionary[key]
        # Loop on the parameters.
        for pi,par in enumerate(params):
            CurvesDataFrameFilt = CurvesDataFrame[par]

            # This case is treated separately because the primary mass distribution is always stored on multiple redshift slices.
            # Since here we plot stationary distributions, we take the first redshift slice.
            if par == 'PrimaryMassDistribution':
                ax[pi].fill_between(CurvesDataFrameFilt['x'], CurvesDataFrameFilt['y'][0][pars['percentiles']['ll']], CurvesDataFrameFilt['y'][0][pars['percentiles']['hh']], color = colors[ki], alpha = 0.1)
                ax[pi].fill_between(CurvesDataFrameFilt['x'], CurvesDataFrameFilt['y'][0][pars['percentiles']['l' ]], CurvesDataFrameFilt['y'][0][pars['percentiles']['h' ]], color = colors[ki], alpha = 0.5)
                ax[pi].plot(        CurvesDataFrameFilt['x'], CurvesDataFrameFilt['y'][0][pars['percentiles']['m' ]],                                                         color = colors[ki], lw    = 0.7)
            else:
                ax[pi].fill_between(CurvesDataFrameFilt['x'], CurvesDataFrameFilt['y'][   pars['percentiles']['ll']], CurvesDataFrameFilt['y'][   pars['percentiles']['hh']], color = colors[ki], alpha = 0.1)
                ax[pi].fill_between(CurvesDataFrameFilt['x'], CurvesDataFrameFilt['y'][   pars['percentiles']['l' ]], CurvesDataFrameFilt['y'][   pars['percentiles']['h' ]], color = colors[ki], alpha = 0.5)
                ax[pi].plot(        CurvesDataFrameFilt['x'], CurvesDataFrameFilt['y'][   pars['percentiles']['m' ]],                                                         color = colors[ki], lw    = 0.7)
            
            # Plot the injected population.
            if par == 'PrimaryMassDistribution':
                if not pars['injected-pop'] == '': ax[pi].plot(m_array, pdf_array, alpha = 0.5, lw = 1.2, color = pars['truth-color'], ls = '--', label = 'Injected population')
            
            ax[pi].set_xlabel(labels[par]['x'])
            ax[pi].set_ylabel(labels[par]['y'])
            if not pars['bounds'] == []: ax[pi].set_xlim(bounds_dict[par])

            # Set the log scale for the primary distributions.
            if 'PrimaryMass' in par:
                ax[pi].set_yscale('log')
                ax[pi].set_ylim(1e-5, 2 )
                if 'DetectorFrame' in par: ax[pi].set_ylim(0.5e-3, 0.1)

    if pars['injected-pop'] == '':
        patch = [mpatches.Patch(facecolor = colors[ci], edgecolor = 'k', alpha = 0.5, label = lp.labels_legend(c)) for ci,c in enumerate(keys)]
    else:
        from matplotlib.lines import Line2D
        patch = [mpatches.Patch(facecolor = colors[2], edgecolor = 'k', alpha = 0.5, label = lp.labels_legend('alpha100-sigma30')), mpatches.Patch(facecolor = colors[1], edgecolor = 'k', alpha = 0.5, label = lp.labels_legend('alpha12-sigma10')), mpatches.Patch(facecolor = colors[0], edgecolor = 'k', alpha = 0.5, label = lp.labels_legend('alpha8-sigma6')), Line2D([0], [0], color = pars['truth-color'], linestyle = '--', linewidth = 1.2, alpha = 0.8, label = lp.labels_legend('injected_pop'))]
    ax[0].legend(handles = patch, loc = 'best', frameon = False)

    if pars['remove-grid']:
        for axx in fig.axes: axx.grid(visible = False)

    plt.tight_layout()
    filename = os.path.join(pars['plots-dir'], 'reconstructed_{name}.pdf'.format(name = pars['stack-mode']))
    plt.savefig(filename, bbox_inches = 'tight', transparent = True)


 
def reconstructed_primary_redshift(pars, CurvesDictionary):
    '''
        Plot reconstructed astropysical and observed primary distributions for redshift dependent populations.
        The file structure should be compatible with the ICAROGW pipeline output.
    '''
    keys = CurvesDictionary.keys()
    if not pars['ordering'] == []:
            if ((set(pars['ordering']) <= set(keys))) and (len(pars['ordering']) == len(keys)): keys = pars['ordering']
            else: raise ValueError('Invalid option for {stack_mode} ordering.'.format(stack_mode = pars['stack-mode']))

    # Check that the distributions are the correct ones.
    primary_distributions = ['PrimaryMassDistribution', 'PrimaryMassDistribution_NoSelectionEffects']
    if not pars['parameters'] == primary_distributions: raise ValueError('Invalid distributions for the redshift dependent populations. Only primary distributions are allowed.')

    if not pars['bounds'] == []: bounds_dict = bounds_dictionary(pars)
    colors = lp.palettes(pars, colormap = True, number_colors = len(keys))
    nz = len(CurvesDictionary[[*CurvesDictionary.keys()][0]]['PrimaryMassDistribution']['y'])

    # Plot the observed distribution in log scale.
    if not pars['obs-primary-nolog']:
        fig, ax = plt.subplots(         nz, 2, figsize = (pars['corner-settings']['figsize'], 2.2 * nz), layout = 'constrained')
    else:
        subplots = [[ni, 'a'] for ni in range(nz)]
        fig, ax = plt.subplot_mosaic(subplots, figsize = (pars['corner-settings']['figsize'], 2.1 * nz), layout = 'constrained')
        
    # Loop on the models.
    for ki,key in enumerate(keys):
        CurvesDataFrame = CurvesDictionary[key]
        # Loop on the parameters.
        for i,par in enumerate(primary_distributions):
            x = CurvesDataFrame[par]['x']
            # Loop on the redshift slices.
            for zi,z in enumerate(CurvesDataFrame[par]['z']):

                zi_inv = nz-1 - zi
                if not pars['obs-primary-nolog']: ax_use = ax[zi_inv][i]
                else:                             ax_use = ax[zi_inv]
                
                if ki == 0: label = '$z={}$'.format(round(z, 2))
                else:       label = ''
                if par == 'PrimaryMassDistribution' or (par == 'PrimaryMassDistribution_NoSelectionEffects' and not pars['obs-primary-nolog']):
                    ax_use.fill_between(x, CurvesDataFrame[par]['y'][zi][pars['percentiles']['ll']],   CurvesDataFrame[par]['y'][zi][pars['percentiles']['hh']],   color = colors[ki], alpha = 0.1)
                    ax_use.fill_between(x, CurvesDataFrame[par]['y'][zi][pars['percentiles']['l' ]],   CurvesDataFrame[par]['y'][zi][pars['percentiles']['h' ]],   color = colors[ki], alpha = 0.5)
                    ax_use.plot(        x, CurvesDataFrame[par]['y'][zi][pars['percentiles']['m' ]],                                                               color = colors[ki], lw = 0.7, label = label)
                else:
                    ax['a'].fill_between(x, CurvesDataFrame[par]['y'][zi][pars['percentiles']['ll']]+z, CurvesDataFrame[par]['y'][zi][pars['percentiles']['hh']]+z, color = colors[ki], alpha = 0.1)
                    ax['a'].fill_between(x, CurvesDataFrame[par]['y'][zi][pars['percentiles']['l' ]]+z, CurvesDataFrame[par]['y'][zi][pars['percentiles']['h' ]]+z, color = colors[ki], alpha = 0.5)
                    ax['a'].plot(        x, CurvesDataFrame[par]['y'][zi][pars['percentiles']['m' ]]+z,                                                             color = colors[ki], lw = 0.7, label = label)

                if not pars['obs-primary-nolog']:
                    if zi_inv == 2: ax[zi_inv][0].set_ylabel('$p(m_1|z)$')
                    else:           ax[zi_inv][0].set_ylabel('')
                    ax[nz-1  ][i].set_xlabel('$m_1\ [M_{\odot}]$')
                else:
                    if zi_inv == 2: ax_use.set_ylabel('$p(m_1|z)$')
                    else:           ax_use.set_ylabel('')
                    ax[nz-1     ].set_xlabel('$m_1\ [M_{\odot}]$')

                ax_use.set_xlim(bounds_dict[par])
                ax_use.set_yscale('log')
                # FIXME: These values should not be hardcoded.
                ax_use.set_ylim(1e-6, 7)
                ax_use.set_yticks([1e-1, 1e-3, 1e-5])
                ax_use.set_xticks([10, 30, 50, 70])
                ax_use.legend(loc = 'best', handletextpad = -2.0, handlelength = 0)  # Revome colors from log plots.
                if not zi_inv == nz-1: plt.setp(ax_use.get_xticklabels(), visible = False)  # Remove x_ticks from log plots.
                if par == 'PrimaryMassDistribution_NoSelectionEffects': ax_use.set_yticks([])

    if pars['obs-primary-nolog']:
        ax['a'].set_xlabel('$m_1\ [M_{\odot}]$')
        ax['a'].set_ylabel('$z$')
        ax['a'].set_xlim(1, 70)
        ax['a'].set_ylim(-0.05, 1.1)
        ax['a'].set_yscale('linear')
        ax['a'].yaxis.tick_right()
        ax['a'].yaxis.set_label_position('right')

    patch = [mpatches.Patch(facecolor = colors[ci], edgecolor = 'k', alpha = 0.5, label = lp.labels_legend(c)) for ci,c in enumerate(keys)]
    if not pars['obs-primary-nolog']: ax[0][1].legend(handles = patch, loc = 'best', frameon = False)
    else:                             ax['a' ].legend(handles = patch, loc = 'best', frameon = False)

    if pars['remove-grid']:
        for axx in fig.axes: axx.grid(visible = False)
    fig.set_constrained_layout_pads(hspace = 0.0, wspace = 0.0, h_pad = 0.0, w_pad = 0.0)

    filename = os.path.join(pars['plots-dir'], 'reconstructed_{name}.pdf'.format(name = pars['stack-mode']))
    plt.savefig(filename, bbox_inches = 'tight', transparent = True)
