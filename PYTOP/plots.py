import matplotlib.pyplot as plt
from corner import corner
import numpy as np, os, pandas as pd
import matplotlib.lines as mlines
from statsmodels.graphics.boxplots import violinplot
from joypy import joyplot
import utils as utils, labels_palettes as lp


def time_array_from_list(input_keys):

    num_keys = len(input_keys)
    # Convert the list into a numpy array and sort it
    tmp = np.empty(num_keys)
    for i, key in enumerate(input_keys):
        tmp[i] = float(key.strip('M'))
    times = np.sort(tmp)

    switch = False
    for key in input_keys:
        if '.' in key: switch = True
    if not switch: times = times.astype(int)

    return times

def sort_times_list(input_keys):

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
        keys[i] = tmp + 'M'

    # IMPROVEME: find a cleaner way to check a generic condition
    try:
        idx = keys.index('-5M')
        keys[idx] = '-05M'
    except: pass
    try:
        idx = keys.index('0M')
        keys[idx] = '00M'
    except: pass
    try:
        idx = keys.index('5M')
        keys[idx] = '05M'
    except: pass

    return keys

def bounds_dictionary(pars):

    bounds_dict = {}
    for i,par in enumerate(pars['parameters']):
        bounds_dict[par] = pars['bounds'][i]
    
    return bounds_dict


def corner_plots(pars, SampDataFrame, PriorDataFrame):

    if not pars['compare'] == '': comp_pars = pd.unique(SampDataFrame[pars['compare']])
    else:                         comp_pars = 'a'

    keys = pd.unique(SampDataFrame[pars['stack-mode']])
    if pars['stack-mode'] == 'time': keys = sort_times_list(keys)
    if not pars['ordering'] == []:
        if ((set(pars['ordering']) <= set(keys))) and (len(pars['ordering']) == len(keys)): keys = pars['ordering']
        else: raise ValueError('Invalid option for {stack_mode} ordering.'.format(stack_mode = pars['stack-mode']))

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

        fig = plt.figure(figsize = pars['corner-settings']['figsize'])
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
                title_kwargs     = {"fontsize":12},
                use_math_text    = True,
                no_fill_contours = True,
                smooth           = pars['corner-settings']['smooth'],
            )
        fig.axes[np.shape(samp)[-1]-1].legend(*fig.axes[0].get_legend_handles_labels(), loc = 'center', frameon = False)
        for axx in fig.axes: axx.grid(visible = True)

        # Plot prior samples if required
        if pars['include-prior']:
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
                        title_kwargs     = {"fontsize":12},
                        use_math_text    = True,
                        no_fill_contours = True,
                        smooth           = pars['corner-settings']['smooth'],
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
                    title_kwargs     = {"fontsize":12},
                    use_math_text    = True,
                    no_fill_contours = True,
                    smooth           = pars['corner-settings']['smooth'],
                )
        
        if not pars['compare'] == '': filename = os.path.join(pars['plots-dir'], 'corner_{name}_{comp}.pdf'.format(name = pars['stack-mode'], comp = comp))
        else:                         filename = os.path.join(pars['plots-dir'], 'corner_{name}.pdf'.format(name = pars['stack-mode']))
        fig.savefig(filename, bbox_inches = 'tight', transparent = True)



def violin_plots(pars, SampDataFrame, PriorDataFrame, EvidenceDataFrame):

    keys = pd.unique(SampDataFrame[pars['stack-mode']])
    if pars['stack-mode'] == 'time': keys = sort_times_list(keys)
    if not pars['ordering'] == []:
        if ((set(pars['ordering']) <= set(keys))) and (len(pars['ordering']) == len(keys)): keys = pars['ordering']
        else: raise ValueError('Invalid option for {stack_mode} ordering.'.format(stack_mode = pars['stack-mode']))
    if not (pars['compare'] == ''):
        comp_pars = pd.unique(SampDataFrame[pars['compare']])
        if not pars['compare-ordering'] == []:
            if ((set(pars['compare-ordering']) <= set(comp_pars))) and (len(pars['compare-ordering']) == len(comp_pars)): comp_pars = pars['compare-ordering']
            else: raise ValueError('Invalid option for {compare} ordering.'.format(compare = pars['compare-ordering']))
    else: comp_pars = 'a'
    if pars['stack-mode'] == 'time': label_x = time_array_from_list(keys)
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
            if not pars['plot-cpnest'] == '':
                if pars['evidence-top']: params = [pars['plot-cpnest']] + params
                else:                    params = params + [pars['plot-cpnest']]

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
            if (not par == 'BF_comparison') and (not par == pars['plot-cpnest']):
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
                    if not pars['bounds'] == []: ax[pi].set_ylim(bounds_dict[par])
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
                    ax[pi].tick_params(labelsize = 8, axis = 'x')
                elif par == pars['plot-cpnest']:
                    EvidenceDataFrame['ordering'] = pd.Categorical(EvidenceDataFrame[pars['stack-mode']], categories = keys, ordered = True)
                    for c,comp in enumerate(comp_pars):
                        if   pars['plot-cpnest'] == 'bayes-factor':
                            value     = EvidenceDataFrame[EvidenceDataFrame[pars['compare']] == comp].sort_values('ordering').lnB
                            value_err = EvidenceDataFrame[EvidenceDataFrame[pars['compare']] == comp].sort_values('ordering').lnZ_error
                        elif pars['plot-cpnest'] == 'information':
                            value     = EvidenceDataFrame[EvidenceDataFrame[pars['compare']] == comp].sort_values('ordering').H
                            value_err = 0
                        elif pars['plot-cpnest'] == 'likelihood':
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
                        ax[pi].scatter(keys, value, s = 50, c = colors[ci], alpha = pars['violin-settings']['alpha'] )
                    ax[pi].set_ylabel(label_evidence)
                    ax[pi].tick_params(labelsize = 8, axis = 'x')
    [l.set_rotation(pars['violin-settings']['rotation']) for l in ax[len(params)-1].get_xticklabels()]
    if pars['stack-mode'] == 'time': plt.xlabel('$Time\ [M_{f}]$')
    if not pars['compare'] == '':
        import matplotlib.patches as mpatches
        patch = [mpatches.Patch(facecolor = colors[ci], edgecolor = 'k', alpha = pars['violin-settings']['alpha'], label = comp_pars[ci]) for ci,c in enumerate(comp_pars)]
        fig.axes[0].legend(handles = patch, loc = 2, frameon = False)
        for axx in fig.axes: axx.grid(visible = True)
    filename = os.path.join(pars['plots-dir'], 'violin_{name}.pdf'.format(name = pars['stack-mode']))
    fig.savefig(filename, bbox_inches = 'tight', transparent = True)



def ridgeline_plots(pars, SampDataFrame, PriorDataFrame):

    import pandas as pd

    keys = pd.unique(SampDataFrame[pars['stack-mode']])
    if pars['stack-mode'] == 'time': keys = sort_times_list(keys)
    if not pars['ordering'] == []:
        if ((set(pars['ordering']) <= set(keys))) and (len(pars['ordering']) == len(keys)): keys = pars['ordering']
        else: raise ValueError('Invalid option for {stack_mode} ordering.'.format(stack_mode = pars['stack-mode']))

    _, labels_dict = lp.labels_parameters(pars['parameters'])
    fig, ax = plt.subplots(len(keys), len(pars['parameters']), figsize = pars['ridgeline-settings']['figsize'])

    for pi,par in enumerate(pars['parameters']):

        if pars['bounds'] == []: range = None
        else:                    range = pars['bounds'][pi]
        if pi == 0: flag = True
        else:       flag = False

        if pars['compare'] == '':
            comp_pars = keys
            colors = lp.palettes(pars, colormap = True, number_colors = 0)
            SampDataFrame['ordering'] = pd.Categorical(SampDataFrame[pars['stack-mode']], categories = keys, ordered = True)
            subset = ax[:,pi]
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
                x_range   = range,
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

            subset = ax[:,pi]
            if pi == 0: flag = True
            else:       flag = False
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
                x_range   = range,
                ax        = subset,
                xlabels   = labels_dict[par]
            )
        ax[len(keys)-1][pi].xaxis.set_visible(True)
        ax[len(keys)-1][pi].set_xlabel(labels_dict[par])

    if pars['compare'] == '': colors = lp.palettes(pars, colormap = False, number_colors = len(comp_pars))
    import matplotlib.patches as mpatches
    patch = [mpatches.Patch(facecolor = colors[ci], edgecolor = 'k', alpha = pars['violin-settings']['alpha'], label = comp_pars[ci]) for ci,c in enumerate(comp_pars)]
    fig.axes[0].legend(handles = patch, loc = 2, frameon = False)
    if pars['stack-mode'] == 'time': ax[round(len(keys)/2)][0].set_ylabel('$Time\ [M_{f}]$')

    filename = os.path.join(pars['plots-dir'], 'ridgeline_{name}.pdf'.format(name = pars['stack-mode']))
    plt.savefig(filename, bbox_inches = 'tight', transparent = True)
