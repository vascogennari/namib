
usage = """

\nWelcome to the namib helper!\n

    # ----- #
    # input #
    # ----- #

        samp-dir           Directory from which input samples are read. The path is relative to the 'samples' namib directory. Default: ''
        output             Directory in which the output is saved. The path is relative to the 'results' namib directory. Default: ''
        screen-output      Flag to deviate the standard output and error to file. Default: 0

        stack-mode         Name of the key used to collect the input files. The value is read directly from the file names. Options: {'event', 'pipeline', 'model', 'submodel', 'time', 'GR_tag'}
        compare            Name of the key used to compare the input files filtered on stack-mode. The value is read directly from the file names. Options: {'event', 'pipeline', 'model', 'submodel', 'time', 'GR_tag'}

        parameters         List of parameters used in the analysis. The parameters name need to be consistent with those of the samples. Default: ['m1', 'm2', 'chi1', 'chi2']
        bounds             List of parameters bounds that are used in the plots. Default: []. Syntax: [[a, b], [c, d], ...]
        ordering           List of the selected parameters ordering. Default: []
        compare-ordering   List of the compared parameters ordering. Default: []
        
        compare-hard       Flag to activate a hard comparison in case only two options are compared. This parameter is used only in violin and ridgeline plots. Default: 0
        evidence           Flag to read the evidence from the samples. Default: 0
        include-prior      Flag to plot the prior. Option currently implemented only in corner plot and not fully implemented. Default: 0
        truths             List of true values of the selected parameters. Currently implemented only with corner-sns=0. Default: []

        modes              List of modes for which the QNMs are computed. This option is used only when QNMs are computed from {Mf, af}. Default: [(2,2)]
        ds-scaling         Flag to convert the damping time in [ms] and scale amplitudes as [1e-21]. The option is used to compare samples from Damped Sinusoids with other models. Default: 0
        qnms-pyRing        Flag to use pyRing fits to compute the QNMs. Default: 1
        remnant-pyRing     Flag to use pyRing fits to compute the remnant parameters. Default: 1

        save-post          Flag to save the imput samples filtered on the selected parameters. They are saved in 'output/reduced_posteriors'. Default: 0
        save-medians       Flag to save the medians of the selected parameters. They are saved in 'output/output_medians'. Default: 0
        downsample         Option to downsample the input samples, taking the corresponding percent value (downsampling=1, 100% of initial samples). Default: 1

    # ----- #
    # plots #
    # ----- #

        corner             Flag to produce corner plot. Default: 0
        violin             Flag to produce violin plot. Default: 0
        ridgeline          Flag to produce ridgeline plot. Default: 0
        TGR-plot           Flag to produce TGR plot. Default: 0
        corner-sns         Flag to produce corner plot with seaborn. Default: 1

        corner-settings    Dictionary for additional corner settings. Options with seaborn: {'figsize', 'figname', 'alpha', 'linewidth'}. Options with corner: {'figsize', 'figname', 'smooth'}. Default: {'figsize': 8, 'figname': 'corner', 'alpha': 0.5, 'smooth': 0, 'linewidth': 1}
        violin-settings    Dictionary for additional violin settings. Options: {'figsize', 'figname', 'alpha', 'rotation', 'pad'}. Default: {'figsize': (15, 25), 'figname': 'violin', 'alpha': 0.5, 'rotation': 0, 'pad': -0.5}
        ridgeline-settings Dictionary for additional ridgeline settings. Options: {'figsize', 'figname', 'alpha', 'overlap', 'fade', 'borderaxespad'}. Default: {'figsize': (20, 10), 'figname': 'ridgeline', 'alpha': 0.5, 'overlap': 0.5, 'fade': 0, 'borderaxespad': 0.5}
        label-sizes        Dictionary to set labels size. Options: {'xtick', 'ytick', 'legend', 'axes'}. Default: {'xtick': 15, 'ytick': 15, 'legend': 17, 'axes': 17}
        palette            Option to set the colors used in the plots. If string is passed colors are read from that colormap, otherwise a list of colors needs to be passed. Syntax: 'cmap_name' or ['#AB3507', '#0771AB', ...]. Default: 'crest'
        
        extra-row          Option to include extra row with additional parameter in violin plot. Available options: ['bayes-factor', 'information', 'likelihood']. Default: ''
        BF-comparison      Flag to compute the Bayes factor between two competing compare options. The option is available only for violin plot. Default: 0
        evidence-top       Flag to insert the additional parameter as the top row. If 0, the parameter is inserted as the bottom row. This option follows both 'plot-cpnest' and 'BF-comparison'. Default: 0
        time-percentiles   Option to include a shaded region in the additional row of violin plot. This option follows both 'plot-cpnest' and 'BF-comparison'. Syntax: [-3.5, 3.0]. Default: []
        automatic-bounds   Flag to automatically set the bounds on the parameters plot from the width of the posteriors. The option is currently implemeted only for the ridgeline plot. Default: 0

        horizontal-legend  Flag to set the legend horizontally in ridgeline plot. Default: 0
        event-name         Option to add string to the top corner of the plot. The option is implemented only in violin and TGR plots. Default: ''
        remove-xticks      Flag to remove x-ticks from violin plot. Default: 0
        remove-legend      Flag to remove legend from ridgeline plot. Default: 0
        fix-dimensions     Flag to save violin plot with fixed proper size. The option is useful to combine multiple violin plots. Default: 0

        single-prior       Option to select the key from which the prior is read. Default: ''
        prior-color        Color used to plot the prior. Default: '#828F61'
        truth-color        Color used to plot the truth values. Default: 'k'

"""