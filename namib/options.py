
usage = """

\nWelcome to the namib helper!\n

    # ----- #
    # input #
    # ----- #

        samp-dir           Directory from which input samples are read. The path is relative to the 'samples' namib directory. Default: ''
        output             Directory in which the output is saved. The path is relative to the 'results' namib directory. Default: ''
        file-path          Option to pass a global path for the input samples that is non relative to the 'samples' namib directory. Default: ''

        stack-mode         Name of the key used to collect the input files. The value is read directly from the file names. Options: {'event', 'pipeline', 'model', 'submodel', 'time', 'GR_tag'}
        compare            Name of the key used to compare the input files filtered on stack-mode. The value is read directly from the file names. Options: {'event', 'pipeline', 'model', 'submodel', 'time', 'GR_tag'}

        parameters         List of parameters used in the analysis. The parameters name need to be consistent with those of the samples. Default: ['m1', 'm2', 'chi1', 'chi2']
        bounds             List of parameters bounds that are used in the plots. Default: []. Syntax: [[a, b], [c, d], ...]
        ordering           List of the selected parameters ordering. Default: []
        compare-ordering   List of the compared parameters ordering. Default: []
        
        compare-hard       Flag to activate a hard comparison in case only two options are compared. This parameter is used only in violin and ridgeline plots. Default: 0
        evidence           Flag to read the evidence from the samples. Default: 0
        include-prior      Flag to plot the prior. Option currently implemented only in corner plot and not fully implemented. Default: 0

        modes              List of modes for which the QNMs are computed. This option is used only when QNMs are computed from {Mf, af}. Default: [(2,2)]
        ds-scaling         Flag to convert the damping time in [ms] and scale amplitudes as [1e-21]. The option is used to compare samples from Damped Sinusoids with other models. Default: 0

        save-post          Flag to save the imput samples filtered on the selected parameters. They are saved in 'output/reduced_posteriors'. Default: 0
        save-medians       Flag to save the medians of the selected parameters. They are saved in 'output/output_medians'. Default: 0

    # ----- #
    # plots #
    # ----- #

        corner             Flag to produce corner plot. Default: 0
        violin             Flag to produce violin plot. Default: 0
        ridgeline          Flag to produce ridgeline plot. Default: 0
        TGR-plot           Flag to produce TGR plot. Default: 0

        corner-settings    Dictionary for additional corner settings. Options: {'figsize', 'smooth'}. Default: {'figsize': (15, 15), 'smooth': 0}
        violin-settings    Dictionary for additional violin settings. Options: {'figsize', 'alpha', 'rotation', 'pad'}. Default: {'figsize': (15, 25), 'alpha': 0.5, 'rotation': 0, 'pad': -0.5},
        ridgeline-settings Dictionary for additional ridgeline settings. Options: {'figsize', 'alpha', 'overlap', 'fade'}. Default: {'figsize': (20, 10), 'alpha': 0.5, 'overlap': 0.5, 'fade': 0},
        label-sizes        Dictionary to set labels size. Options: {'xtick', 'ytick', 'legend', 'axes'}. Default: {'xtick': 15, 'ytick': 15, 'legend': 17, 'axes': 17},
        palette            Option to set the colors used in the plots. If string is passed colors are read from that colormap, otherwise a list of colors needs to be passed. Syntax: 'cmap_name' or ['#AB3507', '#0771AB', ...]. Default: 'crest'
        
        plot-cpnest        Option to include additional parameter in violin plot. Available options: ['bayes-factor', 'information', 'likelihood']. Default: ''
        BF-comparison      Flag to compute the Bayes factor between two competing compare options. The option is available only for violin plot. Default: 0
        evidence-top       Flag to insert the additional parameter as the top row. If 0, the parameter is inserted as the bottom row. This option follows both 'plot-cpnest' and 'BF-comparison'. Default: 0
        time-percentiles   Option to include a shaded region in the additional row of violin plot. This option follows both 'plot-cpnest' and 'BF-comparison'. Syntax: [-3.5, 3.0]. Default: []

        horizontal-legend  Flag to set the legend horizontally in ridgeline plot. Default: 0
        event-name         Option to add string to the top corner of the plot. The option is implemented only in violin and TGR plots. Default: ''
        remove-xticks      Flag to remove x-ticks from violin plot. Default: 0
        remove-legend      Flag to remove legend from ridgeline plot. Default: 0
        fix-dimensions     Flag to save violin plot with fixed proper size. The option is useful to combine multiple violin plots. Default: 0

        single-prior       Option to select the key from which the prior is read. Default: ''
        prior-color        Color used to plot the prior. Default: '#828F61'

"""