import os

ringdown  = 0
cosmology = 1

if  ringdown:
    # ------------------------- #
    input_path  = '/Users/vascogennari/Documents/work/code/python/results/GWTC-3_TEOBPM_TGR/PROD4'
    output_path = '/Users/vascogennari/Documents/work/code/python/namib/samples/GWTC-3_TEOBPM/GWTC-3_TEOBPM_TGR'

    elements = {
        'name'     : '',
        'pipeline' : 'pyRing',
        'model'    : '',
        'submodel' : '',
        'time'     : '',
        'GR_tag'   : 'nGR',
    }

    sampler    = 'raynest'  # Options: [raynest, cpnest]
    deviations = 1
    # ------------------------- #

    print('\nCopying samples.\nFrom:\t{input_path}\nTo:\t{output_path}\n'.format(input_path = input_path, output_path = output_path))

    # Create output Evidences directory
    output_evidence_dir = os.path.join(output_path, 'noise_evidences')
    if not os.path.exists(output_evidence_dir): os.makedirs(output_evidence_dir)

    # Create output SNR directory
    output_SNR_dir      = os.path.join(output_path, 'SNR_samples')
    if not os.path.exists(output_SNR_dir): os.makedirs(output_SNR_dir)

    # Loop on the different runs
    for file in os.listdir(input_path):
        if (not file == '.DS_Store') and (not file == 'TGR'):

            i = 0
            keys = []
            tmp = file.split('_')
            for element in elements:
                if elements[element] == '':
                    keys.append(tmp[i])
                    i += 1
                else:
                    keys.append(elements[element])

            if deviations: keys[5] = tmp[-1]   # Set deviations as nonGR_tag

            filename_tmp = '{}_{}_{}_{}_{}_{}'.format(keys[0], keys[1], keys[2], keys[3], keys[4], keys[5]) # Root filename
            nested_sampler_path = os.path.join(input_path,  file, 'Nested_sampler')                         # Nested Sampler directory path

            # Samples
            filename            = filename_tmp + '.h5'
            input_filename      = os.path.join(nested_sampler_path, '{}.h5'.format(sampler))
            output_filename     = os.path.join(output_path, filename)

            # Evidences
            filename_evidence   = filename_tmp + '_noise.txt'
            input_evidence      = os.path.join(nested_sampler_path, 'Evidence.txt')
            output_evidence     = os.path.join(output_evidence_dir, filename_evidence)

            # SNR
            filename_SNR        = filename_tmp + '_SNR.dat'
            input_SNR           = os.path.join(nested_sampler_path, 'optimal_SNR_TD.dat')
            output_SNR          = os.path.join(output_SNR_dir, filename_SNR)

            os.system('scp -r {input_filename} {output_filename}'.format(input_filename = input_filename, output_filename = output_filename))
            os.system('scp -r {input_evidence} {output_evidence}'.format(input_evidence = input_evidence, output_evidence = output_evidence))
            os.system('scp -r {input_SNR} {output_SNR}'.format(               input_SNR = input_SNR,           output_SNR = output_SNR))

elif cosmology:
 # ------------------------- #
    input_path  = '/Users/vgennari/Documents/work/code/python/icarogw/results/simulations/evolving_population/PROD1/low_SNR/inj-NNN'
    output_path = '/Users/vgennari/Documents/work/code/python/namib/samples/cosmology/simulations_evolving_population_PROD1_inj-NNN'

    elements = {
        'name'     : '',
        'pipeline' : '',
        'model'    : '',
        'submodel' : '',
        'time'     : 'NN',
        'GR_tag'   : '',
    }
    # ------------------------- #

    print('\nCopying samples.\nFrom:\t{input_path}\nTo:\t{output_path}\n'.format(input_path = input_path, output_path = output_path))

    # Create output directory
    if not os.path.exists(output_path): os.makedirs(output_path)

    # Loop on the different runs
    for file in os.listdir(input_path):
        if (not file == '.DS_Store'):

            i = 0
            keys = []
            tmp = file.split('_')
            for element in elements:
                if elements[element] == '':
                    keys.append(tmp[i+1])
                    i += 1
                else:
                    keys.append(elements[element])

            filename_tmp = '{}_{}_{}_{}_{}_{}'.format(keys[0], keys[1], keys[2], keys[3], keys[4], keys[5]) # Root filename
            nested_sampler_path = os.path.join(input_path,  file, 'sampler')                                # Samples directory path

            # Samples
            filename            = filename_tmp + '.json'
            input_filename      = os.path.join(nested_sampler_path, 'label_result.json')
            output_filename     = os.path.join(output_path, filename)

            os.system('scp -r {input_filename} {output_filename}'.format(input_filename = input_filename, output_filename = output_filename))


print('Finished.\n')
