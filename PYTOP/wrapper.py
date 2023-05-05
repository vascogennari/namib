import os

# ------------------------- #
input_path  = '/Users/vascogennari/Documents/work/code/python/results/GW190521B_HMs/PROD0/Kerr-amp'
output_path = '/Users/vascogennari/Documents/work/code/python/PYTOP/samples/GW190521B/PROD0/Kerr_amp'

elements = {
    'name'     : '',
    'pipeline' : 'pyRing',
    'model'    : '',
    'submodel' : '',
    'time'     : '',
    'GR_tag'   : 'GR',
}
# ------------------------- #

print('\nCopying samples.\nFrom:\t{input_path}\nTo:\t{output_path}\n'.format(input_path=input_path, output_path=output_path))

for file in os.listdir(input_path):
    if not file == '.DS_Store':

        i = 0
        keys = []
        tmp = file.split('_')
        for element in elements:
            if elements[element] == '':
                keys.append(tmp[i])
                i += 1
            else:
                keys.append(elements[element])

        filename_tmp = '{}_{}_{}_{}_{}_{}'.format(keys[0], keys[1], keys[2], keys[3], keys[4], keys[5])
        filename          = filename_tmp + '.h5'
        filename_evidence = filename_tmp + '_noise.txt'

        nested_sampler_path = os.path.join(input_path,  file, 'Nested_sampler')
        input_filename  = os.path.join(nested_sampler_path, 'cpnest.h5')
        output_filename = os.path.join(output_path, filename)
        input_evidence  = os.path.join(nested_sampler_path, 'Evidence.txt')
        output_evidence_dir = os.path.join(output_path, 'noise_evidences')
        if not os.path.exists(output_evidence_dir): os.makedirs(output_evidence_dir)
        output_evidence = os.path.join(output_evidence_dir, filename_evidence)

        
        os.system('scp -r {input_filename} {output_filename}'.format(input_filename=input_filename, output_filename=output_filename))
        os.system('scp -r {input_evidence} {output_evidence}'.format(input_evidence=input_evidence, output_evidence=output_evidence))

print('Finished.')
