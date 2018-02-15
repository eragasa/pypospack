import os,shutil

i_iteration = 0
n_ranks = 8
_data_directory = 'data'
_filename_out = os.path.join(
    _data_directory,
    'pyposmat.results.{}.out'.format(i_iteration))

if i_iteration == 0:
    if os.path.isdir(_data_directory):
        shutil.rmtree(_data_directory)
    os.mkdir(_data_directory)

str_list = []
for i_rank in range(n_ranks):
    _filename = os.path.join(
        'rank_{}'.format(i_rank),
        'pyposmat.results.out')
    print('...merging {}'.format(_filename))
    if i_rank == 0:
        with open(_filename,'r') as f:
            lines = f.readlines()
            names = [v for v in lines[0].strip().split(',')]
            types = [v for v in lines[1].strip().split(',')]
            str_list.append(",".join(names))
            str_list.append(",".join(types))
    if i_iteration != 0:
        _filename_kde = os.path.join(
            _data_directory,
            'pyposmat.kde.{}.out'.format(i_iteration))
        with open(_filename_kde,'r') as f:
            lines = f.readlines()
        for k,line in enumerate(lines):
            if k>1:
                line = [v for v in lines[0].strip().split(',')]
                str_list.append(",".join(line[len(names):]))
    with open(_filename,'r') as f:
        lines = f.readlines()
    lines = [line.strip().split(',') for line in lines]
    for k,line in enumerate(lines):
        if k>1:
            line[0] = '{}_{}_{}'.format(i_iteration,i_rank,line[0])
            str_list.append(",".join(line))

with open(_filename_out,'w') as f_out:
    f_out.write("\n".join(str_list))
