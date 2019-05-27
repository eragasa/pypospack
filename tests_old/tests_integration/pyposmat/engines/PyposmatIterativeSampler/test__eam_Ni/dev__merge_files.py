import os
def merge_files(i_iteration):
    n_ranks = 1

    _data_directory = 'data'

    _filename_kde = os.path.join(
        _data_directory,
        'pyposmat.kde.{}.out'.format(i_iteration))
    
    _filename_out = os.path.join(
        _data_directory,
        'pyposmat.results.{}.out'.format(i_iteration))
    # delete directory if directory exists, then create it
    if i_iteration == 0:
        if os.path.isdir(_data_directory):
            shutil.rmtree(_data_directory)
        os.mkdir(_data_directory)

    # gather all the lines into a list, then flush all at once
    str_list = []

    # grab the names, from the rank0 data file
    i_rank = 0
    _rank_filename = os.path.join(
        'rank_{}'.format(i_rank),
        'pyposmat.results.out')
    with open(_rank_filename,'r') as f:
        lines = f.readlines()
    names = [v for v in lines[0].strip().split(',')]
    types = [v for v in lines[1].strip().split(',')]
    str_list.append(",".join(names))
    str_list.append(",".join(types))
    
    # first merge the kde file from the this iteration
    if i_iteration > 0:
        with open(_filename_kde,'r') as f:
            lines = f.readlines()
        for line in lines[2:]:
            line = [v.strip() for v in line.strip().split(',')]
            line = [line[0]] + [float(s) for s in line[1:]]
            line = [s for i,s in enumerate(line) if i<len(names)]
            str_list.append(",".join(line[len(names):]))
        print(str_list)
    # process all the results from this simulations
    for i_rank in range(n_ranks):
        # get the filename of the ith rank
        _rank_filename = os.path.join(
            'rank_{}'.format(i_rank),
            'pyposmat.results.out')
        print('...merging {}'.format(_rank_filename))

        with open(_rank_filename,'r') as f:
            lines = f.readlines()
        lines = [line.strip().split(',') for line in lines]
        for line in lines[2:]:
            line[0] = '{}_{}_{}'.format(i_iteration,i_rank,line[0])
            str_list.append(",".join(line))
    # write results to file
    with open(_filename_out,'w') as f_out:
        f_out.write("\n".join(str_list))

merge_files(1)
