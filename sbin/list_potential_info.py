from pypospack.potential import PotentialObjectMap

_potential_type = 'all'

if _potential_type == 'all':
    map_info = PotentialObjectMap(potential_type=_potential_type)
    
    _max_key_length = 0
    _max_module_length = 0
    _max_class_length = 0

    for k,v in map_info.items():
        _max_key_length = max(_max_key_length,len(k))
        _max_module_length = max(_max_module_length,len(v['module']))
        _max_class_length = max(_max_class_length,len(v['class']))

    _fmt_key = '{:^'+str(_max_key_length)+'}'
    _fmt_module = '{:^'+str(_max_module_length)+'}'
    _fmt_class = '{:^'+str(_max_class_length)+'}'

    # define formats for different row types
    _fmt_header_row = " ".join([_fmt_key,_fmt_module,_fmt_class])
    _fmt_break_row = _fmt_header_row.format(
            _max_key_length*'-',
            _max_module_length*'-',
            _max_class_length*'-')
    _fmt_result_row = _fmt_header_row

    # print header row
    s = _fmt_header_row.format('name','module','class')
    print(s)

    # print break row
    s = _fmt_break_row
    print(s)

    # print result rows
    for k,v in map_info.items():
        s = _fmt_result_row.format(k,v['module'],v['class'])
        print(s)
else:
    _module,_class = PotentialObjectMap(potential_type=_potential_type)

