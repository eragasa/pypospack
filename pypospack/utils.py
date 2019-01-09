import os

def get_pypospack_root_directory():
    _pythonpath_dirs = [v.strip() for v in os.environ["PYTHONPATH"].split(':') if v.endswith('pypospack')]
    _pypospack_root_dir = ""
    for d in _pythonpath_dirs:
        if d.find('pypospack'):
            _pypospack_root_dir = d
    return _pypospack_root_dir


def tail(fname, n_lines, offset=None):
    """ replicates the tail command from unix like operations systems
    
    Args:
        fname (str): filename
        n_lines (int): the number of lines

    Note:
        this is dependent upon the tail command in the unx operating system
        this should actually be rewritten as to be OS agnostic.
    """
    cmd_str = "/usr/bin/tail -n {} {}".format(str(n_lines), fname)
    p = subprocess.Popen(cmd_str, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    lines = stdout.decode('ascii').splitlines()
    return lines
