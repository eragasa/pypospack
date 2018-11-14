import os

def get_pypospack_root_directory():
    python_root_dir = [v.strip() for v in os.environ["PYTHONPATH"].split(':') if v.endswith('pypospack')][0]

    return python_root_dir
