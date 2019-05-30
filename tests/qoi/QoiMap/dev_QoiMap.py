import importlib.util
from pathlib import Path
import os
MODULE_EXTENSIONS = '.py'

class QoiMap(object):

    def __init__(self):
        pass

    def read_qoi_directory():
        pass

def package_contents(package_name):
    spec = importlib.util.find_spec(package_name)
    if spec is None:
        return set()

    pathname = Path(spec.origin).parent
    print(pathname)
    print(type(pathname))
    ret = set()
    with os.scandir(pathname) as entries:
        for entry in entries:
            if entry.name.startswith('__'):
                continue
            current = '.'.join((package_name, entry.name.partition('.')[0]))
            if entry.is_file():
                if entry.name.endswith(MODULE_EXTENSIONS):
                    ret.add(current)
                elif entry.is_dir():
                    ret.add(current)
                    ret |= package_contents(current)

    return ret
if __name__ == "__main__":
    import pkgutil
    package_name = 'pypospack.qois'
    package_path = Path(importlib.util.find_spec(package_name).origin).parent
    module_names = [module_name for _,module_name,_ in pkgutil.iter_modules([str(package_path)])]
    print(module_names)

    import pyclbr
    for module_name in module_names:
        module_info = pyclbr.readmodule('.'.join([package_name,module_name]))
        for item in module_info.values():
            print(module_name,item.name)
    #print(package_contents('pypospack.qois'))
