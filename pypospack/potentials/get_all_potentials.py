import os
import pkgutil
import pypospack.potentials

potential_pkg_path = os.path.dirname(pypospack.potentials.__file__)

for a,b,c in pkgutil.iter_modules([potential_pkg_path]):
    print(a,b,c)
potential_module_names = [name for _, name, _ in pkgutil.iter_modules([potential_pkg_path])]

for mod_name in potential_module_names:
    print(mod_name)
    mod = mod_name
    potential_class_names = dict([(name, cls) for name, cls in mod.__dict__.items() if isinstance(cls, type)])
    print(potential_class_names)

