import pytest

from pypospack.pyposmat.data import PyposmatConfigurationFile

configuration_fn = 'data/pyposmat.config.in'

config = PyposmatConfigurationFile()
config.read(filename=configuration_fn)

print(config.sampling_type)
print(config.configuration['sampling_type'])
