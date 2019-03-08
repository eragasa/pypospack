import os
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.engines import PyposmatFileSampler
#from pypospack.exceptions import PyposmatEngineError

n_smallest = 10

for v in os.environ["PYTHONPATH"].split(":"):
    if v.endswith('pypospack'):
        data_directory = os.path.join(v,"data/MgO_pareto/data")
config_directory = "./data"
output_directory = "./"

config_fn = os.path.join(config_directory,'pyposmat.config.in')
datafile_in_fn = os.path.join(data_directory,'culled_005.out')
datafile_out_fn = os.path.join(output_directory,'qoiplus_005.out')

o_config=PyposmatConfigurationFile()
o_config.read(filename=config_fn)

reference_potential_names = ["LC","BG1","BG2"]

def test__all_resources_exist():
    assert os.path.isdir(config_directory)
    assert os.path.isdir(data_directory)
    assert os.path.isdir(output_directory)

    assert os.path.isfile(config_fn)
    assert os.path.isfile(datafile_in_fn)
    assert os.path.isfile(datafile_out_fn)

def test____init__():

    o_sampler = PyposmatFileSampler(
            config_fn=config_fn,
            data_in_fn=datafile_in_fn,
            data_out_fn=datafile_out_fn)

    assert o_sampler.reference_potentials == reference_potentials

def reference_potential_names_to_string(o_sampler):
    s = 80*'-'+"\n"
    s += "{:^80}\n".format('CHECK REFERENCE POTENTIALS')
    s += 80*'-'+"\n"
    for v in o_sampler.reference_potentials:
        s += "{}.{}\n".format(
                v,
                v in o_sampler.reference_potentials)
    return s

if __name__ == "__main__":
    print ("config_directory:{}".format(config_directory))
    print ("data_directory:{}".format(data_directory))
    print ("output_directory:{}".format(output_directory))

    o_sampler = PyposmatFileSampler(
            config_fn=config_fn,
            data_in_fn=datafile_in_fn,
            data_out_fn=datafile_out_fn)

    print(reference_potential_names_to_string(o_sampler=o_sampler))
