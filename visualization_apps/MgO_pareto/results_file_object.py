import os
import numpy as np
import pandas as pd

class ResultsFileObject(object):

    def __init__(self, fname):
        self.result_file_name = fname
        self._file_object = self.load_data_file()

    @property
    def get_total_dataframe(self):
        return self._file_object['dataframe']

    @property
    def get_param_names(self):
        return self._file_object['param_names']

    @property
    def get_qoi_names(self):
        return self._file_object['qoi_names']

    @property
    def get_err_names(self):
        return self._file_object['err_names']

    def load_data_file(self):
        names = []
        names_types = []
        data = []

        # read in the line
        lines = None
        try:
            with open(self.result_file_name, 'r') as f:
                lines = f.readlines()
        except:
            raise

        for i, line in enumerate(lines):
            line = line.strip()
            if i == 0:
                names = [v.strip() for v in line.split(',')]
            elif i == 1:
                name_types = [v.strip() for v in line.split(',')]
            else:
                data_line = [float(v.strip()) for v in line.split(',')]
                data.append(data_line)

        data = np.array(data)

        assert len(names) == len(name_types)
        # block below organizes data names by type into lists
        param_names = []
        param_key_index = []
        qoi_names = []
        qoi_key_index = []
        err_names = []
        err_key_index = []
        for i, v in enumerate(name_types):
            if v == 'param':
                param_names.append(names[i])
                param_key_index.append(i)
            elif v == 'qoi':
                qoi_names.append(names[i])
                qoi_key_index.append(i)
            elif v == 'err':
                err_names.append(names[i])
                err_key_index.append(i)

        # split array by data type
        param_data = data[:, min(param_key_index):max(param_key_index) + 1]
        qoi_data = data[:, min(qoi_key_index):max(qoi_key_index) + 1]
        err_data = data[:, min(err_key_index):max(err_key_index) + 1]

        # generate pandas dataframes
        param_df = pd.DataFrame(data=param_data, columns=param_names)
        qoi_df = pd.DataFrame(data=qoi_data, columns=qoi_names)
        err_df = pd.DataFrame(data=err_data, columns=err_names)
        total_df = pd.concat(
            [param_df,
             qoi_df,
             err_df], axis=1)

        # make copies to the class for persistence
        param_names = list(param_names)
        qoi_names = list(qoi_names)
        err_names = list(err_names)
        return {'dataframe': total_df, 'param_names': param_names,
                'qoi_names': qoi_names, 'err_names': err_names}


if __name__ == "__main__":
    data_dir = 'data'
    filename = 'culled_009.out'
    result_file = ResultsFileObject(fname=os.path.join(data_dir, filename))
