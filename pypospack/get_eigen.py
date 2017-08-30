import copy
def get_eigenvalues_from_vasp(outcar):
    eigenvalues = []
    eigenvectors = []
    str_regex_eigen_1 = '\d+ f  =\s+\d+.\d+\sTHz\s+\d+.\d+\s2PiTHz\s+\d+.\d+\scm-1\s+\d+.\d+\smeV'
    str_regex_eigen_2 = '\d+ f/i=\s+\d+.\d+\sTHz\s+\d+.\d+\s2PiTHz\s+\d+.\d+\scm-1\s+\d+.\d+\smeV'
    regex_eigen = [re.compile(p) for p in [str_regex_eigen_1, str_regex_eigen_2]]
    file = open(outcar,'r')

    with open(outcar,'r') as f:
        while line in f:
            if(re.match(r' Eigenvectors and eigenvalues of the dynamical matrix',line)):
                is_done_eigen = False
                while not is_done_eigen:
                    line = f.readline()
                    for regex in regex_eigen:
                        if regex.search(line):
                            reResults = re.findall(r'\d+.\d+',line)
                            eigenenergy = reResults[2] # 0 in THz, 1 in 2PiThz, 2 in cm-1, 3 in meV
                            eigenvalues.append(float(eigenenergy))
                            f.readline()                   # skipping header line

                            eigenvector = []
                            for i in range(n_atoms):
                                line = f.readline()
                                reResults = re.findall(r'[-+]?\d*\.\d+|d+',line)
                                
                                eigenvector.append([reResults[0],reResults[1],reResults[2],
                                                    reResults[3],reResults[4],reResults[5]])
                            eigenvectors.append(copy.deepcopy(eigenvector))
                            f.readline()                  # skipping empty line
                    if(re.search(r'Finite differences POTIM',line)):
                        is_done_eigen = 1
                    if(re.search(r'--------------------------------------------------------------------------------------------------------',line)):
                        is_done_eigen = 1
