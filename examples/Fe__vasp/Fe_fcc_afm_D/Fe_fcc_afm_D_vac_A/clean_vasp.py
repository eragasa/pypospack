import os

filenames_delete = [
    'CHG',
    'CHGCAR',
    'CONTCAR',
    'DOSCAR',
    'EIGENVAL',
    'IBZKPT',  
    'job.err',
    'job.out',  
    'OSZICAR',  
    'PCDAT',  
    'REPORT',  
    'vasp.log',  
    'vasprun.xml',  
    'WAVECAR',  
    'XDATCAR'
]

for filename in filenames_delete:
    try:
        os.remove(filename)
        msg = "{} removed.".format(filename)
    except FileNotFoundError as e:
        msg = "{} does not exist.".format(filename)
    except:
        raise
    print(msg)
