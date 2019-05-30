import os,shutil

def delete_rank_directories():
    for k in os.listdir():
        if k.startswith('rank_'):
            shutil.rmtree(k)

def delete_generated_datafiles_in_data_directory(data_dir):
    for k in os.listdir(data_dir):
        fn = os.path.join(data_dir,k)
        if fn.endswith('.py'): pass
        elif fn == 'pyposmat.config.in':
            pass
        elif fn == 'pyposmat.kde.0.in':
            pass
        else:
            if os.path.isfile(fn):
                os.remove(fn)
            else:
                print(k)

if __name__ == "__main__":
    data_dir = 'data'
    delete_rank_directories()
    delete_generated_datafiles_in_data_directory(data_dir)
