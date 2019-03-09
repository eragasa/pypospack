import os
import shutil

def remove_rank_directories(directory=None):
    if directory is None:
        rank_directories = [os.path.join(os.getcwd(),d) for d in os.listdir() if d.startswith('rank_')]
    else:
        rank_directories = [os.path.join(directory,d) for d in os.listdir(directory) if d.startswith('rank_')]

    for d in rank_directories:
        print(d)
        if os.path.isdir(d):
            shutil.rmtree(d)

def remove_data_files(data_directory=None):
    for f in os.listdir(data_directory):
        if f.endswith(".py"):
            print('skip:',f)
        elif f == "pyposmat.kde.0.out":
            print('skip:',f)
        elif f == "pyposmat.config.in":
            print('skip:',f)
        else:
            print('delete:',f)
            if os.path.isfile(os.path.join(data_directory,f)):
                os.remove(os.path.join(data_directory,f))


if __name__ == "__main__":
    remove_rank_directories()
    remove_data_files(data_directory='data')


