import sys, inspect
import pypospack.qoi

class QoiMap(object):

    def __init__(self):
        pass

    def read_qoi_directory():
        pass

if __name__ == "__main__":
    class_members = inspect.getmembers(
        sys.modules[pypospack.qoi],
        inspect.isclass
    )
    for class_ in class_members:
        print(class_)
