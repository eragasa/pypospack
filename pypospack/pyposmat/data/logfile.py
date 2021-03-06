import time

try:
    import fcntl
except ModuleNotFoundError as e:
    import platform
    if platform.system() == 'Windows':
        print("fnctl cannot be loaded on a windows platform")
        pass
    else:
        raise
import errno


class PyposmatLogFile(object):

    def __init__(self, filename=None):

        self.filename = filename
        self.time_sleep = 0.1
        self.msg_format = "{}\n"

    def write(self, s):
        if type(s) is str:
            _s = s
        elif type(s) is list:
            _s = "\n".join(s)
        else:
            _s = s

        with open(self.filename, 'a') as f:

            # attempt to get lock
            while True:
                try:
                    fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    break
                except IOError as e:
                    if e.errno != errno.EAGAIN:
                        raise
                    else:
                        time.sleep(self.time_sleep)

            # write message
            f.write(self.msg_format.format(_s))

            # release the lock
            fcntl.flock(f, fcntl.LOCK_UN)
