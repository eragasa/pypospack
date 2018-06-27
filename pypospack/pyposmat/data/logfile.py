import time
import fcntl
import errno


class PyposmatLogFile(object):

    def __init__(self, filename=None):

        self.filename = filename
        self.time_sleep = 0.1
        self.msg_format = "{}\n"

    def write(self, s):

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
            f.write(self.msg_format.format(s))
            
            # release the lock
            fcntl.flock(f, fcntl.LOCK_UN)
