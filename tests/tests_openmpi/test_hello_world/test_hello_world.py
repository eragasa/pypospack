from subprocess import call
import os

# possibly the greatest hack of a test to ever be written
def test_hello_world():
    call(["sbatch", "runjob_hipergator.sh"])

    while True:
        r0 = os.path.exists("rank_0")
        r1 = os.path.exists("rank_1")
        r2 = os.path.exists("rank_2")
        r3 = os.path.exists("rank_3")
        err = os.path.exists("job.err")
        if all([r0, r1, r2, r3]):
            os.remove("rank_0")
            os.remove("rank_1")
            os.remove("rank_2")
            os.remove("rank_3")
            assert True
            return
        if err:
            os.remove("job.err")
            os.remove("job.out")
            assert False
            return
