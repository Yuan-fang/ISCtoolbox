# import profile

from ISCspace import *

ds = DataSet('../toy_data/sub03.nii.gz', '../toy_data/mask.nii.gz')

conn = Connectivity(ds)

conn.compute()

conn.save('../toy_data/')

isc = Intersubj('../toy_data/hdf5file')

isc.compute()

isc.save('../toy_data/')


# if __name__ == "__main__":
#    profile.run("ds")

class A(object):
    def __init__(self):
        class B:
            p = []

        self.patient = B()
        self.control = B()

    def compute(self):
        self.patient.p = np.array([1, 2, 3])
        self.control.p = 2

    def test(self):
        print self.patient.shape


class A(object):
    def __init__(self, para1, para2, para3):
        value = 'xxx'
        global value
        self.a1 = value
        self.a2 = 'yyy'
        self.a3 = 'zzz'

class B(object):
    def __init__(self):
        self.b1 = value

"test"
print B.b1
