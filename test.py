# import profile

from ISCspace import *

ds = DataSet('../toy_data/sub03.nii.gz', '../toy_data/mask.nii.gz')

conn = Connectivity(ds)

conn.compute()

conn.save('../toy_data/')

isc = Intersubj('../toy_data/hdf5file')

isc.compute()

isc.save('../toy_data/')
