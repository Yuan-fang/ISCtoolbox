import nibabel as nib
import numpy as np
import h5py
import os
from scipy.spatial import distance


class UserDefinedException(Exception):
    """
    Exception defined by user
    """

    def __init__(self, str):
        """

        Parameters
        ----------
        str: a string to indicate the exception

        """
        Exception.__init__(self)
        self._str = str


def load_img(fimg):
    """
    Load Nifti1Image

    Parameters
    ----------
    fimg: an image file

    Returns
    -------
    img: a Nifti1Image object

    """
    # load nifti image with nibabel
    try:
        img = nib.load(fimg)
    except:
        raise UserDefinedException('Wrong Image!')

    return img


def hdfreader(file_name, nrow):
    """

    Parameters
    ----------
    file_name: a hdf5 file directory
    nrow: the Nth row of the file

    Returns
    -------

    """
    print "voxel {:d}".format(nrow)
    with h5py.File(file_name, 'r') as f:
        row = f['data'][nrow, :]

    return row


class DataSet(object):
    def __init__(self, ftarg_img, fnode_img):
        """

        Parameters
        ----------
        ftarg_img: target image file
        fnode_img: node image file

        Returns
        -------

        """
        # load target image
        targ_img = load_img(ftarg_img)
        if len(targ_img.shape) is 4:
            targ = targ_img.get_data()
        else:
            raise UserDefinedException('Target image is not a 4D Nifti volume!')

        # load node image
        node_img = load_img(fnode_img)
        if (len(node_img.shape) is 3) and (node_img.shape == targ_img.shape[:3]):
            node = node_img.get_data()
            # node mask
            nmas = (node != 0)
        else:
            raise UserDefinedException('Node image and target image are not match!')

        # time course
        self.tc = targ[nmas, :]


class Connectivity(object):
    def __init__(self, ds):
        """

        Parameters
        ----------
        ds: DataSet object

        Returns
        -------

        """
        self.ds = ds
        # mat will be assigned in self.compute()
        self.mat = []

    def compute(self):
        """

        Returns
        -------
        self: A Connectivity object

        """
        ds = self.ds
        self.mat = np.corrcoef(ds.tc)

        # fisher-z transformation.
        self.mat = np.arctanh(self.mat)

        return self

    def save(self, out_dir='.'):
        """

        Parameters
        ----------
        out_dir: dir to save the connectivity matrix

        Returns
        -------

        """
        fi = os.path.join(out_dir, 'conn.hdf')
        with h5py.File(fi, 'w') as f:
            f.create_dataset('data', data=self.mat)


class Intersubj(object):
    def __init__(self, base_dir, sessid, fnode_img):
        """

        Parameters
        ----------
        base_dir: directory for hdf5 data
        sessid: a text of subjects
        fnode_img: a mask nifti file

        Returns
        -------

        """
        # subjects list
        try:
            with open(sessid, 'r') as f:
                file_dir = [os.path.join(base_dir, line.strip(), 'conn.hdf') for line in f]
            self.list = file_dir
        except:
            raise UserDefinedException('A text file of hdf file shall be provided!')

        # sample size
        self.nsubs = len(file_dir)

        # voxel size
        node_img = load_img(fnode_img)
        if len(node_img.shape) is 3:
            node = node_img.get_data()
            self.nodes = np.count_nonzero(node)
        else:
            raise UserDefinedException('Node image shall be a 3d mask nifti file!')

        # preallocate memory for later use
        self.mat = np.empty((self.nodes, self.nsubs, self.nsubs))

    def compute(self):
        """

        Returns
        -------
        self: An Intersubj object

        """
        # compute inter-subject pattern correlation
        for node in xrange(self.nodes):
            dn = [hdfreader(file_dir, node) for file_dir in self.list]
            # change to numpy array
            dn = np.array(dn)
            # remove columns corresponding to the diagonal
            dn = np.delete(dn, node, 1)
            self.mat[node, :, :] = distance.cdist(dn, dn, 'correlation')

        # correlation distance to correlation coefficient
        self.mat = 1 - self.mat
        # fisher-z transformation
        self.mat = np.arctanh(self.mat)
        # replace the diagonal in y,z plane to nan
        i, j, k = np.indices(self.mat.shape)
        self.mat[j == k] = np.nan

        return self

    def save(self, out_dir='.'):
        """
        Parameters
        ----------
        out_dir: dir to save the connectivity matrix
        Returns
        -------
        """
        fi = os.path.join(out_dir, 'inter.hdf')
        with h5py.File(fi, 'w') as f:
            f.create_dataset('data', data=self.mat)


class Statistic(object):
    def __init__(self, data_hdf, group):
        """

        Parameters
        ----------
        data_hdf: a hdf file for 3-D array of inter-subject correlation [nx: voxels; ny(z): subjects]
        group: a file or 1-D array indexing group labels: 0 - patients; 1 - controls

        Returns
        -------

        """
        # check if group legal
        try:
            # when group is a 1-d numpy array
            if type(group) is np.ndarray and group.ndim is 1:
                self.grp = group
            # when group is a text file
            else:
                self.grp = np.loadtxt(group)
        except:
            raise UserDefinedException('group shall be a text file or 1-D array indexing group labels!')

        # fetch hdf5 dataset into memory (currently no better method)
        with h5py.File(data_hdf, 'r') as f:
            self.data = f['data'].value

        # number of voxels
        self.nvox = len(self.data)

        # sequence of voxels
        self.svox = np.arange(self.nvox)

        # results container will be later filled
        self.patients_value = []
        self.controls_value = []
        self.patients_relative_value = []

    def compute(self):
        """

        Parameters
        ----------
        Returns
        -------

        """
        # mask for patients and controls.
        dp_mas = (self.grp != 1)
        ct_mas = (self.grp != 0)

        # calculate within groups and between-group subject-based similarity,
        # resulting 2-D array: i, voxel, j, subject
        dp_2d = np.nanmean(self.data[np.ix_(self.svox, dp_mas, dp_mas)], axis=2)
        ct_2d = np.nanmean(self.data[np.ix_(self.svox, ct_mas, ct_mas)], axis=2)
        inter_2d = np.nanmean(self.data[np.ix_(self.svox, dp_mas, ct_mas)], axis=2)

        # individual level map
        patients = tuple(np.transpose(dp_2d))
        self.patients_value = np.stack(patients)

        controls = tuple(np.transpose(ct_2d))
        self.controls_value = np.stack(controls)

        patients_relative = tuple(np.transpose(inter_2d))
        self.patients_relative_value = np.stack(patients_relative)

    def save(self, fnode_img, out_dir='.'):
        """

        Parameters
        ----------
        self
        fnode_img
        out_dir

        Returns
        -------

        """
        node_img = load_img(fnode_img)
        header = node_img.header
        node = node_img.get_data()
        dim = node_img.header.get_data_shape()
        ncoords = np.transpose(np.nonzero(node))

        patients_value = np.zeros((dim[0], dim[1], dim[2], self.patients_value.shape[0]))
        controls_value = np.zeros((dim[0], dim[1], dim[2], self.controls_value.shape[0]))
        patients_relative_value = np.zeros((dim[0], dim[1], dim[2], self.patients_relative_value.shape[0]))

        for i in xrange(self.patients_value.shape[0]):
            patients_value[ncoords[:, 0], ncoords[:, 1], ncoords[:, 2], i] = self.patients_value[i, :]

        for i in xrange(self.controls_value.shape[0]):
            controls_value[ncoords[:, 0], ncoords[:, 1], ncoords[:, 2], i] = self.controls_value[i, :]

        for i in xrange(self.patients_relative_value.shape[0]):
            patients_relative_value[ncoords[:, 0], ncoords[:, 1], ncoords[:, 2], i] = self.patients_relative_value[i, :]

        # remove nan
        patients_value[np.isnan(patients_value)] = 0
        controls_value[np.isnan(controls_value)] = 0
        patients_relative_value[np.isnan(patients_relative_value)] = 0

        # save individual patients maps
        header['cal_max'] = patients_value.max()
        header['cal_min'] = patients_value.min()
        patients_value_img = nib.Nifti1Image(patients_value, None, header)
        nib.save(patients_value_img, os.path.join(out_dir, 'patients_value') + '.nii.gz')

        # save individual controls maps
        header['cal_max'] = controls_value.max()
        header['cal_min'] = controls_value.min()
        controls_value_img = nib.Nifti1Image(controls_value, None, header)
        nib.save(controls_value_img, os.path.join(out_dir, 'controls_value') + '.nii.gz')

        # save individual patients relative controls maps
        header['cal_max'] = patients_relative_value.max()
        header['cal_min'] = patients_relative_value.min()
        patients_relative_value_img = nib.Nifti1Image(patients_relative_value, None, header)
        nib.save(patients_relative_value_img, os.path.join(out_dir, 'patients_relative_value') + '.nii.gz')
