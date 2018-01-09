import nibabel as nib
import numpy as np
import h5py
import os
from scipy.spatial import distance
from scipy import stats
import time


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
    def __init__(self, hdf_list, fnode_img):
        """

        Parameters
        ----------
        hdf_list: a list of subjects
        fnode_img: a mask nifti file

        Returns
        -------

        """
        # subjects list
        try:
            with open(hdf_list, 'r') as f:
                file_dir = [line.strip() for line in f]
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
    def __init__(self, data_hdf, group, method='ttest'):
        """

        Parameters
        ----------
        data_hdf: a hdf file for 3-D array of inter-subject correlation [nx: voxels; ny(z): subjects]
        group: a file or 1-D array indexing group labels: 0 - patients; 1 - controls
        method: string, statistical method: "ttest" or "permutation", if not specified, "ttest" is used.

        Returns
        -------

        """
        # check if method legal
        if method is not 'ttest' and method is not 'permutation':
            raise UserDefinedException("Method shall be 'ttest' or 'permutation'!")

        self.method = method

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

        # nested properties
        class B:
            p = []
            t = []

        self.patient = B()
        self.control = B()
        self.variance_diff = B()
        self.pattern_diff = B()

        # results container will be later filled
        self.tvalue = []
        self.pvalue = []
        self.value = []

    def compute(self, m=None):
        """

        Parameters
        ----------
        m: number of permutations
        Returns
        -------

        """
        t1 = time.time()

        # mask for patients and controls.
        dp_mas = (self.grp != 1)
        ct_mas = (self.grp != 0)

        # calculate within groups and between-group subject-based similarity,
        # resulting 2-D array: i, voxel, j, subject
        dp_2d = np.nanmean(self.data[np.ix_(self.svox, dp_mas, dp_mas)], axis=2)
        ct_2d = np.nanmean(self.data[np.ix_(self.svox, ct_mas, ct_mas)], axis=2)
        inter_2d = np.nanmean(self.data[np.ix_(self.svox, dp_mas, ct_mas)], axis=2)

        # statistical testing
        # 1 sample t test for patient ~= 0
        self.patient.t, self.patient.p = stats.ttest_1samp(dp_2d, 0, axis=1)
        # 1 sample t test for control ~= 0
        self.control.t, self.control.p = stats.ttest_1samp(ct_2d, 0, axis=1)
        # independent two-sample t test for patient - control ~= 0
        self.variance_diff.t, self.variance_diff.p = stats.ttest_ind(dp_2d, ct_2d, axis=1)
        # independent two-sample t test for patient_relative_control - control
        self.pattern_diff.t, self.pattern_diff.p = stats.ttest_ind(inter_2d, ct_2d, axis=1)

        elapse = time.time() - t1

        # permutation test
        if self.method is 'permutation':

            # number of permutations
            try:
                if m is None:
                    m = 1000
                print "{:d} permutations will begin, estimated time cost is {:.2f} hours ...".format(m, elapse*m/3600.0)
                # delay for 5 seconds for users to decide to go on
                time.sleep(5)
            except:
                raise UserDefinedException("m shall be an integer")

            # preallocate memory for later use. i: number of permutations; j: voxels
            t_patient = np.empty((m, self.nvox))
            t_control = np.empty((m, self.nvox))
            t_variance_diff = np.empty((m, self.nvox))
            t_pattern_diff = np.empty((m, self.nvox))

            # loops for correlation. This may be very slow and can be later improved.
            n = 0
            while n < m:
                print "permutation {:d}".format(n+1)
                perms = np.random.permutation(self.grp)
                # mask for patients and controls
                dp_mas_perms = (perms != 1)
                ct_mas_perms = (perms != 0)
                # calculate within groups and between-group subject-based similarity,
                #  resulting 2-D array: i, voxel, j, subject
                dp_2d_perms = np.nanmean(self.data[np.ix_(self.svox, dp_mas_perms, dp_mas_perms)], axis=2)
                ct_2d_perms = np.nanmean(self.data[np.ix_(self.svox, ct_mas_perms, ct_mas_perms)], axis=2)
                inter_2d_perms = np.nanmean(self.data[np.ix_(self.svox, dp_mas_perms, ct_mas_perms)], axis=2)
                # statistics
                t_patient[n, :] = stats.ttest_1samp(dp_2d_perms, 0, axis=1)[0]
                t_control[n, :] = stats.ttest_1samp(ct_2d_perms, 0, axis=1)[0]
                t_variance_diff[n, :] = stats.ttest_ind(dp_2d_perms, ct_2d_perms, axis=1)[0]
                t_pattern_diff[n, :] = stats.ttest_ind(inter_2d_perms, ct_2d_perms, axis=1)[0]
                n += 1

            # statistical testing: 2-tailed
            self.patient.p = np.mean(t_patient > abs(self.patient.t[np.newaxis, :]), axis=0) * 2
            self.control.p = np.mean(t_control > abs(self.control.t[np.newaxis, :]), axis=0) * 2
            self.variance_diff.p = np.mean(t_variance_diff > abs(self.variance_diff.t[np.newaxis, :]), axis=0) * 2
            self.pattern_diff.p = np.mean(t_pattern_diff > abs(self.pattern_diff.t[np.newaxis, :]), axis=0) * 2

        # results combine; i: condition (subjects) j: voxel
        self.tvalue = np.stack((self.patient.t, self.control.t, self.variance_diff.t, self.pattern_diff.t))
        self.pvalue = np.stack((self.patient.p, self.control.p, self.variance_diff.p, self.pattern_diff.p))
        patients = tuple(np.transpose(inter_2d))
        self.value = np.stack(patients)

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

        tvalue = np.zeros((dim[0], dim[1], dim[2], self.tvalue.shape[0]))
        pvalue = np.zeros((dim[0], dim[1], dim[2], self.pvalue.shape[0]))
        value = np.zeros((dim[0], dim[1], dim[2], self.value.shape[0]))

        for i in xrange(self.tvalue.shape[0]):
            # note the "advance indexing"
            tvalue[ncoords[:, 0], ncoords[:, 1], ncoords[:, 2], i] = self.tvalue[i, :]
            pvalue[ncoords[:, 0], ncoords[:, 1], ncoords[:, 2], i] = self.pvalue[i, :]

        for i in xrange(self.value.shape[0]):
            value[ncoords[:, 0], ncoords[:, 1], ncoords[:, 2], i] = self.value[i, :]

        # remove nan
        tvalue[np.isnan(tvalue)] = 0
        pvalue[np.isnan(pvalue)] = 0
        value[np.isnan(value)] = 0

        # save t-value
        header['cal_max'] = tvalue.max()
        header['cal_min'] = tvalue.min()
        tvalue_img = nib.Nifti1Image(tvalue, None, header)
        nib.save(tvalue_img, os.path.join(out_dir, '_'.join(('tvalue', self.method)) + '.nii.gz'))

        # save p-value
        header['cal_max'] = pvalue.max()
        header['cal_min'] = pvalue.min()
        pvalue_img = nib.Nifti1Image(pvalue, None, header)
        nib.save(pvalue_img, os.path.join(out_dir, '_'.join(('pvalue', self.method)) + '.nii.gz'))

        # save individual value
        header['cal_max'] = value.max()
        header['cal_min'] = value.min()
        value_img = nib.Nifti1Image(value, None, header)
        nib.save(value_img, os.path.join(out_dir, '_'.join(('value', self.method)) + '.nii.gz'))
