# ISC
inter-subject similarity/correlation

This is a tool for calculating inter-subject correlation based on the resting-state global brain connectivity pattern. As the required memory for this calculation is far exceeding the affordance of modern home PCs, the excellent hdf5 format is applied to save the data into hard disk and later read them sequentially for computation. This approach trades off between memory usage and computational efficiency. It also requires large space of hard disk (fortunately hard disk is cheap!). For a single subject with 40000 voxels, it typically requires about 13 Gb hard disk storage.

# Usage in ipython console

### create a dataset object
ds = ISCspace.DataSet(data_dir, mask_dir)

### computing global connectivity metric
conn = ISCspace.Connectivity(ds).compute()

### save the connectivity in hdf5 format
conn.save(output_dir)

### create isc object
isc =ISCspace.Intersubj(hdf5_list)

### compute isc statistic (t-test or permutation)
isc.compute()

### save into nii.gz image (t-value, p-value, and individual metric)
isc.save(result_dir)

