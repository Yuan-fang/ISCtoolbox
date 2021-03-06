# Space-based inter-subject similarity/correlation

This is a tool for calculating inter-subject correlation based on the resting-state global brain connectivity pattern. As the required memory for this calculation is far exceeding the affordance of modern home PCs, the excellent [hdf5 format](https://support.hdfgroup.org/HDF5/whatishdf5.html) is applied to save the data into hard disk and later read them sequentially for computation. This approach trades off between memory usage and computational efficiency. It also requires large space of hard disk (fortunately hard disk is cheap!). For a single subject with 40000 voxels, it typically requires about 13 Gb hard disk storage.

# Usage in ipython console

* ***create a dataset object***  
ds = ISCspace.DataSet(data_dir, mask_dir)

* ***computing global connectivity metric***  
conn = ISCspace.Connectivity(ds).compute()

* ***save the connectivity in hdf5 format***  
conn.save(output_dir)

* ***create isc object***  
isc = ISCspace.Intersubj(hdf5_list)

* ***compute isc***  
isc.compute()

* ***save the isc data in hdf5 format***  
isc.save(result_dir)

* ***create stat object***  
stat = ISCspace.Statistic(data_hdf, group)

* ***compute individual-based similarity metric***  
stat.compute()

* ***save results into nii.gz image***  
stat.save(mask, output_dir)
