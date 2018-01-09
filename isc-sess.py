import ISCspace
import os

f = open('xxxx')
subject_list = [line.strip() for line in f]

node = 'xxx'

base_dir = '/nfs/h1/workingshop/zhaoyuanfang/DPrest'
control_dir = 'restnii_OldControl'
DP_dir = 'restnii_DP'
template = '%s/resting/002/sm6_inorm_bp0.01_0.1_csf_wm_gs_mc_confrm_lin_3mm.nii.gz'
save_dir = '/nfs/e4/yuanfang_tmp'

for sub in subject_list:

    if sub[0] is 'D':
        targ = os.path.join(base_dir, DP_dir, template) % sub
    else:
        targ = os.path.join(base_dir, control_dir, template) % sub

    ds = ISCspace.DataSet(targ, node)
    conn = ISCspace.Connectivity(ds).compute()

    out_dir = os.path.join(save_dir, sub)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    conn.save(out_dir)



