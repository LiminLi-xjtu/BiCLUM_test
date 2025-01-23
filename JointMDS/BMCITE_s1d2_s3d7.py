
import os
import torch
import argparse
import numpy as np

from joint_mds import JointMDS
from scipy.spatial.distance import pdist
import utils.scores as scores
from utils.utils import plot_embedding, geodesic_dist
import anndata
import sys

from eva.load_result import *
from eva.evaluation import evaluate
####################################################################
ncomp = 8
neighbor = 50
dataset_name = 'BMCITE_s1d1_s1d2'
batch = 's1d1_s1d2'
dataset_dir = '../../data/' + dataset_name
result_dir = "../results/" + dataset_name
####################################################################

# use_cuda = torch.cuda.is_available()
# device = torch.device('cuda:0' if use_cuda else 'cpu')
# print("Running on", device)


def main(dataset_dir):
    parser = argparse.ArgumentParser(
        description="Joint MDS for solving protein structure alignment",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--components", type=int, default=ncomp, help="number of components"
    )

    args = parser.parse_args()

    np.random.seed(0)
    torch.random.manual_seed(0)

    data1 = anndata.read_h5ad(os.path.join(dataset_dir, 'bm.rna.' + batch + '.h5ad'))
    data2 = anndata.read_h5ad(os.path.join(dataset_dir, 'bm.adt.' + batch + '.h5ad'))
    X1 = data1.X
    X2 = data2.X


    D1 = torch.from_numpy(geodesic_dist(X1, k=neighbor, mode="connectivity", metric="cosine")) # .to(device) # correlation
    D2 = torch.from_numpy(geodesic_dist(X2, k=neighbor, mode="connectivity", metric="cosine")) # .to(device) # correlation
    # dis = np.load(os.path.join(path, 'bm.rna.s1d2_s3d7.npy'), allow_pickle=True).item()
    # D1 = dis['dis1'].to(device)
    # D2 = dis['dis2'].to(device)


    JMDS = JointMDS(
        n_components=args.components,
        alpha=0.3,
        eps=0.1,
        max_iter=500,
        eps_annealing=False,
        dissimilarity="precomputed",
    )
    Z1, Z2, P = JMDS.fit_transform(D1, D2)

    Z1, Z2 = Z1.numpy(), Z2.numpy() #.detach().cpu()

    inte = []
    inte.append(Z1)
    inte.append(Z2)

    JM_inte = dict({"inte": inte})


    adata, anno_rna, anno_other = load_cite(dataset_name, dataset_dir, result_dir, ['JointMDS_new'], args.batch)
    eva_metrics = evaluate(adata, anno_rna, anno_other, args.paired)

    print(eva_metrics)

    path = result_dir
    if not os.path.exists(path):
        os.makedirs(path)

    np.save(os.path.join(path, 'JointMDS.npy'), JM_inte)

if __name__ == "__main__":
    main(dataset_dir)
