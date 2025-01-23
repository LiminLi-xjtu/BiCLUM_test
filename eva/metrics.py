r"""
Performance evaluation metrics
"""

from typing import Tuple
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.spatial
import sklearn.metrics
import sklearn.neighbors
from anndata import AnnData
from scipy.sparse.csgraph import connected_components
from sklearn.neighbors import NearestNeighbors, KNeighborsClassifier

from typehint import RandomState, get_rs



def mean_average_precision(
        x: np.ndarray, y: np.ndarray, neighbor_frac: float = 0.01, **kwargs
) -> float:
    r"""
    Mean average precision

    Parameters
    ----------
    x
        Coordinates
    y
        Cell type labels
    neighbor_frac
        Nearest neighbor fraction
    **kwargs
        Additional keyword arguments are passed to
        :class:`sklearn.neighbors.NearestNeighbors`

    Returns
    -------
    map
        Mean average precision
    """
    k = max(round(y.shape[0] * neighbor_frac), 1)
    nn = sklearn.neighbors.NearestNeighbors(
        n_neighbors=min(y.shape[0], k + 1), **kwargs
    ).fit(x)
    nni = nn.kneighbors(x, return_distance=False)
    match = np.equal(y[nni[:, 1:]], np.expand_dims(y, 1))
    return np.apply_along_axis(_average_precision, 1, match).mean().item()


def _average_precision(match: np.ndarray) -> float:
    if np.any(match):
        cummean = np.cumsum(match) / (np.arange(match.size) + 1)
        return cummean[match].mean().item()
    return 0.0


def normalized_mutual_info(x: np.ndarray, y: np.ndarray, **kwargs) -> float:
    r"""
    Normalized mutual information with true clustering

    Parameters
    ----------
    x
        Coordinates
    y
        Cell type labels
    **kwargs
        Additional keyword arguments are passed to
        :func:`sklearn.metrics.normalized_mutual_info_score`

    Returns
    -------
    nmi
        Normalized mutual information

    Note
    ----
    Follows the definition in `OpenProblems NeurIPS 2021 competition
    <https://openproblems.bio/neurips_docs/about_tasks/task3_joint_embedding/>`__
    """
    x = AnnData(X=x, dtype=x.dtype)
    sc.pp.neighbors(x, n_pcs=0, use_rep="X")
    nmi_list = []
    for res in (np.arange(20) + 1) / 10:
        sc.tl.leiden(x, resolution=res)
        leiden = x.obs["leiden"]
        nmi_list.append(sklearn.metrics.normalized_mutual_info_score(
            y, leiden, **kwargs
        ).item())
    return max(nmi_list)


def avg_silhouette_width(x: np.ndarray, y: np.ndarray, **kwargs) -> float:
    r"""
    Cell type average silhouette width

    Parameters
    ----------
    x
        Coordinates
    y
        Cell type labels
    **kwargs
        Additional keyword arguments are passed to
        :func:`sklearn.metrics.silhouette_score`
 
    Returns
    -------
    asw
        Cell type average silhouette width

    Note
    ----
    Follows the definition in `OpenProblems NeurIPS 2021 competition
    <https://openproblems.bio/neurips_docs/about_tasks/task3_joint_embedding/>`__
    """
    return (sklearn.metrics.silhouette_score(x, y, **kwargs).item() + 1) / 2


def graph_connectivity(
        x: np.ndarray, y: np.ndarray, **kwargs
) -> float:
    r"""
    Graph connectivity

    Parameters
    ----------
    x
        Coordinates
    y
        Cell type labels
    **kwargs
        Additional keyword arguments are passed to
        :func:`scanpy.pp.neighbors`

    Returns
    -------
    conn
        Graph connectivity
    """
    x = AnnData(X=x, dtype=x.dtype)
    sc.pp.neighbors(x, n_pcs=0, use_rep="X", **kwargs)
    conns = []
    for y_ in np.unique(y):
        x_ = x[y == y_]
        _, c = connected_components(
            x_.obsp['connectivities'],
            connection='strong'
        )
        counts = pd.value_counts(c)
        conns.append(counts.max() / counts.sum())
    return np.mean(conns).item()


def seurat_alignment_score(
        x: np.ndarray, y: np.ndarray, neighbor_frac: float = 0.01,
        n_repeats: int = 4, random_state: RandomState = None, **kwargs
) -> float:
    r"""
    Seurat alignment score

    Parameters
    ----------
    x
        Coordinates
    y
        Batch labels
    neighbor_frac
        Nearest neighbor fraction
    n_repeats
        Number of subsampling repeats
    random_state
        Random state
    **kwargs
        Additional keyword arguments are passed to
        :class:`sklearn.neighbors.NearestNeighbors`

    Returns
    -------
    sas
        Seurat alignment score
    """
    rs = get_rs(random_state)
    idx_list = [np.where(y == u)[0] for u in np.unique(y)]
    min_size = min(idx.size for idx in idx_list)
    repeat_scores = []
    for _ in range(n_repeats):
        subsample_idx = np.concatenate([
            rs.choice(idx, min_size, replace=False)
            for idx in idx_list
        ])
        subsample_x = x[subsample_idx]
        subsample_y = y[subsample_idx]
        k = max(round(subsample_idx.size * neighbor_frac), 1)
        nn = sklearn.neighbors.NearestNeighbors(
            n_neighbors=k + 1, **kwargs
        ).fit(subsample_x)
        nni = nn.kneighbors(subsample_x, return_distance=False)
        same_y_hits = (
            subsample_y[nni[:, 1:]] == np.expand_dims(subsample_y, axis=1)
        ).sum(axis=1).mean()
        repeat_score = (k - same_y_hits) * len(idx_list) / (k * (len(idx_list) - 1))
        repeat_scores.append(min(repeat_score, 1))  # score may exceed 1, if same_y_hits is lower than expected by chance
    return np.mean(repeat_scores).item()


def avg_silhouette_width_batch(
        x: np.ndarray, y: np.ndarray, ct: np.ndarray, **kwargs
) -> float:
    r"""
    Batch average silhouette width

    Parameters
    ----------
    x
        Coordinates
    y
        Batch labels
    ct
        Cell type labels
    **kwargs
        Additional keyword arguments are passed to
        :func:`sklearn.metrics.silhouette_samples`

    Returns
    -------
    asw_batch
        Batch average silhouette width

    Note
    ----
    Follows the definition in `OpenProblems NeurIPS 2021 competition
    <https://openproblems.bio/neurips_docs/about_tasks/task3_joint_embedding/>`__
    """
    s_per_ct = []
    for t in np.unique(ct):
        mask = ct == t
        try:
            s = sklearn.metrics.silhouette_samples(x[mask], y[mask], **kwargs)
        except ValueError:  # Too few samples
            s = 0
        s = (1 - np.fabs(s)).mean()
        s_per_ct.append(s)
    return np.mean(s_per_ct).item()


def neighbor_conservation(
        x: np.ndarray, y: np.ndarray, batch: np.ndarray,
        neighbor_frac: float = 0.01, **kwargs
) -> float:
    r"""
    Neighbor conservation score

    Parameters
    ----------
    x
        Cooordinates after integration
    y
        Coordinates before integration
    b
        Batch
    **kwargs
        Additional keyword arguments are passed to
        :class:`sklearn.neighbors.NearestNeighbors`

    Returns
    -------
    nn_cons
        Neighbor conservation score
    """
    nn_cons_per_batch = []
    for b in np.unique(batch):
        mask = batch == b
        x_, y_ = x[mask], y[mask]
        k = max(round(x.shape[0] * neighbor_frac), 1)
        nnx = sklearn.neighbors.NearestNeighbors(
            n_neighbors=min(x_.shape[0], k + 1), **kwargs
        ).fit(x_).kneighbors_graph(x_)
        nny = sklearn.neighbors.NearestNeighbors(
            n_neighbors=min(y_.shape[0], k + 1), **kwargs
        ).fit(y_).kneighbors_graph(y_)
        nnx.setdiag(0)  # Remove self
        nny.setdiag(0)  # Remove self
        n_intersection = nnx.multiply(nny).sum(axis=1).A1
        n_union = (nnx + nny).astype(bool).sum(axis=1).A1
        nn_cons_per_batch.append((n_intersection / n_union).mean())
    return np.mean(nn_cons_per_batch).item()


def foscttm(
        x: np.ndarray, y: np.ndarray, **kwargs
) -> Tuple[np.ndarray, np.ndarray]:
    r"""
    Fraction of samples closer than true match (smaller is better)

    Parameters
    ----------
    x
        Coordinates for samples in modality X
    y
        Coordinates for samples in modality y
    **kwargs
        Additional keyword arguments are passed to
        :func:`scipy.spatial.distance_matrix`

    Returns
    -------
    foscttm_x, foscttm_y
        FOSCTTM for samples in modality X and Y, respectively

    Note
    ----
    Samples in modality X and Y should be paired and given in the same order
    """
    if x.shape != y.shape:
        raise ValueError("Shapes do not match!")
    d = scipy.spatial.distance_matrix(x, y, **kwargs)
    foscttm_x = (d < np.expand_dims(np.diag(d), axis=1)).mean(axis=1)
    foscttm_y = (d < np.expand_dims(np.diag(d), axis=0)).mean(axis=0)
    return foscttm_x, foscttm_y


def transfer_accuracy(domain1, domain2, type1, type2):
    knn = KNeighborsClassifier()
    knn.fit(domain2, type2)
    type1_predict = knn.predict(domain1)
    count = 0
    for label1, label2 in zip(type1_predict, type1):
        if label1 == label2:
            count += 1
    return count / len(type1)

def get_k_neigh_ind(X,k):
    neigh = NearestNeighbors(n_neighbors=k)
    neigh.fit(X)
    neigh_dist, neigh_ind = neigh.kneighbors(X)
    return neigh_dist, neigh_ind

def neigh_overlap(z_rna, z_atac, k = 30):
    dsize = z_rna.shape[0]
    _, neigh_ind = get_k_neigh_ind(np.concatenate((z_rna, z_atac), axis = 0), k = k)
#     print(neigh_ind)
    z1_z2 = ((neigh_ind[:dsize,:] - dsize - np.arange(dsize)[:, None]) == 0)
#     print(z1_z2)
    z2_z1 = (neigh_ind[dsize:,:] - np.arange(dsize)[:, None] == 0)
#     print(z2_z1)
    return 0.5 * (np.sum(z1_z2) + np.sum(z2_z1))/dsize

def all_metrics(combined_latent, combined_cell_type, anno_rna, anno_other, combined_domain, paired=True):

    graph_conn = graph_connectivity(combined_latent, combined_cell_type)
    seurat_align_score = seurat_alignment_score(combined_latent, combined_domain, random_state=0)
    avg_silh_width_batch = avg_silhouette_width_batch(combined_latent, combined_domain,
                                                                    combined_cell_type)

    if paired:
        if len(anno_rna) != len(anno_other):

            z_rna_ = combined_latent[0:len(anno_rna), :]
            z_other_ = combined_latent[len(anno_rna):, :]

            name_other = pd.DataFrame(anno_other.index)
            name_rna = pd.DataFrame(anno_rna.index)
            name_other_ = name_other[0].apply(lambda x: '-'.join(x.split('-')[:-1]))
            name_rna_ = name_rna[0].apply(lambda x: '-'.join(x.split('-')[:-1]))

            name_rna_dict = pd.Series(name_rna_.index, index=name_rna_.values).to_dict()
            match_positions = name_other_.map(name_rna_dict)
            z_rna = z_rna_[match_positions, :]
            z_other = z_other_

            neighborhood_overlap_score = neigh_overlap(z_rna, z_other, k=30)
            omics_mixing = (neighborhood_overlap_score + graph_conn + seurat_align_score + avg_silh_width_batch) / 4

            fosc = np.concatenate(foscttm(z_rna, z_other)).mean().item()

        else:
            z_rna_ = combined_latent[0:len(anno_rna), :]
            z_other_ = combined_latent[len(anno_rna):, :]
            neighborhood_overlap_score = neigh_overlap(z_rna_, z_other_, k=30)
            omics_mixing = (neighborhood_overlap_score + graph_conn + seurat_align_score + avg_silh_width_batch) / 4

            fosc = np.concatenate(foscttm(z_rna_, z_other_)).mean().item()

    else:
        neighborhood_overlap_score = None
        omics_mixing = (graph_conn + seurat_align_score + avg_silh_width_batch) / 3

        fosc = None


    mean_avg_precision = mean_average_precision(combined_latent, combined_cell_type)
    avg_silh_width = avg_silhouette_width(combined_latent, combined_cell_type)
    norm_mutual_info = normalized_mutual_info(combined_latent, combined_cell_type)
    biology_conservation = (mean_avg_precision + avg_silh_width + norm_mutual_info) / 3


    type_rna = combined_cell_type[0:len(anno_rna)]
    type_other = combined_cell_type[len(anno_rna):]
    transfer_acc = transfer_accuracy(z_rna_, z_other_, type_rna, type_other)

    metrics_list = [omics_mixing, biology_conservation, transfer_acc, fosc]


    return metrics_list
