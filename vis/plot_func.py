import os
import pandas as pd
import scanpy as sc
from seaborn import heatmap

import matplotlib.pyplot as plt
from sklearn import metrics
import numpy as np
from munkres import Munkres



# plot metrics
def metrics_plot(eva_metrics, method_colors, vis_dir):

    categories = eva_metrics.index.tolist()
    eva_metrics.columns = ['Omics mixing', 'Cell Type conservation', 'transfer accuracy', 'FOSCTTM']
    method_colors = {key: value for key, value in method_colors.items() if key in categories}


    eva_metrics_overall = eva_metrics.iloc[:, :2]
    columns = eva_metrics_overall.columns
    fig, ax = plt.subplots(figsize=(8, 5))
    for category, x, y in zip(categories, eva_metrics_overall[columns[0]], eva_metrics_overall[columns[1]]):
        ax.scatter(x, y, label=category, color=method_colors[category], s=50)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.autoscale(True)
    handles = [plt.Line2D([0], [0], color=method_colors[cat], marker='o', lw=0) for cat in categories]
    legend = ax.legend(handles=handles, labels=categories, loc='center left', bbox_to_anchor=(1.05, 0.5),
                       title="Methods", frameon=False)
    fig.text(0.36, 0.02, 'Omics mixing', ha='center', fontsize=14)
    fig.text(0.02, 0.5, 'Cell Type conservation', va='center', rotation='vertical', fontsize=14)
    plt.tight_layout(rect=[0.05, 0.05, 0.85, 1])
    plt.show()
    fig.savefig(os.path.join(vis_dir, "metric_overall.png"), dpi=300, bbox_inches='tight')


    eva_metrics_add = eva_metrics.iloc[:, 2:]
    columns = eva_metrics_add.columns
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    for ax, column in zip(axs.flat, columns):
        values = eva_metrics_add[column]
        if np.any(np.isnan(values)):
            ylim_max = 1
        else:
            max_value = values.max()
            import math
            ylim_max = math.ceil(max_value / 0.1) * 0.1
        bars = ax.bar(categories, values, color=[method_colors[cat] for cat in categories])
        ax.set_ylim(0, ylim_max)
        ax.set_xticks([])  # Remove x-axis labels
    handles = [plt.Line2D([0], [0], color=method_colors[cat], lw=4) for cat in categories]
    legend = axs[1].legend(handles=handles, labels=categories, loc='center left', bbox_to_anchor=(1.05, 0.5),
                           title="Methods", frameon=False)
    fig.text(0.02, 0.5, 'Values', va='center', rotation='vertical', fontsize=14)
    plt.tight_layout(rect=[0.05, 0.05, 0.85, 1])
    plt.show()
    fig.savefig(os.path.join(vis_dir, "metric_add.png"), dpi=300, bbox_inches='tight')


def umap_plot_colors(adata, dataset_name, omic_colors, cell_type_colors, vis_dir):

    vis_dir_umap = vis_dir + '/Umap_truth_colors'
    if not os.path.exists(vis_dir_umap):
        os.makedirs(vis_dir_umap)

    obsm_names = adata.obsm_keys()


    for i, method in enumerate(obsm_names):
        sc.pp.neighbors(adata, use_rep=method)
        sc.tl.umap(adata)


        fig, ax = plt.subplots(figsize=(16, 7))
        ax.axis('off')

        inner_ax1 = ax.inset_axes([0.05, 0.1, 0.4, 0.8])
        inner_ax2 = ax.inset_axes([0.55, 0.1, 0.4, 0.8])
        inner_ax1.axis('off')
        inner_ax2.axis('off')

        fig_left = sc.pl.umap(adata, size=80, color='omic_id', palette=omic_colors, title='', ax=inner_ax1, show=False)
        fig_right = sc.pl.umap(adata, size=100, color="cell_type", palette=cell_type_colors, title='', ax=inner_ax2, show=False)


        plt.tight_layout(pad=3.0)
        plt.show()
        fig.savefig(os.path.join(vis_dir_umap, method + ".png"), dpi=300, bbox_inches='tight')
        plt.close(fig)


def paga_plot_colors(adata, obsm_names, cell_type_colors, vis_dir, dataset_name, threshold=0.03):

    vis_dir_paga_trajectory = vis_dir + '/paga_trajectory_colors' + '_' + str(threshold)
    if not os.path.exists(vis_dir_paga_trajectory):
        os.makedirs(vis_dir_paga_trajectory)

    if dataset_name == 'BMCITE_s1d1_s1d2' or dataset_name == 'BMCITE_s1d2_s3d7':
        cell_types = np.array(adata.obs['cell_type2'].array)
        adata.uns['cell_type_colors'] = cell_type_colors
        adata.obs['ground_truth'] = adata.obs['cell_type2'].copy()
    else:
        cell_types = np.array(adata.obs['cell_type'].array)
        adata.uns['cell_type_colors'] = cell_type_colors
        adata.obs['ground_truth'] = adata.obs['cell_type'].copy()

    for method in obsm_names:

        sc.pp.neighbors(adata, use_rep=method)
        sc.tl.umap(adata)
        sc.tl.paga(adata, groups='ground_truth')

        fig, ax = plt.subplots(figsize=(30, 20))
        sc.pl.paga(adata, threshold=threshold, labels=None, show=False, ax=ax, node_size_scale=10, edge_width_scale=2)
        ax.axis('off')
        for artist in ax.get_children():
            if isinstance(artist, plt.Text):
                artist.set_visible(False)
        handles = [plt.Line2D([0], [0], marker='o', color='w', label=cell_type, markersize=10, markerfacecolor=color)
                   for cell_type, color in zip(cell_types, cell_type_colors)]
        legend = ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1, 0.5), title='Cell Types')
        legend.get_frame().set_linewidth(0)
        plt.setp(legend.get_texts(), fontsize=28)
        legend.get_title().set_fontsize(30)

        plt.show()
        plt.savefig(os.path.join(vis_dir_paga_trajectory, method + ".png"), dpi=300, bbox_inches="tight", pad_inches=0)


def calculate_cost_matrix(C, n_clusters):
    cost_matrix = np.zeros((n_clusters, n_clusters))

    # cost_matrix[i,j] will be the cost of assigning cluster i to label j
    for j in range(n_clusters):
        s = np.sum(C[:, j])  # number of examples in cluster i
        for i in range(n_clusters):
            t = C[i, j]
            cost_matrix[j, i] = s - t
    return cost_matrix
def get_cluster_labels_from_indices(indices):
    n_clusters = len(indices)
    clusterLabels = np.zeros(n_clusters)
    for i in range(n_clusters):
        clusterLabels[i] = indices[i][1]
    return clusterLabels
def get_y_preds(y_true, cluster_assignments, n_clusters):
    """Computes the predicted labels, where label assignments now
        correspond to the actual labels in y_true (as estimated by Munkres)

        Args:
            cluster_assignments: array of labels, outputted by kmeans
            y_true:              true labels
            n_clusters:          number of clusters in the dataset

        Returns:
            a tuple containing the accuracy and confusion matrix,
                in that order
    """
    confusion_matrix = metrics.confusion_matrix(y_true, cluster_assignments, labels=None)
    # compute accuracy based on optimal 1:1 assignment of clusters to labels
    cost_matrix = calculate_cost_matrix(confusion_matrix, n_clusters)
    indices = Munkres().compute(cost_matrix)
    kmeans_to_true_cluster_labels = get_cluster_labels_from_indices(indices)

    if np.min(cluster_assignments) != 0:
        cluster_assignments = cluster_assignments - np.min(cluster_assignments)
    y_pred = kmeans_to_true_cluster_labels[cluster_assignments]
    return y_pred
def classification_metric(y_true, y_pred, average='macro', verbose=True, decimals=4):
    """Get classification metric"""
    # confusion matrix
    confusion_matrix = metrics.confusion_matrix(y_true, y_pred)
    # ACC
    accuracy = metrics.accuracy_score(y_true, y_pred)
    accuracy = np.round(accuracy, decimals)

    # precision
    precision = metrics.precision_score(y_true, y_pred, average=average)
    precision = np.round(precision, decimals)

    # recall
    recall = metrics.recall_score(y_true, y_pred, average=average)
    recall = np.round(recall, decimals)

    # F-score
    f_score = metrics.f1_score(y_true, y_pred, average=average)
    f_score = np.round(f_score, decimals)

    return {'accuracy': accuracy, 'precision': precision, 'recall': recall, 'f_measure': f_score}, confusion_matrix

def clustering_metric(y_true, y_pred, n_clusters, verbose=True, decimals=4):
    """Get clustering metric"""
    y_pred_ajusted = get_y_preds(y_true, y_pred, n_clusters)

    classification_metrics, confusion_matrix = classification_metric(y_true, y_pred_ajusted)

    # AMI
    ami = metrics.adjusted_mutual_info_score(y_true, y_pred)
    ami = np.round(ami, decimals)
    # NMI
    nmi = metrics.normalized_mutual_info_score(y_true, y_pred)
    nmi = np.round(nmi, decimals)
    # ARI
    ari = metrics.adjusted_rand_score(y_true, y_pred)
    ari = np.round(ari, decimals)

    return dict({'AMI': ami, 'NMI': nmi, 'ARI': ari}, **classification_metrics), confusion_matrix

def confusion_plot(adata, obsm_names_sel, CM, accuracy, dataset_name, vis_dir, keys = 'confusion_plot'):

    vis_dir_confusion = vis_dir + '/' + keys
    if not os.path.exists(vis_dir_confusion):
        os.makedirs(vis_dir_confusion)

    if dataset_name == 'BMCITE_s1d1_s1d2' or dataset_name == 'BMCITE_s1d2_s3d7':
        ground_truth = adata.obs['cell_type2']
    else:
        ground_truth = adata.obs['cell_type']

    for method in obsm_names_sel:
        cm = CM[method]
        cm = cm / cm.sum(axis=1, keepdims=True)
        fig, ax = plt.subplots(1, 1, figsize=(30, 20))
        labels = np.unique(ground_truth)

        ax = heatmap(cm, ax=ax, xticklabels=labels, yticklabels=labels, cmap='Greens') #'BuGn'
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=40)
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=40)
        ax.set_title(f'{method} (Accuracy: {accuracy[method]:.2f})', fontsize=50)
        cbar = plt.gca().collections[0].colorbar
        cbar.ax.tick_params(labelsize=40)
        plt.show()
        fig.savefig(os.path.join(vis_dir_confusion, method + ".png"), dpi=300, bbox_inches='tight')




