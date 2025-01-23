import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import LabelEncoder
import anndata
from plot_func import *
from scipy.stats import pearsonr
def load_clusters(anno_other, anno_rna):

    all_cell_types = list(set(anno_other.tolist() + anno_rna.tolist()))
    label_encoder = LabelEncoder()
    label_encoder.fit_transform(all_cell_types)
    clu_other = label_encoder.transform(anno_other)
    clu_rna = label_encoder.transform(anno_rna)

    return clu_rna, clu_other

# dataset_name = "BMCITE_s1d1_s1d2"
# dataset_batch = 's1d1_s1d2'
# res_RNA = [1.38, 2.5, 1.4]
# res_Other = [1, 1.1, 1.1]

dataset_name = "BMCITE_s1d2_s3d7"
dataset_batch = 's1d2_s3d7'
res_RNA = [1.2, 2.8, 1.52]
res_Other = [1.1, 1.1, 0.9]
dataset_type = 'RNA_ADT'
paired = True

path = "../"
dataset_dir = path + '../data/'
result_dir = path + "results/" + dataset_name
eva_dir = path + "eva/" + dataset_name
vis_dir = path + 'vis/' + dataset_name
if not os.path.exists(vis_dir):
    os.makedirs(vis_dir)
if not os.path.exists(eva_dir):
    os.makedirs(eva_dir)

RNA_data = anndata.read(os.path.join(dataset_dir, dataset_name + '/bm.rna.' + dataset_batch + '.h5ad'))
Other_data = anndata.read(os.path.join(dataset_dir, dataset_name + '/bm.adt.' + dataset_batch + '.h5ad'))
Other_data.obs['omic_id'] = 'Other'
RNA_data.obs['omic_id'] = 'RNA'

fea_pairs = pd.read_csv(os.path.join(dataset_dir, dataset_name + '/graph.csv'))
other_hvg = Other_data.var_names
rna_hvg_1 = set(RNA_data.var_names[np.where(RNA_data.var['vst.variable'] == 1)])
rna_hvg_2 = set(fea_pairs.values[:, 1])
rna_hvg = list(rna_hvg_1.union(rna_hvg_2))

other_hvg_2 = list(fea_pairs.values[:, 2])

RNA_data = RNA_data[:, rna_hvg]
RNA_data.obsm['encoding_genes'] = RNA_data[:, list(rna_hvg_2)].X
RNA_data.obsm['Modality_RNA'] = RNA_data.X
Other_data = Other_data[:, other_hvg]
Other_data.obsm['encoded_proteins'] = Other_data[:, other_hvg_2].X
Other_data.obsm['Modality_Other'] = Other_data.X
RNA_data.obs['cell_type'], Other_data.obs['cell_type'] = RNA_data.obs['cell_type2'], Other_data.obs['cell_type2']

anno_rna, anno_other = RNA_data.obs['cell_type2'], Other_data.obs['cell_type2']
clu_rna, clu_other = load_clusters(anno_other, anno_rna)
num_other, num_rna = len(clu_other), len(clu_rna)

n_clusters = len(np.unique(clu_rna))


#########################################################################################################

print('correlation between cells')
other_names = fea_pairs.values[:, 2]
rna_names = fea_pairs.values[:, 1]

data_other = Other_data[: ,other_names].X.toarray()
data_rna = RNA_data[:, rna_names].X.toarray()

pearson_cell = np.zeros((len(clu_rna), 1))
for p in range(len(clu_rna)):
    pearson_cell[p, 0], _ = pearsonr(data_other[p,:], data_rna[p,:])
pearson_cell_avg = np.mean(pearson_cell)
pearson_cell_clean = pearson_cell[~np.isnan(pearson_cell)]
print('correlation between gene')
pearson_gene = np.zeros((len(rna_names), 1))
for q in range(len(rna_names)):
    pearson_gene[q, 0], _ = pearsonr(data_other[:, q], data_rna[:, q])
pearson_gene_avg = np.mean(pearson_gene)

data_to_plot = [pearson_cell_clean.flatten(), pearson_gene.flatten()]

# 绘制箱线图
plt.figure(figsize=(8, 6))
plt.boxplot(data_to_plot, labels=["Cell-wise", "Gene-wise"], showmeans=True)

# 添加标题和轴标签
plt.title("Distribution of Pearson Correlation")
plt.ylabel("Pearson Correlation Coefficient")
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.show()

###########################################################################################################
eva_cm = np.load(os.path.join(vis_dir, 'eva_cm_louvain_modality.npy'), allow_pickle=True).item()

def cluster_obsm(adata, res):
    obsm_names = adata.obsm_keys()
    accuracy = {}
    CM = {}
    CluMe = {}
    eval_resolution = {}
    for i, method in enumerate(obsm_names):
        print(method)
        adata_ = anndata.AnnData(adata.obsm[method])
        sc.pp.neighbors(adata_)
        eval_resolution[method] = res[i]
        sc.tl.louvain(adata_, key_added=method+"_louvain", resolution=eval_resolution[method])
        clustering = adata_.obs[method+"_louvain"].values.astype('int')

        Y = clu_rna
        y_preds = get_y_preds(Y, clustering, len(np.unique(Y))).astype('int')
        scores = clustering_metric(Y, y_preds, len(np.unique(Y)))
        accuracy[method] = scores[0]['accuracy']
        CluMe[method] = scores[0]
        CM[method] = scores[1]

    eva_cm = dict({"CluMe": CluMe, "res": eval_resolution, 'CM': CM})
    return eva_cm


def cluster_obsm_kmeans(adata, res):
    obsm_names = adata.obsm_keys()
    accuracy = {}
    CM = {}
    CluMe = {}
    eval_resolution = {}
    for i, method in enumerate(obsm_names):
        kmeans = KMeans(n_clusters=len(np.unique(anno_rna)), random_state=0).fit(
            adata.obsm[method])
        clustering = kmeans.labels_

        Y = clu_rna
        y_preds = get_y_preds(Y, clustering, len(np.unique(Y))).astype('int')
        scores = clustering_metric(Y, y_preds, len(np.unique(Y)))
        accuracy[method] = scores[0]['accuracy']
        CluMe[method] = scores[0]
        CM[method] = scores[1]

    eva_cm = dict({"CluMe": CluMe, "res": eval_resolution, 'CM': CM})
    return eva_cm


eva_cm_RNA = cluster_obsm_kmeans(RNA_data, res_RNA)
eva_cm_Other = cluster_obsm_kmeans(Other_data, res_Other)

eva_cm = dict({'eva_cm_Other': eva_cm_Other, 'eva_cm_RNA':eva_cm_RNA})
np.save(os.path.join(vis_dir, 'eva_cm_kmeans_modality.npy'), eva_cm)


omic_colors = ['#92a5d1', '#feb29b']

cell_type_colors_all = ['#f7b3ac', '#e9212c', '#f5c4db', '#d998b6', '#a61e4d',
          '#bdb5e1', '#a879b8', '#fdeeba', '#f7e16f', '#e59c58',
          '#ff8831', '#c65323', '#958431', '#b0d992', '#82b7a3', '#008d00',
          '#27483e', '#b4d7f0', '#5591dc', '#377483', '#a415bb',
          '#919191', '#5b534d', '#2017a5', '#ab888d']

cell_type_colors = cell_type_colors_all[:n_clusters]

umap_plot_colors(RNA_data, dataset_name, omic_colors, cell_type_colors, vis_dir)
umap_plot_colors(Other_data, dataset_name, omic_colors, cell_type_colors, vis_dir)
######################################################################################################

