
from sklearn.cluster import KMeans
import anndata
from plot_func import *


dataset_name = "PBMC_paired"
dataset_type = 'RNA_ATAC'
GAM_name = 'ArchR'
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

adata = anndata.read_h5ad(result_dir + '/adata.h5ad')
adata_gasm = anndata.read_h5ad(result_dir + '/adata_gasm.h5ad')
annotations = np.load(result_dir + "/annotations.npz", allow_pickle=True)
anno_rna = annotations['anno_rna']
anno_other = annotations['anno_other']
clu_rna = annotations['clu_rna']
clu_other = annotations['clu_other']
anno_rna_gasm = annotations['anno_rna_gasm']
anno_other_gasm = annotations['anno_other_gasm']
clu_rna_gasm = annotations['clu_rna_gasm']
clu_other_gasm = annotations['clu_other_gasm']


obsm_names = adata.obsm_keys()
obsm_names_gasm = adata_gasm.obsm_keys()

n_clusters = len(np.unique(adata.obs['cell_type']))
cell_types = np.unique(adata.obs['cell_type'])



######################################################################################################

omic_colors = ['#92a5d1', '#feb29b']

cell_type_colors_all = ['#f7b3ac', '#e9212c', '#f5c4db', '#d998b6', '#a61e4d',
          '#bdb5e1', '#a879b8', '#f7e16f',
          '#ff8831', '#958431', '#b0d992', '#82b7a3', '#008d00',
          '#27483e', '#b4d7f0', '#5591dc',
          '#919191', '#2017a5', '#ab888d']

cell_type_colors = cell_type_colors_all[:n_clusters]

adata_ = anndata.AnnData(adata_gasm.obsm['raw_data'])
adata_.obsm['raw_legend'] = adata_.X
adata_.obs = adata_gasm.obs
umap_plot_colors(adata_, dataset_name, omic_colors, cell_type_colors, vis_dir)

umap_plot_colors(adata, dataset_name, omic_colors, cell_type_colors, vis_dir)
umap_plot_colors(adata_gasm, dataset_name, omic_colors, cell_type_colors, vis_dir)

######################################################################################################

###########################################################################################

accuracy = {}
CluMe = {}
CM = {}
for i, method in enumerate(obsm_names):
    kmeans = KMeans(n_clusters=len(np.unique(np.hstack((anno_rna, anno_other)))), random_state=0).fit(
        adata.obsm[method])
    clustering = kmeans.labels_
    Y = np.hstack((clu_rna, clu_other))
    y_preds = get_y_preds(Y, clustering, len(np.unique(Y)))
    scores = clustering_metric(Y, clustering, len(np.unique(Y)))
    accuracy[method] = scores[0]['accuracy']
    CM[method] = scores[1]


for i, method in enumerate(obsm_names_gasm):
    kmeans = KMeans(n_clusters=len(np.unique(np.hstack((anno_rna_gasm, anno_other_gasm)))), random_state=0).fit(
        adata_gasm.obsm[method])
    clustering = kmeans.labels_
    Y = np.hstack((clu_rna_gasm, clu_other_gasm))
    y_preds = get_y_preds(Y, clustering, len(np.unique(Y)))
    scores = clustering_metric(Y, clustering, len(np.unique(Y)))
    accuracy[method] = scores[0]['accuracy']
    CluMe[method] = scores[0]
    CM[method] = scores[1]

eva_cm = dict({"CluMe": CluMe, 'CM': CM, 'accuracy': accuracy})

np.save(os.path.join(vis_dir, 'eva_cm_kmeans.npy'), eva_cm)

###########################################################################################

confusion_plot(adata, obsm_names, CM, accuracy, dataset_name, vis_dir, keys = 'confusion_plot')
confusion_plot(adata_gasm, obsm_names_gasm, CM, accuracy, dataset_name, vis_dir, keys='confusion_plot')

###########################################################################################

paga_plot_colors(adata, obsm_names, cell_type_colors, vis_dir, dataset_name)
paga_plot_colors(adata_gasm, obsm_names_gasm, cell_type_colors, vis_dir, dataset_name)

###########################################################################################

eva_metrics = pd.read_csv(eva_dir + '/eva_metrics_' + GAM_name + '.csv', index_col=0, header=0).T

obsm_names_metrics = []
for i, method in enumerate(list(eva_metrics.index)):
    if len(method.split('_'))==2 and method != 'raw_data':
        obsm_names_metrics.append(method.split('_')[0])
    else:
        obsm_names_metrics.append(method)

colors = ["#E41A1C", "#FF7F00", "#DAA520", "#FFFF33", "#66C2A5", "#008000", "#808000", "#2F4F4F", '#1E90FF',
          "#377EB8", "#0000FF", "#DB7093", "#984EA3", "#A65628", "#999999"]

method_colors = {}
for i in range(len(obsm_names_metrics)):
    method_type_tmp = obsm_names_metrics[i]
    color_tmp = colors[i]
    method_colors[method_type_tmp] = color_tmp


metrics_plot(eva_metrics, method_colors, vis_dir)











