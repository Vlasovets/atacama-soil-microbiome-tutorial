import pandas as pd
import numpy as np
import re
import bisect

import seaborn as sns
import networkx as nx

import plotly.express as px
import plotly.figure_factory as ff
import plotly.graph_objs as go

import scipy
import scipy.cluster.hierarchy as sch
from scipy import stats
from scipy.spatial import distance

from sklearn import preprocessing

from math import pi
from itertools import chain

from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
from pyvis.network import Network

from bokeh.io import output_notebook, show, save
from bokeh.models import Range1d, Circle, ColumnDataSource, MultiLine, HoverTool, LabelSet, PointDrawTool
from bokeh.plotting import figure
from bokeh.plotting import from_networkx
from bokeh.palettes import RdBu, Blues8
from bokeh.models import HoverTool, Panel, Tabs, ColorBar, LinearColorMapper
from bokeh.layouts import row

def save_dataframe(data, filename):
    """
    data: pd.Dataframe
    filename: str
    """
    data.to_csv(filename, index=True)
    return data

def load_data(filename):
    """
    filename: .tsv file
    """
    if ".tsv" in filename:
        data = pd.read_csv(filename, sep="\t", index_col=0)
    elif ".csv" in filename:
        data = pd.read_csv(filename, index_col=0)
    return data

def filter_zero_features(data: pd.DataFrame, threshold: float=0.9):
    """
    data: pd.Dataframe
    (p, N) count table where p - features, N -samples
    
    threshold: float
    percentage of zeros will be filtered out
    """
    # Transpose the DataFrame
    data_T = data.transpose()

    # Calculate zero-inflation per ASV across all samples
    zero_perc = (data_T == 0).mean()

    # Filter columns with zero-inflation > 90%
    filtered_cols = zero_perc[zero_perc <= threshold].index
    data_filt = data_T[filtered_cols].transpose()

    # Check if any columns contain only zeros
    any_zero_cols = (data_filt == 0).all().any()

    # Drop columns with all zeros
    zero_cols = data_filt.columns[(data_filt == 0).all()]
    filtered_data = data_filt.drop(zero_cols, axis=1)

    # Check if any columns contain only zeros after filtering
    any_zero_cols_final = (filtered_data == 0).all().any()

    return filtered_data


def filter_zero_samples(data, threshold=0.95):
    """
    data: pd.Dataframe
    (p, N) count table where p - features, N -samples
    
    threshold: float
    percentage of zeros will be filtered out
    """
    # Calculate zero-inflation per sample  across all ASVs
    zero_perc = (data == 0).mean()
    
    # Filter columns with zero-inflation > 95%
    high_zero_cols = zero_perc[zero_perc > threshold].index
    filtered_data = data.drop(columns=high_zero_cols)
    
    return filtered_data


def update_index(filtered_data, processed_taxa):
    # Create a dictionary mapping genus to species
    genus = processed_taxa['genus'].to_dict()

    # Map the index values using the species dictionary
    filtered_data.index = filtered_data.index.map(genus)

    # Fill missing values with 'unknown species'
    filtered_data.index = filtered_data.index.fillna('unknown genus')

    # Check and update index values with length less than 5
    taxa_names = []
    i = 1
    for idx in filtered_data.index:
        if idx == "unknown genus":
            new_idx = f"ASV{i}"
            i += 1
        else:
            new_idx = idx
        taxa_names.append(new_idx)

    filtered_data.index = taxa_names

    return filtered_data


def process_taxonomy(taxa: pd.DataFrame):
    """
    taxa: pd.Dataframe
    Taxonomy information associated with count data.
    """
    taxonomy_levels = {
        "domain": '^d__',
        "phylum": '^p__',
        "class": '^c__',
        "order": '^o__',
        "family": '^f__',
        "genus": '^g__',
        "species": '^s__'
    }

    # Split taxonomic ranks into different columns
    taxa_sep = taxa['Taxon'].str.split(';', expand=True)

    # Rename taxonomic ranks with full names
    taxa_sep.columns = taxonomy_levels.keys()

    # Drop missing species
    taxa_sep = taxa_sep[taxa_sep['species'].notnull()]

    # Remove blank spaces from taxonomic ranks
    taxa_sep = taxa_sep.apply(lambda x: x.str.strip())

    # Subtract "s" from the string names
    taxa_sep['species'] = taxa_sep['genus'] + taxa_sep['species'].str.lstrip('s')

    return taxa_sep


def clean_meta_data(meta_data):
    """
    Clean and preprocess the meta data.
    
    Parameters:
        meta_data (pandas.DataFrame): The meta data to be cleaned.
        
    Returns:
        pandas.DataFrame: The cleaned meta data.
    """
    # Select only numeric features
    meta_data = meta_data.loc[:, meta_data.iloc[0, :] != 'categorical']
    meta_data = meta_data.apply(pd.to_numeric, errors='coerce')

    # Drop QIIME2 header
    meta_data = meta_data.iloc[1:]

    # Fill missing values with zeros
    meta_data = meta_data.fillna(0)
    
    return meta_data


def select_covariates(meta_data, selected_columns):
    """
    Select specific columns as covariates from the meta data.
    
    Parameters:
        meta_data (pandas.DataFrame): The meta data.
        selected_columns (list): The list of columns to be selected as covariates.
        
    Returns:
        pandas.DataFrame: The meta data with selected covariates.
    """
    # Select most interesting columns
    meta_data = meta_data[selected_columns]
    return meta_data


def merge_data(counts, covariates):
    """
    Join dataframes 'counts' and 'covariates' based on their indices.

    Args:
        counts (pandas.DataFrame): The dataframe to be transposed and joined.
        covariates (pandas.DataFrame): The dataframe to be joined.

    Returns:
        pandas.DataFrame: The joined dataframe.

    """
    counts_transposed = counts.transpose()
    df = counts_transposed.merge(covariates, left_index=True, right_index=True)
    return df


def calculate_covariance(data, n_cov: int, method: str=None):
    """
    Calculate the covariance matrix of the data after excluding scaled covariates.

    Args:
        data (pd.DataFrame): The input dataframe containing the data.
        n_cov (int): The number of covariates.
        method (str): The method to use for covariance calculation.

    Returns:
        np.ndarray: The covariance matrix of the data.

    """
    if n_cov is None:
        asv = data
    else:
        asv = data.iloc[:, :-n_cov]
    asv_names = asv.columns
    S = np.cov(asv.T.values, bias=True)
    res = pd.DataFrame(S, columns=asv_names, index=asv_names)
    
    if method == "corr":
        corr = scale_array_by_diagonal(S)
        res = pd.DataFrame(corr, columns=asv_names, index=asv_names)
        
    elif method == "latent":
        clean_types_asv = get_tps(asv)
        print("Estimating latent correlation with datatypes: {0}".format(clean_types_asv))
        lat_cor_asv = latentcor(asv, tps = clean_types_asv, method ='original', use_nearPD=False)
        res = lat_cor_asv["R"]
    
    return res


def scale_meta_data(meta_data):
    """
    Scale the meta data using StandardScaler.
    
    Parameters:
        meta_data (pandas.DataFrame): The meta data to be scaled.
        
    Returns:
        pandas.DataFrame: The scaled meta data.
    """
    scaler = preprocessing.StandardScaler().fit(meta_data)
    meta_scaled = scaler.transform(meta_data)
    meta_scaled = pd.DataFrame(meta_scaled, index=meta_data.index, columns=meta_data.columns)
    return meta_scaled


def create_lambda_mask(counts, p, p_meta):
    """
    Create a lambda matrix mask based on the given parameters.

    Args:
        p_meta (int): Number of rows and columns for the meta block.
        p (int): Number of rows and columns for the ASV block.
        counts (pd.DataFrame): DataFrame used for column names of the mask.

    Returns:
        np.ndarray: Lambda matrix mask.
    """
    shape_meta = (p_meta, p_meta)
    mask = np.zeros(shape_meta)
    mask = mask + 0.01
    asv_block = np.ones((p, p))
    mask[0:p, 0:p] += asv_block - 0.01
    return mask


def create_label_dict(df):
    """
    Creates a dictionary to map column indices to column labels and another 
    dictionary to map column indices to reversed column labels.
    Args:
        df (pd.DataFrame): The input DataFrame.

    Returns:
        tuple: A tuple containing two dictionaries. The first dictionary maps column indices to column labels, 
        and the second dictionary maps column indices to reversed column labels.
    """
    n_labels = len(df.columns)
    labels_dict = dict(zip(range(n_labels), df.columns))
    labels_dict_reversed = dict(zip(range(n_labels),list(labels_dict.values())[::-1]))
    
    return labels_dict, labels_dict_reversed


def _get_order(data: pd.DataFrame, method: str = 'average', metric: str = 'euclidean'):
    """
    Performs hierarchical clustering on the input DataFrame and returns the cluster order.

    Args:
        data (pd.DataFrame): The input DataFrame.
        method (str, optional): The clustering method. Defaults to 'average'.
        metric (str, optional): The distance metric. Defaults to 'euclidean'.

    Returns:
        list: The cluster order as a list of indices.
    """
    grid = sns.clustermap(data, method=method, metric=metric, robust=True)
    plt.close()
    clust_order = grid.dendrogram_row.reordered_ind

    return clust_order


def hierarchical_clustering(data: pd.DataFrame, order: list, n_covariates: int = None):
    """
    Performs hierarchical clustering on the input DataFrame using 
    the provided cluster order and returns the reorganized DataFrame.

    Args:
        data (pd.DataFrame): The input DataFrame.
        order (list): The cluster order as a list of indices.
        n_covariates (int, optional): The number of covariates. Defaults to None.

    Returns:
        pd.DataFrame: The reorganized DataFrame.
    """
    if n_covariates is None:
        re_data = data.iloc[order, order]

    else:
        asv_part = data.iloc[:-n_covariates, :-n_covariates]
        re_asv_part = asv_part.iloc[order, order]
        cov_asv_part = data.iloc[:-n_covariates, -n_covariates:].iloc[order, :]
        cov_part = data.iloc[-n_covariates:, -n_covariates:]

        res = np.block([[re_asv_part.values, cov_asv_part.values],
                        [cov_asv_part.T.values, cov_part.values]])

        labels = list(re_asv_part.columns) + list(cov_part.columns)
        re_data = pd.DataFrame(res, index=labels, columns=labels)

    return re_data


def plot_ordered_heatmap(data: pd.DataFrame, order: list(), n_covariates: int=None,
                 title: str=" ", width: int=1500, height: int=1500, label_size: str="16pt"):
    """
    Plots a heatmap using the input DataFrame and additional parameters.

    Args:
        data (pd.DataFrame): The input DataFrame.
        order: The order as a list of indices.
        n_covariates (int, optional): The number of covariates. Defaults to None.
        title (str, optional): The title of the heatmap. Defaults to " ".
        width (int, optional): The width of the heatmap in pixels. Defaults to 1500.
        height (int, optional): The height of the heatmap in pixels. Defaults to 1500.
        label_size (str, optional): The size of the labels. Defaults to "16pt".

    Returns:
        [type]: The plotted heatmap.
    """
    clust_data = hierarchical_clustering(data, order=order, n_covariates=n_covariates)
    lables_clust, re_labels_clust = create_label_dict(clust_data)

    p = _make_heatmap(data=clust_data, title=title, width=width, height=height, label_size=label_size,
                      labels_dict=lables_clust, labels_dict_reversed=re_labels_clust)
    return p


def create_network_visualization(G_adapt, height: int = 1500, width: int = 1800, show_labels: bool = False, size_degree: bool = False,
                                 scale_edge: int = 2, scale_node: int = 1):
    """
    Creates a network visualization using Pyvis library and returns the Pyvis Network object.

    Parameters:
        G_adapt (networkx.Graph): The graph data loaded from NetworkX.
        height (int): The height of the network visualization in pixels. Default is 1500.
        width (int): The width of the network visualization in pixels. Default is 1800.
        show_labels (bool): Whether to show edge labels. Default is False.
        size_degree (bool): Whether to adjust node sizes based on degree. Default is False.
        scale_edge (int): Scaling factor for edge widths. Default is 2.
        scale_node (int): Scaling factor for node sizes when using size_degree. Default is 1.

    Returns:
        pyvis.network.Network: The Pyvis Network object representing the network visualization.
    """
    # Create a Pyvis network
    net = Network(height=f"{height}px", width=f"{width}px", directed=False, cdn_resources='in_line', notebook=True)

    # Load the graph data from NetworkX
    net.from_nx(G_adapt)

    # Disable physics simulation
    net.toggle_physics(False)

    if size_degree:
        # Calculate the degree of each node
        degrees = dict(G_adapt.degree())

        # Set node sizes based on degree
        for node in net.nodes:
            node_id = node['id']
            if node_id in degrees:
                node['size'] = degrees[node_id] * scale_node  # Adjust the scaling factor as needed

    # Set edge and node styles
    for edge in net.edges:
        edge['width'] = abs(edge['covariance']) * scale_edge
        edge['length'] = 2000

        if show_labels:
            edge['label'] = str(round(edge['covariance'], 2))
        if edge['covariance'] < 0:
            edge['color'] = '#f1ac8b' #red
            if show_labels:
                edge['font'] = {'multi': 'true', 'size': 15, 'color': 'blue', 'face': 'arial', 'align': 'top'}
        else:
            edge['color'] = "#abe4ff" #blue
            if show_labels:
                edge['font'] = {'multi': 'true', 'size': 15, 'color': 'red', 'face': 'arial', 'align': 'top'}

    for node in net.nodes:
        if "ASV" in node['label'] or "g_" in node['label']:
            node['color'] = '#610053'
        else:
            node['color'] = '#FA4665'

        node['font'] = {'multi': 'true', 'size': 20, 'color': 'black', 'face': 'arial', 'align': 'left'}

    return net



def normalize(X):
    """
    transforms to the simplex
    X should be of a pd.DataFrame of form (p,N)
    """
    return X / X.sum(axis=0)


def geometric_mean(x, positive=False):
    """
    calculates the geometric mean of a vector
    """
    assert not np.all(x == 0)

    if positive:
        x = x[x > 0]
    a = np.log(x)
    g = np.exp(a.sum() / len(a))
    return g


def log_transform(X, transformation=str, eps=0.1):
    """
    log transform, scaled with geometric mean
    X should be a pd.DataFrame of form (p,N)
    """
    if transformation == "clr":
        assert not np.any(X.values == 0), "Add pseudo count before using clr"
        g = X.apply(geometric_mean)
        Z = np.log(X / g)
    elif transformation == "mclr":
        g = X.apply(geometric_mean, positive=True)
        X_pos = X[X > 0]
        Z = np.log(X_pos / g)
        Z = Z + abs(np.nanmin(Z.values)) + eps
        Z = Z.fillna(0)
    return Z


def zero_imputation(df: pd.DataFrame, pseudo_count: int = 1):
    X = df.copy()
    original_sum = X.sum(axis=0) # sum per row (sample)
    for col in X.columns:
        X[col].replace(to_replace=0, value=pseudo_count, inplace=True)
    shifted_sum = X.sum(axis=0) # sum per row (sample)
    scaling_parameter = original_sum.div(shifted_sum)
    X = X.mul(scaling_parameter, axis =1) # multiply by column

    return X


def transform_features(X: pd.DataFrame, transformation: str = "clr", pseudo_count: int = 1) -> pd.DataFrame:
    """
    Project compositional data to Euclidean space.

    Parameters
    ----------
    pseudo_count: int, optional
        Add pseudo count, only necessary for transformation = "clr".
    table: biom.Table
        A table with count microbiome data.
    transformation: str
        If 'clr' the data is transformed with center log-ratio method by Aitchison (1982).
        If 'mclr' the data is transformed with modified center log-ratio method by Yoon et al. (2019).

    Returns
    -------
    X: pd.Dataframe
        Count data projected to Euclidean space.

    """
    columns = X.columns

    if transformation == "clr":
        X = zero_imputation(X, pseudo_count=pseudo_count)
        X = normalize(X)
        X = log_transform(X, transformation=transformation)

        return pd.DataFrame(X, columns=columns)

    elif transformation == "mclr":
        X = normalize(X)
        X = log_transform(X, transformation=transformation)

        return pd.DataFrame(X, columns=columns)

    else:
        raise ValueError(
            "Unknown transformation name, use clr and not %r" % transformation
        )
        

def PCA(X, L, inverse=True):
    """
    Perform Principal Component Analysis (PCA) on the given data.

    Parameters:
        X (pandas.DataFrame): Transformed data matrix.
        L (numpy.ndarray): Low-rank matrix.
        inverse (bool, optional): Determines whether to compute inverse PCA or regular PCA. Defaults to True.

    Returns:
        tuple: A tuple containing the following:
            - zu (numpy.ndarray): The projected data matrix.
            - loadings (numpy.ndarray): The loadings matrix.
            - sig (numpy.ndarray): The rounded singular values.

    """
    sig, V = np.linalg.eigh(L)

    # sort eigenvalues in descending order
    sig = sig[::-1]
    V = V[:, ::-1]

    ind = np.argwhere(sig > 1e-9)

    if inverse:
        loadings = V[:, ind] @ np.diag(np.sqrt(1 / sig[ind]))
    else:
        loadings = V[:, ind] @ np.diag(np.sqrt(sig[ind]))

    # compute the projection
    zu = X.values @ loadings

    return zu, loadings, np.round(sig[ind].squeeze(), 3)



def scale_array_by_diagonal(X, d = None):
    """
    scales a 2d-array X with 1/sqrt(d), i.e.
    
    X_ij/sqrt(d_i*d_j)
    in matrix notation: W^-1 @ X @ W^-1 with W^2 = diag(d)
    
    if d = None, use square root diagonal, i.e. W^2 = diag(X)
    see (2.4) in https://fan.princeton.edu/papers/09/Covariance.pdf
    """
    assert len(X.shape) == 2
    if d is None:
        d = np.diag(X)
    else:
        assert len(d) == X.shape[0]
        
    scale = np.tile(np.sqrt(d),(X.shape[0],1))
    scale = scale.T * scale
    
    return X/scale



### Visualization

def plot_covariates(unscaled_covariates: pd.DataFrame, scaled_covariates: pd.DataFrame):
    for col in scaled_covariates.columns:
        fig, axes = plt.subplots(1, 2, figsize=(15, 7), sharey=True)

        unscaled_covariates[col].hist(ax=axes[0], color="#610053")
        axes[0].set_xlabel('Unscaled', fontsize=24)
        axes[0].set_ylabel('Samples', fontsize=24)
        axes[0].grid(False)
        axes[0].annotate('Color: #610053', xy=(0.5, 0.9), xycoords='axes fraction', fontsize=12, ha='center', color='white')

        scaled_covariates[col].hist(ax=axes[1], color="#FA4665")
        axes[1].set_xlabel('Scaled', fontsize=24)
        axes[1].grid(False)

        titles = {
            'toc': 'Total organic carbon (μg/g)',
            'ec': 'Electric conductivity (S/m in SI)',
            'average-soil-relative-humidity': 'Average soil humidity (%)',
            'elevation': 'Elevation (m)',
            'average-soil-temperature': 'Average soil temperature (t°)',
            'ph': 'pH'
        }
        subtitle = titles.get(col, col)

        plt.suptitle(subtitle, fontsize=36)
        plt.tight_layout()
        #plt.savefig("plots/{0}.pdf".format(col), format='pdf');
        # plt.close()

def plotly_heatmap(z, x, y, x_label: str, y_label: str, title: str=None, zmin: int=None, zmax: int=None,
            height: int=1200, width: int=1200, colorscale: str='RdPu'):
    # Create a Plotly heatmap using the correlation matrix
    heatmap = go.Heatmap(z=z, x=x, y=y, colorscale=colorscale, zmin = zmin, zmax = zmax)
    
    # Modify the colorbar font size
    heatmap.colorbar = dict(
        title='counts',
        titleside='top',
        titlefont=dict(
            size=24  # Set the desired font size
        ),
        tickfont=dict(
            size=24  # Set the desired font size
        )
    )
    # Create a layout for the heatmap
    layout = go.Layout(title=title, xaxis=dict(title=x_label), yaxis=dict(title=y_label), 
                       height=height, width=width, xaxis_tickangle=45, 
                       xaxis_title_font=dict(size=54), yaxis_title_font=dict(size=54)
    )
        
    # Create a figure object and add the heatmap to it
    fig = go.Figure(data=[heatmap], layout=layout)
    
    
    
    return fig

def _make_heatmap(data: pd.DataFrame(), title: str = None, labels_dict: dict = None,
                  labels_dict_reversed: dict = None,
                  width: int = 1500, height: int = 1500, label_size: str = "5pt",
                  title_size: str = "24pt", not_low_rank: bool = True):
    nlabels = len(labels_dict)
    shifted_labels_dict = {k + 0.5: v for k, v in labels_dict.items()}
    shifted_labels_dict_reversed = {k + 0.5: v for k, v in labels_dict_reversed.items()}

    df = data.iloc[::-1]  # rotate matrix 90 degrees
    df = pd.DataFrame(df.stack(), columns=['covariance']).reset_index()
    df.columns = ["taxa_y", "taxa_x", "covariance"]
    df = df.replace({"taxa_x": labels_dict, "taxa_y": labels_dict})

    color_list, colors = _get_colors(df=df)
    # min_value = df['covariance'].min()
    # max_value = df['covariance'].max()
    # mapper = LinearColorMapper(palette=colors, low=min_value, high=max_value)
    mapper = LinearColorMapper(palette=colors, low=-1, high=1)
    color_bar = ColorBar(color_mapper=mapper, location=(0, 0))

    bottom, top, left, right = _get_bounds(nlabels=nlabels)

    source = ColumnDataSource(
        dict(top=top, bottom=bottom, left=left, right=right, color_list=color_list,
             taxa_x=df['taxa_x'], taxa_y=df['taxa_y'], covariance=df['covariance']))

    bokeh_tools = ["save, zoom_in, zoom_out, wheel_zoom, box_zoom, crosshair, reset, hover"]

    p = figure(plot_width=width, plot_height=height, x_range=(0, nlabels), y_range=(0, nlabels),
               title=title, title_location='above', x_axis_location="below",
               tools=bokeh_tools, toolbar_location='left')

    p.quad(top="top", bottom="bottom", left="left", right="right", line_color='white',
           color="color_list", source=source)
    p.xaxis.major_label_orientation = pi / 4
    p.yaxis.major_label_orientation = "horizontal"
    p.xaxis.major_label_text_font_size = label_size
    p.yaxis.major_label_text_font_size = label_size
    p.title.text_font_size = title_size
    p.add_layout(color_bar, 'right')
    p.toolbar.autohide = True
    p.xaxis.ticker = [x + 0.5 for x in
                      list(range(0, nlabels))]  ### shift label position to the center
    p.yaxis.ticker = [x + 0.5 for x in list(range(0, nlabels))]
    p.xaxis.major_label_overrides = shifted_labels_dict
    p.yaxis.major_label_overrides = shifted_labels_dict_reversed

    hover = p.select(dict(type=HoverTool))
    hover.tooltips = [
        ("taxa_x", "@taxa_x"),
        ("taxa_y", "@taxa_y"),
        ("covariance", "@covariance"),
    ]

    return p


def plot_network(G, title, width, height, node_size=None, amplify_x=10):
    #Establish which categories will appear when hovering over each node
    HOVER_TOOLTIPS = [("Character", "@index")]
    hover = HoverTool(tooltips=[('','@index')])
    tools = ["save, zoom_in, zoom_out, wheel_zoom, box_zoom, crosshair, reset, hover, pan"]

    #Create a plot — set dimensions, toolbar, and title
    plot = figure(tooltips = HOVER_TOOLTIPS, plot_width=width, plot_height=height,
                  tools=tools, active_scroll='wheel_zoom',
                x_range=Range1d(-10.1, 10.1), 
                  y_range=Range1d(-10.1, 10.1), title=title)
    
    
    
    color_map = ["#88CCEE" if "ASV" in j else "#DDCC77" for j in G.nodes()] 
    nx.set_node_attributes(G, {j: {'color': color_map[i]} for i, j in enumerate(G.nodes())})

    if node_size is not None:
        n_degrees = {k: 15*v for k,v in G.degree()} 
        nx.set_node_attributes(G, n_degrees, 'node_size')
        node_size = 'node_size'
    else:
        node_size = 40

    network_graph = from_networkx(G, nx.spring_layout, scale=10, center=(0, 0))


    #Set node size and color
    network_graph.node_renderer.glyph = Circle(size=node_size,  fill_color="color")
    
    #Set edge width and color green - positive, red - negative
    network_graph.edge_renderer.data_source.data["line_width"] = [G.get_edge_data(a,b)['covariance']*amplify_x for a, b in G.edges()]  ### amplify edges strengh
    network_graph.edge_renderer.data_source.data["line_color"] = ["#117733" if G.get_edge_data(a, b)['covariance'] >= 0 else "#CC6677" for a, b in G.edges()]
    network_graph.edge_renderer.glyph.line_width = {'field': 'line_width'} 
    network_graph.edge_renderer.glyph.line_color = {'field': 'line_color'}

    #Add network graph to the plot
    plot.renderers.append(network_graph)
    
    x, y = zip(*network_graph.layout_provider.graph_layout.values())
    node_labels = list(G.nodes)
    source = ColumnDataSource({'x': x, 'y': y, 'asv': [node_labels[i] for i in range(len(x))]})
    labels = LabelSet(x='x', y='y', text='asv', x_offset=30, y_offset=-15, source=source, render_mode='canvas', text_font_size='12pt')

    plot.renderers.append(labels)    

    return plot


def create_graph(corr_matrix: pd.DataFrame(), threshold: float):
    #take the upper part only
    upper = np.triu(np.ones(corr_matrix.shape)).astype(bool)
    df = corr_matrix.where(upper)
    df = pd.DataFrame(corr_matrix.stack(), columns=['covariance']).reset_index()
    df.columns = ["source", "target", "covariance"]
    
    #remove diagonal entries
    #df = df[df['covariance'] <= threshold]
    df = df[abs(df['covariance']) >= threshold]
    #remove diagonal entries
    df = df[df['source'] != df['target']]
    #remove zero entries
    df = df[df['covariance'] != 0]
    
    #build graph
    G = nx.from_pandas_edgelist(df, edge_attr="covariance")
    
    return G


def add_labels(df):
    i = 1
    for col in df.columns:
        # length of ASVs identifier
        if len(col) == 32:
            asv_name = "ASV_{0}".format(i)
            id_dict[asv_name] = col
            df.rename(columns={col: asv_name}, inplace=True)

            i += 1
    return df


def _get_bounds(nlabels: int):
    bottom = list(chain.from_iterable([[ii] * nlabels for ii in range(nlabels)]))
    top = list(chain.from_iterable([[ii + 1] * nlabels for ii in range(nlabels)]))
    left = list(chain.from_iterable([list(range(nlabels)) for ii in range(nlabels)]))
    right = list(chain.from_iterable([list(range(1, nlabels + 1)) for ii in range(nlabels)]))

    return bottom, top, left, right


def _get_colors(df: pd.DataFrame()):
    rdbu = plt.get_cmap('RdBu')
    cmap = ListedColormap(rdbu(np.arange(256)))
    
    # Create a list of hex color codes from the colormap
    colors = [cmap(i)[:3] for i in range(256)]
    colors = ['#' + ''.join([format(int(c * 255), '02x') for c in color]) for color in colors]
    colors = colors[::-1]  # red - positive, blue - negative

    ccorr = np.arange(-1, 1, 1 / (len(colors) / 2))
    color_list = []
    for value in df.covariance.values:
        ind = bisect.bisect_left(ccorr, value)  # smart array insertion
        if ind == 0:  # avoid ind == -1 on the next step
            ind = ind + 1
        color_list.append(colors[ind - 1])
    return color_list, colors


def create_label_dict(df):
    n_labels = len(df.columns)
    labels_dict = dict(zip(range(n_labels), df.columns))
    labels_dict_reversed = dict(zip(range(n_labels),list(labels_dict.values())[::-1]))
    
    return labels_dict, labels_dict_reversed


def project_covariates(transformed_counts=pd.DataFrame(), raw_counts = pd.DataFrame(), metadata=pd.DataFrame(), L=np.ndarray, y=str, PC=0):
    """
    Perform covariate projection and create a scatter plot using PCA results.

    Parameters:
        transformed_counts (pandas.DataFrame, optional): Transformed count data. Default is an empty DataFrame.
        raw_counts (pandas.DataFrame, optional): Raw count data. Default is an empty DataFrame.
        metadata (pandas.DataFrame, optional): Metadata associated with the samples. Default is an empty DataFrame.
        L (numpy.ndarray): Eigenvalues matrix.
        y (str): Name of the variable to plot on the y-axis.
        PC (int): Index of the principal component to plot on the x-axis. Default is 0.

    Returns:
        bokeh.layouts.row: A row layout containing the scatter plot and color bar.

    """
    r = np.linalg.matrix_rank(L)
    proj, loadings, eigv = PCA(transformed_counts, L, inverse=True)

    eigv_sum = np.sum(eigv)
    var_exp = [(value / eigv_sum) for value in sorted(eigv, reverse=True)]

    counts_sum = raw_counts.sum(axis=0)
    depth = pd.DataFrame(data=counts_sum, columns=["sequencing depth"])
    metadata = depth.join(metadata)

    pc_columns = list('PC{0} ({1}%)'.format(i+1, str(100 * var_exp[i])[:4]) for i in range(0, r))
    df_proj = pd.DataFrame(proj, columns=pc_columns, index=Z_mclr.index)
    df = df_proj.join(metadata)
    
    varName1 = 'PC{0} ({1}%)'.format(PC+1, str(100 * var_exp[PC])[:4])
    varName2 = y
    # varName2 = 'PC{0} ({1}%)'.format(PC+2, str(100 * var_exp[1])[:4])
    df['x'] = df[varName1]
    df['y'] = df[varName2]

    source = ColumnDataSource(df)

    p0 = figure(tools='save, zoom_in, zoom_out, wheel_zoom, box_zoom, reset', plot_width=800, plot_height=800,
                active_scroll="wheel_zoom",
                x_axis_label=varName1, y_axis_label=varName2,
                tooltips=[(varName1, "@" + varName1),
                          (varName2, "@" + varName2)
                          ],
                title=varName1 + " vs " + varName2)
    
    
    
    rdbu = plt.get_cmap('RdPu_r')
    cmap = ListedColormap(rdbu(np.arange(256)))
    # Create a list of hex color codes from the colormap
    colors = [cmap(i)[:3] for i in range(256)]
    colors = ['#' + ''.join([format(int(c * 255), '02x') for c in color]) for color in colors]
    colors = colors[::-1]  # red - positive, blue - negative
    exp_cmap = LinearColorMapper(palette=colors, low=depth.values.min(), high=depth.values.max())
    
    #exp_cmap = LinearColorMapper(palette=Blues8[::-1], low=min(df['sequencing depth'].values), high=max(df['sequencing depth'].values))
    p0.circle('x', 'y', source=source, size=15, line_color=None, fill_color={"field": "sequencing depth", "transform": exp_cmap}, fill_alpha=0.3)

    color_bar_plot = figure(title='sequencing depth', title_location="right",
                            height=500, width=150, toolbar_location=None, min_border=0,
                            outline_line_color=None)

    bar = ColorBar(color_mapper=exp_cmap, location=(1, 1))
    #bar = ColorBar(color_mapper=exp_cmap, location=(1, 1))

    color_bar_plot.add_layout(bar, 'right')
    color_bar_plot.title.align = "center"
    color_bar_plot.title.text_font_size = '12pt'

    layout = row(p0, color_bar_plot)

    return layout


def scater_plot(x, y, width=800, height=600, size=3, color='#610053'):
    """
    Create a scatter plot using the given x and y data.

    Parameters:
        x (pandas.Series or numpy.ndarray): Data for the x-axis.
        y (pandas.Series or numpy.ndarray): Data for the y-axis.
        width (int, optional): Width of the plot. Default is 800.
        height (int, optional): Height of the plot. Default is 600.
        size (int, optional): Size of the data points. Default is 3.
        color (str, optional): Color of the data points. Default is '#610053'.

    Returns:
        bokeh.plotting.figure.Figure: Scatter plot figure.

    """
    bokeh_tools = ["save, zoom_in, zoom_out, wheel_zoom, box_zoom, crosshair, reset, hover"]
    p = figure(plot_width=width, plot_height=height, tools=bokeh_tools, toolbar_location='left')

    source = ColumnDataSource({'x': x, 'y': y})

    p.circle("x", "y", size=3*size, source=source, line_color=None, fill_color=color)

    p.xaxis.axis_label = x.name
    p.yaxis.axis_label = y.name
    
    p.xgrid.grid_line_color = None  # Remove x-axis grid lines
    p.ygrid.grid_line_color = None  # Remove y-axis grid lines
    
    return p



def cluster_corr(corr_array, inplace=False):
    """
    Rearranges the correlation matrix, corr_array, so that groups of highly 
    correlated variables are next to eachother 
    
    Parameters
    ----------
    corr_array : pandas.DataFrame or numpy.ndarray
        a NxN correlation matrix 
        
    Returns
    -------
    pandas.DataFrame or numpy.ndarray
        a NxN correlation matrix with the columns and rows rearranged
    """
    pairwise_distances = sch.distance.pdist(corr_array)
    linkage = sch.linkage(pairwise_distances, method='complete')
    cluster_distance_threshold = pairwise_distances.max()/2
    idx_to_cluster_array = sch.fcluster(linkage, cluster_distance_threshold, 
                                        criterion='distance')
    idx = np.argsort(idx_to_cluster_array)
    
    if not inplace:
        corr_array = corr_array.copy()
    
    if isinstance(corr_array, pd.DataFrame):
        return corr_array.iloc[idx, :].T.iloc[idx, :]
    return corr_array[idx, :][:, idx]


def plot_heatmap(cov, precision, mask, low_rank=None, low=False):
    
    meta_ticks = np.array(cov.columns[-14:]) # TO DO: change according to the data
    bug_ticks = np.arange(len(cov.columns[:-14]))
    ticks = np.hstack((bug_ticks, meta_ticks))

    left = 0.125  # the left side of the subplots of the figure
    right = 0.9  # the right side of the subplots of the figure
    bottom = 0.1  # the bottom of the subplots of the figure
    top = 0.9  # the top of the subplots of the figure
    wspace = -0.6  # the amount of width reserved for blank space between subplots,
    hspace = 0.5  # the amount of height reserved for white space between subplots,
    fontsize = 56
    cmap = "coolwarm"
    vmin = -0.5
    vmax = 0.5
    linewidth = .5
    square = True
    cbar = False

    if low:
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(90, 35))

            plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

            ax1.get_shared_y_axes().join(ax2)
            ax3.get_shared_y_axes().join(ax4)

            g1 = sns.heatmap(cov, cmap=cmap, cbar=cbar, ax=ax1, vmin=vmin, vmax=vmax, linewidth=linewidth, square=square,
                             xticklabels=ticks, yticklabels=ticks)
            g1.set_ylabel('')
            g1.set_xlabel('Covariance', fontsize=fontsize)

            g2 = sns.heatmap(precision, cmap=cmap, cbar=cbar, ax=ax2, vmin=vmin, vmax=vmax, linewidth=linewidth, square=square,
                             xticklabels=ticks, yticklabels=ticks)
            g2.set_ylabel('')
            g2.set_xlabel('Inverse covariance', fontsize=fontsize)
            g2.set_yticks([])

            g3 = sns.heatmap(low_rank, cmap=cmap, ax=ax3, cbar=cbar, vmin=vmin, vmax=vmax, linewidth=linewidth, square=square,
                             xticklabels=ticks, yticklabels=ticks)
            g3.set_ylabel('')
            g3.set_xlabel('Low-rank solution', fontsize=fontsize)
            g3.set_yticks([])

            g4 = sns.heatmap(mask, cmap=cmap, ax=ax4, cbar=cbar, vmin=vmin, vmax=vmax, linewidth=linewidth, square=square,
                             xticklabels=ticks, yticklabels=ticks)
            g4.set_ylabel('')
            g4.set_xlabel('Mask', fontsize=fontsize)
            g4.set_yticks([])
    else:

        wspace = 0.5  # the amount of width reserved for blank space between subplots,
        hspace = 0.5

        fig, (ax1, ax2 ,ax3) = plt.subplots(1, 3, figsize=(90, 35))

        plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

        ax1.get_shared_y_axes().join(ax2, ax3)

        g1 = sns.heatmap(cov, cmap=cmap, cbar=cbar, ax=ax1, vmin=vmin, vmax=vmax, linewidth=linewidth, square=square,
                         xticklabels=ticks, yticklabels=ticks)
        g1.set_ylabel('')
        g1.set_xlabel('Covariance', fontsize=fontsize)

        g2 = sns.heatmap(precision, cmap=cmap, cbar=cbar, ax=ax2, vmin=vmin, vmax=vmax, linewidth=linewidth,
                         square=square,
                         xticklabels=ticks, yticklabels=ticks)
        g2.set_ylabel('')
        g2.set_xlabel('Inverse covariance', fontsize=fontsize)
        g2.set_yticks([])

        g3 = sns.heatmap(mask, cmap=cmap, ax=ax3, cbar=cbar, vmin=vmin, vmax=vmax, linewidth=linewidth, square=square,
                         xticklabels=ticks, yticklabels=ticks)
        g3.set_ylabel('')
        g3.set_xlabel('Mask', fontsize=fontsize)
        g3.set_yticks([])

    return fig