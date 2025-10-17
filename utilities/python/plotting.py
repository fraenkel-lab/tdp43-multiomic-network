import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr, pearsonr
from adjustText import adjust_text
from collections import defaultdict


def facet_scatterplot(
    df,
    x_col,
    y_col,
    hue_col=None,
    col_col=None,
    row_col=None,
    # Facet-grid sharing
    sharex=False,
    sharey=False,
    # Optionally fix x and y limits across all facets
    xlim=None,
    ylim=None,
    # Correlation options
    correlation=None,    # 'spearman', 'pearson', or None
    corr_location=(0.05, 0.95),  # (x, y) in Axes coords
    point_size = 60,
    # Highlighting sets of genes -> { 'Label': (set_of_genes, 'color'), ... }
    highlight_dict=None,
    highlight_marker_size=70,
    # Labeling a subset of genes
    label_genes=None,       # A set/list of gene names to label
    label_offset=(0.05, 0.05),
    label_fontsize=10,
    # Aesthetic toggles
    dotted_axes=False,
    dotted_axes_color = 'black',
    dotted_grid_lines = False,
    remove_facet_titles=True,
    annotate_num_points=False,
    # Additional arguments to sns.relplot
    height=4,
    aspect=1.2,
    alpha=0.7,
    palette=None,
    hue_order=None,
    # Additional facet_kws beyond sharex/sharey
    extra_facet_kws=None,
    # If there's only one facet and no hue, optionally set a single color
    single_facet_color=None,
    # Whether to use adjustText to minimize overlapping labels
    use_adjust_text=False,
    # Resolution of the figure in dots-per-inch
    dpi=100,
    **kwargs
):
    """
    Create a faceted scatter plot with optional correlation annotation,
    highlighting specific genes, labeling selected points, and
    optional dotted reference lines or shared axes.
    
    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame containing the plotting data.
    x_col : str
        Column name for x-axis values.
    y_col : str
        Column name for y-axis values.
    hue_col : str, optional
        Column name for hue encoding. (Default: None)
    col_col : str, optional
        Column name to facet across columns. (Default: None)
    row_col : str, optional
        Column name to facet across rows. (Default: None)
    sharex, sharey : bool, optional
        Whether to share x and y axes across facets. (Default: False)
    xlim, ylim : tuple of floats, optional
        If given, applies these axis limits to every facet. (Default: None)
    correlation : {'spearman', 'pearson', None}, optional
        If not None, compute correlation for each facet
        and annotate on the plot. (Default: None)
    corr_location : tuple of float
        Axes coordinates for correlation text (Default: (0.05, 0.95)).
    highlight_dict : dict, optional
        Dictionary specifying gene sets to highlight with custom borders.
        Example: {
            'TFs': (set_of_tfs, 'black'),
            'RBPs': (set_of_rbps, 'red')
        }
    highlight_marker_size : int
        Size of the highlight markers. (Default: 70)
    label_genes : set or list, optional
        Genes to label with text on the scatter plot. (Default: None)
    label_offset : (float, float)
        Offset in data coordinates for labeling. (Default: (0.05, 0.05))
    label_fontsize : int
        Font size for the gene labels. (Default: 10)
    dotted_axes : bool
        If True, draws dotted grey lines at x=0 and y=0. (Default: False)
    remove_facet_titles : bool
        Whether to remove the default facet titles. (Default: True)
    annotate_num_points : bool
        If True, annotate top-right corner with N=<count> in each facet. (Default: False)
    height, aspect : float
        Control the size of each subplot. (Defaults: 4, 1.2)
    alpha : float
        Opacity for the main scatter points. (Default: 0.7)
    palette : str or list, optional
        Color palette for the hue dimension. (Default: None)
    hue_order : list, optional
        Specific order for the hue categories. (Default: None)
    extra_facet_kws : dict, optional
        Additional facet_kws to pass to sns.relplot. (Default: None)
    single_facet_color : str, optional
        If row_col, col_col, and hue_col are None, color all points with this single color.
    use_adjust_text : bool
        If True, use the adjustText package to avoid overlapping labels. (Default: False)
    dpi : int
        Figure resolution in dots-per-inch. (Default: 100)
    **kwargs
        Additional keyword arguments passed directly to sns.relplot().

    Returns
    -------
    g : sns.FacetGrid
        The resulting FacetGrid object for further customization.
    """

    # We will create a temporary rc_context to control figure dpi
    # This ensures that relplot uses your specified dpi.
    with plt.rc_context({'figure.dpi': dpi}):

        # Prepare facet_kws
        facet_kws = {'sharex': sharex, 'sharey': sharey}
        if extra_facet_kws:
            facet_kws.update(extra_facet_kws)

        # If single facet is requested (no row/col/hue),
        # and single_facet_color is set, force the scatter color.
        if (col_col is None and row_col is None and hue_col is None and single_facet_color):
            # We'll pass color to the scatterplot. Also remove palette to avoid conflicts.
            kwargs['color'] = single_facet_color
            palette_for_plot = None
        else:
            palette_for_plot = palette

        # Create a relplot -> FacetGrid
        g = sns.relplot(
            data=df,
            x=x_col,
            y=y_col,
            hue=hue_col,
            col=col_col,
            row=row_col,
            alpha=alpha,
            height=height,
            aspect=aspect,
            palette=palette_for_plot,
            hue_order=hue_order,
            facet_kws=facet_kws,
            s = point_size,
            **kwargs
        )

        # If there's no actual faceting, g.ax is the single Axes
        if row_col is None and col_col is None:
            ax_dict = {'single': g.ax}
        else:
            ax_dict = g.axes_dict

        # We'll collect the Text objects for each facet if using adjustText
        texts_dict = defaultdict(list)

        # Loop over each facet Axes
        for facet_values, ax in ax_dict.items():

            # 1) Determine the facet data
            if row_col is None and col_col is None:
                # No faceting => entire DataFrame is used
                facet_data = df
            elif row_col is not None and col_col is None:
                # Facet by row only
                row_val = facet_values
                facet_data = df[df[row_col] == row_val]
            elif row_col is None and col_col is not None:
                # Facet by col only
                col_val = facet_values
                facet_data = df[df[col_col] == col_val]
            else:
                # Facet by both row and col
                row_val, col_val = facet_values
                facet_data = df[(df[row_col] == row_val) & (df[col_col] == col_val)]

            # 2) Remove facet titles if requested
            if remove_facet_titles:
                ax.set_title('')

            # 3) Optional dotted reference lines at x=0, y=0
            if dotted_axes:
                ax.axvline(0, color=dotted_axes_color, linestyle=':', linewidth=1)
                ax.axhline(0, color=dotted_axes_color, linestyle=':', linewidth=1)
                
            if dotted_grid_lines:
                # Dotted grid lines
                ax.grid(True, linestyle=':', color='gray', linewidth=0.5)

            # 4) Enforce global xlim/ylim if provided
            if xlim is not None:
                ax.set_xlim(xlim)
            if ylim is not None:
                ax.set_ylim(ylim)

            # 5) Correlation annotation
            if correlation is not None and len(facet_data) > 1:
                xvals = facet_data[x_col]
                yvals = facet_data[y_col]
                if correlation.lower() == 'spearman':
                    corr_val, p_val = spearmanr(xvals, yvals)
                elif correlation.lower() == 'pearson':
                    corr_val, p_val = pearsonr(xvals, yvals)
                else:
                    corr_val, p_val = None, None  # or raise an error

                if corr_val is not None:
                    ax.text(
                        corr_location[0],
                        corr_location[1],
                        f"{correlation.title()} r={corr_val:.2f}\nP={p_val:.1e}",
                        transform=ax.transAxes,
                        ha='left',
                        va='top',
                        fontsize=14,
                        color='black'
                    )

            # 6) Annotate number of points in the top-right corner
            if annotate_num_points:
                num_points = len(facet_data)
                ax.text(
                    0.95, 0.95,
                    f"N={num_points}",
                    transform=ax.transAxes,
                    ha='right',
                    va='top',
                    fontsize=14,
                    color='black',
                    weight='bold'
                )

            # 7) Highlight sets of genes with borders
            if highlight_dict:
                for label, (gene_set, border_color) in highlight_dict.items():
                    sub_highlight = facet_data[facet_data['gene_name'].isin(gene_set)]
                    if not sub_highlight.empty:
                        ax.scatter(
                            x=sub_highlight[x_col],
                            y=sub_highlight[y_col],
                            facecolors='none',
                            edgecolors=border_color,
                            linewidths=1.2,
                            s=highlight_marker_size,
                            label=label
                        )

            # 8) Label specific genes
            if label_genes:
                sub_labels = facet_data[facet_data['gene_name'].isin(label_genes)]
                for _, row_ in sub_labels.iterrows():
                    txt = ax.text(
                        row_[x_col] + label_offset[0],
                        row_[y_col] + label_offset[1],
                        row_['gene_name'],
                        fontsize=label_fontsize,
                        color='black',
                        ha='left',
                        va='bottom'
                    )
                    # If using adjustText, store this Text object
                    if use_adjust_text:
                        texts_dict[facet_values].append(txt)

        # Set overall axis labels
        g.set_axis_labels(x_col, y_col)
        plt.tight_layout()

        # Optionally call adjustText on each facet to minimize label overlap
        # must be called after everything else
        if use_adjust_text:
            for facet_values, ax in ax_dict.items():
                
                # Determine facet valeus again
                # 1) Determine the facet data
                if row_col is None and col_col is None:
                    # No faceting => entire DataFrame is used
                    facet_data = df
                elif row_col is not None and col_col is None:
                    # Facet by row only
                    row_val = facet_values
                    facet_data = df[df[row_col] == row_val]
                elif row_col is None and col_col is not None:
                    # Facet by col only
                    col_val = facet_values
                    facet_data = df[df[col_col] == col_val]
                else:
                    # Facet by both row and col
                    row_val, col_val = facet_values
                    facet_data = df[(df[row_col] == row_val) & (df[col_col] == col_val)]
                
                xvals = facet_data[x_col]
                yvals = facet_data[y_col]
                if texts_dict[facet_values]:  # only if we actually have text labels
                    adjust_text(
                        texts_dict[facet_values],
                        x = xvals,
                        y = yvals,
                        ax=ax,
                        # You can add arrowprops if you want lines connecting labels
                        arrowprops=dict(arrowstyle="->", color='black'),
                        lw = 0.7,
                        expand=(1.5, 2),
                        avoid_self=True,
                    )

    return g