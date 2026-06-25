# uplot.py


import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
from matplotlib.ticker import ScalarFormatter

from udf import long2wide
from umath import format_equation, format_pvalue
from umath import stat_LR as stat_LR1



def ax_format_wrapper(func):
    """Wrapper for formatting ax.
    func: function name to generate plot.
    label_font_size (float): font size for x & y tick and axis labels.
    remove_spline (bool): whether to remove the top and right splines.
    sci_notation (dict): whether to use scientific notations.
        xaxis (bool)
        yaxis (bool)
    """
    def wrapper(*args, **kwargs):
        ax = func(*args, **kwargs)
        kws = kwargs.get('ax_kws')
        if kws is None:
            kws = dict()       # to allow some default settings below.
            
        label_font_size = kws.get('label_font_size')
        if label_font_size is not None:
            ax = ax_set_label_fontsize(ax, label_font_size)
            
        xlabel_pos = kws.get('xlabel_pos')
        if xlabel_pos is not None:
            ax.xaxis.set_label_coords(xlabel_pos[0], xlabel_pos[1])

        remove_spline = kws.get('remove_spline')
        if remove_spline is None:
            remove_spline = True
        if remove_spline:
            ax = ax_remove_splines(ax)
        
        # Use below link to move the sci-notation.
        # https://stackoverflow.com/questions/38331332/move-scientific-notation-exponential-to-right-side-of-y-axis-in-matplotlib
        sci_notation = kws.get('sci_notation')
        if sci_notation is not None:
            ax = ax_set_sci_notation(
                ax, 
                xaxis = sci_notation.get('xaxis', False),
                yaxis = sci_notation.get('yaxis', False),
                xkws = sci_notation.get('xkws'),
                ykws = sci_notation.get('ykws')
            )
        return(ax)
    return wrapper    
    
    

def ax_remove_splines(ax):
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    return(ax)



def ax_set_label_fontsize(ax, font_size):
    for axis_name in ('x', 'y'):
        axis = ax.xaxis if axis_name == 'x' else ax.yaxis
        axis.label.set_size(font_size)       # axis label.
        ax.tick_params(                      # tick label.
            axis = axis_name, which = 'major', labelsize = font_size)
    return(ax)



def ax_set_sci_notation(
    ax, 
    xaxis = False, yaxis = False,
    xkws = None, ykws = None
):
    for axis_name, use_sci, kws in zip(('x', 'y'), (xaxis, yaxis), (xkws, ykws)):
        if use_sci:
            axis = ax.xaxis if axis_name == 'x' else ax.yaxis
            useMathText = True
            fontsize = 9
            pos = None
            if kws:
                if 'useMathText' in kws:
                    useMathText = kws['useMathText']
                if 'fontsize' in kws:
                    fontsize = kws['fontsize']
                if 'position' in kws:
                    pos = kws['position']
            ax.ticklabel_format(axis = axis_name, style = 'sci', 
                                useOffset = False, useMathText = useMathText)
            axis.major.formatter.set_powerlimits((0, 0))
            axis.get_offset_text().set_fontsize(fontsize)
            if pos is not None:
                axis.get_offset_text().set_position(pos)
    return(ax)



def savefig_wrapper(func):
    return savefig_wrapper_v2(func)


def savefig_wrapper_v2(func):
    """Wrapper for saving specific plot.
    func: function name to generate plot.
    fname (None or str): path to the saved figure.
    figsize (None or tuple of (float, float)): width and height in inches.
    dpi (int).
    """
    def wrapper(*args, **kwargs):
        kws = kwargs.get('savefig_kws')
        if kws is None:
            g = func(*args, **kwargs)
            plt.tight_layout()
            return(g)

        figsize = kws.get('figsize')
        fname = kws.get('fname')
        dpi = kws.get('dpi', 300)
        
        fig = plt.gcf()
        if figsize is not None:
            fig.set_size_inches(*figsize, forward = True)
            
        g = func(*args, **kwargs)
        plt.tight_layout()
        
        fig = plt.gcf()
        if figsize is not None:
            fig.set_size_inches(*figsize, forward = True)

        fig.savefig(fname, dpi = dpi)
        return(g)
    return wrapper



### Figures


@savefig_wrapper
@ax_format_wrapper
def box_plot(
    x, y, ylabel, data = None, log_scale = False, 
    colors = None, showfliers = True, rotation = 45, fill = False,
    ylims = None,
    ax_kws = None,
    savefig_kws = None
):
    """Box plot.
    x, y (array or str): array-like values or column names in `data`.
    ylabel (str): label of y axis.
    data (None or DataFrame).
    log_scale (bool): whether to transform `y` with log10(i+1).
    colors, showfliers, fill: pass to sns.boxplot().
    rotation (int): rotation of x axis. 
    """
    ax = None
    if data is not None:
        x = data[x]
        y = data[y]
    if log_scale:
        y = np.log10(y + 1)
    ax = sns.boxplot(
        x = x, y = y, 
        fill = fill, palette = colors, hue = x, legend = False,
        showfliers = showfliers
    )
    if ylims is not None:
        ax.set_ylim(*ylims)
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels(
        ax.get_xticklabels(), 
        rotation = rotation, rotation_mode = "anchor",
        ha = 'right'
    )
    ax.set_xlabel(None)
    ax.set_ylabel(ylabel)
    return(ax)



def stat_LR(
    data,
    metric = "mean", group = "X_name", gvars = None, 
    log_scale = False
):
    """Statistics for linear regression.
    gvars (None or DataFrame).
    """
    df = long2wide(data, columns = group, values = metric)
    if log_scale:
        df = df.map(lambda x: np.log10(x + 1))
        
    if gvars is None:
        groups = data[group].unique()
        cx, cy = [], []
        n = len(groups)
        for i in range(n):
            for j in range(i):
                cx.append(groups[j])
                cy.append(groups[i])
        gvars = pd.DataFrame(data = dict(x = cx, y = cy))
    else:
        assert 'x' in gvars.columns
        assert 'y' in gvars.columns
        
    res = gvars[['x', 'y']].copy()
    res['a'] = np.nan
    res['b'] = np.nan
    res['r2'] = np.nan
    res['p_value'] = np.nan

    for i in range(res.shape[0]):
        cx = res['x'].loc[i]
        cy = res['y'].loc[i]
        a, b, r2, p_value = stat_LR1(
            x = df[cx].to_numpy(), 
            y = df[cy].to_numpy()
        )
        res.loc[i, 'a'] = a
        res.loc[i, 'b'] = b
        res.loc[i, 'r2'] = r2
        res.loc[i, 'p_value'] = p_value
    return(res)



@savefig_wrapper
@ax_format_wrapper
def scatter_plot_LR(
    x, y, data, 
    lr_kws,
    pos_r, pos_fit,
    point_size, 
    title = None,
    xlabel = None, ylabel = None,
    log_scale = False,
    fit_line = True, fit_kws = None,
    ref_line = False, ref_kws = None,
    xlims_line = None,
    ylims = None,
    ax = None,
    ax_kws = None,
    savefig_kws = None
):
    """Scatter plot for linear regression.
    x, y (str): column names.
    data (DataFrame): wide-format.
    lr_kws (dict): linear regression statistics, such as `a`, `b`,
        `r2`, and `p_value`.
    pos_r (tuple of float): position of annotation text for regression R2 
        and p value.
    pos_fit (tuple of float): position of annotation text for regression 
        equation.
    point_size (float): scatter point size.
    fit_kws (None or dict): keywords for fitted line.
        Note the 'color' keyward also applies to text of `pos_fit`.
    ref_kws (None or dict): keywords for reference line `y = x`.
    xlims_line (None or tuple of float): range of xlims for fitted or 
        reference line.
    ax (None or ax): the ax object. None if use `plt.gca()`.
    """
    data = data.copy()
    cx, cy = x, y
    x, y = data[cx], data[cy]
    
    if log_scale:
        x = np.log10(x + 1)
        y = np.log10(y + 1)

    if ax is None:
        ax = plt.gca()
    ax.plot(x, y, marker = '.', markersize = point_size, linestyle = 'none')
    
    for key in ('a', 'b', 'r2', 'p_value'):
        assert key in lr_kws
    a = lr_kws['a']
    b = lr_kws['b']
    r2 = lr_kws['r2']
    p_value = lr_kws['p_value']    
    
    x_vals = None
    if xlims_line is None:
        x_vals = np.linspace(max(ax.get_xlim()[0], np.min(x)), 
                             min(ax.get_xlim()[1], np.max(x)), 100)
    else:
        x_vals = np.linspace(xlims_line[0], xlims_line[1], 100)

    if fit_line:
        if fit_kws is None:
            fit_kws = dict()
        color_fit = fit_kws.get('color', 'gray')
        ax.plot(
            x_vals, a*x_vals+b, 
            linestyle = fit_kws.get('linestyle', '--'), 
            color = color_fit, 
            linewidth = fit_kws.get('linewidth', 1.0)
        )
        equation_str = format_equation(a, b)
        ax.annotate(equation_str, xy = pos_fit, xycoords = ax.transAxes,
                    color = color_fit)

    if ref_line:
        if ref_kws is None:
            ref_kws = dict()
        ax.plot(
            x_vals, x_vals, 
            linestyle = ref_kws.get('linestyle', '--'), 
            color = ref_kws.get('color', 'black'), 
            linewidth = ref_kws.get('linewidth', 1.0)
        )
    
    if ylims is not None:
        ax.set_ylim(*ylims)

    pv_str = format_pvalue(p_value) 
    ax.annotate(r"R$^2$=%.2f, %s" % (r2, pv_str), xy = pos_r, 
                xycoords = ax.transAxes)

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    return(ax)



@savefig_wrapper
def mulax_LR(
    data, st_lr, 
    pos_r, pos_fit, 
    columns = None,
    xlabel = None, ylabel = None,
    log_scale = False,
    point_size = 25,
    fit_line = True, fit_kws = None,
    ref_line = False, ref_kws = None,
    ylims = None,
    xlabel_pos = None, ylabel_pos = None,
    title_font_size = 10,
    nrows = 1, ncols = 1,
    ax_kws = None,
    savefig_kws = None
):
    """Multiple axes for linear regressions.
    """
    for key in ('x', 'y', 'a', 'b', 'r2', 'p_value'):
        assert key in st_lr.columns
    
    if columns is None:
        columns = np.asarray(data.columns)
    df = data[columns].copy()
    if log_scale:
        df = df.map(lambda x: np.log10(x + 1))
        
    gvars = st_lr[['x', 'y']].copy()
    gvars["xy_group"] = gvars['y'] + ' ~ ' + gvars['x']
    
    g = sns.FacetGrid(gvars, col = 'xy_group', col_wrap = ncols)
    
    x_min, x_max = None, None
    cols = np.union1d(gvars['x'], gvars['y'])
    for c in cols:
        if x_min is None:
            x_min = df[c].min()
        elif df[c].min() < x_min:
            x_min = df[c].min()
        if x_max is None:
            x_max = df[c].max()
        elif df[c].max() > x_max:
            x_max = df[c].max()
    xlims = None
    if x_min is not None and x_max is not None:
        xlims = (x_min, x_max)
    
    k = 0
    for i in range(gvars.shape[0]):
        ax = g.axes[k]
        cx = gvars['x'].loc[i]
        cy = gvars['y'].loc[i]
        d = st_lr.loc[(st_lr['x'] == cx) & (st_lr['y'] == cy)]
        assert d.shape[0] == 1
        lr_kws = dict(
            a = d['a'].iloc[0],
            b = d['b'].iloc[0],
            r2 = d['r2'].iloc[0],
            p_value = d['p_value'].iloc[0]
        )
        scatter_plot_LR(
            x = cx, y = cy, data = df,
            lr_kws = lr_kws,
            pos_r = pos_r, pos_fit = pos_fit,  
            point_size = point_size,
            fit_line = fit_line, fit_kws = fit_kws,
            ref_line = ref_line, ref_kws = ref_kws,
            xlims_line = xlims,
            ylims = ylims,
            ax = ax, ax_kws = ax_kws
        )
        ax.set_title(cy + ' ~ ' + cx, fontsize = title_font_size)
        k += 1


    g.set_xlabels(label = xlabel, clear_inner = True)
    if xlabel_pos is not None:
        for ax in g.axes:
            if ax.get_xlabel():
                ax.xaxis.set_label_coords(xlabel_pos[0], xlabel_pos[1])
        
    g.set_ylabels(label = ylabel, clear_inner = True)
    if ylabel_pos is not None:
        for ax in g.axes:
            if ax.get_ylabel():
                ax.yaxis.set_label_coords(ylabel_pos[0], ylabel_pos[1])

    for ax in g.axes:
        for s in ["top", "right"]:
            ax.spines[s].set_visible(True)
    return(g)



@savefig_wrapper
@ax_format_wrapper
def stack_plot_perc(
    x, y, data, hue, 
    ylabel, legend_title, 
    font_size, bbox_to_anchor, colors,
    rotation = 45,
    min_v = 0,
    ax_kws = None,
    savefig_kws = None
):
    """Stacked bar plot for percentage of each group.
    x, y (str): column names in `data`. Columns containing x and y axis values.
    data (DataFrame).
    hue (str): column name in `data`. Column for the names of stacked groups.
    ylabel (str): label of y axis.
    legend_title (str): legend title. Legends are groups in the `hue` column.
    rotation (int): rotation of x axis.
    """
    # Note a quick implementation is using seaborn histplot
    # e.g., https://stackoverflow.com/questions/69846902/how-to-plot-stacked-100-bar-plot-with-seaborn-for-categorical-data
    df = long2wide(data, index = x, columns = hue, values = y)
    df[x] = df.index
    hue_groups = data[hue].unique()
    cumsum = np.array([0.0] * df.shape[0])
    assert len(colors) == len(hue_groups)
    for i, h in enumerate(hue_groups):
        if i == 0:
            plt.bar(df[x], df[h], label = h, color = colors[i])
        else:
            ax = plt.gca()
            ax.bar(df[x], df[h], label = h, bottom = cumsum, color = colors[i])
        cumsum += df[h]
        
    for c in ax.containers:
        labels = [str(round(v.get_height(), 1)) + "%"    \
                  if v.get_height() > min_v else '' for v in c]
        ax.bar_label(c, label_type = 'center', labels = labels, 
                     fontsize = font_size)

    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels(
        ax.get_xticklabels(), 
        rotation = rotation, rotation_mode = "anchor",
        ha = 'right'
    )
    ax.set_ylabel(ylabel)
    ax.set_xlabel("")

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(
        handles = handles[::-1],
        labels = labels[::-1],
        title = legend_title,
        loc = 'best',
        bbox_to_anchor = bbox_to_anchor,
        title_fontsize = font_size + 1,
        fontsize = font_size + 1,
        alignment = 'left'
    )
    return(ax)



@savefig_wrapper
@ax_format_wrapper
def violin_plot(
    x, y, ylabel, data = None, log_scale = False, 
    colors = None, rotation = 45, inner = None, cut = 2, fill = False,
    ylims = None,
    ax_kws = None,
    savefig_kws = None
):
    """Violin plot.
    x, y (array or str): array-like values or column names in `data`.
    ylabel (str): label of y axis.
    data (None or DataFrame).
    log_scale (bool): whether to transform `y` with log10(i+1).
    colors, inner, cut, fill: pass to sns.violinplot().
    rotation (int): rotation of x axis. 
    """
    ax = None
    if data is not None:
        x = data[x]
        y = data[y]
    if log_scale:
        y = np.log10(y + 1)
    ax = sns.violinplot(
        x = x, y = y, 
        fill = fill, palette = colors, hue = x, legend = False,
        inner = inner, cut = cut
    )
    if ylims is not None:
        ax.set_ylim(*ylims)
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels(
        ax.get_xticklabels(), 
        rotation = rotation, rotation_mode = "anchor",
        ha = 'right'
    )
    ax.set_xlabel(None)
    ax.set_ylabel(ylabel)
    return(ax)
