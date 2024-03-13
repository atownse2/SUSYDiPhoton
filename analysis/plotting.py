import os

import numpy as np
import hist as bh

import matplotlib as mpl
import matplotlib.pyplot as plt
import mplhep as hep

import analysis.histograms as hgm


hep.style.use(hep.style.CMS)

fig_size = (10, 10)
cms_font_size = 20

top_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
plot_dir = f"{top_dir}/outputs/plots"
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)

SEABORN_PALETTES = dict(
    deep=["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3",
          "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD"],
    deep6=["#4C72B0", "#55A868", "#C44E52",
           "#8172B3", "#CCB974", "#64B5CD"],
    muted=["#4878D0", "#EE854A", "#6ACC64", "#D65F5F", "#956CB4",
           "#8C613C", "#DC7EC0", "#797979", "#D5BB67", "#82C6E2"],
    muted6=["#4878D0", "#6ACC64", "#D65F5F",
            "#956CB4", "#D5BB67", "#82C6E2"],
    pastel=["#A1C9F4", "#FFB482", "#8DE5A1", "#FF9F9B", "#D0BBFF",
            "#DEBB9B", "#FAB0E4", "#CFCFCF", "#FFFEA3", "#B9F2F0"],
    pastel6=["#A1C9F4", "#8DE5A1", "#FF9F9B",
             "#D0BBFF", "#FFFEA3", "#B9F2F0"],
    bright=["#023EFF", "#FF7C00", "#1AC938", "#E8000B", "#8B2BE2",
            "#9F4800", "#F14CC1", "#A3A3A3", "#FFC400", "#00D7FF"],
    bright6=["#023EFF", "#1AC938", "#E8000B",
             "#8B2BE2", "#FFC400", "#00D7FF"],
    dark=["#001C7F", "#B1400D", "#12711C", "#8C0800", "#591E71",
          "#592F0D", "#A23582", "#3C3C3C", "#B8850A", "#006374"],
    dark6=["#001C7F", "#12711C", "#8C0800",
           "#591E71", "#B8850A", "#006374"],
    colorblind=["#0173B2", "#DE8F05", "#029E73", "#D55E00", "#CC78BC",
                "#CA9161", "#FBAFE4", "#949494", "#ECE133", "#56B4E9"],
    colorblind6=["#0173B2", "#029E73", "#D55E00",
                 "#CC78BC", "#ECE133", "#56B4E9"]
)

# plt.rcParams['axes.prop_cycle'] = plt.cycler(color=SEABORN_PALETTES['colorblind6'])

def add_version_tag(fig, tag):
    fig.text(0.07, 0.01, tag, ha='left', va='bottom', fontsize=16)

def get_error_bars(hist):

    def calc_garwood(k):
        import scipy.stats as stats
        alpha = 1 - 0.682689492137086
        lower = stats.chi2.ppf(alpha/2, 2*k)/2
        upper = stats.chi2.ppf(1-alpha/2, 2*(k+1))/2

        return k-lower, upper-k
    pass

def plot_stack(hists, labels, title=None,
               setlogy=False):
    fig, ax = plt.subplots()
    fig.set_size_inches(*fig_size)
    if setlogy:
        ax.set_yscale('log')
    hep.histplot(hists, ax=ax, histtype='fill', stack=True, label=labels)
    ax.set_title(title, fontsize=12)
    ax.legend(loc='upper left', frameon=False)
    hep.cms.label(ax=ax, year='2018')

def plot1d(
    h,
    yerr=None,
    flow='none',
    subtitle=None, subtitle_coord=(0.5, 0.9),
    cms_label='WIP', is_data=False, cms_year='2018', do_cms_label=True,
    ylabel=None, setlogy=False, 
    xlim=None,
    fig=None, ax=None, legend=None, label=None,
    density=False,
    save_as=None,
    ):

    if fig is None or ax is None:
        fig, ax = plt.subplots()
        fig.set_size_inches(*fig_size)
    if setlogy:
        ax.set_yscale('log')
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if xlim is not None:
        ax.set_xlim(xlim)
    else:
        xlim=get_1D_axis_lims(h)
        ax.set_xlim(xlim)
    if subtitle is not None:
        ax.text(*subtitle_coord, subtitle, ha='center', va='bottom', transform=ax.transAxes, fontsize=16)
    hep.histplot(
        h,
        flow=flow,
        ax=ax,
        label=label,
        density=density, 
        color='black' )
    # ax.errorbar(h.axes[0].centers, h.values(), yerr=np.sqrt(h.values()), fmt='o', markersize=4, linestyle='None', label=label)
    if do_cms_label:
        hep.cms.label(ax=ax, year=cms_year, label=cms_label, data=is_data, fontsize=cms_font_size)

    if legend is not None:
        ax.legend(loc=legend, fontsize=16)
    plt.tight_layout()
    if save_as is not None:
        plt.savefig(save_as)
    return fig, ax

def plot2d(h, title=None, fig_size=(10,10),
           fig=None, ax=None,
            cmap='binary', alpha=1, cbar=True,
            xticks=None, yticks=None,
            display_bin_values=False,
            xlim=None, ylim=None,
            setlogy=False, setlogx=False, setlogz=False,
            signal_point=None,
            is_data=False, cms_label='Work in Progress',
            year=None,
            save_as=None):
    
    if isinstance(h, bh.Hist):
        if 'dataset' in h.axes.name:
            h = h[{'dataset':sum}]

    if fig is None or ax is None:
        fig, ax = plt.subplots()
        fig.set_size_inches(*fig_size)
    ax.set_title(title)
    if setlogy:
        ax.set_yscale('log')
    if setlogx:
        ax.set_xscale('log')
    if xticks is not None:
        ax.set_xticks(xticks)
    if yticks is not None:
        ax.set_yticks(yticks)

    if display_bin_values:
        alpha=0.5
        H, xedges, yedges = h.to_numpy()
        for i in range(len(xedges)-1):
            for j in range(len(yedges)-1):
                ax.text((xedges[i]+xedges[i+1])/2, (yedges[j]+yedges[j+1])/2, np.round(H[i,j], decimals=2), 
                        ha='center', va='center', color='black')
    hep.hist2dplot(
        h,
        ax=ax,
        cmap=cmap,
        flow=None,
        norm=mpl.colors.LogNorm() if setlogz else None,
        alpha=alpha,
        cbar=cbar)
    
    hep.cms.label(
        ax=ax,
        year=year,
        label=cms_label,
        data=is_data,
        fontsize=cms_font_size)

    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    if signal_point is not None:
        m_Bkk, m_Radion = signal_point
        ax.text(0.05, 0.95, f'$m_{{Bkk}} = {m_Bkk}$ GeV\n$m_{{Radion}} = {m_Radion}$ GeV',
                ha='left', va='top', transform=ax.transAxes, fontsize=12)
        ax.set_xlim(m_Bkk-50, m_Bkk+50)
        ax.set_ylim(m_Radion-2, m_Radion+2)

    fig.tight_layout()
    if save_as is not None:
        plt.savefig(save_as)

def get_2D_axis_lims(hist):
    x_hist = hist[:,sum]
    y_hist = hist[sum,:]

    x_lim = get_1D_axis_lims(x_hist)
    y_lim = get_1D_axis_lims(y_hist)

    return x_lim, y_lim

def get_1D_axis_lims(hists):

    x_min = np.inf
    x_max = -np.inf

    if isinstance(hists, bh.Hist):
        hists = [hists]

    for hist in hists:
        bin_contents = hist.values()
        # Mask to ignore zero bins
        mask = bin_contents != 0

        # Apply the mask and find the minimum non-zero bin
        bin_indicies = np.arange(len(bin_contents))

        min_non_zero_bin_index = bin_indicies[mask][0]
        max_non_zero_bin_index = bin_indicies[mask][-1] + 3

        if max_non_zero_bin_index > len(bin_contents):
            max_non_zero_bin_index = len(bin_contents)

        x_min = min(x_min, hist.axes[0].edges[min_non_zero_bin_index])
        x_max = max(x_max, hist.axes[0].edges[max_non_zero_bin_index])
        

    xlim = (x_min, x_max)
    return xlim

def plot1ds(hists, labels, sub_title=None,
            cms_label='Work in Progress', is_data=False,
            setlogy=False, ylabel=None,
            legend='best',
            xlim=None, xticks=None, 
            density=False,
            save_as=None):
    fig, ax = plt.subplots()

    xlim=get_1D_axis_lims(hists)

    for hist, label in zip(hists,labels):
        plot1d(hist, label=label,
               sub_title=sub_title, cms_label=cms_label, 
               is_data=is_data, 
               setlogy=setlogy, 
               fig=fig, ax=ax,
               legend=legend,
               xlim=xlim,
               density=density)

def plot_residuals(hist, model, labels=None, title=None, setlogy=False, ylabel=None, xlim=None, xticks=None, save_as=None):
    # Plot histograms with residuals underneath
    fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]}, sharex=True)

    ax1.set_title(title)
    if setlogy:
        ax1.set_yscale('log')
    if ylabel is not None:
        ax1.set_ylabel(ylabel)

    x = np.linspace(0,4000,10000)

    hep.histplot(hist, ax=ax1, histtype='fill', stack=True)
    ax1.plot(x, model(x))

    # Residual
    hist_vals, bin_edges = hist.to_numpy()
    bin_centers = (bin_edges[1:] + bin_edges[:-1])/2

    model_vals = model(bin_centers)
    residuals = model_vals - hist_vals
    standardized_residuals = residuals / np.sqrt(hist_vals)

    ax2.errorbar(bin_centers, standardized_residuals, fmt='o', color='black', markersize=2, linestyle='None')
    # hep.histplot(hists, ax=ax2, histtype='errorbar', label=labels, stack=True, density=True)

    hep.cms.label(ax=ax1, year='2018')
    ax1.legend(loc='upper left', frameon=False)
    ax1.set_xlim(xlim)
    ax2.set_xlabel(None)
    ax2.set_ylabel('Residual')
    ax2.axhline(1, color='black', linestyle='--')
    if save_as is not None:
        plt.savefig(save_as)

def plot1d_ratio(h1,h2, ax, label=None, color='black'):
    
    bin_edges = h1.axes[0].edges
    bins = [(bin_edges[i], bin_edges[i+1]) for i in range(len(bin_edges)-1)]

    y1 = h1.values()
    y2 = h2.values()
    v1 = h1.variances()
    v2 = h2.variances()

    zero_mask = y2 > 0
    y1 = y1[zero_mask]
    y2 = y2[zero_mask]
    v1 = v1[zero_mask]
    v2 = v2[zero_mask]
    bins = [ b for i, b in enumerate(bins) if zero_mask[i]]
    bin_edges = [lower for lower, upper in bins] + [bins[-1][1]]

    ratio = y1/y2
    ratio_variance = v1*(1/y2)**2 + v2*(2*y1/y2**2)**2

    hep.histplot((ratio, bin_edges), yerr=np.sqrt(ratio_variance), ax=ax, label=label)

def plot_analysis_bins(
        h,
        label: str,
        yerr=None,
        fig=None, axs=None,
        ylabel="Events",
        is_data=False,
        print_bdt_labels=True,
        year=None,
        legend=True,
        cms_label='Work in Progress',
        save_as=None):

    bdt_labels = ['Low BDT', 'Med BDT', 'High BDT']

    if fig is None and axs is None:
        fig, axs = plt.subplots(1,3, sharey=True, figsize=(15,5))
    elif fig is None or axs is None:
        raise ValueError('Must provide both fig and axs or neither.')

    for i, ax in enumerate(axs):
        if isinstance(h, bh.Hist):
            h1d = h[{'bdt':i}]
            values, edges = h1d.to_numpy(flow=True)
            yerr = np.sqrt(h1d.variances(flow=True))
        elif isinstance(h, hgm.SimpleHistogram):
            values = h.values[i]
            edges = h.bins[1]
            yerr = np.sqrt(h.variances[i])
        else:
            raise ValueError('h must be a hist or SimpleHistogram')
            
        edges[-1] = edges[-2] + (edges[-2] - edges[-3])*1.5
        
        hep.histplot(
            (values,edges),
            yerr=yerr,
            ax=ax,
            label=label)
        
        if print_bdt_labels:
            ax.text(
                0.05, 0.92,
                bdt_labels[i],
                transform=ax.transAxes,
                ha='left', va='top', 
                fontsize=18)

        bin_centers = (edges[1:] + edges[:-1])/2
        ax.set_xticks(bin_centers)
        edges = [round(edge, 0) for edge in edges]
        xticklabels = [f"[{edges[i]}, {edges[i+1]})" for i in range(len(edges)-2)] + [f"[{edges[-2]}, Inf)"]
        ax.set_xticklabels(
            xticklabels,
            rotation=-35,
            ha='left',
            fontsize=18)

        ax.set_xlim(edges[0], edges[-1])
    
    if legend:
        axs[0].legend(
            fontsize=18,
            loc='upper right',
            bbox_to_anchor=(0.9, 0.7))

    axs[0].set_ylabel(ylabel)

    # Create a new full-figure axes, but make it invisible
    fig_label = fig.add_subplot(111, frame_on=False)
    fig_label.set_xlabel(
        'Hard MET $p_T$ [GeV]',
        fontsize=18,
        labelpad=100)
    
    hep.cms.label(
        ax=fig_label,
        year=year,
        label=cms_label,
        data=is_data,
        fontsize=cms_font_size)

    fig_label.set_xticks([])
    fig_label.set_yticks([])

    fig.subplots_adjust(wspace=0)
    
    if save_as is not None:
        plt.savefig(save_as)
    
    return fig, axs

