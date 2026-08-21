"""
Module for generating azimuth vs elevation summary heatmaps
for multipath and C/N₀ (SNR) across all constellations and signals.

Produces rectangular 2-D heatmaps with:
  - X-axis: Azimuth  (0–360°)
  - Y-axis: Elevation (0–90°)
  - Color:  Mean multipath [m]  or  mean C/N₀ [dB-Hz]

These give a single panoramic overview of antenna performance,
making it easy to spot reflection zones and signal obstructions.

Requested in GitHub issue #54 / #55.

Made by: Per Helge Aarnes
E-mail: per.helge.aarnes@gmail.com
"""

import os
import logging
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

logger = logging.getLogger(__name__)

_NAME2CODE = {'GPS': 'G', 'GLONASS': 'R', 'Galileo': 'E', 'BeiDou': 'C'}


def _collect_multipath(analysisResults, system_filter=None):
    """Collect azimuth, elevation and |multipath| values.

    Parameters
    ----------
    analysisResults : dict
    system_filter : str or None
        If given (e.g. ``'GPS'``), collect data only for that system.
        If ``None``, collect across all systems.
    """
    all_az, all_el, all_mp = [], [], []
    systems = [system_filter] if system_filter else analysisResults.get('GNSSsystems', [])
    for system in systems:
        sys_code = _NAME2CODE.get(system)
        if sys_code is None:
            continue
        try:
            sat_az = analysisResults['Sat_position'][sys_code]['azimuth']
            sat_el = analysisResults['Sat_position'][sys_code]['elevation']
        except KeyError:
            logger.warning("SkyPlotSummary: No satellite position data for %s, skipping.", system)
            continue

        for band in analysisResults[system].get('Bands', []):
            codes = analysisResults[system][band].get('Codes', [])
            codes = [c for c in codes if c and c != []]
            for code in codes:
                if code not in analysisResults[system][band]:
                    continue
                mp = analysisResults[system][band][code].get('multipath_range1')
                if mp is None:
                    continue
                all_az.append(sat_az.ravel())
                all_el.append(sat_el.ravel())
                all_mp.append(np.abs(mp).ravel())

    if not all_az:
        return None, None, None
    return np.concatenate(all_az), np.concatenate(all_el), np.concatenate(all_mp)


def _collect_snr(analysisResults, system_filter=None):
    """Collect azimuth, elevation and SNR values.

    Parameters
    ----------
    analysisResults : dict
    system_filter : str or None
        If given (e.g. ``'GPS'``), collect data only for that system.
        If ``None``, collect across all systems.
    """
    all_az, all_el, all_snr = [], [], []
    systems = [system_filter] if system_filter else analysisResults.get('GNSSsystems', [])
    for system in systems:
        sys_code = _NAME2CODE.get(system)
        if sys_code is None:
            continue
        try:
            sat_az = analysisResults['Sat_position'][sys_code]['azimuth']
            sat_el = analysisResults['Sat_position'][sys_code]['elevation']
        except KeyError:
            continue

        snr_dict = analysisResults[system].get('SNR', {})
        for code, snr_data in snr_dict.items():
            snr = np.array(snr_data, dtype=float)
            snr[snr == 0] = np.nan
            all_az.append(sat_az.ravel())
            all_el.append(sat_el.ravel())
            all_snr.append(snr.ravel())

    if not all_az:
        return None, None, None
    return np.concatenate(all_az), np.concatenate(all_el), np.concatenate(all_snr)


def _bin_data(az, el, values, az_bin_size=5, el_bin_size=2):
    """Bin scattered (az, el, value) data into a 2-D grid of mean values.

    Returns
    -------
    az_edges : 1-D array, azimuth bin edges  (degrees)
    el_edges : 1-D array, elevation bin edges (degrees)
    grid     : 2-D array [n_el_bins, n_az_bins], mean value per cell (NaN where no data)
    counts   : 2-D array, number of observations per cell
    """
    valid = np.isfinite(az) & np.isfinite(el) & np.isfinite(values)
    az, el, values = az[valid], el[valid], values[valid]

    az_edges = np.arange(0, 360 + az_bin_size, az_bin_size)
    el_edges = np.arange(0, 90 + el_bin_size, el_bin_size)
    n_az = len(az_edges) - 1
    n_el = len(el_edges) - 1

    sum_grid = np.zeros((n_el, n_az))
    cnt_grid = np.zeros((n_el, n_az), dtype=int)

    az_idx = np.clip(((az) / az_bin_size).astype(int), 0, n_az - 1)
    el_idx = np.clip(((el) / el_bin_size).astype(int), 0, n_el - 1)

    np.add.at(sum_grid, (el_idx, az_idx), values)
    np.add.at(cnt_grid, (el_idx, az_idx), 1)

    with np.errstate(invalid='ignore'):
        grid = np.where(cnt_grid > 0, sum_grid / cnt_grid, np.nan)

    return az_edges, el_edges, grid, cnt_grid


def _plot_heatmap(az_edges, el_edges, grid, title, cbar_label, cmap, graph_dir,
                  filename, vmin=None, vmax=None, use_latex=False):
    """Render and save a single rectangular heatmap."""
    prev_usetex = plt.rcParams.get('text.usetex', False)
    try:
        plt.rcParams['text.usetex'] = use_latex

        fig, ax = plt.subplots(figsize=(16, 6), dpi=200)

        if vmin is None:
            vmin = np.nanmin(grid) if np.any(np.isfinite(grid)) else 0
        if vmax is None:
            vmax = np.nanmax(grid) if np.any(np.isfinite(grid)) else 1

        norm = Normalize(vmin=vmin, vmax=vmax)

        # pcolormesh wants (n_el+1, n_az+1) edges
        im = ax.pcolormesh(az_edges, el_edges, grid, cmap=cmap, norm=norm, shading='flat')

        cbar = fig.colorbar(im, ax=ax, orientation='vertical', shrink=0.85, pad=0.02, aspect=25)
        cbar.set_label(cbar_label, fontsize=14)
        cbar.ax.tick_params(labelsize=12)

        ax.set_xlabel('Azimuth [deg]', fontsize=14)
        ax.set_ylabel('Elevation [deg]', fontsize=14)
        ax.set_title(title, fontsize=16, pad=10)
        ax.set_xlim(0, 360)
        ax.set_ylim(0, 90)
        ax.set_xticks(np.arange(0, 361, 30))
        ax.set_yticks(np.arange(0, 91, 10))
        ax.tick_params(labelsize=12)
        ax.grid(True, alpha=0.3, linewidth=0.5)

        # Mark empty cells with light grey
        masked = np.ma.masked_where(~np.isnan(grid), grid)
        ax.pcolormesh(az_edges, el_edges, masked, cmap=matplotlib.colors.ListedColormap(['#f0f0f0']), shading='flat')

        fig.tight_layout()
        fig.savefig(os.path.join(graph_dir, filename), dpi=300, bbox_inches='tight')
        plt.close(fig)
    finally:
        plt.rcParams['text.usetex'] = prev_usetex


def _generate_multipath_heatmap(analysisResults, graph_dir, az_bin_size, el_bin_size,
                                use_latex, system_filter=None):
    """Generate a multipath heatmap for one system or all combined."""
    az, el, mp = _collect_multipath(analysisResults, system_filter=system_filter)
    if az is None or len(az) == 0:
        return False
    az_edges, el_edges, grid, _ = _bin_data(az, el, mp, az_bin_size, el_bin_size)
    vmax = round(2 * np.nanmean(grid[np.isfinite(grid)]), 2) if np.any(np.isfinite(grid)) else 1.0
    if system_filter:
        title = f'Mean Multipath vs Azimuth and Elevation ({system_filter})'
        filename = f'Summary_Multipath_AzEl_{system_filter}.png'
    else:
        title = 'Mean Multipath vs Azimuth and Elevation (all constellations / signals)'
        filename = 'Summary_Multipath_AzEl.png'
    _plot_heatmap(az_edges, el_edges, grid, title=title,
                  cbar_label='Mean |Multipath| [m]', cmap='cividis',
                  graph_dir=graph_dir, filename=filename,
                  vmin=0, vmax=vmax, use_latex=use_latex)
    return True


def _generate_snr_heatmap(analysisResults, graph_dir, az_bin_size, el_bin_size,
                           use_latex, system_filter=None):
    """Generate a C/N₀ heatmap for one system or all combined."""
    az, el, snr = _collect_snr(analysisResults, system_filter=system_filter)
    if az is None or len(az) == 0:
        return False
    az_edges, el_edges, grid, _ = _bin_data(az, el, snr, az_bin_size, el_bin_size)
    if system_filter:
        title = f'Mean C/N$_0$ vs Azimuth and Elevation ({system_filter})'
        filename = f'Summary_SNR_AzEl_{system_filter}.png'
    else:
        title = r'Mean C/N$_0$ vs Azimuth and Elevation (all constellations / signals)'
        filename = 'Summary_SNR_AzEl.png'
    _plot_heatmap(az_edges, el_edges, grid, title=title,
                  cbar_label=r'Mean C/N$_0$ [dB-Hz]', cmap='jet',
                  graph_dir=graph_dir, filename=filename,
                  use_latex=use_latex)
    return True


def make_skyplot_summary(analysisResults, graph_dir, az_bin_size=5, el_bin_size=2, use_latex=False):
    """Generate the azimuth-vs-elevation summary heatmaps.

    Creates PNG files in *graph_dir*:
      - ``Summary_Multipath_AzEl.png``               (combined, all systems)
      - ``Summary_Multipath_AzEl_<system>.png``       (one per system)
      - ``Summary_SNR_AzEl.png``                      (combined, all systems)
      - ``Summary_SNR_AzEl_<system>.png``             (one per system)

    Parameters
    ----------
    analysisResults : dict
        The result dictionary returned by ``GNSS_MultipathAnalysis()``.
    graph_dir : str
        Directory where the plots are saved.
    az_bin_size : float
        Azimuth bin width in degrees (default 5).
    el_bin_size : float
        Elevation bin width in degrees (default 2).
    use_latex : bool
        Whether to use LaTeX rendering for text (default False).
    """
    os.makedirs(graph_dir, exist_ok=True)
    systems = analysisResults.get('GNSSsystems', [])

    # ── Multipath heatmaps ────────────────────────────────────────────────
    # Combined
    if not _generate_multipath_heatmap(analysisResults, graph_dir, az_bin_size, el_bin_size, use_latex):
        logger.warning("SkyPlotSummary: No multipath data available for summary plot.")
    # Per system
    for system in systems:
        _generate_multipath_heatmap(analysisResults, graph_dir, az_bin_size, el_bin_size, use_latex, system_filter=system)

    # ── SNR / C/N₀ heatmaps ──────────────────────────────────────────────
    # Combined
    if not _generate_snr_heatmap(analysisResults, graph_dir, az_bin_size, el_bin_size, use_latex):
        logger.warning("SkyPlotSummary: No SNR data available for summary plot.")
    # Per system
    for system in systems:
        _generate_snr_heatmap(analysisResults, graph_dir, az_bin_size, el_bin_size, use_latex, system_filter=system)
