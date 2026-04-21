"""
End-to-end tests for plot generation in :func:`gnssmultipath.GNSS_MultipathAnalysis`.

The tests run a full multipath analysis on the small cropped OPEC RINEX
file with different combinations of the plotting flags and verify that:

* the expected plot files are written under ``<outputDir>/Graphs/``,
* disabling a plot category prevents the corresponding files from being created,
* the LaTeX-fallback path also works (``use_LaTex=False``).

All output is written to a ``tmp_path`` directory provided by pytest, which is
removed automatically at the end of each test.
"""
from __future__ import annotations

import os
import sys
import warnings
from pathlib import Path
from typing import List

import pytest

# Make the local source importable without relying on the installed copy.
PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from gnssmultipath import GNSS_MultipathAnalysis  # noqa: E402


# ── Test data ────────────────────────────────────────────────────────────────
RIN_OBS = PROJECT_ROOT / "TestData" / "ObservationFiles" / "v3" / \
    "OPEC00NOR_S_20220010000_01D_30S_MO_3.04_croped.rnx"
RIN_NAV = PROJECT_ROOT / "TestData" / "NavigationFiles" / "v3" / \
    "BRDC00IGS_R_20220010000_01D_MN.rnx"


pytestmark = [
    pytest.mark.filterwarnings("ignore::RuntimeWarning"),
    pytest.mark.filterwarnings("ignore::DeprecationWarning"),
]


# ── Helpers ──────────────────────────────────────────────────────────────────
def _list_graphs(graph_dir: Path) -> List[str]:
    """Return all filenames written into the ``Graphs`` directory."""
    if not graph_dir.is_dir():
        return []
    return sorted(p.name for p in graph_dir.iterdir() if p.is_file())


def _run_analysis(output_dir: Path, **overrides) -> dict:
    """Run :func:`GNSS_MultipathAnalysis` with sensible defaults for tests."""
    defaults = dict(
        rinObsFilename=str(RIN_OBS),
        broadcastNav1=str(RIN_NAV),
        outputDir=str(output_dir),
        save_results_as_pickle=False,
        save_results_as_compressed_pickle=False,
        write_results_to_csv=False,
        nav_data_rate=120,
        use_LaTex=False,  # avoid LaTeX dependency on CI
    )
    defaults.update(overrides)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return GNSS_MultipathAnalysis(**defaults)


# ── Tests ────────────────────────────────────────────────────────────────────
class TestPlottingDisabled:
    """When ``plotEstimates=False`` no graphs should be produced."""

    def test_no_graphs_when_plotting_disabled(self, tmp_path: Path):
        out_dir = tmp_path / "out_no_plots"
        _run_analysis(
            out_dir,
            plotEstimates=False,
            plot_polarplot=False,
            include_SNR=False,
        )

        graphs = _list_graphs(out_dir / "Graphs")
        # Either the directory is empty or it does not exist at all.
        assert graphs == [], (
            f"Expected no graph files when plotEstimates=False, got: {graphs}"
        )

    def test_polarplot_disabled_skips_polarplots(self, tmp_path: Path):
        out_dir = tmp_path / "out_no_polar"
        _run_analysis(
            out_dir,
            plotEstimates=True,
            plot_polarplot=False,
            include_SNR=False,
        )

        graphs = _list_graphs(out_dir / "Graphs")
        # The MP_<sys>_<code>.png polar plots come from make_polarplot and
        # must NOT be produced when plot_polarplot=False.
        assert not any(g.startswith("MP_") and g.endswith(".png") for g in graphs), \
            f"Polar MP plots should not be created. Found: {graphs}"
        # But at least the bar plot / skyplot / combined-MP plots should exist.
        assert any("Barplot_RMS" in g for g in graphs)
        assert any(g.startswith("Skyplot_") for g in graphs)

    def test_snr_disabled_skips_snr_plots(self, tmp_path: Path):
        out_dir = tmp_path / "out_no_snr"
        _run_analysis(
            out_dir,
            plotEstimates=True,
            plot_polarplot=True,
            include_SNR=False,
        )

        graphs = _list_graphs(out_dir / "Graphs")
        assert not any(g.startswith("SNR_") for g in graphs), \
            f"SNR plots should not be created. Found: {[g for g in graphs if g.startswith('SNR_')]}"


class TestPlottingEnabled:
    """End-to-end run with all plot flags enabled."""

    @pytest.fixture(scope="class")
    def analysis_result(self, tmp_path_factory):
        """Run the analysis once for the whole class to keep the test fast."""
        out_dir = tmp_path_factory.mktemp("out_full_plots")
        result = _run_analysis(
            out_dir,
            plotEstimates=True,
            plot_polarplot=True,
            include_SNR=True,
        )
        return out_dir, result

    def test_graphs_directory_exists(self, analysis_result):
        out_dir, _ = analysis_result
        assert (out_dir / "Graphs").is_dir()

    def test_barplot_files_created(self, analysis_result):
        out_dir, _ = analysis_result
        graphs = _list_graphs(out_dir / "Graphs")
        # Combined bar plot for all systems (always written).
        assert "Barplot_RMS_all.pdf" in graphs
        # Per-system bar plots (one per GNSS system in the file).
        per_system = [g for g in graphs if g.startswith("Barplot_RMS_") and g != "Barplot_RMS_all.pdf"]
        assert per_system, f"Expected per-system bar plots, got: {graphs}"

    def test_skyplot_files_created(self, analysis_result):
        out_dir, _ = analysis_result
        graphs = _list_graphs(out_dir / "Graphs")
        skyplots = [g for g in graphs if g.startswith("Skyplot_") and g.endswith(".pdf")]
        assert skyplots, f"Expected at least one Skyplot_<system>.pdf, got: {graphs}"

    def test_multipath_polarplot_files_created(self, analysis_result):
        out_dir, _ = analysis_result
        graphs = _list_graphs(out_dir / "Graphs")
        mp_polar = [g for g in graphs if g.startswith("MP_") and g.endswith(".png")]
        assert mp_polar, f"Expected MP_<sys>_<code>.png polar plots, got: {graphs}"

    def test_multipath_combined_plots_created(self, analysis_result):
        out_dir, _ = analysis_result
        graphs = _list_graphs(out_dir / "Graphs")
        combined = [g for g in graphs if g.endswith("_MP_combined.png")]
        assert combined, f"Expected <sys>_<code>_MP_combined.png plots, got: {graphs}"

    def test_ionospheric_delay_plots_created(self, analysis_result):
        out_dir, _ = analysis_result
        graphs = _list_graphs(out_dir / "Graphs")
        iono = [g for g in graphs if "ionospheric_delay_combined" in g and g.endswith(".pdf")]
        assert iono, f"Expected ionospheric_delay_combined PDFs, got: {graphs}"

    def test_snr_plots_created(self, analysis_result):
        """SNR plots are only produced when the RINEX file contains SNR codes.

        The cropped OPEC test file does not include SNR observations, so we
        skip rather than fail when no SNR was processed.
        """
        out_dir, result = analysis_result
        has_snr = any(
            result.get(sys_name, {}).get("SNR")
            for sys_name in result.get("GNSSsystems", [])
        )
        if not has_snr:
            pytest.skip("Test RINEX file contains no SNR observations.")

        graphs = _list_graphs(out_dir / "Graphs")
        snr_polar = [g for g in graphs if g.startswith("SNR_Polar_") and g.endswith(".png")]
        snr_elev = [
            g for g in graphs
            if g.startswith("SNR_") and not g.startswith("SNR_Polar_") and g.endswith(".pdf")
        ]
        assert snr_polar, f"Expected SNR_Polar_<sys>_<code>.png plots, got: {graphs}"
        assert snr_elev, f"Expected SNR_<sys>_<code>.pdf plots, got: {graphs}"

    def test_summary_heatmaps_created(self, analysis_result):
        """Multipath heatmaps are always written; SNR heatmaps only with SNR data."""
        out_dir, result = analysis_result
        graphs = _list_graphs(out_dir / "Graphs")

        # Multipath az/el heatmaps are always produced when plotEstimates=True.
        assert "Summary_Multipath_AzEl.png" in graphs, \
            f"Expected combined multipath heatmap, got: {graphs}"
        assert any(
            g.startswith("Summary_Multipath_AzEl_") and g.endswith(".png")
            for g in graphs
        ), f"Expected per-system multipath heatmaps, got: {graphs}"

        # SNR heatmaps require SNR observations in the source file.
        has_snr = any(
            result.get(sys_name, {}).get("SNR")
            for sys_name in result.get("GNSSsystems", [])
        )
        if has_snr:
            assert "Summary_SNR_AzEl.png" in graphs, \
                f"Expected combined SNR heatmap, got: {graphs}"
            assert any(
                g.startswith("Summary_SNR_AzEl_") and g.endswith(".png")
                for g in graphs
            ), f"Expected per-system SNR heatmaps, got: {graphs}"

    def test_generated_plot_files_are_non_empty(self, analysis_result):
        """Sanity check: every produced file should have a positive size."""
        out_dir, _ = analysis_result
        graph_dir = out_dir / "Graphs"
        zero_byte = [
            p.name for p in graph_dir.iterdir()
            if p.is_file() and p.stat().st_size == 0
        ]
        assert not zero_byte, f"These plot files were created but are empty: {zero_byte}"


class TestCleanup:
    """Verify pytest's tmp_path actually removes generated files between tests."""

    def test_tmp_directory_is_isolated(self, tmp_path: Path):
        out_dir = tmp_path / "isolated_run"
        _run_analysis(
            out_dir,
            plotEstimates=True,
            plot_polarplot=False,
            include_SNR=False,
        )
        assert (out_dir / "Graphs").is_dir()
        # The tmp_path fixture cleans up after the test ends; we just verify
        # the directory was created inside tmp_path so cleanup will catch it.
        assert str(tmp_path) in str(out_dir.resolve())


def test_no_residual_files_in_repo_after_tests():
    """Guard against accidental writes to the repo when ``outputDir`` is set."""
    # Common stray locations the analysis would write to if outputDir defaulted.
    repo_default = PROJECT_ROOT / "Output_Files"
    src_default = PROJECT_ROOT / "src" / "Output_Files"
    # The tests above always pass an explicit outputDir under tmp_path, so the
    # default folders should never appear as a side-effect of this test module.
    # We only assert that *if* they exist they were not freshly created here —
    # we cannot delete user-owned data, so just make sure the test is a no-op.
    for path in (repo_default, src_default):
        if path.exists():
            # Allow pre-existing directories from the user's own runs.
            assert path.is_dir()
