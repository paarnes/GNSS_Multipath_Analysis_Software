"""
Tests for the data-rate (decimation) features of the RINEX readers.

Covers:
  * Observation reader (`readRinexObs`) ``desired_data_rate`` argument
    for both RINEX 3.04 and RINEX 2.11 files.
  * Internal helper :func:`_decimate_rinex_obs_data`.
  * Navigation reader (`RinexNav.read_nav`) ``data_rate`` argument.
"""
import os
import sys
import numpy as np
import pytest

project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.append(os.path.join(project_path, 'src'))

from gnssmultipath.readers.readRinexObs import (
    readRinexObs,
    _decimate_rinex_obs_data,
    RinexObsData,
)
from gnssmultipath.readers.RinexNav import RinexNav


# ── Observation-reader decimation (RINEX 3.04, native 30 s) ──────────────────


class TestObsDecimationRinex304:
    """Down-sampling a 30 s RINEX 3.04 file."""

    @pytest.fixture(scope="class")
    def native(self, rinex304_obs_file):
        return readRinexObs(rinex304_obs_file)

    def test_native_interval_is_30s(self, native):
        assert native.tInterval == pytest.approx(30.0, abs=0.5)

    @pytest.mark.parametrize("desired,expected_interval,expected_stride", [
        (60.0, 60.0, 2),
        (90.0, 90.0, 3),
        (300.0, 300.0, 10),
    ])
    def test_decimation_strides(self, rinex304_obs_file, native,
                                 desired, expected_interval, expected_stride):
        decimated = readRinexObs(rinex304_obs_file, desired_data_rate=desired)

        # Interval and epoch count
        assert decimated.tInterval == pytest.approx(expected_interval, abs=0.5)
        expected_n = (native.nepochs + expected_stride - 1) // expected_stride
        assert decimated.nepochs == expected_n

        # First epoch must match native's first epoch
        assert np.allclose(decimated.time_epochs[0], native.time_epochs[0])

        # Every retained epoch is exactly stride positions away in the native
        # array (i.e. index i in decimated == index i*stride in native).
        for i, row in enumerate(decimated.time_epochs):
            assert np.allclose(row, native.time_epochs[i * expected_stride])

    def test_no_decimation_when_rate_smaller(self, rinex304_obs_file, native):
        """Asking for a finer rate than native must keep every epoch."""
        same = readRinexObs(rinex304_obs_file, desired_data_rate=10.0)
        assert same.nepochs == native.nepochs
        assert same.tInterval == pytest.approx(native.tInterval, abs=0.5)

    def test_no_decimation_when_rate_equal(self, rinex304_obs_file, native):
        same = readRinexObs(rinex304_obs_file, desired_data_rate=native.tInterval)
        assert same.nepochs == native.nepochs

    def test_default_keeps_all_epochs(self, rinex304_obs_file, native):
        default = readRinexObs(rinex304_obs_file)
        assert default.nepochs == native.nepochs

    def test_decimated_obs_dict_keys_are_contiguous(self, rinex304_obs_file):
        decimated = readRinexObs(rinex304_obs_file, desired_data_rate=300.0)
        for sys_char, epoch_dict in decimated.GNSS_obs.items():
            keys = sorted(epoch_dict.keys())
            assert keys == list(range(1, decimated.nepochs + 1)), (
                f"System {sys_char!r} epoch keys not contiguous 1..N: {keys[:5]}..."
            )

    def test_decimated_observations_match_native(self, rinex304_obs_file, native):
        """A decimated observation must be byte-identical to the corresponding native epoch."""
        stride = 4
        decimated = readRinexObs(rinex304_obs_file,
                                 desired_data_rate=native.tInterval * stride)
        for sys_char, native_dict in native.GNSS_obs.items():
            dec_dict = decimated.GNSS_obs.get(sys_char, {})
            for new_e, old_e in zip(range(1, decimated.nepochs + 1),
                                    range(1, native.nepochs + 1, stride)):
                assert np.array_equal(dec_dict[new_e], native_dict[old_e]), (
                    f"System {sys_char} epoch mismatch (new={new_e}, old={old_e})"
                )

    def test_gnss_svs_rows_are_decimated(self, rinex304_obs_file, native):
        decimated = readRinexObs(rinex304_obs_file, desired_data_rate=60.0)
        for sys_char, native_arr in native.GNSS_SVs.items():
            dec_arr = decimated.GNSS_SVs[sys_char]
            assert dec_arr.shape[0] == decimated.nepochs
            assert dec_arr.shape[1] == native_arr.shape[1]
            # Stride 2: row i in decimated == row i*2 in native
            for i in range(decimated.nepochs):
                assert np.array_equal(dec_arr[i], native_arr[i * 2])

    def test_tlastobs_updated(self, rinex304_obs_file):
        decimated = readRinexObs(rinex304_obs_file, desired_data_rate=300.0)
        # tLastObs must correspond to the last retained epoch (within 1 s)
        last_tow = decimated.time_epochs[-1, 1]
        last_sec = float(np.asarray(decimated.tLastObs[5]).item())
        assert last_sec == pytest.approx(last_tow % 60, abs=1.0)

    def test_negative_rate_raises(self, rinex304_obs_file):
        with pytest.raises(ValueError):
            readRinexObs(rinex304_obs_file, desired_data_rate=-1)

    def test_zero_rate_raises(self, rinex304_obs_file):
        with pytest.raises(ValueError):
            readRinexObs(rinex304_obs_file, desired_data_rate=0)


# ── Observation-reader decimation (RINEX 2.11) ───────────────────────────────


class TestObsDecimationRinex211:
    """Down-sampling a RINEX 2.11 file."""

    @pytest.fixture(scope="class")
    def native(self, rinex211_obs_file):
        return readRinexObs(rinex211_obs_file)

    def test_decimation_halves_epochs(self, rinex211_obs_file, native):
        stride = 2
        decimated = readRinexObs(rinex211_obs_file,
                                 desired_data_rate=native.tInterval * stride)
        expected_n = (native.nepochs + stride - 1) // stride
        assert decimated.nepochs == expected_n
        assert decimated.tInterval == pytest.approx(native.tInterval * stride, abs=0.5)

    def test_decimation_observations_match_native(self, rinex211_obs_file, native):
        stride = 3
        decimated = readRinexObs(rinex211_obs_file,
                                 desired_data_rate=native.tInterval * stride)
        for sys_char, native_dict in native.GNSS_obs.items():
            if not native_dict:
                continue
            dec_dict = decimated.GNSS_obs.get(sys_char, {})
            for new_e, old_e in zip(range(1, decimated.nepochs + 1),
                                    range(1, native.nepochs + 1, stride)):
                assert np.array_equal(dec_dict[new_e], native_dict[old_e])


# ── Direct unit tests for the helper ─────────────────────────────────────────


def _make_synthetic_rinex_obs_data(nepochs=10, t_interval=1.0):
    """Build a tiny synthetic RinexObsData object for unit testing."""
    GNSS_obs = {'G': {i: np.full((37, 2), float(i)) for i in range(1, nepochs + 1)}}
    GNSS_LLI = {'G': {i: np.zeros((37, 2)) for i in range(1, nepochs + 1)}}
    GNSS_SS  = {'G': {i: np.zeros((37, 2)) for i in range(1, nepochs + 1)}}
    GNSS_SVs = {'G': np.arange(nepochs * 5).reshape(nepochs, 5).astype(np.int16)}
    time_epochs = np.column_stack([
        np.full(nepochs, 2200),
        np.arange(nepochs) * t_interval,
    ]).astype(float)
    return RinexObsData(
        GNSS_obs=GNSS_obs,
        GNSS_LLI=GNSS_LLI,
        GNSS_SS=GNSS_SS,
        GNSS_SVs=GNSS_SVs,
        time_epochs=time_epochs,
        nepochs=nepochs,
        GNSSsystems={1: 'G'},
        obsCodes={1: {'G': ['C1C', 'L1C']}},
        tInterval=float(t_interval),
        tFirstObs=np.array([[2024], [1], [1], [0], [0], [0]], dtype=float),
        tLastObs=np.array([[2024], [1], [1], [0], [0], [float(nepochs - 1) * t_interval]], dtype=float),
    )


class TestDecimateHelper:
    """Direct tests for :func:`_decimate_rinex_obs_data`."""

    def test_keeps_first_epoch(self):
        data = _make_synthetic_rinex_obs_data(nepochs=20, t_interval=1.0)
        out = _decimate_rinex_obs_data(data, desired_data_rate=5.0)
        assert out.nepochs == 4  # epochs 0,5,10,15
        assert out.tInterval == pytest.approx(5.0)
        assert np.array_equal(out.GNSS_obs['G'][1], data.GNSS_obs['G'][1])
        assert np.array_equal(out.GNSS_obs['G'][2], data.GNSS_obs['G'][6])
        assert np.array_equal(out.GNSS_obs['G'][4], data.GNSS_obs['G'][16])

    def test_returns_same_when_stride_one(self):
        data = _make_synthetic_rinex_obs_data(nepochs=5, t_interval=10.0)
        out = _decimate_rinex_obs_data(data, desired_data_rate=10.0)
        # Stride 1 => same object returned (no copy)
        assert out is data

    def test_returns_same_when_finer_rate(self):
        data = _make_synthetic_rinex_obs_data(nepochs=5, t_interval=30.0)
        out = _decimate_rinex_obs_data(data, desired_data_rate=10.0)
        assert out is data

    def test_negative_rate_raises(self):
        data = _make_synthetic_rinex_obs_data()
        with pytest.raises(ValueError):
            _decimate_rinex_obs_data(data, desired_data_rate=-1.0)

    def test_zero_nepochs_returns_unchanged(self):
        data = _make_synthetic_rinex_obs_data(nepochs=0, t_interval=1.0)
        out = _decimate_rinex_obs_data(data, desired_data_rate=10.0)
        assert out.nepochs == 0

    def test_nan_interval_returns_unchanged(self):
        data = _make_synthetic_rinex_obs_data()
        data.tInterval = float('nan')
        out = _decimate_rinex_obs_data(data, desired_data_rate=10.0)
        assert out is data

    def test_time_epochs_decimated(self):
        data = _make_synthetic_rinex_obs_data(nepochs=12, t_interval=1.0)
        out = _decimate_rinex_obs_data(data, desired_data_rate=4.0)
        assert out.time_epochs.shape == (3, 2)  # epochs 0,4,8
        assert out.time_epochs[0, 1] == 0
        assert out.time_epochs[1, 1] == 4
        assert out.time_epochs[2, 1] == 8

    def test_gnss_svs_decimated(self):
        data = _make_synthetic_rinex_obs_data(nepochs=10, t_interval=1.0)
        out = _decimate_rinex_obs_data(data, desired_data_rate=2.0)
        assert out.GNSS_SVs['G'].shape == (5, 5)
        assert np.array_equal(out.GNSS_SVs['G'][0], data.GNSS_SVs['G'][0])
        assert np.array_equal(out.GNSS_SVs['G'][1], data.GNSS_SVs['G'][2])


# ── Navigation-reader data_rate ──────────────────────────────────────────────


class TestNavDataRate:
    """`RinexNav.read_nav` ``data_rate`` argument (interval in minutes)."""

    @pytest.fixture(scope="class")
    def all_records(self, nav_mixed_file):
        return RinexNav.read_nav(nav_mixed_file, data_rate=0)

    def test_zero_disables_filtering(self, all_records):
        assert all_records.nepochs > 0

    def test_higher_rate_reduces_records(self, nav_mixed_file, all_records):
        coarse = RinexNav.read_nav(nav_mixed_file, data_rate=120)
        assert coarse.nepochs <= all_records.nepochs

    def test_rate_is_per_satellite(self, nav_mixed_file):
        """Each retained satellite must be spaced by at least *data_rate* minutes."""
        from datetime import datetime
        interval_min = 60
        nav = RinexNav.read_nav(nav_mixed_file, data_rate=interval_min)
        eph = nav.ephemerides
        if hasattr(eph, "values"):
            eph = eph.values
        # group epochs per PRN
        per_sat = {}
        for row in eph:
            prn = str(row[0])
            epoch = datetime(int(float(row[1])), int(float(row[2])), int(float(row[3])),
                             int(float(row[4])), int(float(row[5])), int(float(row[6])))
            per_sat.setdefault(prn, []).append(epoch)
        for prn, epochs in per_sat.items():
            epochs.sort()
            for prev, curr in zip(epochs, epochs[1:]):
                gap_min = (curr - prev).total_seconds() / 60.0
                # Allow tiny rounding tolerance
                assert gap_min >= interval_min - 1e-6, (
                    f"PRN {prn}: gap {gap_min} min < {interval_min} min"
                )

    def test_default_data_rate(self, nav_mixed_file, all_records):
        """Default data_rate=30 should keep <= all-records count."""
        default = RinexNav.read_nav(nav_mixed_file)
        assert 0 < default.nepochs <= all_records.nepochs
