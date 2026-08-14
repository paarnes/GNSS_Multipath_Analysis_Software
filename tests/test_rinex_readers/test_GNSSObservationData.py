"""
Tests for GNSSObservationData — the pythonic observation data accessor.
"""
import sys
import os
import pytest
import numpy as np

project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.append(os.path.join(project_path, 'src'))

from gnssmultipath.readers.readRinexObs import readRinexObs
from gnssmultipath.readers.GNSSObservationData import (
    GNSSObservationData, SystemObservations, _CodeAccessor,
)


class TestGNSSObservationDataFromSynthetic:
    """Unit tests with hand-crafted synthetic data (no I/O)."""

    @pytest.fixture(autouse=True, scope="class")
    def _build_store(self):
        cls = type(self)
        n_epochs, n_sat, n_codes = 5, 4, 3
        rng = np.random.default_rng(42)

        # Build epoch dicts (1-based keys like readRinexObs produces)
        gps_obs = {ep: rng.random((n_sat, n_codes)) for ep in range(1, n_epochs + 1)}
        gal_obs = {ep: rng.random((n_sat, 2)) for ep in range(1, n_epochs + 1)}
        gps_lli = {ep: np.zeros((n_sat, n_codes)) for ep in range(1, n_epochs + 1)}
        gps_ss  = {ep: rng.random((n_sat, n_codes)) for ep in range(1, n_epochs + 1)}

        cls.gnss_obs = {'G': gps_obs, 'E': gal_obs}
        cls.gnss_lli = {'G': gps_lli, 'E': np.nan}
        cls.gnss_ss  = {'G': gps_ss,  'E': np.nan}
        cls.obs_codes = {
            1: {'G': ['C1C', 'L1C', 'S1C']},
            2: {'E': ['C1X', 'L1X']},
        }
        cls.gnss_systems = {1: 'G', 2: 'E'}

        cls.store = GNSSObservationData(
            cls.gnss_obs, cls.obs_codes, cls.gnss_systems,
            cls.gnss_lli, cls.gnss_ss,
        )
        cls.n_epochs = n_epochs
        cls.n_sat = n_sat

    # ── Store-level tests ────────────────────────────────────────────────

    def test_systems(self):
        assert self.store.systems == ['G', 'E']

    def test_len(self):
        assert len(self.store) == 2

    def test_contains(self):
        assert 'G' in self.store
        assert 'E' in self.store
        assert 'R' not in self.store

    def test_iter(self):
        items = list(self.store)
        assert len(items) == 2
        assert all(isinstance(s, SystemObservations) for s in items)

    def test_getitem_invalid_raises(self):
        with pytest.raises(KeyError, match="System 'R' not available"):
            self.store['R']

    def test_gps_property(self):
        assert self.store.gps.system_code == 'G'

    def test_galileo_property(self):
        assert self.store.galileo.system_code == 'E'

    def test_glonass_property_raises(self):
        with pytest.raises(KeyError):
            self.store.glonass

    def test_beidou_property_raises(self):
        with pytest.raises(KeyError):
            self.store.beidou

    # ── SystemObservations tests (GPS) ───────────────────────────────────

    def test_gps_codes(self):
        assert self.store.gps.codes == ['C1C', 'L1C', 'S1C']

    def test_gps_pseudorange_codes(self):
        assert self.store.gps.pseudorange_codes == ['C1C']

    def test_gps_phase_codes(self):
        assert self.store.gps.phase_codes == ['L1C']

    def test_gps_snr_codes(self):
        assert self.store.gps.snr_codes == ['S1C']

    def test_gps_doppler_codes(self):
        assert self.store.gps.doppler_codes == []

    def test_gps_bands(self):
        assert self.store.gps.bands == ['1']

    def test_gps_n_epochs(self):
        assert self.store.gps.n_epochs == self.n_epochs

    def test_gps_n_satellites(self):
        assert self.store.gps.n_satellites == self.n_sat

    def test_gps_system_name(self):
        assert self.store.gps.system_name == 'GPS'

    def test_gal_system_name(self):
        assert self.store.galileo.system_name == 'Galileo'

    # ── Code-level data access ───────────────────────────────────────────

    def test_getitem_shape(self):
        c1c = self.store.gps['C1C']
        assert c1c.shape == (self.n_epochs, self.n_sat)

    def test_getitem_values_match_raw(self):
        c1c = self.store.gps['C1C']
        # Must equal column 0 of the stacked raw data
        stacked = np.stack(list(self.gnss_obs['G'].values()))
        np.testing.assert_array_equal(c1c, stacked[:, :, 0])

    def test_getitem_second_code(self):
        l1c = self.store.gps['L1C']
        stacked = np.stack(list(self.gnss_obs['G'].values()))
        np.testing.assert_array_equal(l1c, stacked[:, :, 1])

    def test_getitem_invalid_code_raises(self):
        with pytest.raises(KeyError, match="'C5X' not available"):
            self.store.gps['C5X']

    def test_contains_code(self):
        assert 'C1C' in self.store.gps
        assert 'L1C' in self.store.gps
        assert 'C5X' not in self.store.gps

    def test_data_property_shape(self):
        data = self.store.gps.data
        assert data.shape == (self.n_epochs, self.n_sat, 3)

    def test_band_filter(self):
        band1 = self.store.gps.band(1)
        assert sorted(band1.keys()) == ['C1C', 'L1C', 'S1C']
        assert band1['C1C'].shape == (self.n_epochs, self.n_sat)

    def test_band_filter_empty(self):
        band5 = self.store.gps.band(5)
        assert band5 == {}

    def test_by_type(self):
        pseudoranges = self.store.gps.by_type('C')
        assert list(pseudoranges.keys()) == ['C1C']
        assert pseudoranges['C1C'].shape == (self.n_epochs, self.n_sat)

    # ── LLI / SS sub-accessors ───────────────────────────────────────────

    def test_lli_accessor(self):
        lli_l1c = self.store.gps.lli['L1C']
        assert lli_l1c.shape == (self.n_epochs, self.n_sat)
        np.testing.assert_array_equal(lli_l1c, 0)  # all zeros

    def test_ss_accessor(self):
        ss_s1c = self.store.gps.ss['S1C']
        assert ss_s1c.shape == (self.n_epochs, self.n_sat)

    def test_galileo_lli_unavailable(self):
        with pytest.raises(AttributeError, match="LLI data not available"):
            self.store.galileo.lli

    def test_galileo_ss_unavailable(self):
        with pytest.raises(AttributeError, match="Signal-Strength data not available"):
            self.store.galileo.ss

    # ── repr ─────────────────────────────────────────────────────────────

    def test_system_repr(self):
        r = repr(self.store.gps)
        assert "GPS" in r
        assert "C1C" in r

    def test_store_repr(self):
        r = repr(self.store)
        assert "GNSSObservationData" in r


class TestGNSSObservationDataFromRealData:
    """Integration tests using a real RINEX 3.04 observation file."""

    @pytest.fixture(autouse=True, scope="class")
    def _read_obs(self, rinex304_obs_file):
        cls = type(self)
        cls.rinex = readRinexObs(rinex304_obs_file)
        cls.store = cls.rinex.observations

    def test_store_created(self):
        assert isinstance(self.store, GNSSObservationData)

    def test_systems_present(self):
        systems = self.store.systems
        assert 'G' in systems
        assert 'E' in systems

    def test_gps_codes_not_empty(self):
        assert len(self.store.gps.codes) > 0

    def test_gps_has_c1c(self):
        assert 'C1C' in self.store.gps

    def test_gps_c1c_shape(self):
        c1c = self.store.gps['C1C']
        assert c1c.ndim == 2
        assert c1c.shape[0] == self.rinex.nepochs

    def test_gps_c1c_matches_raw(self):
        """Verify the store provides the same data as the raw dict."""
        gps_codes = self.rinex.obsCodes[
            next(k for k, v in self.rinex.GNSSsystems.items() if v == 'G')
        ]['G']
        c1c_col = gps_codes.index('C1C')
        raw_stacked = np.stack(list(self.rinex.GNSS_obs['G'].values()))
        np.testing.assert_array_equal(
            self.store.gps['C1C'],
            raw_stacked[:, :, c1c_col],
        )

    def test_gps_pseudorange_codes_start_with_C(self):
        for code in self.store.gps.pseudorange_codes:
            assert code[0] == 'C'

    def test_gps_phase_codes_start_with_L(self):
        for code in self.store.gps.phase_codes:
            assert code[0] == 'L'

    def test_galileo_band_1(self):
        if 'E' in self.store:
            band1 = self.store.galileo.band(1)
            for code in band1:
                assert code[1] == '1'

    def test_n_epochs_matches(self):
        assert self.store.gps.n_epochs == self.rinex.nepochs

    def test_observations_property_caches(self):
        """Accessing .observations twice returns the same object."""
        obs1 = self.rinex.observations
        obs2 = self.rinex.observations
        assert obs1 is obs2

    def test_lli_available_if_read(self):
        """LLI data should be accessible if the file was read with LLI."""
        if isinstance(self.rinex.GNSS_LLI, dict) and 'G' in self.rinex.GNSS_LLI:
            gps = self.store.gps
            if gps.phase_codes:
                lli = gps.lli[gps.phase_codes[0]]
                assert lli.ndim == 2


class TestSelectAndSignals:
    """Filtering, signal metadata and carrier frequencies on real data."""

    @pytest.fixture(autouse=True, scope="class")
    def _read_obs(self, rinex304_obs_file):
        cls = type(self)
        cls.rinex = readRinexObs(rinex304_obs_file)
        cls.store = cls.rinex.observations

    # ── select() ─────────────────────────────────────────────────────────

    def test_select_by_type(self):
        codes = self.store.gps.select(obs_type='L')
        assert codes and all(c[0] == 'L' for c in codes)

    def test_select_by_band(self):
        codes = self.store.gps.select(band=1)
        assert codes and all(c[1] == '1' for c in codes)

    def test_select_by_attribute(self):
        codes = self.store.gps.select(attribute='C')
        assert all(c[2] == 'C' for c in codes)

    def test_select_combined(self):
        codes = self.store.gps.select(obs_type='C', band=1)
        assert all(c[0] == 'C' and c[1] == '1' for c in codes)

    def test_select_no_match_returns_empty(self):
        assert self.store.gps.select(band=9) == []

    def test_select_pseudorange_matches_property(self):
        gps = self.store.gps
        assert gps.select(obs_type='C') == gps.pseudorange_codes

    def test_store_select_across_systems(self):
        result = self.store.select(obs_type='L', band=1)
        assert 'G' in result
        for sys_code, codes in result.items():
            assert all(c[0] == 'L' and c[1] == '1' for c in codes)

    def test_store_select_single_system(self):
        result = self.store.select(system='G', band=1)
        assert list(result) == ['G']

    def test_store_select_omits_systems_without_match(self):
        result = self.store.select(band=9)
        assert result == {}

    # ── signal metadata ──────────────────────────────────────────────────

    def test_signal_returns_obscode(self):
        sig = self.store.gps.signal('C1C')
        assert sig.code == 'C1C'
        assert sig.system == 'G'
        assert sig.band == 1
        assert sig.is_pseudorange

    def test_signal_is_cached(self):
        gps = self.store.gps
        assert gps.signal('C1C') is gps.signal('C1C')

    def test_signals_covers_all_codes(self):
        gps = self.store.gps
        assert [s.code for s in gps.signals] == gps.codes

    def test_signal_unknown_code_raises(self):
        with pytest.raises(KeyError, match="not available"):
            self.store.gps.signal('C9Z')

    # ── frequencies ──────────────────────────────────────────────────────

    def test_gps_frequency(self):
        assert self.store.gps.frequency('C1C') == pytest.approx(1575.42e6)

    def test_gps_wavelength(self):
        gps = self.store.gps
        assert gps.wavelength('C1C') == pytest.approx(
            299792458.0 / gps.frequency('C1C'))

    def test_galileo_frequency(self):
        assert self.store.galileo.frequency('C5X') == pytest.approx(1176.45e6)

    def test_glonass_frequency_per_prn(self):
        glo = self.store.glonass
        prn = glo.prns[0]
        channel = glo.glonass_channel(prn)
        assert channel is not None
        assert glo.frequency('C1C', prn=prn) == pytest.approx(
            1602.0e6 + channel * 562500.0)

    def test_glonass_frequency_array(self):
        glo = self.store.glonass
        freqs = glo.frequency('C1C')
        assert freqs.shape == (glo.n_satellites,)
        assert np.isfinite(freqs[glo.prns[0]])

    def test_glonass_frequency_unknown_prn_raises(self):
        glo = self.store.glonass
        unknown = next(p for p in range(1, glo.n_satellites)
                       if glo.glonass_channel(p) is None)
        with pytest.raises(KeyError, match="No GLONASS channel number"):
            glo.frequency('C1C', prn=unknown)


class TestDataRetrieval:
    """get(), per-satellite views, PRNs and epoch times on real data."""

    @pytest.fixture(autouse=True, scope="class")
    def _read_obs(self, rinex304_obs_file):
        cls = type(self)
        cls.rinex = readRinexObs(rinex304_obs_file)
        cls.store = cls.rinex.observations

    # ── get() ────────────────────────────────────────────────────────────

    def test_get_masks_zeros(self):
        gps = self.store.gps
        raw = gps['C1C']
        masked = gps.get('C1C')
        assert np.array_equal(raw == 0, np.isnan(masked))

    def test_get_keeps_observed_values(self):
        gps = self.store.gps
        raw = gps['C1C']
        masked = gps.get('C1C')
        observed = raw != 0
        np.testing.assert_array_equal(raw[observed], masked[observed])

    def test_get_without_masking_matches_raw(self):
        gps = self.store.gps
        np.testing.assert_array_equal(gps.get('C1C', mask_missing=False), gps['C1C'])

    def test_get_returns_independent_copy(self):
        gps = self.store.gps
        before = gps['C1C'][0, 1]
        gps.get('C1C')[0, 1] = -1.0
        assert gps['C1C'][0, 1] == before

    # ── column caching ───────────────────────────────────────────────────

    def test_column_is_cached(self):
        gps = self.store.gps
        assert gps['C1C'] is gps['C1C']

    def test_clear_cache_releases_column(self):
        gps = self.store.gps
        first = gps['C1C']
        gps._obs.clear_cache()
        assert gps['C1C'] is not first

    def test_data_matches_columns(self):
        gal = self.store.galileo
        cube = gal.data
        for idx, code in enumerate(gal.codes):
            np.testing.assert_array_equal(cube[:, :, idx], gal[code])

    # ── PRNs ─────────────────────────────────────────────────────────────

    def test_prns_are_sorted_and_positive(self):
        prns = self.store.gps.prns
        assert prns == sorted(prns)
        assert all(p > 0 for p in prns)

    def test_prns_match_gnss_svs(self):
        expected = {int(p) for p in np.unique(self.rinex.GNSS_SVs['G'][:, 1:]) if p > 0}
        assert set(self.store.gps.prns) == expected

    # ── satellite view ───────────────────────────────────────────────────

    def test_sat_returns_1d_slice(self):
        gps = self.store.gps
        prn = gps.prns[0]
        np.testing.assert_array_equal(gps.sat(prn)['C1C'], gps['C1C'][:, prn])

    def test_sat_id_and_prn(self):
        sat = self.store.gps.sat(5)
        assert sat.prn == 5
        assert sat.sv_id == 'G05'
        assert sat.system_name == 'GPS'

    def test_sat_get_masks(self):
        gps = self.store.gps
        prn = gps.prns[0]
        np.testing.assert_array_equal(
            gps.sat(prn).get('C1C'), gps.get('C1C')[:, prn])

    def test_sat_frequency(self):
        assert self.store.gps.sat(1).frequency('C1C') == pytest.approx(1575.42e6)

    def test_sat_lli_slice(self):
        gps = self.store.gps
        prn = gps.prns[0]
        code = gps.phase_codes[0]
        np.testing.assert_array_equal(
            gps.sat(prn).lli[code], gps.lli[code][:, prn])

    def test_sat_out_of_range_raises(self):
        with pytest.raises(KeyError, match="out of range"):
            self.store.gps.sat(999)

    def test_sat_contains_code(self):
        assert 'C1C' in self.store.gps.sat(1)

    # ── epoch times ──────────────────────────────────────────────────────

    def test_store_time_epochs(self):
        assert self.store.time_epochs.shape == (self.rinex.nepochs, 2)

    def test_store_datetimes_length(self):
        assert len(self.store.datetimes) == self.rinex.nepochs

    def test_datetimes_are_increasing(self):
        dt = self.store.datetimes
        assert np.all(np.diff(dt).astype('timedelta64[s]').astype(int) > 0)

    def test_system_epoch_times_shared(self):
        np.testing.assert_array_equal(
            self.store.gps.epoch_times, self.store.time_epochs)

    def test_store_metadata(self):
        assert self.store.n_epochs == self.rinex.nepochs
        assert self.store.interval == self.rinex.tInterval
        assert self.store.approx_position is not None

    # ── summary ──────────────────────────────────────────────────────────

    def test_system_summary_mentions_codes(self):
        text = self.store.gps.summary()
        assert 'GPS' in text
        assert 'C1C' in text

    def test_store_summary_covers_all_systems(self):
        text = self.store.summary()
        for sys_code in self.store.systems:
            assert self.store[sys_code].system_name in text


class TestEpochObservations:
    """Single-epoch views: all signals of one epoch."""

    @pytest.fixture(autouse=True, scope="class")
    def _read_obs(self, rinex304_obs_file):
        cls = type(self)
        cls.rinex = readRinexObs(rinex304_obs_file)
        cls.store = cls.rinex.observations

    # ── Indexing ─────────────────────────────────────────────────────────

    def test_index_and_number(self):
        ep = self.store.gps.epoch(34)
        assert ep.index == 34
        assert ep.number == 35

    def test_negative_index_is_normalised(self):
        gps = self.store.gps
        ep = gps.epoch(-1)
        assert ep.index == gps.n_epochs - 1

    def test_out_of_range_raises(self):
        gps = self.store.gps
        with pytest.raises(IndexError, match="out of range"):
            gps.epoch(gps.n_epochs)

    def test_negative_out_of_range_raises(self):
        gps = self.store.gps
        with pytest.raises(IndexError, match="out of range"):
            gps.epoch(-gps.n_epochs - 1)

    # ── Data access ──────────────────────────────────────────────────────

    def test_getitem_matches_column_slice(self):
        gps = self.store.gps
        np.testing.assert_array_equal(gps.epoch(34)['C1C'], gps['C1C'][34])

    def test_matrix_shape(self):
        gps = self.store.gps
        assert gps.epoch(0).matrix.shape == (gps.n_satellites, len(gps.codes))

    def test_matrix_is_the_stored_block(self):
        """The epoch matrix is the reader's own array, not a copy."""
        gps = self.store.gps
        assert gps.epoch(0).matrix is self.rinex.GNSS_obs['G'][1]

    def test_matrix_matches_all_codes(self):
        gps = self.store.gps
        ep = gps.epoch(10)
        for idx, code in enumerate(gps.codes):
            np.testing.assert_array_equal(ep.matrix[:, idx], ep[code])

    def test_get_masks_missing(self):
        gps = self.store.gps
        ep = gps.epoch(34)
        raw, masked = ep['C1C'], ep.get('C1C')
        assert np.array_equal(raw == 0, np.isnan(masked))

    def test_get_returns_independent_copy(self):
        gps = self.store.gps
        ep = gps.epoch(34)
        before = ep['C1C'][1]
        ep.get('C1C')[1] = -1.0
        assert ep['C1C'][1] == before

    def test_unknown_code_raises(self):
        with pytest.raises(KeyError, match="not available"):
            self.store.gps.epoch(0)['C9Z']

    def test_contains_code(self):
        assert 'C1C' in self.store.gps.epoch(0)

    # ── Satellites ───────────────────────────────────────────────────────

    def test_sat_returns_all_signals(self):
        gps = self.store.gps
        ep = gps.epoch(34)
        prn = ep.prns[0]
        signals = ep.sat(prn)
        assert set(signals) == set(gps.codes)
        assert signals['C1C'] == gps['C1C'][34, prn]

    def test_prns_are_observed_in_that_epoch(self):
        gps = self.store.gps
        ep = gps.epoch(34)
        assert ep.prns
        assert set(ep.prns).issubset(set(gps.prns))
        for prn in ep.prns:
            assert np.any(ep.matrix[prn] != 0)

    def test_prns_match_gnss_svs(self):
        gps = self.store.gps
        svs = self.rinex.GNSS_SVs['G']
        expected = {int(p) for p in svs[34, 1:svs[34, 0] + 1]}
        assert set(gps.epoch(34).prns) == expected

    # ── Time ─────────────────────────────────────────────────────────────

    def test_time_matches_time_epochs(self):
        ep = self.store.gps.epoch(34)
        np.testing.assert_array_equal(np.asarray(ep.time),
                                      self.rinex.time_epochs[34])

    def test_datetime_matches_store(self):
        assert self.store.gps.epoch(34).datetime == self.store.datetimes[34]

    # ── LLI / SS ─────────────────────────────────────────────────────────

    def test_lli_slice(self):
        gps = self.store.gps
        code = gps.phase_codes[0]
        np.testing.assert_array_equal(gps.epoch(34).lli[code], gps.lli[code][34])

    def test_ss_slice(self):
        gps = self.store.gps
        code = gps.codes[0]
        np.testing.assert_array_equal(gps.epoch(34).ss[code], gps.ss[code][34])

    # ── Export ───────────────────────────────────────────────────────────

    def test_dataframe_is_satellites_by_signals(self):
        gps = self.store.gps
        ep = gps.epoch(34)
        frame = ep.to_dataframe()
        assert list(frame.columns) == ['prn'] + gps.codes
        assert len(frame) == len(ep.prns)
        assert frame.index.name == 'sv'

    def test_dataframe_values_match(self):
        gps = self.store.gps
        ep = gps.epoch(34)
        frame = ep.to_dataframe(codes=['C1C'])
        prn = int(frame['prn'].iloc[0])
        assert frame['C1C'].iloc[0] == pytest.approx(gps['C1C'][34, prn])

    def test_dataframe_subset(self):
        gps = self.store.gps
        ep = gps.epoch(34)
        frame = ep.to_dataframe(codes=['C1C', 'L1C'], prns=ep.prns[:2])
        assert list(frame.columns) == ['prn', 'C1C', 'L1C']
        assert len(frame) == 2

    def test_empty_selection_returns_empty_frame(self):
        assert self.store.gps.epoch(0).to_dataframe(prns=[]).empty

    def test_repr(self):
        text = repr(self.store.gps.epoch(34))
        assert 'epoch 35' in text
        assert 'satellites' in text


class TestDataFrameExport:

    @pytest.fixture(autouse=True, scope="class")
    def _read_obs(self, rinex304_obs_file):
        cls = type(self)
        cls.rinex = readRinexObs(rinex304_obs_file)
        cls.store = cls.rinex.observations

    def test_columns(self):
        gps = self.store.gps
        frame = gps.to_dataframe(codes=['C1C'], prns=gps.prns[:1])
        assert list(frame.columns) == ['epoch', 'datetime', 'sv', 'prn', 'code', 'value']

    def test_values_match_arrays(self):
        gps = self.store.gps
        prn = gps.prns[0]
        frame = gps.to_dataframe(codes=['C1C'], prns=[prn])
        expected = gps.get('C1C')[:, prn]
        np.testing.assert_allclose(frame['value'].to_numpy(),
                                   expected[np.isfinite(expected)])

    def test_dropna_removes_missing(self):
        gps = self.store.gps
        prn = gps.prns[0]
        kept = gps.to_dataframe(codes=['C1C'], prns=[prn])
        full = gps.to_dataframe(codes=['C1C'], prns=[prn], dropna=False)
        assert len(full) == gps.n_epochs
        assert len(kept) <= len(full)

    def test_sv_identifier(self):
        gps = self.store.gps
        frame = gps.to_dataframe(codes=['C1C'], prns=[gps.prns[0]])
        assert frame['sv'].iloc[0] == f"G{gps.prns[0]:02d}"

    def test_multiple_codes(self):
        gps = self.store.gps
        frame = gps.to_dataframe(codes=['C1C', 'L1C'], prns=gps.prns[:2])
        assert set(frame['code']) == {'C1C', 'L1C'}

    def test_include_lli_and_ss(self):
        gps = self.store.gps
        frame = gps.to_dataframe(codes=['L1C'], prns=gps.prns[:1],
                                 include_lli=True, include_ss=True)
        assert 'lli' in frame.columns
        assert 'ss' in frame.columns

    def test_satellite_dataframe(self):
        gps = self.store.gps
        prn = gps.prns[0]
        frame = gps.sat(prn).to_dataframe(codes=['C1C'])
        assert set(frame['prn']) == {prn}

    def test_store_dataframe_has_system_column(self):
        frame = self.store.to_dataframe(systems=['G', 'E'], codes=['C1C', 'C1X'])
        assert frame.columns[0] == 'system'
        assert set(frame['system']) == {'G', 'E'}

    def test_empty_selection_returns_empty_frame(self):
        gps = self.store.gps
        assert gps.to_dataframe(codes=['C1C'], prns=[]).empty


class TestRinex211Compatibility:
    """Two-character RINEX 2 observation codes."""

    @pytest.fixture(autouse=True, scope="class")
    def _read_obs(self, rinex211_obs_file):
        cls = type(self)
        cls.rinex = readRinexObs(rinex211_obs_file)
        cls.store = cls.rinex.observations

    def test_codes_are_two_characters(self):
        for code in self.store.gps.codes:
            assert len(code) == 2

    def test_p_codes_counted_as_pseudorange(self):
        gps = self.store.gps
        p_codes = [c for c in gps.codes if c[0] == 'P']
        if p_codes:
            assert set(p_codes).issubset(set(gps.pseudorange_codes))

    def test_signal_parses_two_char_code(self):
        gps = self.store.gps
        sig = gps.signal(gps.codes[0])
        assert sig.attribute == ''

    def test_select_by_empty_attribute(self):
        gps = self.store.gps
        assert gps.select(attribute='') == gps.codes

    def test_frequency_available(self):
        gps = self.store.gps
        code = next(c for c in gps.codes if c[1] in '125')
        assert gps.frequency(code) > 1e9

