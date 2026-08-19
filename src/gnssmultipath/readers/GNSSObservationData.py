"""
Accessors for RINEX observation data.

Wraps the raw dict-of-dicts produced by ``readRinexObs`` into typed,
attribute-accessible objects.  Observation codes follow the RINEX
convention, and both variants are supported.

**RINEX 3/4 — three characters** ``<type><band><attribute>``::

    Character 1 – observation type
        C  pseudorange          (code)
        L  carrier phase        (cycles)
        D  Doppler              (Hz)
        S  signal strength      (dB-Hz)

    Character 2 – frequency band
        1  L1 / E1 / B1C / G1   (1575.42 MHz)
        2  L2 / B1I / G2        (1227.60 / 1561.098 / 1246 MHz)
        3  G3                   (1202.025 MHz)
        4  G1a                  (1600.995 MHz)
        5  L5 / E5a / B2a       (1176.45 MHz)
        6  E6 / B3 / G2a        (1278.75 / 1268.52 / 1248.06 MHz)
        7  E5b / B2b            (1207.14 MHz)
        8  E5(a+b) / B2(a+b)    (1191.795 MHz)

    Character 3 – tracking mode / attribute
        C  C/A code, civilian
        P  P code
        W  Z-tracking (codeless)
        X  Pilot + data combined
        I  Data channel (Galileo I)
        Q  Pilot channel (Galileo Q)
        …  (see RINEX 3.05 / 4.02 specification for full list)

**RINEX 2 — two characters** ``<type><band>``, e.g. ``C1``, ``P2``, ``L1``.
The ``P`` type is a P-code pseudorange and is treated as such.


Usage
-----
::

    from gnssmultipath import readRinexObs

    rinex = readRinexObs("observation.rnx")
    obs   = rinex.observations          # GNSSObservationData

    obs.summary()                       # overview table of systems and codes
    obs.systems                         # ['G', 'R', 'E', 'C']
    obs.select(obs_type='L', band=5)    # {'G': ['L5Q'], 'E': ['L5X'], ...}

    # Per-system accessor (property or bracket)
    gps = obs.gps                       # SystemObservations
    gal = obs['E']                      # SystemObservations

    # List observation codes
    gps.codes                           # ['C1C', 'L1C', 'S1C', 'C5X', ...]
    gps.pseudorange_codes               # ['C1C', 'C5X']
    gps.phase_codes                     # ['L1C', 'L5X']
    gps.snr_codes                       # ['S1C', 'S5X']
    gps.doppler_codes                   # ['D1C']
    gps.select(obs_type='C', band=1)    # ['C1C']

    # Signal metadata
    sig = gps.signal('L1C')             # ObsCode
    sig.band, sig.attribute             # 1, 'C'
    sig.frequency(), sig.wavelength()   # 1575420000.0, 0.19029...

    # Extract a specific observation  →  ndarray [epochs, max_sat]
    gps['C1C']                          # raw values, missing = 0.0
    gps.get('C1C')                      # independent copy, missing = NaN

    # Per-satellite view  →  ndarray [epochs]
    g05 = gps.sat(5)
    g05['C1C']
    g05.sv_id                           # 'G05'

    # Per-epoch view  →  ndarray [max_sat], i.e. all signals of one epoch
    ep = gps.epoch(34)                  # 0-based; epoch(-1) is the last one
    ep['C1C']                           # one signal, all satellites
    ep.matrix                           # [max_sat, n_codes] block, no copy
    ep.sat(23)                          # {'C1C': ..., 'L1C': ..., ...}
    ep.to_dataframe()                   # satellites x signals table

    # Carrier frequencies (GLONASS FDMA is channel-aware)
    gps.frequency('L1C')                # 1575420000.0
    gps.wavelength('L1C')               # 0.19029...
    obs.glonass.frequency('L1C', prn=3) # per-satellite FDMA frequency

    # Filter by frequency band
    gps.band(1)                         # dict {'C1C': arr, 'L1C': arr, 'S1C': arr}
    gps.band(5)                         # dict {'C5X': arr, 'L5X': arr, 'S5X': arr}

    # Check if a code exists
    'C1C' in gps                        # True

    # Number of epochs / satellites
    gps.n_epochs                        # int
    gps.n_satellites                    # int (max PRN + 1)
    gps.prns                            # observed PRNs, e.g. [1, 2, 3, ...]

    # Epoch time stamps
    obs.time_epochs                     # [[gps_week, tow], ...]
    obs.datetimes                       # ndarray of datetime64 (GPS time scale)

    # Tidy pandas export
    gps.to_dataframe(codes=['C1C', 'L1C'], prns=[5, 12])

    # LLI and signal-strength sub-accessors (same interface)
    gps.lli['L1C']                      # ndarray [epochs, max_sat]
    gps.ss['S1C']                       # ndarray [epochs, max_sat]
"""

import numpy as np

from gnssmultipath.constants import (
    BAND_DESCRIPTIONS,
    GLONASS_FDMA,
    OBS_TYPE_NAMES,
    PSEUDORANGE_TYPES,
    SPEED_OF_LIGHT,
    SYSTEM_NAMES,
    band_label,
    carrier_frequency,
)
from gnssmultipath.readers.ObsCode import ObsCode


# ── Observation-type helpers ──────────────────────────────────────────────────

_OBS_TYPE_NAMES = OBS_TYPE_NAMES
_BAND_DESCRIPTIONS = BAND_DESCRIPTIONS
_SYSTEM_NAMES = SYSTEM_NAMES

_GPS_EPOCH = np.datetime64('1980-01-06T00:00:00', 'ns')
_SECONDS_PER_WEEK = 604800.0


def _epochs_to_datetime64(time_epochs):
    """Convert a ``[epochs, 2]`` (gps_week, tow) array to ``datetime64[ns]``.

    Nanoseconds are used because the RINEX epoch field holds 7 decimals
    (0.1 us), which a microsecond resolution would round away.
    """
    if time_epochs is None:
        return None
    arr = np.asarray(time_epochs, dtype=float)
    if arr.ndim != 2 or arr.shape[1] < 2 or arr.shape[0] == 0:
        return None
    # The tow is split into whole and fractional seconds before the week offset
    # (~1e9 s) is added, otherwise float64 rounding would swallow the sub-second
    # part of epochs such as 13:22:14.0001055.
    whole_tow = np.floor(arr[:, 1])
    seconds = (arr[:, 0] * _SECONDS_PER_WEEK + whole_tow).astype('int64')
    nano = np.round((arr[:, 1] - whole_tow) * 1e9).astype('int64')
    return (_GPS_EPOCH + seconds.astype('timedelta64[s]')
            + nano.astype('timedelta64[ns]'))


# ── Code-indexed sub-accessor (shared by obs, LLI, SS) ───────────────────────

class _CodeAccessor:
    """Indexes per-epoch matrices by observation code string.

    Columns are extracted lazily and cached, so requesting a single code
    never materialises the full ``[epochs, sats, codes]`` cube.
    """

    __slots__ = ('_epoch_dict', '_codes', '_stacked', '_columns', '_epoch_keys')

    def __init__(self, epoch_dict, codes):
        self._epoch_dict = epoch_dict if isinstance(epoch_dict, dict) else {}
        self._codes = list(codes)
        self._stacked = None
        self._columns = {}
        self._epoch_keys = None

    @property
    def n_epochs(self):
        return len(self._epoch_dict)

    @property
    def n_satellites(self):
        for arr in self._epoch_dict.values():
            return arr.shape[0]
        return 0

    def _check_code(self, code):
        if code not in self._codes:
            raise KeyError(
                f"Observation code '{code}' not available. "
                f"Available codes: {self._codes}"
            )

    def _column(self, code):
        """Return the cached 2-D array ``[epochs, max_sat]`` for *code*."""
        self._check_code(code)
        cached = self._columns.get(code)
        if cached is not None:
            return cached

        idx = self._codes.index(code)
        epochs = list(self._epoch_dict.values())
        if not epochs:
            out = np.empty((0, 0))
        else:
            first = epochs[0]
            out = np.empty((len(epochs), first.shape[0]), dtype=first.dtype)
            for row, arr in enumerate(epochs):
                out[row] = arr[:, idx]

        self._columns[code] = out
        return out

    def _epoch_matrix(self, index):
        """Return the stored matrix ``[max_sat, n_codes]`` for one epoch.

        This is the layout the readers produce, so no copy is made.
        *index* is 0-based and supports negative indexing.
        """
        if self._epoch_keys is None:
            self._epoch_keys = list(self._epoch_dict)
        n_epochs = len(self._epoch_keys)
        if not -n_epochs <= index < n_epochs:
            raise IndexError(
                f"Epoch index {index} is out of range (file has {n_epochs} epochs)."
            )
        return self._epoch_dict[self._epoch_keys[index]]

    def _stack(self):
        """Return the full 3-D array ``[epochs, max_sat, n_codes]``.

        Materialises every observation code at once and is therefore
        memory-heavy for long, high-rate files.  Prefer ``[code]``.
        """
        if self._stacked is None:
            if self._epoch_dict:
                self._stacked = np.stack(list(self._epoch_dict.values()))
            else:
                self._stacked = np.empty((0, 0, len(self._codes)))
        return self._stacked

    def clear_cache(self):
        """Drop all cached arrays to release memory."""
        self._columns.clear()
        self._stacked = None
        self._epoch_keys = None

    def __getitem__(self, code):
        """Return 2-D array [epochs, max_sat] for the given obs code."""
        return self._column(code)

    def __contains__(self, code):
        return code in self._codes

    def __repr__(self):
        return f"_CodeAccessor(codes={self._codes})"


class _SatelliteCodeView:
    """Per-satellite slice of a :class:`_CodeAccessor` (used for LLI / SS)."""

    __slots__ = ('_accessor', '_prn')

    def __init__(self, accessor, prn):
        self._accessor = accessor
        self._prn = prn

    def __getitem__(self, code):
        return self._accessor[code][:, self._prn]

    def __contains__(self, code):
        return code in self._accessor

    def __repr__(self):
        return f"_SatelliteCodeView(prn={self._prn})"


class _EpochCodeView:
    """Single-epoch slice of a :class:`_CodeAccessor` (used for LLI / SS)."""

    __slots__ = ('_accessor', '_index')

    def __init__(self, accessor, index):
        self._accessor = accessor
        self._index = index

    def __getitem__(self, code):
        self._accessor._check_code(code)
        col = self._accessor._codes.index(code)
        return self._accessor._epoch_matrix(self._index)[:, col]

    def __contains__(self, code):
        return code in self._accessor

    def __repr__(self):
        return f"_EpochCodeView(epoch_index={self._index})"


# ── Per-epoch observations ───────────────────────────────────────────────────

class EpochObservations:
    """All signals of one epoch for a single GNSS system.

    Returned by :meth:`SystemObservations.epoch`.  Where
    ``sys_obs['C1C']`` slices one signal across all epochs, this slices one
    epoch across all signals.  Arrays returned by ``[code]`` are 1-D over
    PRN, and :attr:`matrix` is the ``[max_sat, n_codes]`` block exactly as
    stored by the reader, so no copying takes place.
    """

    __slots__ = ('_system', '_index')

    def __init__(self, system, index):
        self._system = system
        n_epochs = system.n_epochs
        if not -n_epochs <= index < n_epochs:
            raise IndexError(
                f"Epoch index {index} is out of range "
                f"({system.system_name} has {n_epochs} epochs)."
            )
        self._index = index if index >= 0 else n_epochs + index

    # ── Identity ─────────────────────────────────────────────────────────

    @property
    def index(self):
        """0-based epoch index, matching the row index of ``sys_obs[code]``."""
        return self._index

    @property
    def number(self):
        """1-based epoch number, as used by the raw ``GNSS_obs`` dicts."""
        return self._index + 1

    @property
    def system_code(self):
        """Single-character GNSS system code."""
        return self._system.system_code

    @property
    def system_name(self):
        """Full GNSS system name."""
        return self._system.system_name

    @property
    def codes(self):
        """All available observation codes for this system."""
        return self._system.codes

    @property
    def time(self):
        """``(gps_week, time_of_week)`` for this epoch, or None."""
        times = self._system.epoch_times
        if times is None:
            return None
        return tuple(np.asarray(times)[self._index])

    @property
    def datetime(self):
        """Epoch time stamp as ``datetime64[ns]`` (GPS time scale), or None."""
        stamps = self._system.datetimes
        return None if stamps is None else stamps[self._index]

    # ── Data access ──────────────────────────────────────────────────────

    @property
    def matrix(self):
        """The stored ``[max_sat, n_codes]`` block for this epoch (no copy)."""
        return self._system._obs._epoch_matrix(self._index)

    def __getitem__(self, code):
        """``epoch['C1C']`` → ndarray [max_sat] of raw values (missing = 0.0)."""
        self._system._obs._check_code(code)
        col = self._system.codes.index(code)
        return self.matrix[:, col]

    def __contains__(self, code):
        return code in self._system

    def get(self, code, mask_missing=True):
        """Return a copy of this epoch's values for *code*, optionally NaN-masked."""
        values = self[code].astype(np.float64, copy=True)
        if mask_missing:
            values[values == 0] = np.nan
        return values

    def sat(self, prn):
        """Return all signal values for one satellite as ``{code: value}``."""
        row = self.matrix[int(prn)]
        return {code: float(row[i]) for i, code in enumerate(self.codes)}

    @property
    def prns(self):
        """PRNs with at least one observation in this epoch."""
        observed = np.any(self.matrix != 0, axis=1)
        return [int(p) for p in np.flatnonzero(observed) if p > 0]

    # ── Sub-accessors ────────────────────────────────────────────────────

    @property
    def lli(self):
        """Loss-of-Lock Indicator accessor (``epoch.lli['L1C']``)."""
        return _EpochCodeView(self._system.lli, self._index)

    @property
    def ss(self):
        """Signal-Strength accessor (``epoch.ss['S1C']``)."""
        return _EpochCodeView(self._system.ss, self._index)

    # ── Export ───────────────────────────────────────────────────────────

    def to_dataframe(self, codes=None, prns=None, mask_missing=True, dropna=True):
        """Return this epoch as a table of satellites (rows) by signals (columns).

        Parameters
        ----------
        codes : list[str], optional
            Observation codes to include.  Defaults to all codes.
        prns : list[int], optional
            Satellite PRNs to include.  Defaults to the satellites observed
            in this epoch.
        mask_missing : bool
            Replace the ``0.0`` fill value with ``NaN``.
        dropna : bool
            Drop satellites without any observation in this epoch.
        """
        import pandas as pd

        codes = list(codes) if codes is not None else self.codes
        prn_list = np.asarray(prns if prns is not None else self.prns, dtype=int)

        if prn_list.size == 0 or not codes:
            return pd.DataFrame(columns=codes, index=pd.Index([], name='sv'))

        cols = [self._system.codes.index(c) for c in codes]
        values = self.matrix[np.ix_(prn_list, cols)].astype(np.float64, copy=True)
        if mask_missing:
            values[values == 0] = np.nan

        frame = pd.DataFrame(
            values, columns=codes,
            index=pd.Index([f"{self.system_code}{p:02d}" for p in prn_list], name='sv'),
        )
        frame.insert(0, 'prn', prn_list)
        if dropna:
            frame = frame.dropna(subset=codes, how='all')
        return frame

    def __repr__(self):
        stamp = self.datetime
        when = '' if stamp is None else f", {stamp}"
        return (f"EpochObservations({self.system_code}, epoch {self.number}"
                f"{when}, {len(self.prns)} satellites)")


# ── Per-satellite observations ───────────────────────────────────────────────

class SatelliteObservations:
    """Observations of a single satellite, indexed by observation code.

    Returned by :meth:`SystemObservations.sat`.  All arrays are 1-D over
    epochs.
    """

    __slots__ = ('_system', '_prn')

    def __init__(self, system, prn):
        self._system = system
        self._prn = prn

    # ── Identity ─────────────────────────────────────────────────────────

    @property
    def prn(self):
        """Satellite PRN / slot number."""
        return self._prn

    @property
    def sv_id(self):
        """RINEX satellite identifier, e.g. ``'G05'``."""
        return f"{self._system.system_code}{self._prn:02d}"

    @property
    def system_code(self):
        """Single-character GNSS system code."""
        return self._system.system_code

    @property
    def system_name(self):
        """Full GNSS system name."""
        return self._system.system_name

    @property
    def codes(self):
        """All available observation codes for this satellite's system."""
        return self._system.codes

    # ── Data access ──────────────────────────────────────────────────────

    def __getitem__(self, code):
        """``sat['C1C']`` → ndarray [epochs] of raw values (missing = 0.0)."""
        return self._system[code][:, self._prn]

    def __contains__(self, code):
        return code in self._system

    def get(self, code, mask_missing=True):
        """Return a copy of the observations, optionally NaN-masked."""
        return self._system.get(code, mask_missing=mask_missing)[:, self._prn]

    def frequency(self, code):
        """Carrier frequency in Hz for *code* on this satellite."""
        return self._system.frequency(code, prn=self._prn)

    def wavelength(self, code):
        """Carrier wavelength in metres for *code* on this satellite."""
        return self._system.wavelength(code, prn=self._prn)

    @property
    def lli(self):
        """Loss-of-Lock Indicator accessor (``sat.lli['L1C']``)."""
        return _SatelliteCodeView(self._system.lli, self._prn)

    @property
    def ss(self):
        """Signal-Strength accessor (``sat.ss['S1C']``)."""
        return _SatelliteCodeView(self._system.ss, self._prn)

    def to_dataframe(self, codes=None, mask_missing=True, dropna=True):
        """Return a tidy ``DataFrame`` for this satellite."""
        return self._system.to_dataframe(
            codes=codes, prns=[self._prn],
            mask_missing=mask_missing, dropna=dropna,
        )

    def __repr__(self):
        return f"SatelliteObservations({self.sv_id}, codes={self.codes})"


# ── Per-system observations ──────────────────────────────────────────────────

class SystemObservations:
    """Observations for a single GNSS system, indexed by obs code.

    Parameters
    ----------
    epoch_dict : dict[int, ndarray]
        Per-epoch observation matrices ``{epoch: array[max_sat, n_codes]}``.
    obs_codes : list[str]
        Ordered observation codes matching the column axis.
    lli_epoch_dict : dict[int, ndarray] or None
        Per-epoch Loss-of-Lock Indicator matrices (same shape).
    ss_epoch_dict : dict[int, ndarray] or None
        Per-epoch Signal-Strength matrices (same shape).
    system_code : str
        Single-character GNSS system identifier ('G', 'R', 'E', 'C').
    svs : ndarray, optional
        ``GNSS_SVs`` matrix ``[epochs, max_sat]`` where column 0 holds the
        satellite count and the remaining columns the observed PRNs.
    time_epochs : ndarray, optional
        ``[epochs, 2]`` array of (GPS week, time-of-week).
    glo_slot2channel : dict, optional
        ``{glonass_slot: channel_number}`` needed for GLONASS frequencies.
    """

    def __init__(self, epoch_dict, obs_codes, lli_epoch_dict=None,
                 ss_epoch_dict=None, system_code='', *,
                 svs=None, time_epochs=None, glo_slot2channel=None):
        self._obs = _CodeAccessor(epoch_dict, obs_codes)
        self._obs_codes = list(obs_codes)
        self._system_code = system_code
        self._svs = svs
        self._time_epochs = time_epochs
        self._glo_slot2channel = glo_slot2channel if isinstance(glo_slot2channel, dict) else {}
        self._signal_cache = {}
        self._prns = None

        self._lli = (_CodeAccessor(lli_epoch_dict, obs_codes)
                     if isinstance(lli_epoch_dict, dict) else None)
        self._ss = (_CodeAccessor(ss_epoch_dict, obs_codes)
                    if isinstance(ss_epoch_dict, dict) else None)

    # ── Code listing ─────────────────────────────────────────────────────

    @property
    def codes(self):
        """All available observation codes for this system."""
        return list(self._obs_codes)

    @property
    def pseudorange_codes(self):
        """Codes starting with 'C', or RINEX 2 'P' (pseudorange observations)."""
        return [c for c in self._obs_codes if c[0] in PSEUDORANGE_TYPES]

    @property
    def phase_codes(self):
        """Codes starting with 'L' (carrier-phase observations)."""
        return [c for c in self._obs_codes if c[0] == 'L']

    @property
    def doppler_codes(self):
        """Codes starting with 'D' (Doppler observations)."""
        return [c for c in self._obs_codes if c[0] == 'D']

    @property
    def snr_codes(self):
        """Codes starting with 'S' (signal-strength / SNR observations)."""
        return [c for c in self._obs_codes if c[0] == 'S']

    @property
    def bands(self):
        """Set of frequency-band digits present in the observation codes."""
        return sorted({c[1] for c in self._obs_codes})

    def select(self, obs_type=None, band=None, attribute=None):
        """Return the observation codes matching all given criteria.

        Parameters
        ----------
        obs_type : str, optional
            ``'C'`` (pseudorange, also matches RINEX 2 ``'P'``), ``'L'``
            (phase), ``'D'`` (Doppler) or ``'S'`` (SNR).
        band : int or str, optional
            Frequency-band digit, e.g. ``1``, ``2`` or ``5``.
        attribute : str, optional
            Tracking attribute, e.g. ``'C'``, ``'W'``, ``'X'``, ``'Q'``.
            Use ``''`` to select RINEX 2 codes, which carry no attribute.

        Examples
        --------
        ::

            gps.select(obs_type='L', band=1)     # ['L1C', 'L1W']
            gps.select(band=5)                   # ['C5Q', 'L5Q', 'S5Q']
        """
        result = list(self._obs_codes)

        if obs_type is not None:
            wanted = obs_type.upper()
            types = PSEUDORANGE_TYPES if wanted in PSEUDORANGE_TYPES else {wanted}
            result = [c for c in result if c[0] in types]

        if band is not None:
            digit = str(band)
            result = [c for c in result if c[1] == digit]

        if attribute is not None:
            wanted_attr = attribute.upper()
            result = [c for c in result if (c[2] if len(c) > 2 else '') == wanted_attr]

        return result

    # ── Signal metadata ──────────────────────────────────────────────────

    def signal(self, code):
        """Return the parsed :class:`ObsCode` for *code*."""
        self._obs._check_code(code)
        cached = self._signal_cache.get(code)
        if cached is None:
            cached = ObsCode.parse(code, self._system_code)
            self._signal_cache[code] = cached
        return cached

    @property
    def signals(self):
        """All observation codes as parsed :class:`ObsCode` objects."""
        return [self.signal(c) for c in self._obs_codes]

    # ── Carrier frequencies ──────────────────────────────────────────────

    def glonass_channel(self, prn):
        """Return the GLONASS FDMA channel number (k) for *prn*, or None."""
        return self._glo_slot2channel.get(int(prn))

    def frequency(self, code, prn=None):
        """Carrier frequency in Hz for *code*.

        For the GLONASS FDMA bands (1 and 2) the frequency is satellite
        specific: pass *prn* to get a scalar, or omit it to get an array
        indexed by PRN, holding ``NaN`` where the channel number is unknown.
        """
        signal = self.signal(code)

        if self._system_code == 'R' and signal.band in GLONASS_FDMA:
            if prn is not None:
                channel = self.glonass_channel(prn)
                if channel is None:
                    raise KeyError(
                        f"No GLONASS channel number known for PRN {prn}. The channel "
                        "numbers come from the 'GLONASS SLOT / FRQ #' RINEX header record."
                    )
                return carrier_frequency('R', signal.band, channel)

            freqs = np.full(self.n_satellites, np.nan)
            for slot, channel in self._glo_slot2channel.items():
                slot = int(slot)
                if 0 <= slot < freqs.size:
                    freqs[slot] = carrier_frequency('R', signal.band, channel)
            return freqs

        return signal.frequency()

    def wavelength(self, code, prn=None):
        """Carrier wavelength in metres for *code*. See :meth:`frequency`."""
        return SPEED_OF_LIGHT / self.frequency(code, prn=prn)

    # ── Shape info ───────────────────────────────────────────────────────

    @property
    def n_epochs(self):
        """Number of observation epochs."""
        return self._obs.n_epochs

    @property
    def n_satellites(self):
        """Size of the satellite (PRN) axis (includes unused row 0)."""
        return self._obs.n_satellites

    @property
    def system_code(self):
        """Single-character GNSS system code ('G', 'R', 'E', 'C')."""
        return self._system_code

    @property
    def system_name(self):
        """Full GNSS system name."""
        return SYSTEM_NAMES.get(self._system_code, self._system_code)

    @property
    def prns(self):
        """Sorted PRNs that have observations in this file."""
        if self._prns is None:
            self._prns = self._find_prns()
        return list(self._prns)

    def _find_prns(self):
        svs = self._svs
        if isinstance(svs, np.ndarray) and svs.ndim == 2 and svs.shape[1] > 1:
            return [int(p) for p in np.unique(svs[:, 1:]) if p > 0]

        n_sat = self.n_satellites
        if n_sat == 0:
            return []
        observed = np.zeros(n_sat, dtype=bool)
        for arr in self._obs._epoch_dict.values():
            observed |= np.any(arr != 0, axis=1)
        return [int(p) for p in np.flatnonzero(observed) if p > 0]

    # ── Epoch times ──────────────────────────────────────────────────────

    @property
    def epoch_times(self):
        """``[epochs, 2]`` array of (GPS week, time-of-week), or None."""
        return self._time_epochs

    @property
    def datetimes(self):
        """Epoch time stamps as ``datetime64[ns]`` in the GPS time scale."""
        return _epochs_to_datetime64(self._time_epochs)

    # ── Data access ──────────────────────────────────────────────────────

    def __getitem__(self, code):
        """``sys_obs['C1C']`` → ndarray [epochs, max_sat] of raw values.

        Missing observations are ``0.0``.  The returned array is cached and
        shared — use :meth:`get` for an independent, NaN-masked copy.
        """
        return self._obs[code]

    def __contains__(self, code):
        return code in self._obs

    def get(self, code, mask_missing=True):
        """Return an independent copy of the observations for *code*.

        Parameters
        ----------
        code : str
            Observation code, e.g. ``'C1C'``.
        mask_missing : bool
            When True (default) missing observations, stored as ``0.0`` by
            the RINEX readers, are replaced by ``NaN``.
        """
        values = self._obs[code].astype(np.float64, copy=True)
        if mask_missing:
            values[values == 0] = np.nan
        return values

    def sat(self, prn):
        """Return a :class:`SatelliteObservations` view for *prn*."""
        prn = int(prn)
        n_sat = self.n_satellites
        if not 0 <= prn < n_sat:
            raise KeyError(
                f"PRN {prn} is out of range for {self.system_name} "
                f"(valid range 1-{max(n_sat - 1, 0)})."
            )
        return SatelliteObservations(self, prn)

    def epoch(self, index):
        """Return an :class:`EpochObservations` view for a single epoch.

        *index* is 0-based and matches the row index of ``sys_obs[code]``;
        negative indices count from the end (``epoch(-1)`` is the last
        epoch).  Use :attr:`EpochObservations.number` for the 1-based epoch
        number used by the raw ``GNSS_obs`` dicts.
        """
        return EpochObservations(self, int(index))

    @property
    def data(self):
        """Full 3-D observation array [epochs, max_sat, n_codes].

        Materialises every code at once; for long, high-rate files prefer
        ``sys_obs['C1C']`` or :meth:`get`, which build one code at a time.
        """
        return self._obs._stack()

    def band(self, band_num):
        """Return a dict of ``{code: array}`` for a specific frequency band.

        Parameters
        ----------
        band_num : int or str
            Frequency-band digit (1, 2, 5, 6, 7, 8).
        """
        return {c: self._obs[c] for c in self.select(band=band_num)}

    def by_type(self, obs_type):
        """Return a dict of ``{code: array}`` for a specific observation type.

        Parameters
        ----------
        obs_type : str
            Single character: 'C' (pseudorange), 'L' (phase),
            'D' (Doppler), or 'S' (SNR).
        """
        return {c: self._obs[c] for c in self.select(obs_type=obs_type)}

    # ── LLI / SS sub-accessors ───────────────────────────────────────────

    @property
    def lli(self):
        """Loss-of-Lock Indicator accessor (same ``[code]`` interface)."""
        if self._lli is None:
            raise AttributeError(
                "LLI data not available (readLLI was disabled or data missing)."
            )
        return self._lli

    @property
    def ss(self):
        """Signal-Strength accessor (same ``[code]`` interface)."""
        if self._ss is None:
            raise AttributeError(
                "Signal-Strength data not available (readSS was disabled or data missing)."
            )
        return self._ss

    # ── Export ──────────────────────────────────────────────────────────

    def to_dataframe(self, codes=None, prns=None, mask_missing=True, dropna=True,
                     include_lli=False, include_ss=False):
        """Return the observations in tidy (long) form as a ``DataFrame``.

        Columns: ``epoch``, ``datetime`` (when epoch times are available),
        ``sv``, ``prn``, ``code``, ``value``, plus ``lli`` / ``ss`` when
        requested.

        Parameters
        ----------
        codes : list[str], optional
            Observation codes to include.  Defaults to all codes.
        prns : list[int], optional
            Satellite PRNs to include.  Defaults to the observed satellites.
        mask_missing : bool
            Replace the ``0.0`` fill value with ``NaN``.
        dropna : bool
            Drop rows without an observation.  Recommended, as the full cross
            product of epochs, satellites and codes is large.
        include_lli, include_ss : bool
            Add the Loss-of-Lock / Signal-Strength indicator columns.
        """
        import pandas as pd

        codes = list(codes) if codes is not None else self.codes
        prn_list = np.asarray(prns if prns is not None else self.prns, dtype=int)
        n_epochs = self.n_epochs

        if n_epochs == 0 or prn_list.size == 0 or not codes:
            return pd.DataFrame(columns=['epoch', 'sv', 'prn', 'code', 'value'])

        epoch_idx = np.repeat(np.arange(1, n_epochs + 1), prn_list.size)
        prn_idx = np.tile(prn_list, n_epochs)
        sv_ids = [f"{self._system_code}{p:02d}" for p in prn_idx]
        datetimes = self.datetimes

        frames = []
        for code in codes:
            block = {
                'epoch': epoch_idx,
                'sv': sv_ids,
                'prn': prn_idx,
                'code': code,
                'value': self.get(code, mask_missing=mask_missing)[:, prn_list].ravel(),
            }
            if datetimes is not None:
                block['datetime'] = np.repeat(datetimes, prn_list.size)
            if include_lli and self._lli is not None:
                block['lli'] = self._lli[code][:, prn_list].ravel()
            if include_ss and self._ss is not None:
                block['ss'] = self._ss[code][:, prn_list].ravel()
            frames.append(pd.DataFrame(block))

        frame = pd.concat(frames, ignore_index=True)
        if dropna:
            frame = frame.dropna(subset=['value']).reset_index(drop=True)

        ordered = ['epoch', 'datetime', 'sv', 'prn', 'code', 'value', 'lli', 'ss']
        return frame[[c for c in ordered if c in frame.columns]]

    # ── Representation ────────────────────────────────────────────────

    def summary(self):
        """Return a human-readable overview of the available signals."""
        lines = [f"{self.system_name} ({self._system_code}) — "
                 f"{self.n_epochs} epochs, {len(self.prns)} satellites"]
        for band_digit in self.bands:
            lines.append(f"  Band {band_digit}: {band_label(self._system_code, band_digit)}")
            lines.append(f"    {', '.join(self.select(band=band_digit))}")
        return "\n".join(lines)

    def __repr__(self):
        return (f"SystemObservations({self.system_name}, "
                f"codes={self._obs_codes}, "
                f"epochs={self.n_epochs}, sats={self.n_satellites})")


# ── Top-level store ──────────────────────────────────────────────────────────

class GNSSObservationData:
    """Top-level observation container, keyed by GNSS system code.

    Constructed automatically from the raw dicts returned by
    ``readRinexObs``.  Access per-system data via properties
    (``.gps``, ``.glonass``, ``.galileo``, ``.beidou``) or bracket
    notation (``store['G']``).
    """

    def __init__(self, gnss_obs, obs_codes, gnss_systems, gnss_lli=None, gnss_ss=None,
                 *, time_epochs=None, gnss_svs=None, glo_slot2channel=None,
                 t_interval=None, marker_name='', approx_position=None):
        self._systems = {}
        self._time_epochs = time_epochs
        self._t_interval = t_interval
        self._marker_name = marker_name
        self._approx_position = approx_position

        if not isinstance(glo_slot2channel, dict):
            glo_slot2channel = {}

        # Flatten obsCodes: {1: {'G': [...]}, 2: {'R': [...]}} → {'G': [...], 'R': [...]}
        flat_codes = {}
        for sys_idx, sys_code in gnss_systems.items():
            if sys_idx in obs_codes and sys_code in obs_codes[sys_idx]:
                flat_codes[sys_code] = obs_codes[sys_idx][sys_code]

        for sys_idx, sys_code in gnss_systems.items():
            epoch_dict = gnss_obs.get(sys_code, {})
            codes = flat_codes.get(sys_code, [])

            # Skip systems with no data
            if not epoch_dict or not isinstance(epoch_dict, dict):
                continue

            lli_dict = None
            if isinstance(gnss_lli, dict):
                lli_val = gnss_lli.get(sys_code)
                if isinstance(lli_val, dict) and lli_val:
                    lli_dict = lli_val

            ss_dict = None
            if isinstance(gnss_ss, dict):
                ss_val = gnss_ss.get(sys_code)
                if isinstance(ss_val, dict) and ss_val:
                    ss_dict = ss_val

            svs = gnss_svs.get(sys_code) if isinstance(gnss_svs, dict) else None

            self._systems[sys_code] = SystemObservations(
                epoch_dict, codes, lli_dict, ss_dict, sys_code,
                svs=svs, time_epochs=time_epochs,
                glo_slot2channel=glo_slot2channel if sys_code == 'R' else None,
            )

    # ── System listing ───────────────────────────────────────────────────

    @property
    def systems(self):
        """List of available system codes (e.g. ``['G', 'E', 'C']``)."""
        return list(self._systems.keys())

    def select(self, system=None, obs_type=None, band=None, attribute=None):
        """Return ``{system_code: [codes]}`` matching the given criteria.

        Systems without a matching code are omitted.

        Examples
        --------
        ::

            obs.select(obs_type='L', band=5)   # {'G': ['L5Q'], 'E': ['L5X']}
            obs.select(system='G', band=1)     # {'G': ['C1C', 'L1C', 'S1C']}
        """
        if system is None:
            wanted = self.systems
        else:
            wanted = [system] if isinstance(system, str) else list(system)

        result = {}
        for sys_code in wanted:
            codes = self[sys_code].select(
                obs_type=obs_type, band=band, attribute=attribute)
            if codes:
                result[sys_code] = codes
        return result

    # ── Metadata ─────────────────────────────────────────────────────

    @property
    def time_epochs(self):
        """``[epochs, 2]`` array of (GPS week, time-of-week), or None."""
        return self._time_epochs

    @property
    def datetimes(self):
        """Epoch time stamps as ``datetime64[ns]`` in the GPS time scale."""
        return _epochs_to_datetime64(self._time_epochs)

    @property
    def n_epochs(self):
        """Number of observation epochs."""
        if self._time_epochs is not None:
            return int(np.asarray(self._time_epochs).shape[0])
        for sys_obs in self._systems.values():
            return sys_obs.n_epochs
        return 0

    @property
    def interval(self):
        """Observation interval in seconds, or None if unknown."""
        return self._t_interval

    @property
    def marker_name(self):
        """Station / marker name from the RINEX header."""
        return self._marker_name

    @property
    def approx_position(self):
        """Approximate receiver position [X, Y, Z] from the RINEX header."""
        return self._approx_position

    # ── Bracket access ───────────────────────────────────────────────

    def __getitem__(self, sys_code):
        """``store['G']`` → SystemObservations for GPS."""
        if sys_code not in self._systems:
            available = list(self._systems.keys())
            names = [SYSTEM_NAMES.get(s, s) for s in available]
            raise KeyError(
                f"System '{sys_code}' not available. "
                f"Available: {dict(zip(available, names))}"
            )
        return self._systems[sys_code]

    def __contains__(self, sys_code):
        return sys_code in self._systems

    def __iter__(self):
        return iter(self._systems.values())

    def __len__(self):
        return len(self._systems)

    # ── Named properties ─────────────────────────────────────────────────

    @property
    def gps(self):
        """GPS observations."""
        return self['G']

    @property
    def glonass(self):
        """GLONASS observations."""
        return self['R']

    @property
    def galileo(self):
        """Galileo observations."""
        return self['E']

    @property
    def beidou(self):
        """BeiDou observations."""
        return self['C']

    # ── Export ──────────────────────────────────────────────────────────

    def to_dataframe(self, systems=None, codes=None, prns=None,
                     mask_missing=True, dropna=True,
                     include_lli=False, include_ss=False):
        """Return several systems as one tidy ``DataFrame``.

        Adds a ``system`` column to the per-system frame produced by
        :meth:`SystemObservations.to_dataframe`.
        """
        import pandas as pd

        if systems is None:
            wanted = self.systems
        else:
            wanted = [systems] if isinstance(systems, str) else list(systems)

        frames = []
        for sys_code in wanted:
            sys_obs = self[sys_code]
            sys_codes = [c for c in codes if c in sys_obs] if codes is not None else None
            if codes is not None and not sys_codes:
                continue
            frame = sys_obs.to_dataframe(
                codes=sys_codes, prns=prns, mask_missing=mask_missing,
                dropna=dropna, include_lli=include_lli, include_ss=include_ss,
            )
            frame.insert(0, 'system', sys_code)
            frames.append(frame)

        if not frames:
            return pd.DataFrame(columns=['system', 'epoch', 'sv', 'prn', 'code', 'value'])
        return pd.concat(frames, ignore_index=True)

    # ── Representation ────────────────────────────────────────────────

    def summary(self):
        """Return a human-readable overview of all systems and signals."""
        header = f"GNSS observation data — {self.n_epochs} epochs"
        if self._marker_name:
            header += f", marker '{self._marker_name}'"
        if self._t_interval:
            header += f", interval {self._t_interval} s"
        parts = [header, ""]
        parts.extend(sys_obs.summary() for sys_obs in self._systems.values())
        return "\n".join(parts)

    def __repr__(self):
        parts = [repr(s) for s in self._systems.values()]
        return "GNSSObservationData([\n  " + ",\n  ".join(parts) + "\n])"
