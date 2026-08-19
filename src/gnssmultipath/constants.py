"""
Physical and GNSS signal constants shared across the package.

This module is the single source of truth for carrier frequencies, wavelengths,
the speed of light and the per-constellation Earth model parameters.
Frequencies follow the RINEX 3/4 band-digit convention, i.e. the second
character of an observation code such as ``C1C`` (band 1), ``L5Q`` (band 5) or
``C2I`` (band 2).

References
----------
Teunissen, P.J.G. and Montenbruck, O. (eds.), *Springer Handbook of Global
Navigation Satellite Systems*, Springer, 2017.  Table 3.4 "Physical parameters
of GNSS almanac and ephemeris models", p. 80 — source of :data:`GM` and
:data:`EARTH_ROTATION_RATE`.

Usage
-----
::

    from gnssmultipath.constants import carrier_frequency, wavelength

    carrier_frequency('G', 1)              # 1575420000.0  (GPS L1)
    carrier_frequency('E', 5)              # 1176450000.0  (Galileo E5a)
    carrier_frequency('R', 1, -4)          # GLONASS G1, FDMA channel k = -4
    wavelength('G', 2)                     # 0.2442102134...

    from gnssmultipath.constants import earth_gravitational_constant, earth_rotation_rate

    earth_gravitational_constant('C')      # 3.986004418e14  (BeiDou SIS ICD)
    earth_rotation_rate('C')               # 7.292115e-05

Made by: Per Helge Aarnes
E-mail: per.helge.aarnes@gmail.com
"""

from typing import Dict, Mapping, Optional, Tuple

import numpy as np


# ── Physical constants ────────────────────────────────────────────────────────

SPEED_OF_LIGHT = 299792458.0
"""Speed of light in vacuum [m/s] (IAU / CODATA exact value)."""


# ── GNSS system identifiers ───────────────────────────────────────────────────

SYSTEM_NAMES: Dict[str, str] = {
    'G': 'GPS',
    'R': 'GLONASS',
    'E': 'Galileo',
    'C': 'BeiDou',
}

SYSTEM_CODES: Dict[str, str] = {name: code for code, name in SYSTEM_NAMES.items()}


# ── Observation type identifiers (first character of an obs code) ─────────────

OBS_TYPE_NAMES: Dict[str, str] = {
    'C': 'pseudorange',
    'P': 'pseudorange',   # RINEX 2 P-code pseudorange
    'L': 'phase',
    'D': 'doppler',
    'S': 'snr',
    'I': 'iono_phase',
    'X': 'channel',
}

PSEUDORANGE_TYPES = frozenset({'C', 'P'})


# ── Frequency band identifiers (second character of an obs code) ──────────────

BAND_DESCRIPTIONS: Dict[str, str] = {
    '1': 'L1/E1/B1C/G1 (1575.42 MHz)',
    '2': 'L2/B1I/G2 (1227.60 / 1561.098 / 1246 MHz)',
    '3': 'G3 (1202.025 MHz)',
    '4': 'G1a (1600.995 MHz)',
    '5': 'L5/E5a/B2a (1176.45 MHz)',
    '6': 'E6/B3/G2a (1278.75 / 1268.52 / 1248.06 MHz)',
    '7': 'E5b/B2b (1207.14 MHz)',
    '8': 'E5(a+b)/B2(a+b) (1191.795 MHz)',
    '9': 'S-band (IRNSS)',
}

SYSTEM_BANDS: Dict[str, Tuple[int, ...]] = {
    'G': (1, 2, 5),
    'R': (1, 2, 3, 4, 6),
    'E': (1, 5, 6, 7, 8),
    'C': (1, 2, 5, 6, 7, 8),
}


# ── Carrier frequencies [Hz] ──────────────────────────────────────────────────

CARRIER_FREQUENCIES: Dict[str, Dict[int, float]] = {
    'G': {
        1: 1.575420e9,   # L1
        2: 1.227600e9,   # L2
        5: 1.176450e9,   # L5
    },
    'R': {
        # Bands 1 and 2 are FDMA and depend on the satellite channel number;
        # the values below are the k = 0 centre frequencies (see GLONASS_FDMA).
        1: 1.602000e9,   # G1
        2: 1.246000e9,   # G2
        3: 1.202025e9,   # G3  (CDMA)
        4: 1.600995e9,   # G1a (CDMA)
        6: 1.248060e9,   # G2a (CDMA)
    },
    'E': {
        1: 1.575420e9,   # E1
        5: 1.176450e9,   # E5a
        6: 1.278750e9,   # E6
        7: 1.207140e9,   # E5b
        8: 1.191795e9,   # E5(a+b)
    },
    'C': {
        1: 1.575420e9,   # B1C
        2: 1.561098e9,   # B1I
        5: 1.176450e9,   # B2a
        6: 1.268520e9,   # B3
        7: 1.207140e9,   # B2b
        8: 1.191795e9,   # B2(a+b)
    },
}

GLONASS_FDMA: Dict[int, Tuple[float, float]] = {
    1: (1.602000e9, 562500.0),   # G1: 1602 MHz + k * 562.5 kHz
    2: (1.246000e9, 437500.0),   # G2: 1246 MHz + k * 437.5 kHz
}
"""GLONASS FDMA bands as ``{band: (centre_frequency_hz, channel_step_hz)}``."""

BAND_SIGNAL_NAMES: Dict[str, Dict[int, str]] = {
    'G': {1: 'L1', 2: 'L2', 5: 'L5'},
    'R': {1: 'G1', 2: 'G2', 3: 'G3', 4: 'G1a', 6: 'G2a'},
    'E': {1: 'E1', 5: 'E5a', 6: 'E6', 7: 'E5b', 8: 'E5(a+b)'},
    'C': {1: 'B1C', 2: 'B1I', 5: 'B2a', 6: 'B3', 7: 'B2b', 8: 'B2(a+b)'},
}
"""Conventional signal name for each system and RINEX band digit."""

MAX_GLONASS_SLOT = 36
"""Highest GLONASS slot (PRN) number handled by the software."""


# ── Earth model parameters ───────────────────────────────────────────────

# Source: Teunissen, P.J.G. and Montenbruck, O. (eds.), "Springer Handbook of
# Global Navigation Satellite Systems", Springer, 2017.
# Table 3.4 "Physical parameters of GNSS almanac and ephemeris models", p. 80.
#
# Each constellation's Interface Control Document defines its own values, and a
# receiver must reconstruct the orbit with the same constants the control
# segment used to fit the broadcast elements.  Mixing them up shifts the
# satellite along-track: using the GPS rotation rate for BeiDou is worth up to
# ~25 m (MEO) at the end of a BDT week, because of the -omega*toe term.

GM: Dict[str, float] = {
    'G': 398600.5e9,      # IS-GPS-200
    'R': 398600.4418e9,   # GLONASS ICD (PZ-90)
    'E': 398600.4418e9,   # Galileo OS SIS ICD
    'C': 398600.4418e9,   # BeiDou SIS ICD
    'J': 398600.5e9,      # QZSS IS-QZSS (follows GPS)
}
"""Earth gravitational constant ``GM`` [m^3/s^2] per GNSS."""

EARTH_ROTATION_RATE: Dict[str, float] = {
    'G': 7.2921151467e-5,
    'R': 7.292115e-5,
    'E': 7.2921151467e-5,
    'C': 7.2921150e-5,
    'J': 7.2921151467e-5,
}
"""Earth rotation rate ``omega`` [rad/s] per GNSS."""

J2 = 1.0826257e-3
"""Second zonal harmonic coefficient, as used by the GLONASS equations of motion."""

PZ90_SEMI_MAJOR_AXIS = 6378136.0
"""Semi-major axis of the PZ-90 ellipsoid [m]."""

WGS84_SEMI_MAJOR_AXIS = 6378137.0
"""Semi-major axis of the WGS 84 ellipsoid [m]."""

WGS84_SEMI_MINOR_AXIS = 6356752.314245
"""Semi-minor axis of the WGS 84 ellipsoid [m]."""


# ── Public helpers ────────────────────────────────────────────────────────

def earth_gravitational_constant(system: str) -> float:
    """Return ``GM`` [m^3/s^2] for *system*.

    Unknown systems fall back to the GPS value.
    """
    return GM.get(_normalise_system(system), GM['G'])


def earth_rotation_rate(system: str) -> float:
    """Return the Earth rotation rate [rad/s] for *system*.

    Unknown systems fall back to the GPS value.
    """
    return EARTH_ROTATION_RATE.get(_normalise_system(system), EARTH_ROTATION_RATE['G'])


def carrier_frequency(system: str, band: int, glonass_channel: Optional[int] = None) -> float:
    """Return the carrier frequency in Hz for *system* and *band*.

    Parameters
    ----------
    system : str
        Single-character GNSS system code (``'G'``, ``'R'``, ``'E'``, ``'C'``).
    band : int or str
        RINEX frequency band digit, e.g. ``1``, ``2`` or ``5``.
    glonass_channel : int, optional
        GLONASS FDMA channel number (k), required for GLONASS bands 1 and 2.

    Raises
    ------
    KeyError
        If *system* is unknown or *band* is not defined for that system.
    ValueError
        If a GLONASS FDMA band is requested without a channel number.
    """
    sys_code = _normalise_system(system)
    band_num = int(band)

    try:
        bands = CARRIER_FREQUENCIES[sys_code]
    except KeyError:
        raise KeyError(
            f"Unknown GNSS system '{system}'. "
            f"Supported systems: {sorted(CARRIER_FREQUENCIES)}"
        ) from None

    if band_num not in bands:
        raise KeyError(
            f"Band {band_num} is not defined for {SYSTEM_NAMES[sys_code]}. "
            f"Available bands: {sorted(bands)}"
        )

    if sys_code == 'R' and band_num in GLONASS_FDMA:
        if glonass_channel is None:
            raise ValueError(
                f"GLONASS band {band_num} is FDMA — a channel number (k) is required. "
                "The channel numbers are available in the "
                "'GLONASS SLOT / FRQ #' RINEX header record (GLO_Slot2ChannelMap)."
            )
        base, step = GLONASS_FDMA[band_num]
        return base + int(glonass_channel) * step

    return bands[band_num]


def wavelength(system: str, band: int, glonass_channel: Optional[int] = None) -> float:
    """Return the carrier wavelength in metres for *system* and *band*.

    Accepts the same arguments as :func:`carrier_frequency`.
    """
    return SPEED_OF_LIGHT / carrier_frequency(system, band, glonass_channel)


def band_label(system: str, band: int) -> str:
    """Return a short, system-specific description of a frequency band.

    Unlike :data:`BAND_DESCRIPTIONS`, which lists every constellation using
    a band digit, this gives the frequency of *that* system only::

        band_label('G', 1)   # 'L1 (1575.42 MHz)'
        band_label('R', 1)   # 'G1 (1602.000 + k x 0.5625 MHz)'
    """
    sys_code = _normalise_system(system)
    band_num = int(band)
    name = BAND_SIGNAL_NAMES.get(sys_code, {}).get(band_num, f'Band {band_num}')

    if sys_code == 'R' and band_num in GLONASS_FDMA:
        base, step = GLONASS_FDMA[band_num]
        return f"{name} ({_mhz(base)} + k x {step / 1e6:.4f} MHz)"

    freq = CARRIER_FREQUENCIES.get(sys_code, {}).get(band_num)
    if freq is None:
        return name
    return f"{name} ({_mhz(freq)} MHz)"


def _mhz(frequency_hz: float) -> str:
    """Format a frequency in MHz without trailing zeros."""
    return f"{frequency_hz / 1e6:.6f}".rstrip('0').rstrip('.')


def build_frequency_overview(gnss_systems: Mapping[int, str],
                             glo_slot2channel: Optional[Mapping[int, int]] = None,
                             max_glo_id: int = MAX_GLONASS_SLOT) -> Dict[int, np.ndarray]:
    """Build the per-system carrier frequency tables used by the analysis pipeline.

    Parameters
    ----------
    gnss_systems : mapping
        ``{system_index: system_code}`` as produced by ``readRinexObs``,
        e.g. ``{1: 'G', 2: 'R', 3: 'E', 4: 'C'}``.
    glo_slot2channel : mapping, optional
        ``{glonass_slot: channel_number}``.  Mandatory when GLONASS is present.
    max_glo_id : int
        Highest GLONASS slot number to expand the FDMA table for.

    Returns
    -------
    dict
        ``{system_index: ndarray}`` where each array is indexed by
        ``band - 1`` on the first axis.  Non-GLONASS systems get shape
        ``(9, 1)``; GLONASS gets shape ``(9, max_glo_id + 1)`` holding the
        per-satellite FDMA frequencies (column index = slot number).

    Raises
    ------
    ValueError
        If GLONASS is present but no channel mapping was supplied.
    """
    overview: Dict[int, np.ndarray] = {}

    for sys_idx, sys_code in gnss_systems.items():
        sys_code = _normalise_system(sys_code)
        if sys_code == 'R':
            overview[sys_idx] = _glonass_frequency_table(glo_slot2channel, max_glo_id)
        else:
            table = np.full((9, 1), np.nan)
            for band, freq in CARRIER_FREQUENCIES[sys_code].items():
                table[band - 1, 0] = freq
            overview[sys_idx] = table

    return overview


def _glonass_frequency_table(glo_slot2channel: Optional[Mapping[int, int]],
                             max_glo_id: int) -> np.ndarray:
    """Return the GLONASS ``(9, max_glo_id + 1)`` per-slot frequency table."""
    if not glo_slot2channel:
        raise ValueError(
            "ERROR! GLONASS k-numbers do not exist. This is mandatory to be able to "
            "compute GLONASS carrier frequencies. Please add GLONASS SLOT / FRQ to the "
            "RINEX header or use a RINEX navigation file instead of SP3."
        )

    table = np.full((9, max_glo_id + 1), np.nan)
    bands = CARRIER_FREQUENCIES['R']

    for band, freq in bands.items():
        row = band - 1
        base, step = GLONASS_FDMA.get(band, (freq, 0.0))
        for slot, channel in glo_slot2channel.items():
            slot = int(slot)
            if 0 <= slot <= max_glo_id:
                table[row, slot] = base + int(channel) * step

    return table


def _normalise_system(system: str) -> str:
    """Accept either a system code ('G') or a full name ('GPS')."""
    if system in SYSTEM_NAMES:
        return system
    if system in SYSTEM_CODES:
        return SYSTEM_CODES[system]
    return system
