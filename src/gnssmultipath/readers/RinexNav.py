"""
Comprehensive RINEX navigation file reader supporting versions 2, 3, and 4.

Reads broadcast ephemerides for GPS, GLONASS, Galileo, and BeiDou from RINEX
navigation files and stores them in a structured ``RinexNavData`` dataclass.

Supported message types per system:

  * **GPS**     -- LNAV  (v2/v3/v4)
  * **GLONASS** -- FDMA  (v2/v3/v4)
  * **Galileo** -- INAV, FNAV, IFNV  (v3/v4)
  * **BeiDou**  -- D1, D2, D1D2  (v3/v4)

Usage::

    from gnssmultipath.readers.RinexNav import RinexNav
    nav_data = RinexNav.read_nav("path/to/nav.rnx")

Made by: Per Helge Aarnes
E-mail: per.helge.aarnes@gmail.com
"""

from __future__ import annotations

import re
import warnings
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from typing import List, Optional, Union

import numpy as np
from numpy import ndarray
from pandas import DataFrame
from tqdm import tqdm

warnings.filterwarnings("ignore")

# Number of broadcast-orbit data lines per system (excluding the epoch line).
_ORBIT_LINES: dict[str, int] = {"G": 7, "R": 3, "E": 7, "C": 7}

# Body-line count for RINEX v4 EPH records (epoch line + orbit lines).
_RINEX4_BODY_LEN: dict[str, int] = {"G": 8, "R": 5, "E": 8, "C": 8}

# Supported nav-message families for v4 filtering.
_SUPPORTED_V4_MESSAGES: dict[str, set[str]] = {
    "G": {"LNAV"},
    "R": {"FDMA"},
    "E": {"INAV", "FNAV", "IFNV"},
    "C": {"D1", "D2", "D1D2"},
}

# Total columns in the output ephemerides array.
_N_COLS = 36

# Regex that replaces Fortran-style 'D'/'d' exponent notation with 'E'.
_FORTRAN_EXP = re.compile(r"(?<=[0-9.])[Dd](?=[+-]?\d)")

class BroadcastColumns:
    """Column-name definitions for RINEX broadcast ephemeris parameters.

    Each GNSS system defines a tuple of column names that matches the order
    in the ``(n, 36)`` ephemeris array produced by :meth:`RinexNav.read_nav`.
    Inherit this class to gain named access to per-system DataFrames.

    Column layout
    -------------
    * ``0``      -- PRN string (e.g. ``'G01'``, ``'R24'``)
    * ``1-6``    -- Epoch: Year, Month, Day, Hour, Minute, Second
    * ``7-9``    -- Clock parameters (system-specific names)
    * ``10-35``  -- Broadcast-orbit parameters (system-specific)

    For **GLONASS**, only columns 0-21 carry meaningful values; the remaining
    columns are zero-padded to maintain uniform shape.
    """

    _COMMON = ("PRN", "Year", "Month", "Day", "Hour", "Minute", "Second")
    _CLOCK = ("af0", "af1", "af2")
    _CLOCK_GLO = ("TauN", "GammaN", "MessageFrameTime")

    GPS: tuple[str, ...] = _COMMON + _CLOCK + (
        "IODE", "Crs", "Delta_n", "M0",
        "Cuc", "e", "Cus", "sqrt_A",
        "Toe", "Cic", "OMEGA0", "Cis",
        "i0", "Crc", "omega", "OMEGA_DOT",
        "IDOT", "L2_codes", "GPS_week", "L2P_flag",
        "SV_accuracy", "SV_health", "TGD", "IODC",
        "Transmission_time", "Fit_interval",
    )

    GLONASS: tuple[str, ...] = _COMMON + _CLOCK_GLO + (
        "X", "X_vel", "X_acc", "health",
        "Y", "Y_vel", "Y_acc", "freq_num",
        "Z", "Z_vel", "Z_acc", "age",
    )

    GALILEO: tuple[str, ...] = _COMMON + _CLOCK + (
        "IODnav", "Crs", "Delta_n", "M0",
        "Cuc", "e", "Cus", "sqrt_A",
        "Toe", "Cic", "OMEGA0", "Cis",
        "i0", "Crc", "omega", "OMEGA_DOT",
        "IDOT", "Data_sources", "GAL_week", "Spare_1",
        "SISA", "SV_health", "BGD_E5a_E1", "BGD_E5b_E1",
        "Transmission_time", "Spare_2",
    )

    BEIDOU: tuple[str, ...] = _COMMON + _CLOCK + (
        "AODE", "Crs", "Delta_n", "M0",
        "Cuc", "e", "Cus", "sqrt_A",
        "Toe", "Cic", "OMEGA0", "Cis",
        "i0", "Crc", "omega", "OMEGA_DOT",
        "IDOT", "Spare_1", "BDT_week", "Spare_2",
        "SV_accuracy", "SatH1", "TGD1", "TGD2",
        "Transmission_time", "AODC",
    )

    SYSTEM_COLUMNS: dict[str, tuple[str, ...]] = {
        "G": GPS,
        "R": GLONASS,
        "E": GALILEO,
        "C": BEIDOU,
    }

@dataclass
class RinexNavData(BroadcastColumns):
    """Container for parsed RINEX navigation data.

    Inherits column-name definitions from :class:`BroadcastColumns`.

    Attributes
    ----------
    ephemerides : ndarray | DataFrame
        Ephemeris matrix of shape ``(n, 36)``.  Column layout:

        * ``0``     -- PRN string (e.g. ``'G01'``, ``'R24'``)
        * ``1-6``   -- Epoch: year, month, day, hour, minute, second
        * ``7-9``   -- Clock: af0 (bias), af1 (drift), af2 (drift-rate)
        * ``10-35`` -- Broadcast-orbit parameters (system-specific)

        For **GLONASS**, only columns 10-21 carry meaningful values (position,
        velocity, acceleration, health, frequency-number, age-of-data); the
        remaining columns are zero-padded to maintain uniform shape.
    header : list[str]
        Raw header lines (up to and including ``END OF HEADER``).
    nepochs : int
        Number of ephemeris records.
    glonass_fcn : dict[int, int] | None
        Mapping of GLONASS slot number to frequency-channel number (FCN).
        ``None`` when no GLONASS data is present.
    """

    ephemerides: Union[ndarray, DataFrame]
    header: list = field(default_factory=list)
    nepochs: int = 0
    glonass_fcn: Optional[dict] = None


    def _system_df(self, system: str) -> DataFrame:
        """Filter ephemerides for one GNSS system and return a named DataFrame.

        Parameters
        ----------
        system : str
            Single-character system identifier (``'G'``, ``'R'``, ``'E'``, ``'C'``).
        """
        cols = self.SYSTEM_COLUMNS[system]
        eph = self.ephemerides

        raw = eph.values if isinstance(eph, DataFrame) else eph
        if raw.size == 0:
            return DataFrame(columns=cols)

        prns = raw[:, 0].astype(str)
        mask = np.char.startswith(prns, system)
        if not np.any(mask):
            return DataFrame(columns=cols)

        n = len(cols)
        subset = DataFrame(raw[mask][:, :n].copy(), columns=cols)
        float_cols = list(cols[1:])
        subset[float_cols] = subset[float_cols].astype(float)
        return subset.reset_index(drop=True)

    @property
    def gps(self) -> DataFrame:
        """GPS broadcast ephemerides as a :class:`~pandas.DataFrame` with named columns."""
        return self._system_df("G")

    @property
    def glonass(self) -> DataFrame:
        """GLONASS broadcast ephemerides as a :class:`~pandas.DataFrame` with named columns."""
        return self._system_df("R")

    @property
    def galileo(self) -> DataFrame:
        """Galileo broadcast ephemerides as a :class:`~pandas.DataFrame` with named columns."""
        return self._system_df("E")

    @property
    def beidou(self) -> DataFrame:
        """BeiDou broadcast ephemerides as a :class:`~pandas.DataFrame` with named columns."""
        return self._system_df("C")

def _parse_version(first_line: str) -> int:
    """Return the major RINEX version from the first header line."""
    try:
        return int(float(first_line[:9].strip()))
    except ValueError as exc:
        raise ValueError(f"Cannot parse RINEX version from: {first_line.rstrip()!r}") from exc


def _fortran_to_float_str(text: str) -> str:
    """Replace Fortran D/d exponent notation with E."""
    return _FORTRAN_EXP.sub("E", text)


def _parse_v3_data_line(line: str) -> list[str]:
    """Parse a fixed-width (4 x 19-char) v3/v4 broadcast-orbit data line.

    Each orbit line has format ``3X,4D19.12`` -- four 19-character fields
    starting at column 4.  Falls back to whitespace splitting for
    short/truncated lines.
    """
    line = _fortran_to_float_str(line.rstrip())
    fields: list[str] = []
    starts = (4, 23, 42, 61)
    for s in starts:
        chunk = line[s : s + 19].strip() if len(line) > s else ""
        if chunk:
            fields.append(chunk)
    if not fields:
        fields = line.split()
    return fields


def _parse_v3_epoch_line(line: str) -> list:
    """Parse a v3/v4 epoch line into ``[PRN, Y, M, D, H, Min, Sec, af0, af1, af2]``."""
    line = _fortran_to_float_str(line.rstrip())
    prn = line[:3].strip()
    year   = line[4:8].strip()
    month  = line[9:11].strip()
    day    = line[12:14].strip()
    hour   = line[15:17].strip()
    minute = line[18:20].strip()
    second = line[21:23].strip()
    af0 = line[23:42].strip()
    af1 = line[42:61].strip()
    af2 = line[61:80].strip() if len(line) > 61 else "0"
    return [prn, year, month, day, hour, minute, int(float(second)), af0, af1, af2]


def _parse_v2_epoch_line(line: str, sys_char: str = "G") -> list:
    """Parse a RINEX v2 epoch line into v3-compatible tokens.

    Works for both GPS (`sys_char='G'`) and GLONASS (`sys_char='R'`) v2
    navigation files.  The slot/PRN integer in the first field is
    rewritten as ``f"{sys_char}{int:02d}"``.
    """
    line = _fortran_to_float_str(line.rstrip())
    # Ensure space between seconds field and af0 (position 22).
    if len(line) > 22 and line[22] not in (" ", "\t"):
        line = line[:22] + " " + line[22:]
    # Insert spaces after exponent+digits where fields are stuck together.
    line = re.sub(r"([eE][+-]?\d\d)(-)", r"\1 \2", line)
    parts = line.split()
    if not parts:
        return []
    prn = f"{sys_char}{int(parts[0]):02d}"
    year_2d = int(parts[1])
    year = str(2000 + year_2d) if year_2d < 80 else str(1900 + year_2d)
    second = int(float(parts[6]))
    return [prn, year, parts[2], parts[3], parts[4], parts[5], second] + parts[7:10]


def _parse_v2_data_line(line: str) -> list[str]:
    """Parse a RINEX v2 broadcast-orbit data line (``3X,4D19.12``)."""
    line = _fortran_to_float_str(line.rstrip())
    line = re.sub(r"([eE][+-]?\d\d)(-)", r"\1 \2", line)
    return line.split()

def _read_header(filename: str) -> tuple[int, list[str]]:
    """Read header lines and detect RINEX version.

    Returns
    -------
    version : int
        Major RINEX version (2, 3, or 4).
    header : list[str]
        All header lines including ``END OF HEADER``.
    """
    header: list[str] = []
    with open(filename, "r") as f:
        for line in f:
            stripped = line.rstrip()
            header.append(stripped)
            if "END OF HEADER" in stripped:
                break
    if not header:
        raise ValueError(f"Empty or unreadable navigation file: {filename}")
    return _parse_version(header[0]), header


def _read_body_lines(filename: str) -> list[str]:
    """Read all lines after the header."""
    lines: list[str] = []
    past_header = False
    with open(filename, "r") as f:
        for line in f:
            if past_header:
                lines.append(line)
            elif "END OF HEADER" in line:
                past_header = True
    return lines

_V2_EPOCH_RE = re.compile(r"^\s*\d{1,2}\s+\d{2}\s+\d{1,2}\s+\d{1,2}")
_V3_EPOCH_RE = re.compile(r"^[GREC]\d{2}")


def _extract_v2_blocks(body: list[str], n_lines: int = 8) -> list[list[str]]:
    """Extract fixed-size ephemeris blocks from RINEX v2 body lines.

    ``n_lines`` is 8 for GPS (epoch + 7 broadcast-orbit lines) and 4 for
    GLONASS (epoch + 3 broadcast-orbit lines).
    """
    blocks: list[list[str]] = []
    idx = 0
    n = len(body)
    while idx < n:
        if _V2_EPOCH_RE.match(body[idx]):
            end = min(idx + n_lines, n)
            block = body[idx:end]
            if len(block) == n_lines:
                blocks.append(block)
            idx = end
        else:
            idx += 1
    return blocks


def _extract_v3_blocks(body: list[str], desired_systems: set[str]) -> list[list[str]]:
    """Extract ephemeris blocks from RINEX v3 body lines."""
    blocks: list[list[str]] = []
    idx = 0
    n = len(body)
    while idx < n:
        if _V3_EPOCH_RE.match(body[idx]):
            sys_char = body[idx][0]
            n_orbit = _ORBIT_LINES.get(sys_char, 7)
            end = idx + 1 + n_orbit
            if sys_char in desired_systems and end <= n:
                blocks.append(body[idx:end])
            idx = max(end, idx + 1)
        else:
            idx += 1
    return blocks


def _extract_v4_blocks(all_lines: list[str], desired_systems: set[str]) -> list[list[str]]:
    """Extract ephemeris blocks from RINEX v4 ``> EPH`` records.

    Parameters
    ----------
    all_lines : list[str]
        All file lines (including header).
    desired_systems : set[str]
        System letters to include (e.g. ``{'G', 'R', 'E', 'C'}``).

    Returns
    -------
    list[list[str]]
        Body lines of each matching EPH record.
    """
    header_end = 0
    for i, line in enumerate(all_lines):
        if "END OF HEADER" in line:
            header_end = i + 1
            break

    blocks: list[list[str]] = []
    idx = header_end
    n = len(all_lines)
    while idx < n:
        line = all_lines[idx]
        if not line.startswith("> EPH"):
            idx += 1
            continue

        tokens = line.split()
        if len(tokens) < 4:
            idx += 1
            continue

        sat_id = tokens[2]
        sys_char = sat_id[0]
        msg_type = tokens[3].upper()

        # Find the next record header.
        end = idx + 1
        while end < n and not all_lines[end].startswith(">"):
            end += 1

        body = all_lines[idx + 1 : end]
        expected = _RINEX4_BODY_LEN.get(sys_char, _ORBIT_LINES.get(sys_char, 0) + 1)

        if (
            sys_char in desired_systems
            and msg_type in _SUPPORTED_V4_MESSAGES.get(sys_char, set())
        ):
            if len(body) != expected:
                raise ValueError(
                    f"Unexpected number of data lines for {sat_id} ({msg_type}): "
                    f"expected {expected}, got {len(body)}."
                )
            blocks.append(body)

        idx = end

    return blocks

def _parse_v2_block(block: list[str], sys_char: str = "G") -> ndarray:
    """Parse one v2 ephemeris block into a ``(1, 36)`` object ndarray.

    GPS blocks are 8 lines (epoch + 7 orbit lines).  GLONASS blocks are
    4 lines (epoch + 3 orbit lines) and are zero-padded to 36 columns so
    that GLONASS rows interoperate with the GPS layout downstream
    (column 17 = frequency-channel number, matching the v3 layout).
    """
    tokens = _parse_v2_epoch_line(block[0], sys_char=sys_char)
    if not tokens:
        return np.empty((0, _N_COLS), dtype=object)
    row: list = list(tokens)
    for line in block[1:]:
        row.extend(_parse_v2_data_line(line))
    pad = "0" if sys_char == "R" else "nan"
    row = row[:_N_COLS]
    while len(row) < _N_COLS:
        row.append(pad)
    return np.array(row, dtype=object).reshape(1, _N_COLS)


def _parse_v3_block(block: list[str], is_glonass: bool) -> ndarray:
    """Parse one v3/v4 ephemeris block into a ``(1, 36)`` object ndarray.

    GLONASS blocks (fewer orbit parameters) are zero-padded to 36 columns.
    Non-GLONASS blocks with missing trailing fields get ``'nan'``.
    """
    tokens = _parse_v3_epoch_line(block[0])
    row: list = list(tokens)
    for line in block[1:]:
        row.extend(_parse_v3_data_line(line))
    pad = "0" if is_glonass else "nan"
    row = row[:_N_COLS]
    while len(row) < _N_COLS:
        row.append(pad)
    return np.array(row, dtype=object).reshape(1, _N_COLS)

def _filter_on_time(blocks: list[list[str]], interval_minutes: float) -> list[list[str]]:
    """Keep at most one ephemeris per satellite per *interval_minutes*.

    Tracks each satellite independently so that interleaved multi-satellite
    data is filtered correctly regardless of ordering.

    Parameters
    ----------
    blocks : list[list[str]]
        Ephemeris blocks (each is a list of raw lines).
    interval_minutes : float
        Minimum time gap between retained records of the same satellite.
        ``0`` disables filtering.
    """
    if not blocks or interval_minutes <= 0:
        return blocks

    delta = timedelta(minutes=interval_minutes)
    filtered: list[list[str]] = []
    last_epoch: dict[str, datetime] = {}

    for block in blocks:
        parts = block[0].split()
        sat = parts[0]
        epoch = datetime(*map(int, parts[1:6]))
        prev = last_epoch.get(sat)
        if prev is None or (epoch - prev) >= delta:
            filtered.append(block)
            last_epoch[sat] = epoch

    return filtered

def _extract_glonass_fcn(data: ndarray) -> Optional[dict[int, int]]:
    """Extract GLONASS frequency-channel numbers from the ephemeris array.

    Returns ``{slot_number: fcn}`` or ``None`` when no GLONASS rows exist.
    """
    prns = data[:, 0].astype(str)
    glo_mask = np.char.startswith(prns, "R")
    if not np.any(glo_mask):
        return None
    glo_data = data[glo_mask]
    _, idx = np.unique(glo_data[:, 0], return_index=True)
    unique = glo_data[idx]
    return {int(p[1:]): int(float(f)) for p, f in zip(unique[:, 0], unique[:, 17])}

class RinexNav:
    """Unified RINEX navigation-file reader (v2 / v3 / v4).

    Only ``@staticmethod`` and ``@classmethod`` methods -- no instance state.
    Primary entry point: :meth:`read_nav`.

    Examples
    --------
    >>> from gnssmultipath.readers.RinexNav import RinexNav
    >>> nav = RinexNav.read_nav("BRDC00IGS_R_20220010000_01D_MN.rnx")
    >>> nav.ephemerides.shape
    (1234, 36)
    """

    @staticmethod
    def read_nav(
        filename: str,
        desired_GNSS: Optional[List[str]] = None,
        dataframe: bool = False,
        data_rate: float = 30,
    ) -> RinexNavData:
        """Read a RINEX navigation file (auto-detecting version).

        Parameters
        ----------
        filename : str
            Path to the RINEX navigation file.
        desired_GNSS : list[str] | None
            System letters to include (e.g. ``['G', 'R', 'E', 'C']``).
            Defaults to all four global systems.
        dataframe : bool
            If ``True``, return ephemerides as a :class:`~pandas.DataFrame`.
        data_rate : float
            Minimum interval in minutes between retained ephemerides of the
            same satellite.  ``0`` keeps all records.
        """
        if desired_GNSS is None:
            desired_GNSS = ["G", "R", "E", "C"]
        version, header = _read_header(filename)
        if version == 2:
            return _read_v2(filename, header, dataframe=dataframe,
                            desired_systems=set(desired_GNSS))
        return _read_v3v4(
            filename, header, version,
            desired_systems=set(desired_GNSS),
            dataframe=dataframe,
            data_rate=data_rate,
        )

def _v2_system_from_header(header: list[str]) -> str:
    """Detect the GNSS system encoded by a RINEX v2 nav file from its header.

    Returns ``'R'`` for GLONASS nav files (header contains
    ``GLONASS NAV DATA``) and ``'G'`` otherwise (default GPS).
    """
    for line in header:
        upper = line.upper()
        if "RINEX VERSION / TYPE" in upper:
            if "GLONASS" in upper:
                return "R"
            # Galileo / BeiDou v2 nav are not standard; fall through to GPS.
            return "G"
    return "G"


def _read_v2(filename: str, header: list[str], dataframe: bool = False,
             desired_systems: Optional[set[str]] = None) -> RinexNavData:
    """Read a RINEX v2 navigation file (GPS or GLONASS)."""
    print("Reading broadcast ephemeris from RINEX-navigation file.....")
    sys_char = _v2_system_from_header(header)
    if desired_systems is not None and sys_char not in desired_systems:
        empty = np.empty((0, _N_COLS), dtype=object)
        return RinexNavData(ephemerides=empty, header=header, nepochs=0,
                            glonass_fcn=None)

    n_block_lines = 4 if sys_char == "R" else 8
    body = _read_body_lines(filename)
    blocks = _extract_v2_blocks(body, n_lines=n_block_lines)
    if not blocks:
        return RinexNavData(ephemerides=np.empty((0, _N_COLS), dtype=object),
                            header=header, nepochs=0, glonass_fcn=None)

    rows = [_parse_v2_block(b, sys_char=sys_char) for b in blocks if b]
    rows = [r for r in rows if r.size > 0]
    if not rows:
        return RinexNavData(ephemerides=np.empty((0, _N_COLS), dtype=object),
                            header=header, nepochs=0, glonass_fcn=None)

    data = np.vstack(rows).astype(str)
    glo_fcn = _extract_glonass_fcn(data) if sys_char == "R" else None
    n_eph = len(data)
    if dataframe in (True, "yes", "YES"):
        data = DataFrame(data)
    return RinexNavData(ephemerides=data, header=header, nepochs=n_eph,
                        glonass_fcn=glo_fcn)


def _read_v3v4(
    filename: str,
    header: list[str],
    version: int,
    desired_systems: set[str],
    dataframe: bool = False,
    data_rate: float = 30,
) -> RinexNavData:
    """Read a RINEX v3 or v4 navigation file."""
    for s in desired_systems:
        if s not in {"G", "R", "E", "C"}:
            raise ValueError(f"Invalid GNSS system '{s}'. Must be one of 'G', 'R', 'E', 'C'.")

    if version >= 4:
        with open(filename, "r") as f:
            all_lines = f.readlines()
        raw_blocks = _extract_v4_blocks(all_lines, desired_systems)
    else:
        body = _read_body_lines(filename)
        raw_blocks = _extract_v3_blocks(body, desired_systems)

    blocks = _filter_on_time(raw_blocks, data_rate)
    if not blocks:
        empty = np.empty((0, _N_COLS), dtype=object)
        if dataframe:
            empty = DataFrame(empty)
        return RinexNavData(ephemerides=empty, header=header, nepochs=0, glonass_fcn=None)

    bar_fmt = "{desc}: {percentage:3.0f}%|{bar}| ({n_fmt}/{total_fmt})"
    n_blocks = len(blocks)
    n_steps = min(100, n_blocks)
    step = max(1, n_blocks // n_steps)

    rows: list[ndarray] = []
    with tqdm(total=n_steps, desc="Rinex navigation file is being read", bar_format=bar_fmt) as pbar:
        for i, block in enumerate(blocks):
            sys_char = block[0].strip()[0]
            rows.append(_parse_v3_block(block, is_glonass=(sys_char == "R")))
            if (i + 1) % step == 0:
                pbar.update(1)
        pbar.update(n_steps - pbar.n)

    data = np.vstack(rows).astype(str)
    glo_fcn = _extract_glonass_fcn(data)
    n_eph = len(data)

    if dataframe:
        data = DataFrame(data)
        float_cols = list(data.columns[1:])
        data[float_cols] = data[float_cols].astype(float)

    return RinexNavData(ephemerides=data, header=header, nepochs=n_eph, glonass_fcn=glo_fcn)
