"""
Parsed RINEX observation code (signal) metadata.

An :class:`ObsCode` turns an opaque observation code string such as ``'C1C'``
or ``'P2'`` into a typed object that knows its observation type, frequency
band, tracking attribute and carrier frequency.

Both RINEX conventions are supported:

* **RINEX 3/4** — three characters, ``<type><band><attribute>``, e.g. ``L5Q``.
* **RINEX 2**   — two characters, ``<type><band>``, e.g. ``P2`` or ``L1``.
  The tracking attribute is then an empty string, and the RINEX 2 ``'P'``
  type is recognised as a pseudorange.

Usage
-----
::

    from gnssmultipath.readers.ObsCode import ObsCode

    sig = ObsCode.parse('L5Q', 'E')
    sig.obs_type          # 'L'
    sig.band              # 5
    sig.attribute         # 'Q'
    sig.type_name         # 'phase'
    sig.is_phase          # True
    sig.frequency()       # 1176450000.0
    sig.wavelength()      # 0.2548...

    glo = ObsCode.parse('C1C', 'R')
    glo.frequency(glonass_channel=-4)   # FDMA-aware

Made by: Per Helge Aarnes
E-mail: per.helge.aarnes@gmail.com
"""

from dataclasses import dataclass
from typing import Optional

from gnssmultipath.constants import (
    OBS_TYPE_NAMES,
    PSEUDORANGE_TYPES,
    SYSTEM_NAMES,
    band_label,
    carrier_frequency,
    wavelength,
)


@dataclass(frozen=True, slots=True)
class ObsCode:
    """A parsed RINEX observation code.

    Attributes
    ----------
    code : str
        The original code string, e.g. ``'C1C'`` or ``'P2'``.
    system : str
        Single-character GNSS system code (``'G'``, ``'R'``, ``'E'``, ``'C'``).
    obs_type : str
        First character of the code: ``'C'``/``'P'`` (pseudorange), ``'L'``
        (phase), ``'D'`` (Doppler), ``'S'`` (SNR).
    band : int
        RINEX frequency band digit.
    attribute : str
        Tracking-mode attribute (RINEX 3/4 only), ``''`` for RINEX 2 codes.
    """

    code: str
    system: str
    obs_type: str
    band: int
    attribute: str

    # ── Construction ─────────────────────────────────────────────────────

    @classmethod
    def parse(cls, code: str, system: str) -> 'ObsCode':
        """Parse *code* for *system*.

        Raises
        ------
        ValueError
            If the code is not a valid 2- or 3-character RINEX observation code.
        """
        if not isinstance(code, str):
            raise ValueError(f"Observation code must be a string, got {type(code).__name__}")

        text = code.strip()
        if len(text) not in (2, 3):
            raise ValueError(
                f"'{code}' is not a valid RINEX observation code — expected 2 characters "
                "(RINEX 2, e.g. 'P2') or 3 characters (RINEX 3/4, e.g. 'C1C')."
            )

        obs_type = text[0].upper()
        if obs_type not in OBS_TYPE_NAMES:
            raise ValueError(
                f"'{code}' has an unknown observation type '{obs_type}'. "
                f"Valid types: {sorted(OBS_TYPE_NAMES)}"
            )

        if not text[1].isdigit():
            raise ValueError(f"'{code}' has a non-numeric frequency band '{text[1]}'.")

        return cls(
            code=text,
            system=system,
            obs_type=obs_type,
            band=int(text[1]),
            attribute=text[2].upper() if len(text) == 3 else '',
        )

    # ── Type classification ──────────────────────────────────────────────

    @property
    def type_name(self) -> str:
        """Human-readable observation type, e.g. ``'pseudorange'``."""
        return OBS_TYPE_NAMES.get(self.obs_type, 'unknown')

    @property
    def is_pseudorange(self) -> bool:
        """True for ``'C'`` and RINEX 2 ``'P'`` codes."""
        return self.obs_type in PSEUDORANGE_TYPES

    @property
    def is_phase(self) -> bool:
        """True for carrier-phase (``'L'``) codes."""
        return self.obs_type == 'L'

    @property
    def is_doppler(self) -> bool:
        """True for Doppler (``'D'``) codes."""
        return self.obs_type == 'D'

    @property
    def is_snr(self) -> bool:
        """True for signal-strength (``'S'``) codes."""
        return self.obs_type == 'S'

    # ── Descriptive metadata ─────────────────────────────────────────────

    @property
    def system_name(self) -> str:
        """Full GNSS system name, e.g. ``'Galileo'``."""
        return SYSTEM_NAMES.get(self.system, self.system)

    @property
    def band_description(self) -> str:
        """Signal name and frequency for this system, e.g. ``'E5a (1176.45 MHz)'``."""
        return band_label(self.system, self.band)

    # ── Frequency ────────────────────────────────────────────────────────

    def frequency(self, glonass_channel: Optional[int] = None) -> float:
        """Carrier frequency in Hz.

        *glonass_channel* is required for GLONASS FDMA bands 1 and 2.
        """
        return carrier_frequency(self.system, self.band, glonass_channel)

    def wavelength(self, glonass_channel: Optional[int] = None) -> float:
        """Carrier wavelength in metres.

        *glonass_channel* is required for GLONASS FDMA bands 1 and 2.
        """
        return wavelength(self.system, self.band, glonass_channel)

    # ── String interop ───────────────────────────────────────────────────

    def __str__(self) -> str:
        return self.code

    def __eq__(self, other) -> bool:
        """Compare against another ObsCode or directly against a code string."""
        if isinstance(other, str):
            return self.code == other
        if isinstance(other, ObsCode):
            return (self.code, self.system) == (other.code, other.system)
        return NotImplemented

    def __hash__(self) -> int:
        return hash((self.system, self.code))

    def __repr__(self) -> str:
        return (f"ObsCode('{self.code}', system='{self.system}', "
                f"type='{self.type_name}', band={self.band}, "
                f"attribute='{self.attribute}')")
