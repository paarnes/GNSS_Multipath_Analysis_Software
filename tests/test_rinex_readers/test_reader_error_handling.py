"""
Tests for error handling and header edge cases in the RINEX observation readers.
"""
import sys
import os
import textwrap
import pytest
import numpy as np
from pathlib import Path

project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.append(os.path.join(project_path, 'src'))

from gnssmultipath import readRinexObs
from gnssmultipath.readers.readRinexObs import rinexReadObsFileHeader304

NAV_FILE = os.path.join(project_path, "TestData", "NavigationFiles", "v3",
                        "BRDC00IGS_R_20220010000_01D_MN.rnx")

_OBS_BODY = [
    "> 2022 01 01 00 00  0.0000000  0  1",
    "G01  23000000.000 7 120000000.000 7      1000.000       45.000",
]


def _write(tmp_path, name, lines):
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def _header(obs_lines, marker="SENTINEL", version_line=None):
    default_version = ("     3.04           OBSERVATION DATA    M                   "
                       "RINEX VERSION / TYPE")
    return [
        version_line or default_version,
        *obs_lines,
        marker.ljust(60) + "MARKER NAME",
        "  2022     1     1     0     0    0.0000000     GPS         TIME OF FIRST OBS",
        "".ljust(60) + "END OF HEADER",
    ]


class TestWrongFileType:
    """A non-observation file must give an actionable error.

    Regression test: the header reader used to `return` None on this path while
    the caller unpacked 23 values, so the user got
    ``TypeError: cannot unpack non-iterable NoneType object``.
    """

    def test_navigation_file_raises_value_error(self):
        with pytest.raises(ValueError, match="not a RINEX observation file"):
            readRinexObs(NAV_FILE)

    def test_error_names_the_file(self):
        with pytest.raises(ValueError, match="BRDC00IGS"):
            readRinexObs(NAV_FILE)

    def test_error_is_not_a_typeerror(self):
        with pytest.raises(Exception) as excinfo:
            readRinexObs(NAV_FILE)
        assert not isinstance(excinfo.value, TypeError)

    def test_unsupported_system_letter_raises(self, tmp_path):
        version_line = ("     3.04           OBSERVATION DATA    J                   "
                        "RINEX VERSION / TYPE")
        lines = _header(["G    1 C1C                                                  "
                         "SYS / # / OBS TYPES"], version_line=version_line) + _OBS_BODY
        path = _write(tmp_path, "qzss.rnx", lines)
        with pytest.raises(ValueError, match="invalid satellite system type"):
            readRinexObs(path)

    def test_non_string_filename_raises_type_error(self):
        with pytest.raises(TypeError, match="Must be a string or os.PathLike"):
            rinexReadObsFileHeader304(123, 1, 1, ['G'], ['C'], [1])

    def test_pathlib_path_is_accepted(self, tmp_path):
        lines = _header(["G    1 C1C                                                  "
                         "SYS / # / OBS TYPES"]) + _OBS_BODY
        path = Path(_write(tmp_path, "pathlike.rnx", lines))
        assert readRinexObs(path).observations.systems == ['G']


class TestNoObservationRecords:
    """A header-only file must say so instead of failing later."""

    def test_header_only_file_raises(self, tmp_path):
        lines = _header(["G    1 C1C                                                  "
                         "SYS / # / OBS TYPES"])
        path = _write(tmp_path, "headeronly.rnx", lines)
        with pytest.raises(ValueError, match="No epoch records found"):
            readRinexObs(path)

    def test_empty_file_raises(self, tmp_path):
        path = tmp_path / "empty.rnx"
        path.write_text("")
        with pytest.raises(ValueError, match="empty"):
            readRinexObs(str(path))


class TestObsCodeContinuationLines:
    """`SYS / # / OBS TYPES` continuation handling.

    Regression test: the continuation was read whenever the code count hit a
    multiple of 13, so a system with exactly 13*n codes consumed the following
    header record (silently losing e.g. TIME OF FIRST OBS or MARKER NAME).
    """

    @staticmethod
    def _sys_lines(codes):
        first = ("G  %3d " % len(codes)) + " ".join(codes[:13])
        lines = [first.ljust(60) + "SYS / # / OBS TYPES"]
        for start in range(13, len(codes), 13):
            cont = "       " + " ".join(codes[start:start + 13])
            lines.append(cont.ljust(60) + "SYS / # / OBS TYPES")
        return lines

    ALL_CODES = ["C1C", "L1C", "D1C", "S1C", "C1W", "L1W", "D1W", "S1W",
                 "C2W", "L2W", "D2W", "S2W", "C2L", "L2L", "D2L", "S2L",
                 "C5Q", "L5Q", "D5Q", "S5Q", "C1L", "L1L", "D1L", "S1L",
                 "C6C", "L6C", "D6C", "S6C", "C7Q", "L7Q", "D7Q", "S7Q",
                 "C8Q", "L8Q", "D8Q", "S8Q", "C3I", "L3I", "D3I"]

    @pytest.mark.parametrize("n_codes", [1, 12, 13, 14, 25, 26, 27, 39])
    def test_all_codes_are_parsed(self, tmp_path, n_codes):
        codes = self.ALL_CODES[:n_codes]
        lines = _header(self._sys_lines(codes)) + _OBS_BODY
        path = _write(tmp_path, f"codes{n_codes}.rnx", lines)

        data = readRinexObs(path)
        assert data.obsCodes[1]['G'] == codes

    @pytest.mark.parametrize("n_codes", [13, 26, 39])
    def test_following_header_record_is_not_swallowed(self, tmp_path, n_codes):
        """Exact multiples of 13 are the case that used to break."""
        codes = self.ALL_CODES[:n_codes]
        lines = _header(self._sys_lines(codes), marker="SENTINEL") + _OBS_BODY
        path = _write(tmp_path, f"marker{n_codes}.rnx", lines)

        data = readRinexObs(path)
        assert data.markerName == "SENTINEL"

    @pytest.mark.parametrize("n_codes", [13, 26])
    def test_time_of_first_obs_survives(self, tmp_path, n_codes):
        codes = self.ALL_CODES[:n_codes]
        lines = _header(self._sys_lines(codes)) + _OBS_BODY
        path = _write(tmp_path, f"tfirst{n_codes}.rnx", lines)

        data = readRinexObs(path)
        assert int(data.tFirstObs[0][0]) == 2022


class TestMarkerName:
    """markerName must hold the value only, not the RINEX label."""

    def test_label_is_stripped(self, tmp_path):
        lines = _header(["G    1 C1C                                                  "
                         "SYS / # / OBS TYPES"], marker="OPEC") + _OBS_BODY
        path = _write(tmp_path, "marker.rnx", lines)
        data = readRinexObs(path)
        assert data.markerName == "OPEC"
        assert "MARKER NAME" not in data.markerName

    def test_blank_marker_gives_empty_string(self, tmp_path):
        lines = _header(["G    1 C1C                                                  "
                         "SYS / # / OBS TYPES"], marker="") + _OBS_BODY
        path = _write(tmp_path, "blankmarker.rnx", lines)
        data = readRinexObs(path)
        assert data.markerName == ""


class TestInvalidArguments:
    """Invalid flags must raise rather than return None."""

    @pytest.fixture(scope="class")
    @classmethod
    def obs_file(cls):
        return os.path.join(project_path, "TestData", "ObservationFiles", "v3",
                            "OPEC00NOR_S_20220010000_01D_30S_MO_3.04_croped.rnx")

    def test_bad_read_ss(self, obs_file):
        with pytest.raises(ValueError, match="readSS must be either 1 or 0"):
            readRinexObs(obs_file, readSS=5)

    def test_bad_read_lli(self, obs_file):
        with pytest.raises(ValueError, match="readLLI must be either 1 or 0"):
            readRinexObs(obs_file, readLLI=5)


if __name__ == '__main__':
    pytest.main()
