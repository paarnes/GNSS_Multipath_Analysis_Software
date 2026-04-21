"""
Shared configuration for plot tests.

We force the non-interactive ``Agg`` matplotlib backend before importing the
package so the test suite can run on CI / headless workstations and so that no
interactive figure windows are spawned during the tests.
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg", force=True)
