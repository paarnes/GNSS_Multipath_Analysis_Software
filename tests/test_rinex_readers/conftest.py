"""
Shared configuration for the whole test suite.

The package no longer selects a matplotlib backend on import, so the test suite picks the
non-interactive ``Agg`` backend itself. That keeps the tests runnable on CI and headless
workstations, and stops figure windows from being spawned no matter which module a test
imports first.
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg", force=True)
