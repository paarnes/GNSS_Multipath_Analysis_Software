"""
Matplotlib backend selection for the plotting routines.

The analysis writes its figures straight to file and never calls ``plt.show()``, so a batch
run has nothing to gain from a GUI backend and may not even have a display available.
``Agg`` is therefore selected for batch runs, but only when the analysis is started, never
on import: a library that takes the backend away from its host stops the host's own figures
from rendering, which is exactly what happens to a notebook that imports this package.

Made by: Per Helge Aarnes
E-mail: per.helge.aarnes@gmail.com
"""

import sys

import matplotlib


def use_non_interactive_backend() -> bool:
    """Select the non-interactive ``Agg`` backend unless the host already chose one.

    Return:
    ------
    ``True`` when the backend was switched to ``Agg``, ``False`` when the host's backend
    was left untouched.
    """
    if host_controls_the_backend():
        return False

    matplotlib.use('Agg')
    return True


def host_controls_the_backend() -> bool:
    """
    ``True`` when the code runs inside an interactive session that manages its own backend.

    A Jupyter kernel installs the ``inline`` backend, and switching away from it would send
    the user's plots to a file-only backend for the rest of the session.
    """
    return 'ipykernel' in sys.modules or 'matplotlib_inline' in sys.modules
