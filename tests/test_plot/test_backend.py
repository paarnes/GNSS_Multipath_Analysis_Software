"""
Tests for the matplotlib backend handling.

Importing the package must not take the backend away from the host application, but a
batch analysis should still end up on the non-interactive Agg backend.
"""
import os
import subprocess
import sys

project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
SRC = os.path.join(project_path, 'src')


def _run(code: str) -> str:
    """Run a snippet in a clean interpreter and return its stdout."""
    env = dict(os.environ, PYTHONPATH=SRC)
    result = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True, env=env)
    assert result.returncode == 0, result.stderr
    return result.stdout.strip()


def test_import_does_not_change_the_backend():
    """A library must leave the host application's backend alone."""
    backend = _run(
        "import matplotlib;"
        "matplotlib.use('Agg');"
        "before = matplotlib.get_backend();"
        "import gnssmultipath;"
        "import gnssmultipath.plot.make_polarplot;"
        "import gnssmultipath.plot.plotResults;"
        "import gnssmultipath.plot.SkyPlotSummary;"
        "print(before == matplotlib.get_backend())"
    )
    assert backend == 'True'


def test_import_keeps_a_non_agg_backend():
    """The template backend stands in for a host backend that must survive the import."""
    backend = _run(
        "import matplotlib;"
        "matplotlib.use('template');"
        "import gnssmultipath.plot.make_polarplot;"
        "print(matplotlib.get_backend())"
    )
    assert backend == 'template'


def test_batch_run_selects_agg():
    """Outside an interactive session the analysis switches to Agg itself."""
    backend = _run(
        "import matplotlib;"
        "matplotlib.use('template');"
        "from gnssmultipath.plot.backend import use_non_interactive_backend;"
        "switched = use_non_interactive_backend();"
        "print(f'{switched} {matplotlib.get_backend().lower()}')"
    )
    assert backend == 'True agg'


def test_interactive_session_keeps_its_backend():
    """An imported ipykernel marks a session that manages its own backend."""
    backend = _run(
        "import sys, types;"
        "sys.modules['ipykernel'] = types.ModuleType('ipykernel');"
        "import matplotlib;"
        "matplotlib.use('template');"
        "from gnssmultipath.plot.backend import use_non_interactive_backend;"
        "switched = use_non_interactive_backend();"
        "print(f'{switched} {matplotlib.get_backend()}')"
    )
    assert backend == 'False template'
