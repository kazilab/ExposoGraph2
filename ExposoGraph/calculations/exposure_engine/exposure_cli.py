"""CLI wrapper for :mod:`ExposoGraph.exposure_engine`.

Use ``python -m ExposoGraph.exposure_cli`` to avoid the ``runpy`` warning that
can appear when executing ``ExposoGraph.exposure_engine`` directly while the
package re-exports that module from ``ExposoGraph.__init__``.
"""

from __future__ import annotations

from .exposure_engine import main

if __name__ == "__main__":
    raise SystemExit(main())
