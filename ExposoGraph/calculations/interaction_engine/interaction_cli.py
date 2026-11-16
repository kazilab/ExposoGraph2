"""CLI wrapper for :mod:`ExposoGraph.interaction_engine`.

Use ``python -m ExposoGraph.interaction_cli`` to avoid the ``runpy`` warning
that can appear when executing ``ExposoGraph.interaction_engine`` directly
while the package re-exports that module from ``ExposoGraph.__init__``.
"""

from __future__ import annotations

from .interaction_engine import main

if __name__ == "__main__":
    raise SystemExit(main())
