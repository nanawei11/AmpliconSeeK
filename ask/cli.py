from __future__ import annotations

import sys
from pathlib import Path


def _ensure_legacy_import_path() -> None:
    package_dir = Path(__file__).resolve().parent
    text = str(package_dir)
    if text not in sys.path:
        sys.path.insert(0, text)


def ask_main() -> None:
    _ensure_legacy_import_path()
    from . import ask_cmd

    ask_cmd.main()


def ecdna_search_main() -> None:
    _ensure_legacy_import_path()
    from . import ecDNA_search

    ecDNA_search.main()
