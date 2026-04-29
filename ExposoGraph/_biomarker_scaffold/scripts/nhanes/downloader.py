"""Download and read NHANES XPT files with local caching."""

from __future__ import annotations

from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

import pandas as pd


def download_file(url: str, cache_dir: str | Path = "data/nhanes/raw", force: bool = False) -> Path:
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    filename = url.rstrip("/").split("/")[-1]
    out = cache_dir / filename
    if out.exists() and out.stat().st_size > 0 and not force:
        return out
    try:
        with urlopen(url, timeout=60) as response:
            content = response.read()
    except HTTPError as exc:
        raise RuntimeError(f"NHANES download failed with HTTP {exc.code}: {url}") from exc
    except URLError as exc:
        raise RuntimeError(f"NHANES download failed: {url}: {exc}") from exc
    if not content:
        raise RuntimeError(f"NHANES download returned an empty file: {url}")
    out.write_bytes(content)
    if out.stat().st_size == 0:
        raise RuntimeError(f"Cached NHANES file is empty: {out}")
    return out


def read_xpt(url: str, cache_dir: str | Path = "data/nhanes/raw", force: bool = False) -> pd.DataFrame:
    path = download_file(url, cache_dir=cache_dir, force=force)
    try:
        return pd.read_sas(path, format="xport")
    except Exception as exc:  # pandas can raise several parser exceptions
        raise RuntimeError(f"Failed to read NHANES XPT file {path}: {exc}") from exc
