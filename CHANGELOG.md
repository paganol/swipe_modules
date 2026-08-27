# HEAD

- Add `.gitignore` for build/dist artifacts, notebook checkpoints, and editor/OS files [#9]
- Fix the README's Usage snippet, which called `np.deg2rad(...)` without importing numpy
  (`NameError` if copy-pasted as-is), and bump its `conda create ... python=3.8` example,
  which predated `pyproject.toml`'s `requires-python = ">=3.10,<3.14"` [#10]
- Add `examples/test_swipe_raster_fort_sumner.ipynb`: a `SwipeRasterScanningStrategy`
  example from Fort Sumner, NM (the NASA long-duration balloon launch site) that generates
  an input sky (CMB + foregrounds), observes it, bins it back into a map, and overlays the
  Sun, Moon, outer planets, the Galactic plane, and the Crab Nebula on the scan coverage [#11]

# Version 0.2.1 (2026-08-27)

- Bump `litebird_sim` requirement to `>=0.17.0` and unpin `astropy` (was hard-pinned to
  `6.1.7`, incompatible with the numpy 2.x that litebird_sim 0.17 pulls in) [#5]
- Fix a pointing bug where both scanning strategies double-counted precession in the
  sidereal-time rotation: `_ct_jd_to_lst_rad` already includes secular precession via its
  own polynomial, but an extra `_equinox_precession_rad` term was being added on top. This
  was off by 19-33 arcmin for 2023-2040 -- over 10x the ~2 arcmin SWIPE beam tolerance the
  code targets. Verified against `pyerfa`'s IAU2006 reference implementation
  (`obl06`/`p06e`/`gmst06`) [#6]
- Fix `SwipeSpinScanningStrategy`'s balloon-trajectory mode, which broke when `balloon_time`
  was passed as a plain list of `Time` objects -- the type its own constructor annotation
  advertised as valid -- rather than a vectorized `Time` array [#6]
- Remove unused dead-code precession helpers (`_chiA_rad`, `_zitaA_rad`, `_zetaA_rad`,
  `_thetaA_rad`) from `common.py`; they had no callers [#6]
- Add a test suite (`tests/`) covering the astronomy formulas in `common.py` (cross-checked
  against `pyerfa`), both scanning strategies, and the raster sawtooth waveform [#6]
- Add CI (GitHub Actions: `ruff` + `pytest` on Python 3.11/3.13) [#6]

Note: this repo has stale `v0.1.0`/`v0.2.0` git tags from December 2021/January 2022,
predating the `pyproject.toml` version scheme this changelog tracks (reset to `0.1.0`
independently of those tags). This release is tagged `v0.2.1` to avoid colliding with the
old `v0.2.0` tag.

# Version 0.1.0

Baseline packaging version prior to this changelog. See the git history for details.
