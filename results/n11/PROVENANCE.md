# n = 11 run (archived from version 1)

Produced by `verify_m11.py` (version 1 of this repository) on a dual-core Intel Core i5
(1.6 GHz, 8 GB RAM) with 8 worker processes: 268,715 s wall-clock, 961,275 s total worker CPU.

* `certificates/subtree_K.json` — per-root-box certificates (same format as version 2).
* `datasets/surviving_boxes.json` — all 7,385 surviving boxes as [lower corner, upper corner].
* `datasets/execution_telemetry.json`, `execution.log` — the original summary and log.

Audit with `python -m woods_verify summarize --n 11`. Re-run any subtree with version 2 and
compare with `python -m woods_verify compare` — the certificates must agree exactly.
