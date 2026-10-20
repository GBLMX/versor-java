#!/bin/sh
# Build and run the Java ga-benchmark driver, then validate its JSON.
# Optional: set GA_BENCHMARK_DIR to a ga-benchmark checkout to run gabm.py.
set -eu
cd "$(dirname "$0")"

OUT=${1:-GABenchmark_versor-java.json}
case "$OUT" in
  /*) ;;
  *) OUT="$(pwd)/$OUT" ;;
esac

./bench.sh gabm "$OUT"
python3 - "$OUT" <<'PY'
import json, sys
path = sys.argv[1]
with open(path) as f:
    data = json.load(f)
benchmarks = data.get("benchmarks", [])
if not benchmarks:
    raise SystemExit("no benchmarks in " + path)
required = {"name", "iterations", "cpu_time", "time_unit", "aggregate_name"}
for i, bench in enumerate(benchmarks):
    missing = required - bench.keys()
    if missing:
        raise SystemExit("benchmark %d missing %s" % (i, sorted(missing)))
    if "Solution:versor-java" not in bench["name"]:
        raise SystemExit("benchmark %d has no versor-java solution tag" % i)
print("validated %d benchmarks in %s" % (len(benchmarks), path))
PY

if [ -n "${GA_BENCHMARK_DIR:-}" ]; then
    TMP=$(mktemp -d)
    trap 'rm -rf "$TMP"' EXIT
    cp "$GA_BENCHMARK_DIR/gabm.py" "$TMP/gabm.py"
    python3 - "$TMP/gabm.py" <<'PY'
import pathlib, sys
path = pathlib.Path(sys.argv[1])
src = path.read_text()
needle = "'Versor': {'description': 'Versor'}"
if "versor-java" not in src:
    if needle not in src:
        raise SystemExit("unsupported gabm.py: cannot register versor-java")
    src = src.replace(needle, needle + ",\n        'versor-java': {'description': 'versor-java'}", 1)
    src = src.replace("np.float_", "np.float64")
    path.write_text(src)
PY
    (cd "$TMP" && python3 gabm.py -i "$(dirname "$OUT")" -o "$GA_BENCHMARK_DIR/results")
fi
