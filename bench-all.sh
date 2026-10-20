#!/bin/sh
# Full test matrix for the current tree. Writes human-readable logs to stdout.
set -e
cd "$(dirname "$0")"
export JAVA=${JAVA:-$HOME/.jdks/jbrsdk-25.0.3/bin/java}
echo "== build/tests =="
sh build.sh
echo "== sweep/batch/stress/jit/joml =="
./bench.sh all
echo "== sparse =="
./bench.sh sparse
echo "== support =="
./bench.sh support
echo "== robust =="
./bench.sh robust
echo "== render =="
./bench.sh render 1000000 1000
echo "== batchscale =="
./bench.sh batchscale
echo "== fuzz =="
./bench.sh fuzz
echo "== multicore =="
./bench.sh multicore
echo "== jmh dense =="
./bench.sh jmh JmhDense -prof gc
echo "== jmh fast 1t =="
./bench.sh jmh JmhFast -t 1
echo "== jmh fast 4t =="
./bench.sh jmh JmhFast -t 4
echo "== jit c2 =="
./bench.sh jit > /tmp/jit-c2.log
echo "== jit c1 =="
BENCH_JAVA_OPTS="-XX:TieredStopAtLevel=1" ./bench.sh jit > /tmp/jit-c1.log
echo "== ga-benchmark json =="
./bench.sh gabm /tmp/GABenchmark_versor-java.json
