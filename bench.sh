#!/bin/sh
# Standalone performance harness (pure math core). Optional JOML comparison.
#   ./bench.sh [all|sweep|batch|stress|jit|sparse|render]
#   ./bench.sh jmh <JMH args>
set -e
cd "$(dirname "$0")"

MODE=${1:-all}
[ $# -gt 0 ] && shift
JAVA=${JAVA:-java}
mkdir -p out-bench
JOML_JAR=${JOML_JAR:-$(ls "$HOME"/.local/share/PrismLauncher/libraries/org/joml/joml/*/joml-*.jar \
    "$HOME"/.m2/repository/org/joml/joml/*/joml-*.jar lib/joml*.jar 2>/dev/null | sort -V | tail -1)}

case "$MODE" in
  sparse)
    javac --release 17 -d out-bench $(find src/main/java -name '*.java') \
        src/bench/java/vsr/GaBench.java src/bench/java/vsr/SparseBench.java
    $JAVA ${BENCH_JAVA_OPTS:-} -Xmx6g -cp out-bench vsr.SparseBench
    exit 0 ;;
  render)
    if [ -z "$JOML_JAR" ] || [ ! -f "$JOML_JAR" ]; then
        echo "RenderBench needs JOML (set JOML_JAR)." >&2; exit 1
    fi
    javac --release 17 -Xlint:all,-classfile -cp "$JOML_JAR" -d out-bench \
        $(find src/main/java src/blaze3d/java src/bench/java -name '*.java')
    $JAVA ${BENCH_JAVA_OPTS:-} -Xmx4g -cp "out-bench:$JOML_JAR" vsr.RenderBench "$@"
    exit 0 ;;
  jmh)
    exec ./bench-jmh.sh "$@" ;;
  robust)
    javac --release 17 -d out-bench $(find src/main/java -name '*.java') \
        src/bench/java/vsr/GaBench.java src/bench/java/vsr/SparseRobustBench.java
    $JAVA ${BENCH_JAVA_OPTS:-} -Xmx6g -cp out-bench vsr.SparseRobustBench
    exit 0 ;;
  batchscale)
    javac --release 17 -d out-bench $(find src/main/java -name '*.java') \
        src/bench/java/vsr/GaBench.java src/bench/java/vsr/BatchScalingBench.java
    $JAVA ${BENCH_JAVA_OPTS:-} -Xmx10g -cp out-bench vsr.BatchScalingBench
    exit 0 ;;
  fuzz|gabm)
    if [ -z "$JOML_JAR" ] || [ ! -f "$JOML_JAR" ]; then
        echo "$MODE needs JOML for compile (set JOML_JAR)." >&2; exit 1
    fi
    javac --release 17 -cp "$JOML_JAR" -d out-bench $(find src/main/java -name '*.java') \
        src/bench/java/vsr/GaBench.java src/bench/java/vsr/NumericalFuzz.java \
        src/bench/java/vsr/GaBenchmarkJson.java
    if [ "$MODE" = fuzz ]; then
        $JAVA ${BENCH_JAVA_OPTS:-} -cp "out-bench:$JOML_JAR" vsr.NumericalFuzz
    else
        $JAVA ${BENCH_JAVA_OPTS:-} -cp "out-bench:$JOML_JAR" vsr.GaBenchmarkJson "$@"
    fi
    exit 0 ;;
esac

javac --release 17 -d out-bench $(find src/main/java -name '*.java') src/bench/java/vsr/GaBench.java
$JAVA ${BENCH_JAVA_OPTS:-} -Xmx4g -cp out-bench vsr.GaBench "$MODE"

if [ -n "$JOML_JAR" ] && [ -f "$JOML_JAR" ]; then
    echo "JOML: $JOML_JAR"
    javac --release 17 -Xlint:all,-classfile -cp "$JOML_JAR" -d out-bench \
        $(find src/main/java src/blaze3d/java src/bench/java -name '*.java')
    $JAVA ${BENCH_JAVA_OPTS:-} -Xmx2g -cp "out-bench:$JOML_JAR" vsr.JomlCompare
fi
