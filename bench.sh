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
    $JAVA ${BENCH_JAVA_OPTS:-} -Xmx4g -cp "out-bench:$JOML_JAR" vsr.RenderBench
    exit 0 ;;
  jmh)
    exec ./bench-jmh.sh "$@" ;;
esac

javac --release 17 -d out-bench $(find src/main/java -name '*.java') src/bench/java/vsr/GaBench.java
$JAVA ${BENCH_JAVA_OPTS:-} -Xmx4g -cp out-bench vsr.GaBench "$MODE"

if [ -n "$JOML_JAR" ] && [ -f "$JOML_JAR" ]; then
    echo "JOML: $JOML_JAR"
    javac --release 17 -Xlint:all,-classfile -cp "$JOML_JAR" -d out-bench \
        $(find src/main/java src/blaze3d/java src/bench/java -name '*.java')
    $JAVA ${BENCH_JAVA_OPTS:-} -Xmx2g -cp "out-bench:$JOML_JAR" vsr.JomlCompare
fi
