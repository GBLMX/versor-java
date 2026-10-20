#!/bin/sh
# Standalone performance harness (pure math core). Optional JOML comparison.
set -e
cd "$(dirname "$0")"

MODE=${1:-all}
JAVA=${JAVA:-java}
rm -rf out-bench
mkdir -p out-bench
javac --release 17 -d out-bench $(find src/main/java -name '*.java') src/bench/java/vsr/GaBench.java
$JAVA ${BENCH_JAVA_OPTS:-} -Xmx4g -cp out-bench vsr.GaBench "$MODE"

JOML_JAR=${JOML_JAR:-$(ls "$HOME"/.local/share/PrismLauncher/libraries/org/joml/joml/*/joml-*.jar \
    "$HOME"/.m2/repository/org/joml/joml/*/joml-*.jar lib/joml*.jar 2>/dev/null | sort -V | tail -1)}
if [ -n "$JOML_JAR" ] && [ -f "$JOML_JAR" ]; then
    echo "JOML: $JOML_JAR"
    javac --release 17 -Xlint:all,-classfile -cp "$JOML_JAR" -d out-bench \
        $(find src/main/java src/blaze3d/java src/bench/java -name '*.java')
    $JAVA ${BENCH_JAVA_OPTS:-} -Xmx2g -cp "out-bench:$JOML_JAR" vsr.JomlCompare
fi
