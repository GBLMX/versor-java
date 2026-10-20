#!/bin/sh
# JMH benchmark runner. Downloads JMH 1.37 jars into bench-jmh/lib (gitignored).
set -e
cd "$(dirname "$0")"
mkdir -p bench-jmh/lib out-jmh
fetch() {
    name=$1; path=$2
    [ -f "bench-jmh/lib/$name" ] || curl -fsSL -o "bench-jmh/lib/$name" \
        "https://repo1.maven.org/maven2/$path/$name"
}
fetch jmh-core-1.37.jar org/openjdk/jmh/jmh-core/1.37
fetch jmh-generator-annprocess-1.37.jar org/openjdk/jmh/jmh-generator-annprocess/1.37
fetch jopt-simple-5.0.4.jar net/sf/jopt-simple/jopt-simple/5.0.4
fetch commons-math3-3.6.1.jar org/apache/commons/commons-math3/3.6.1

CP="bench-jmh/lib/jmh-core-1.37.jar:bench-jmh/lib/jopt-simple-5.0.4.jar:bench-jmh/lib/commons-math3-3.6.1.jar"
javac --release 17 -cp "$CP" -processorpath "$CP:bench-jmh/lib/jmh-generator-annprocess-1.37.jar" \
    -d out-jmh $(find src/main/java src/jmh/java -name '*.java')
${JAVA:-java} ${BENCH_JAVA_OPTS:-} -Xmx4g -cp "out-jmh:$CP" org.openjdk.jmh.Main "$@"
