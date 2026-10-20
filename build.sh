#!/bin/sh
set -e
cd "$(dirname "$0")"
rm -rf out out-test
mkdir -p out out-test
javac --release 17 -Xlint:all -d out $(find src/main/java -name '*.java')
javac --release 17 -cp out -d out-test $(find src/test/java -name '*.java')
java -cp out:out-test vsr.GaTest
