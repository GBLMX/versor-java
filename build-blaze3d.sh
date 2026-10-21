#!/bin/sh
# Builds the optional Blaze3D (JOML) adapter. The core Gradle build stays
# dependency-free.
set -e
cd "$(dirname "$0")"

if [ -z "$JOML_JAR" ]; then
    JOML_JAR=$(ls "$HOME"/.local/share/PrismLauncher/libraries/org/joml/joml/*/joml-*.jar \
                  "$HOME"/.m2/repository/org/joml/joml/*/joml-*.jar \
                  lib/joml*.jar 2>/dev/null | sort -V | tail -1)
fi
if [ -z "$JOML_JAR" ] || [ ! -f "$JOML_JAR" ]; then
    echo "JOML jar not found. Set JOML_JAR=/path/to/joml.jar" >&2
    exit 1
fi
echo "JOML: $JOML_JAR"

rm -rf out-blaze3d
mkdir -p out-blaze3d
javac --release 17 -Xlint:all,-classfile -cp "$JOML_JAR" -d out-blaze3d \
    $(find src/main/java src/blaze3d/java -name '*.java')
javac --release 17 -Xlint:all,-classfile -cp "out-blaze3d:$JOML_JAR" -d out-blaze3d \
    $(find src/blaze3dTest/java -name '*.java')
java -cp "out-blaze3d:$JOML_JAR" vsr.blaze3d.Blaze3DTest
