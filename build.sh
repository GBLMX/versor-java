#!/bin/sh
# Builds the library, runs the self-contained checks, and writes release jars.
set -e
cd "$(dirname "$0")"

VERSION=1.0.0
rm -rf out out-test out-javadoc dist
mkdir -p out out-test dist

javac --release 17 -Xlint:all -d out $(find src/main/java/vsr -name '*.java')
javac --release 17 -cp out -d out-test $(find src/test/java -name '*.java')
java -cp out:out-test vsr.GaTest

cat > dist/MANIFEST.MF <<MANIFEST
Manifest-Version: 1.0
Automatic-Module-Name: vsr
Implementation-Title: versor-java
Implementation-Version: $VERSION
MANIFEST

mkdir -p out/META-INF
cp LICENSE NOTICE out/META-INF/
jar --create --file "dist/versor-java-$VERSION.jar" --manifest dist/MANIFEST.MF -C out .
jar --create --file "dist/versor-java-$VERSION-sources.jar" -C src/main/java .
javadoc -quiet -Xdoclint:none -d out-javadoc -sourcepath src/main/java vsr
jar --create --file "dist/versor-java-$VERSION-javadoc.jar" -C out-javadoc .
rm dist/MANIFEST.MF

ls -lh dist
