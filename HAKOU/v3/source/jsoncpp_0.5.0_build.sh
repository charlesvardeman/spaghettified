# http://sourceforge.net/projects/jsoncpp/files/jsoncpp/0.5.0/jsoncpp-src-0.5.0.tar.gz/download

set -e

tar -xzf jsoncpp-src-0.5.0.tar.gz
cd jsoncpp-src-0.5.0
scons platform=linux-gcc

GCCVER=$(gcc -dumpversion)

cp libs/linux-gcc-${GCCVER}/libjson_linux-gcc-${GCCVER}_libmt.so ../libs/libjson.so
cp libs/linux-gcc-${GCCVER}/libjson_linux-gcc-${GCCVER}_libmt.a ../libs/libjson.a
cp -r include/json ..

cd ..
rm -r jsoncpp-src-0.5.0