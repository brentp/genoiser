set -euo pipefail

NIM_VERSION=devel
NIMBLE_VERSION=devel

base=$(pwd)

git clone -b $NIM_VERSION --depth 1 git://github.com/nim-lang/nim nim-$NIM_VERSION/
cd nim-$NIM_VERSION
git clone -b master --depth 1 git://github.com/nim-lang/csources csources/

cd csources
sh build.sh
cd ..
rm -rf csources
bin/nim c koch
./koch boot -d:release

export PATH=$PATH:$base/nim-$NIM_VERSION/bin/:$PATH:$base/nimble/src

cd $base

git clone -b master --depth 1 git://github.com/nim-lang/nimble.git
cd nimble
nim c src/nimble
cp ./src/nimble /usr/bin/

cd $base
git clone --depth 1 git://github.com/brentp/hts-nim.git
cd hts-nim
grep -v requires hts.nimble > k.nimble && mv k.nimble hts.nimble
nimble install -y

cd $base
git clone --depth 1 git://github.com/brentp/kexpr-nim.git
cd kexpr-nim
nimble install -y

cd $base
git clone --depth 1 git://github.com/docopt/docopt.nim.git
cd docopt.nim
nimble install -y

cd $base
git clone --depth 1 git://github.com/brentp/genoiser.git
cd genoiser
nim c -d:release --passC:-flto --passL:-s --threads:on src/genoiser.nim
cp ./src/genoiser /io
