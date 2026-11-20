#!/bin/bash
# Install Geant4 11.4.3 + GDCM 2.6.8 + OpenTOPAS v4.3.1 in your HOME, without
# root access and without graphics (Qt), following the official
# OpenTOPAS_quickStart_Debian.md adapted for an HPC cluster.
#
# Usage (two stages, because compute nodes often have no internet):
#   bash scripts/install_opentopas.sh download   # on the LOGIN node (internet)
#   bash scripts/install_opentopas.sh build      # on the login node or in a batch job
#   bash scripts/install_opentopas.sh envfile    # (re)write ~/Applications/topas_env.sh only
#   bash scripts/install_opentopas.sh test       # quick check that topas runs
#
# Before running: load a compiler with C++17 support and CMake >= 3.16, e.g.
#   module load gcc cmake        (names differ per cluster: `module avail gcc cmake`)
# Optional: APPS=/some/other/dir JOBS=8 bash scripts/install_opentopas.sh build
#
# Using a Geant4 that is ALREADY INSTALLED on the cluster:
#   make `geant4-config` available (module load / source .../bin/geant4.sh), then
#   USE_SYSTEM_G4=1 bash scripts/install_opentopas.sh download
#   USE_SYSTEM_G4=1 bash scripts/install_opentopas.sh build
# The OpenTOPAS release matching that Geant4 version is chosen automatically and
# Geant4 itself is neither downloaded nor built.
set -euo pipefail

APPS=${APPS:-$HOME/Applications}
JOBS=${JOBS:-8}
G4V=11.4.3
TOPASV=v4.3.1
G4=$APPS/GEANT4
TOPAS=$APPS/TOPAS
GDCM=$APPS/GDCM
DATASETS="G4NDL.4.7.1 G4EMLOW.8.8 G4PhotonEvaporation.6.1.2 G4RadioactiveDecay.6.1.2
G4PARTICLEXS.4.2 G4ABLA.3.3 G4INCL.1.3 G4ENSDFSTATE.3.0 G4CHANNELING.2.0 G4NUDEXLIB.1.0
G4URRPT.1.1 G4PII.1.3 G4RealSurface.2.2 G4SAIDDATA.2.0 G4TENDL.1.4 LEND_GNDS2.0_ENDF.BVII.1"

stage=${1:-}
say() { echo -e "\n=== $* ==="; }

USE_SYSTEM_G4=${USE_SYSTEM_G4:-0}
if [ "$USE_SYSTEM_G4" = 1 ]; then
  command -v geant4-config >/dev/null || { echo "USE_SYSTEM_G4=1 but geant4-config is not on PATH (load the Geant4 module first)"; exit 1; }
  SYS_G4_VERSION=$(geant4-config --version)
  SYS_G4_PREFIX=$(geant4-config --prefix)
  case "$SYS_G4_VERSION" in
    11.4.*) TOPASV=v4.3.1 ;;
    11.3.*) TOPASV=v4.2.3 ;;
    11.1.*) TOPASV=v4.1.1 ;;
    *) echo "Cluster Geant4 $SYS_G4_VERSION has no matching OpenTOPAS release (need 11.1, 11.3 or 11.4)."
       echo "Run without USE_SYSTEM_G4 to build Geant4 $G4V yourself."; exit 1 ;;
  esac
  # data directory = parent folder of the installed data sets
  SYS_G4_DATA=$(geant4-config --datasets 2>/dev/null | awk '$1=="G4EMLOW"{print $3}' | xargs -r dirname)
  echo "Using cluster Geant4 $SYS_G4_VERSION at $SYS_G4_PREFIX -> OpenTOPAS $TOPASV"
  echo "Geant4 data directory: ${SYS_G4_DATA:-NOT FOUND}"
fi

download() {
  mkdir -p "$TOPAS"
  if [ "$USE_SYSTEM_G4" = 1 ]; then
    download_topas; return
  fi
  mkdir -p "$G4/G4DATA"
  say "Geant4 $G4V source"
  cd "$G4"
  [ -d geant4-v$G4V ] || { wget -nc https://gitlab.cern.ch/geant4/geant4/-/archive/v$G4V/geant4-v$G4V.tar.gz \
                           && tar -zxf geant4-v$G4V.tar.gz; }
  say "Geant4 data sets (several GB, takes a while)"
  cd "$G4/G4DATA"
  for d in $DATASETS; do
    wget -nc -q --show-progress https://cern.ch/geant4-data/datasets/$d.tar.gz
    tar -zxf $d.tar.gz
  done
  download_topas
}

download_topas() {
  say "OpenTOPAS $TOPASV source"
  cd "$TOPAS"
  [ -d OpenTOPAS ] || git clone https://github.com/OpenTOPAS/OpenTOPAS.git
  (cd OpenTOPAS && git fetch --tags -q && git checkout -q $TOPASV)
  say "Download stage finished. Now run:  bash $0 build"
}

libdir() {  # print <prefix>/lib64 if it exists, else <prefix>/lib (never fails)
  if [ -d "$1/lib64" ]; then echo "$1/lib64"; else echo "$1/lib"; fi
}

build() {
  command -v cmake >/dev/null || { echo "cmake not found: module load cmake"; exit 1; }
  echo "cmake: $(cmake --version | head -1)"; echo "c++:   $(c++ --version | head -1)"
  [ -d "$TOPAS/OpenTOPAS" ] || { echo "Run the download stage first"; exit 1; }
  if [ "$USE_SYSTEM_G4" = 1 ]; then
    G4PREFIX=$SYS_G4_PREFIX
    G4DATA=$SYS_G4_DATA
    [ -n "$G4DATA" ] && [ -d "$G4DATA" ] || { echo "Cluster Geant4 data sets not found; ask HPC support where they are"; exit 1; }
  else
    [ -d "$G4/geant4-v$G4V" ] || { echo "Run the download stage first"; exit 1; }
    say "Build Geant4 (30-90 min)"
    cd "$G4"; rm -rf geant4-build geant4-install; mkdir geant4-build geant4-install; cd geant4-build
    cmake ../geant4-v$G4V -DCMAKE_INSTALL_PREFIX=../geant4-install \
          -DGEANT4_INSTALL_DATA=OFF -DGEANT4_BUILD_MULTITHREADED=ON \
          -DGEANT4_USE_QT=OFF -DGEANT4_USE_OPENGL_X11=OFF -DGEANT4_USE_SYSTEM_EXPAT=OFF \
          -DCMAKE_BUILD_TYPE=Release
    make -j"$JOBS" install
    G4PREFIX=$G4/geant4-install
    G4DATA=$G4/G4DATA
  fi

  say "Build GDCM 2.6.8 (bundled with OpenTOPAS)"
  mkdir -p "$GDCM"; cd "$GDCM"
  [ -d gdcm-2.6.8 ] || tar -zxf "$TOPAS/OpenTOPAS/gdcm-2.6.8.tar.gz"
  rm -rf gdcm-build gdcm-install; mkdir gdcm-build gdcm-install; cd gdcm-build
  cmake ../gdcm-2.6.8 -DGDCM_BUILD_SHARED_LIBS=ON -DGDCM_BUILD_DOCBOOK_MANPAGES:BOOL=OFF \
        -DCMAKE_INSTALL_PREFIX=../gdcm-install -DCMAKE_BUILD_TYPE=Release
  make -j"$JOBS" install

  say "Build OpenTOPAS $TOPASV"
  cd "$TOPAS"; rm -rf OpenTOPAS-build OpenTOPAS-install; mkdir OpenTOPAS-build OpenTOPAS-install; cd OpenTOPAS-build
  export Geant4_DIR=$G4PREFIX
  export GDCM_DIR=$GDCM/gdcm-install
  cmake ../OpenTOPAS -DCMAKE_INSTALL_PREFIX=../OpenTOPAS-install -DCMAKE_BUILD_TYPE=Release
  make -j"$JOBS" install
  write_env
}

write_env() {
  # G4PREFIX/G4DATA come from build(); when run on its own ("envfile"), work them out
  if [ -z "${G4PREFIX:-}" ]; then
    if [ "$USE_SYSTEM_G4" = 1 ]; then G4PREFIX=$SYS_G4_PREFIX; G4DATA=$SYS_G4_DATA
    else G4PREFIX=$G4/geant4-install; G4DATA=$G4/G4DATA; fi
  fi
  [ -x "$TOPAS/OpenTOPAS-install/bin/topas" ] || { echo "topas executable not found: build first"; exit 1; }
  say "Environment file: $APPS/topas_env.sh"
  G4LIB=$(libdir "$G4PREFIX")
  TLIB=$(libdir "$TOPAS/OpenTOPAS-install")
  GLIB=$(libdir "$GDCM/gdcm-install")
  cat > "$APPS/topas_env.sh" <<EOF
# source this file (also in SLURM scripts) to use OpenTOPAS
export TOPAS_G4_DATA_DIR=$G4DATA
export LD_LIBRARY_PATH=$TLIB:$G4LIB:$GLIB:\${LD_LIBRARY_PATH:-}
export PATH=$TOPAS/OpenTOPAS-install/bin:\$PATH
EOF
  cat "$APPS/topas_env.sh"
  say "Build finished. Next:  bash $0 test"
}

test_topas() {
  source "$APPS/topas_env.sh"
  command -v topas >/dev/null || { echo "topas not on PATH"; exit 1; }
  d=$(mktemp -d); cd "$d"
  cat > onebox.txt <<'EOF'
s:Ge/World/Material = "G4_AIR"
d:Ge/World/HLX = 50 cm
d:Ge/World/HLY = 50 cm
d:Ge/World/HLZ = 50 cm
s:Ge/Box/Type = "TsBox"
s:Ge/Box/Parent = "World"
s:Ge/Box/Material = "G4_WATER"
d:Ge/Box/HLX = 10 cm
d:Ge/Box/HLY = 10 cm
d:Ge/Box/HLZ = 10 cm
i:Ge/Box/ZBins = 20
sv:Ph/Default/Modules = 1 "g4em-standard_opt4"
s:So/Beam/Type = "Beam"
s:So/Beam/Component = "BeamPosition"
s:So/Beam/BeamParticle = "gamma"
d:So/Beam/BeamEnergy = 6 MeV
s:So/Beam/BeamPositionDistribution = "None"
s:So/Beam/BeamAngularDistribution = "None"
d:Ge/BeamPosition/TransZ = 40 cm
i:So/Beam/NumberOfHistoriesInRun = 10000
s:Sc/Dose/Quantity = "DoseToMedium"
s:Sc/Dose/Component = "Box"
s:Sc/Dose/OutputType = "csv"
s:Sc/Dose/OutputFile = "Dose"
s:Sc/Dose/IfOutputFileAlreadyExists = "Overwrite"
b:Gr/Enable = "False"
b:Ts/PauseBeforeQuit = "False"
EOF
  topas onebox.txt > topas.log 2>&1 || { tail -30 topas.log; echo "TOPAS FAILED, log: $d/topas.log"; exit 1; }
  if [ -s Dose.csv ]; then
    echo "SUCCESS: topas ran 10,000 photons in water. Output: $d/Dose.csv"
    grep -i "elapsed\|cpu" topas.log | tail -2 || true
  else
    echo "topas ran but produced no Dose.csv; see $d/topas.log"; exit 1
  fi
}

case "$stage" in
  download) download ;;
  build) build ;;
  envfile) write_env ;;
  test) test_topas ;;
  *) echo "usage: bash $0 {download|build|envfile|test}"; exit 1 ;;
esac
