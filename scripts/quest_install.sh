#!/bin/bash
# quest_install.sh
#
# One-shot installer for bisulfitehic on Northwestern's Quest cluster.
#
# What it does:
#   1. Backs up any existing install at $TARGET (rename to ${TARGET}_old_<date>)
#   2. Clones https://github.com/epifluidlab/bisulfiteHic.git to $TARGET
#   3. Loads Quest modules and exports JAVA_HOME
#   4. Runs the repo's install.sh (jbwa native build + mvn assembly:single)
#   5. Verifies the NOMe-HiC v1.0 patches (XG/XR tags, chimera-aware TLEN reset)
#      are present in the built fat jar
#   6. Prints a summary + rollback hint
#
# Usage on Quest:
#   scp this file from your laptop to Quest:
#     scp scripts/quest_install.sh quest:~/
#   then on a Quest login node (NOT a compute node — clone + maven need internet):
#     bash ~/quest_install.sh
#
# To re-run after the first install, you must `rm -rf $TARGET` or pick a new
# $TARGET. The script refuses to overwrite an existing $TARGET to avoid clobber.

set -euo pipefail

# ─── CONFIG ──────────────────────────────────────────────────────────────────
#
# Adjust paths and module names for your environment. Defaults match the
# Yoshii lab layout on Quest.

TARGET="${TARGET:-/gpfs/projects/b1198/epifluidlab/yoshii/software/bisulfitehic}"
REPO_URL="${REPO_URL:-https://github.com/epifluidlab/bisulfiteHic.git}"
# Pin to a specific commit or branch (default: master HEAD)
REPO_REF="${REPO_REF:-master}"

# Quest module names — adjust if Quest's catalogue changes
MOD_JAVA8="${MOD_JAVA8:-java/jdk1.8.0_191}"
MOD_MAVEN="${MOD_MAVEN:-maven/3.6.3}"
MOD_GCC="${MOD_GCC:-gcc/4.8.3}"
MOD_SAMTOOLS="${MOD_SAMTOOLS:-samtools/1.10.1}"

# Whether to also write a versioned copy of the built jar (NOMe-HiC v1.0).
# Set TAG_JAR=0 to skip.
TAG_JAR="${TAG_JAR:-1}"
JAR_VERSION_TAG="${JAR_VERSION_TAG:-NOMeHiC_v1.0}"

# ─── HELPERS ─────────────────────────────────────────────────────────────────

log()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
fail() { echo "[ERROR] $*" >&2; exit 1; }

# ─── PRE-FLIGHT ──────────────────────────────────────────────────────────────

log "Pre-flight checks..."

PARENT="$(dirname "$TARGET")"
[ -d "$PARENT" ]    || fail "Parent dir does not exist: $PARENT"
[ -w "$PARENT" ]    || fail "Parent dir not writable: $PARENT (check group / fs permissions)"

# Don't clobber a fresh-but-incomplete install.
if [ -d "$TARGET" ] && [ ! -f "$TARGET/install.sh" ]; then
    fail "$TARGET exists but doesn't look like a bisulfitehic checkout (no install.sh). Refusing to touch it. Inspect manually."
fi

# Backup naming
BACKUP_SUFFIX="_old_$(date +%Y%m%d)"
BACKUP_PATH="${TARGET}${BACKUP_SUFFIX}"
if [ -d "$BACKUP_PATH" ]; then
    fail "Backup target already exists from earlier today: $BACKUP_PATH. Remove it or rename it, then re-run."
fi

# Check git / scp inputs
command -v git >/dev/null    || fail "git not on PATH"

# ─── 1. BACKUP OLD INSTALL ───────────────────────────────────────────────────

if [ -d "$TARGET" ]; then
    log "Backing up existing install: $TARGET → $BACKUP_PATH"
    mv "$TARGET" "$BACKUP_PATH"
else
    log "No existing install at $TARGET — fresh install."
fi

# ─── 2. CLONE FROM GITHUB ────────────────────────────────────────────────────

log "Cloning $REPO_URL @ $REPO_REF → $TARGET"
git clone "$REPO_URL" "$TARGET"
cd "$TARGET"
if [ "$REPO_REF" != "master" ]; then
    git checkout "$REPO_REF"
fi
CURRENT_REF=$(git rev-parse --short HEAD)
log "Cloned. HEAD = $CURRENT_REF"

# ─── 3. LOAD MODULES + JAVA_HOME ─────────────────────────────────────────────

log "Loading modules..."
# `module` is a shell function on Quest, not an executable — source the LMOD
# initializer if `module` isn't already in this shell.
if ! command -v module >/dev/null 2>&1; then
    if [ -f /etc/profile.d/modules.sh ]; then
        # shellcheck disable=SC1091
        . /etc/profile.d/modules.sh
    fi
fi
command -v module >/dev/null 2>&1 || fail "module command not available. Load LMOD manually and re-run."

module purge       || true
module load "$MOD_JAVA8"
module load "$MOD_MAVEN"
module load "$MOD_GCC"
module load "$MOD_SAMTOOLS" || log "(samtools module load failed — only needed at runtime, not build. Continuing.)"

# jbwa's Makefile needs JAVA_HOME pointing at a JDK with jni.h
JAVA_BIN=$(command -v java)
[ -n "$JAVA_BIN" ] || fail "java not on PATH after module load"
export JAVA_HOME=$(dirname "$(dirname "$(readlink -f "$JAVA_BIN")")")
[ -f "$JAVA_HOME/include/jni.h" ] || fail "JAVA_HOME=$JAVA_HOME has no include/jni.h — wrong JDK module."
log "JAVA_HOME=$JAVA_HOME"
java -version 2>&1 | sed 's/^/  /'

# ─── 4. BUILD ────────────────────────────────────────────────────────────────

log "Running install.sh (jbwa make + mvn assembly:single)..."
bash install.sh

JAR="$TARGET/target/bisulfitehic-0.38-jar-with-dependencies.jar"
[ -f "$JAR" ] || fail "Build finished but expected jar not found: $JAR"
log "Built jar: $JAR ($(du -h "$JAR" | awk '{print $1}'))"

# ─── 5. VERIFY PATCHES ───────────────────────────────────────────────────────

log "Verifying NOMe-HiC v1.0 patches are baked into the jar..."

# 5a. isLikelyChimericLigation method present in Bhmem bytecode
CHIM=$(javap -p -cp "$JAR" main.java.edu.mit.compbio.bisulfitehic.mapping.Bhmem 2>/dev/null \
       | grep -c 'isLikelyChimericLigation' || true)
if [ "$CHIM" -lt 1 ]; then
    fail "Bhmem.class is missing isLikelyChimericLigation — chimera-aware TLEN reset patch NOT in this build."
fi
log "  ✓ isLikelyChimericLigation present ($CHIM bytecode refs)"

# 5b. XG/XR setAttribute strings present in Bhmem class constant pool
unzip -p "$JAR" main/java/edu/mit/compbio/bisulfitehic/mapping/Bhmem.class \
    | strings | grep -qE '^XG$' \
    || fail "Bhmem.class is missing XG constant — XG/XR tag patch NOT in this build."
unzip -p "$JAR" main/java/edu/mit/compbio/bisulfitehic/mapping/Bhmem.class \
    | strings | grep -qE '^XR$' \
    || fail "Bhmem.class is missing XR constant — XG/XR tag patch NOT in this build."
log "  ✓ XG and XR constants present"

# 5c. ccInference slim jar present and referenced by pom.xml
SLIM="$TARGET/lib/ccInference-0.19-slim.jar"
[ -f "$SLIM" ] || fail "Slim ccInference jar missing: $SLIM"
SLIM_SIZE=$(du -k "$SLIM" | awk '{print $1}')
[ "$SLIM_SIZE" -lt 1024 ] || fail "ccInference jar at $SLIM is unexpectedly large (${SLIM_SIZE} KB) — expected <1 MB slim jar."
log "  ✓ ccInference slim jar present (${SLIM_SIZE} KB)"

# ─── 6. OPTIONAL — TAG THE JAR WITH A VERSION LABEL ──────────────────────────

if [ "$TAG_JAR" = "1" ]; then
    TAGGED="$TARGET/target/bisulfitehic-0.38-${JAR_VERSION_TAG}-jar-with-dependencies.jar"
    cp "$JAR" "$TAGGED"
    log "Tagged jar: $TAGGED"
fi

# ─── 7. SUMMARY ──────────────────────────────────────────────────────────────

cat <<EOF

═══════════════════════════════════════════════════════════════════════════════
INSTALL SUCCESS

  Install path : $TARGET
  Commit       : $CURRENT_REF
  Built jar    : $JAR
$([ "$TAG_JAR" = "1" ] && echo "  Tagged jar   : $TAGGED")
  Backup of previous install: ${BACKUP_PATH}$([ -d "$BACKUP_PATH" ] || echo " (none — fresh install)")

NEXT STEPS

  1. Update your BISHIC env var (install.sh already appended to ~/.bash_profile,
     but if your current shell still points at the old install path, open a
     fresh shell or run:  export BISHIC=$TARGET

  2. To roll back to the previous install (if anything looks wrong):
       rm -rf $TARGET && mv $BACKUP_PATH $TARGET

  3. Run a test alignment against a known-good fastq pair before deleting the
     backup. See scripts/cluster_run_test.sh in bissnp_upgrade for the test
     pipeline that compares against the Fix-D baseline.

═══════════════════════════════════════════════════════════════════════════════
EOF
