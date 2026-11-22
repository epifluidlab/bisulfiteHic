# scripts/

Helper scripts for deploying and testing bisulfitehic.

| File | Purpose |
|---|---|
| `quest_install.sh` | One-shot installer for Northwestern's Quest cluster. Backs up any existing install, clones fresh from GitHub, loads Quest modules, runs `install.sh`, and verifies the NOMe-HiC v1.0 patches (XG/XR tags + chimera-aware TLEN reset) are present in the built fat jar. It then runs `test/regression/run_tests.sh` against the built jar (`RUN_TESTS=0` to skip). Config via env vars at the top — see comments in the script. Generated 2026-05-11; updated 2026-11-21 (current Quest Maven module, JDK 8 forced first on `PATH` for the jbwa build, module check before the old install is moved, working XG/XR check, default `TARGET` = `.../software/bisulfiteHic`). Tested end-to-end on Quest on 2026-11-21 by installing into a scratch `TARGET`. |
