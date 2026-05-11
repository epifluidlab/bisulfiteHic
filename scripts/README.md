# scripts/

Helper scripts for deploying and testing bisulfitehic.

| File | Purpose |
|---|---|
| `quest_install.sh` | One-shot installer for Northwestern's Quest cluster. Backs up any existing install, clones fresh from GitHub, loads Quest modules, runs `install.sh`, and verifies the NOMe-HiC v1.0 patches (XG/XR tags + chimera-aware TLEN reset) are present in the built fat jar. Config via env vars at the top — see comments in the script. Generated 2026-05-11. |
