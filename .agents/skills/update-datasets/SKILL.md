---
name: update-datasets
description: Build and publish Bactopia's version-pinned datasets to Cloudflare R2. Currently implements the AMRFinder+ database. It verifies the amrfinderplus module is at the latest bioconda version, rebuilds amrfinderplus.tar.gz inside the module's pinned container, and (after confirmation) uploads it via rclone to datasets/v<version>/amrfinderplus.tar.gz. Use when asked to update datasets, rebuild the amrfinderplus database, publish a dataset bundle, refresh the version-pinned datasets, or prepare datasets for a release.
---

# Update Datasets

Build and publish Bactopia's **version-pinned** datasets to Cloudflare R2.

`bactopia datasets` downloads bundles from `https://datasets.bactopia.com/datasets/`.
Most bundles are shared across releases, but a few are pinned to the pipeline
version. The only version-pinned dataset today is **AMRFinder+**
(`conf/params.config`: `amrfinderplus_url = ".../datasets/v${params.bactopia_version}/amrfinderplus.tar.gz"`),
so this skill implements that path concretely; the name stays generic for future
version-pinned datasets.

The publish version is `versions.yml`'s `bactopia` value, and the R2 layout is
`datasets/v<version>/amrfinderplus.tar.gz`. The tool version that **builds** the
database must equal the version the pipeline **pins** (an outdated binary can
produce a database the pinned binary cannot load), so the skill hard-gates on
the `amrfinderplus` module being current before building. This skill is the sole
owner of the build recipe — the old `modules/amrfinderplus/update/` module that
carried it was removed.

## Steps

1. **Currency gate (hard stop).** Confirm the module tool is the latest bioconda
   release before building:
   ```
   bash .agents/skills/update-module/scripts/run-bactopia-update.sh \
       --bactopia-path /home/rpetit3/repos/bactopia/bactopia --module amrfinderplus --json --silent
   ```
   Parse the entry whose `tool == "ncbi-amrfinderplus"`.
   - `needs_update == true` -> **STOP**. Tell the user to run `/update-module`
     (bumps `modules/amrfinderplus/run/module.config`), then re-run this skill.
     Do not build against a stale version.
   - `latest_version == null` (API failure) -> **STOP** and report the failure.
   - `needs_update == false` -> proceed.

2. **Build.** Rebuild the database in the module's pinned container:
   ```
   bash .agents/skills/update-datasets/scripts/build-amrfinderplus-db.sh \
       --bactopia-path /home/rpetit3/repos/bactopia/bactopia
   ```
   Add `--runtime singularity` if the user asks or docker is unavailable. This
   downloads the latest NCBI database and produces the tarball; expect it to take
   several minutes and hundreds of MB. The script reads the container image from
   `modules/amrfinderplus/run/module.config`, so version bumps flow automatically.

3. **Report the build.** Show the script's summary fields (tarball path, sha256,
   tool version, database version). Read the publish version:
   ```
   awk '/^bactopia:/{print $2}' /home/rpetit3/repos/bactopia/bactopia/versions.yml
   ```
   The R2 key is `datasets/v<version>/amrfinderplus.tar.gz`.

4. **Upload (gated behind explicit confirmation).** The rclone destination base
   is `bactopia-r2:bactopia` (remote `bactopia-r2`, bucket `bactopia`);
   `$BACTOPIA_R2_DEST` overrides it if set. The public URL
   `https://datasets.bactopia.com/datasets/...` maps to
   `bactopia-r2:bactopia/datasets/...`, so the full key is
   `<base>/datasets/v<version>/amrfinderplus.tar.gz`.
   - Verify reachability: `rclone lsd bactopia-r2:bactopia`. On failure, report
     it and fall back to printing the manual command below — do not upload.
   - Show the exact command and **ask for confirmation** before running it:
     ```
     rclone copyto "<tarball>" "bactopia-r2:bactopia/datasets/v<version>/amrfinderplus.tar.gz" --s3-no-check-bucket --progress
     ```
   - `--s3-no-check-bucket` is **required**: R2 API tokens cannot `CreateBucket`,
     which rclone otherwise attempts before the first upload and fails with a
     `403 AccessDenied`. The bucket already exists, so skip the check.
   - Only on an explicit yes, run it. `copyto` overwrites an existing key, so
     re-publishing a version is idempotent.

5. **Verify publication.** After upload, confirm the public URL resolves:
   ```
   curl -s -o /dev/null -w '%{http_code}' -I "https://datasets.bactopia.com/datasets/v<version>/amrfinderplus.tar.gz"
   ```
   `200` -> published (this is exactly `/release-checklist` check 15's sentinel).
   Report the HTTP code.

## Notes

- The currency gate reuses `/update-module`'s `run-bactopia-update.sh`; the build
  image is read from `modules/amrfinderplus/run/module.config`, so a version bump
  there is picked up without editing this skill.
- **Generic by design:** to add a future version-pinned dataset, add a sibling
  build script and a step block — the R2 layout, publish-version logic, and
  upload gating are shared.

### Sibling skills

- `/update-module` — the currency gate's remediation (step 1).
- `/release-checklist` — check 15 probes this exact bundle and recommends this
  skill when it 404s.
