# TranslonScorer release process

Operationally mirrors RiboMetric: single-source semantic version bumped by
`bump2version` (commit + tag), CI on every push/PR, and tag-triggered PyPI
publication via **Trusted Publishing** (OIDC, no stored token).

---

## 1. Versioning model

- **Scheme:** Semantic Versioning `MAJOR.MINOR.PATCH`; tags are `vX.Y.Z`.
- **Sources of truth:** `TranslonScorer/__init__.py` (`__version__`) and the
  `[project]` `version` in `pyproject.toml`. They are kept in lockstep by
  `bump2version` (`.bumpversion.cfg`). (RiboMetric keeps version only in
  `__init__.py`; TranslonScorer's `__init__` imports the CLI, so attr-based
  dynamic versioning is unsafe — hence both files are bumped together.)
- **Bumping** (commits and tags automatically):
  ```bash
  make bump-patch      # 0.1.1 -> 0.1.2   (bugfixes)
  make bump-minor      # 0.1.1 -> 0.2.0   (features, backwards compatible)
  make bump-major      # 0.1.1 -> 1.0.0   (breaking)
  ```
  Verify first with `bump2version --dry-run --verbose patch`.

  Every `make` target runs through one interpreter, `$(PYTHON)` (default
  `python3`), and invokes tools as `$(PYTHON) -m <tool>` rather than as bare
  executables — a `pip install --user` puts them in a scripts dir that is
  usually not on `PATH`, which is what used to make these fail with "command
  not found". Point it at a venv when you have one:
  ```bash
  make bump-minor PYTHON=.venv/bin/python
  python3 -m bumpversion --dry-run --verbose minor    # check before bumping
  ```

## 2. CI — the gate (`.github/workflows/ci.yml`)

There is **one** definition of "does this pass". It runs on push/PR to `main`
and `dev`, and `release.yml` *calls* it (`workflow_call`) so a tagged release
re-runs the identical gate instead of a drifting copy. Adding a job here
automatically enforces it at release time.

- **lint** — `ruff check` + `black --check`. Fatal. (Previously defined in the
  `Makefile` but never run by CI, so it sat red unnoticed.)
- **typecheck** — `mypy` (`mypy.ini`). Report-only while alpha; drop the
  `|| true` to make it fail-on-error before 1.0 (RiboMetric fails on warnings).
- **test** — matrix Python 3.11 + 3.12; installs samtools + htslib system libs
  (for pysam / pyBigWig); `pip install -e '.[full]'`; `pytest -n auto`;
  coverage + Codecov upload on 3.12.
- **build** — `python -m build` → sdist + wheel, `twine check`, **and an
  assertion that both filenames start with `translonscorer-`**. That last check
  exists because it is precisely what broke every release before v0.3.1 (§8):
  `twine check` passes a non-normalised wheel, and PyPI rejects it with a 400
  only at upload, long after a green build. The artifact is uploaded as `dist`
  and is what `release.yml` publishes — the bits that ship are the bits CI
  checked.
- **image** — on pushes to `main` only, gated on lint+test: builds and pushes
  `ghcr.io/<owner>/translonscorer:main` and `:sha-<short>`.

Locally, `make preflight` mirrors the **build** job exactly (build + twine
check + filename assertion).

## 3. Release (`.github/workflows/release.yml`)

The single release path. Pushing a `vX.Y.Z` tag (which `bump2version` creates)
runs, each step gated on the previous:

1. **gate** — calls `ci.yml` in full, on the tag.
2. **pypi** — downloads the `dist` artifact the gate built and publishes with
   `pypa/gh-action-pypi-publish` from the `pypi` GitHub environment using OIDC
   Trusted Publishing (no API token in the repo).
3. **image** — pushes `ghcr.io/<owner>/translonscorer:vX.Y.Z`, `:X.Y` and
   `:latest`, then smoke-tests the pushed digest with `--help`.
4. **release** — creates a GitHub Release with the distributions attached and
   generated notes.

> **Gotcha — reusable-workflow permissions.** A called workflow's jobs may not
> request more permission than the *calling* job grants, and GitHub validates
> that when the run starts, **before any job-level `if` is evaluated**. `ci.yml`'s
> `image` job declares `packages: write`, so the `gate` job here must grant it
> too — even though that job is gated to main-branch pushes and can never run on
> a tag. Getting this wrong fails the entire run as `startup_failure` with zero
> jobs and the message "This run likely failed because of a workflow file issue",
> which names neither the job nor the permission. This bit v0.3.1's first tag.

### Container tag meanings

| tag | written by | means |
|---|---|---|
| `:latest` | `release.yml` | the most recent release |
| `:vX.Y.Z`, `:X.Y` | `release.yml` | that release |
| `:main`, `:sha-<short>` | `ci.yml` | tip of `main`, lint+test green |

`:latest` used to mean "tip of main" — two separate workflows
(`docker-latest.yml`, `publish-translonscorer.yml`) both pushed it on every
main push, one of them with **no test gate**, racing each other for the same
tag. Both are deleted; nothing but `release.yml` writes `:latest` now.

**Consumers to be aware of:** `ensembl-genes-nf/pipelines/riboseq/modules/translonscorer.nf`
and its `-feature` twin default to `ghcr.io/jackcurragh/translonscorer:latest`.
Under the new meaning they follow releases rather than every main commit —
which is what you want for reproducible annotation runs. Pin
`params.translonscorer_container` to `:vX.Y.Z` for a fixed run, or to `:main`
to keep the old always-newest behaviour.

## 4. One-time setup — **done, and proven working**

1. ✅ **PyPI Trusted Publisher (pending flow).** Registered 2026-08-09 at the
   account level (https://pypi.org/manage/account/publishing/ → "Add a new
   pending publisher"): project `translonscorer`, owner `JackCurragh`, repo
   `TranslonScorer`, workflow `release.yml`, environment `pypi`. On the first
   successful publish PyPI creates the project and converts this into a normal
   publisher at https://pypi.org/manage/project/translonscorer/.

   Two fields are easy to get wrong and both fail as an opaque "not authorised":
   **"Workflow name" wants the filename** (`release.yml`, not the `name:` inside
   the YAML), and the environment name must match `environment: pypi` exactly.
2. ✅ **GitHub Environment `pypi`** (Settings → Environments). Created.
   Optionally add a required reviewer to make each publish pause for approval.
3. **Codecov:** enable the repo (no token needed for public repos).
4. Confirm the repo's default branches are `main`/`dev` (CI triggers).

**Steps 1 and 2 are not merely configured — they are verified.** The v0.3.0 run
(32828505888) got as far as `Uploading distributions to
https://upload.pypi.org/legacy/` and failed with a **400 on the filename**, not
a 403 on auth. The OIDC token exchange succeeded. The only thing that has ever
blocked a publish is §8.

**TestPyPI rehearsal is not a useful pre-flight for this path** and there is
deliberately no target for it: it uses a stored token, exercising neither OIDC
nor the pending-publisher binding (TestPyPI needs its own separate publisher).
`make preflight` is the real rehearsal — it reproduces CI's `build` job.

## 5. Release runbook

Releases are cut from **`main`**. (`v0.3.0` was tagged on `dev`; that is no
longer the practice. `release.yml` fires on any `v*` tag regardless of branch,
so the discipline is the runbook's, not the workflow's.)

```bash
# 0. Clean tree, on main, level with origin
git switch main && git pull && git status

# 1. Move [Unreleased] items under the new version
$EDITOR CHANGELOG.md

# 2. Pre-flight — the same checks CI will run
make lint test preflight PYTHON=.venv/bin/python

# 3. Bump (commits + creates the vX.Y.Z tag) and push
make bump-minor PYTHON=.venv/bin/python      # or bump-patch / bump-major
git push origin main --follow-tags

# 4. Fast-forward dev so the branches do not diverge
git branch -f dev main && git push origin dev

# 5. The tag triggers release.yml: gate -> PyPI -> GHCR -> GitHub Release.
#    Verify: https://pypi.org/project/translonscorer/
#            ghcr.io/jackcurragh/translonscorer:vX.Y.Z and :latest
#            the GitHub Release page
```

## 6. Gap closed vs. RiboMetric

| capability | RiboMetric | TranslonScorer (now) |
|---|---|---|
| Semantic version, single bump command | bump2version | ✅ `.bumpversion.cfg`, `make bump-*` |
| CI tests (py3.11/3.12) + coverage | ci.yml + Codecov | ✅ ci.yml |
| Lint enforced in CI | ruff/black | ✅ ci.yml `lint` job (fatal) |
| Type checking | mypy (fail) | ✅ mypy (report-only → tighten to fail) |
| Tag → PyPI Trusted Publishing | release.yml | ✅ release.yml |
| Changelog / citation | CHANGELOG.md / CITATION.cff | ✅ added |
| Make targets | Makefile | ✅ Makefile |
| Container image | GH Actions | ✅ ci.yml (`:main`) + release.yml (`:latest`) |
| GitHub Release with artifacts | — | ✅ release.yml |

## 7. Remaining before 1.0

- Make `mypy` fatal (drop `|| true` in `ci.yml`). 229 errors outstanding,
  overwhelmingly `X | None` passed into non-optional parameters.
- Un-skip the contract tests in `tests/test_alignments.py` and
  `tests/test_counts.py`. Both carry a module-level `skipif` on the local
  `data/global_partitioned` cohort, so all 29 skip in CI — including the ones
  their own comments call "contract tests (fast, no I/O beyond schema
  inspection)", which need no cohort data. The provenance backbone (VersionKey,
  schemas) currently has zero CI coverage.
- Decide whether `requires-python` needs an upper bound.

## 8. Post-mortem: why v0.2.0 and v0.3.0 never published

Three release runs, three failures, no package on PyPI:

| run | tag | date | outcome |
|---|---|---|---|
| 28224587497 | v0.2.0 | 2026-06-26 | failed |
| 31315960329 | v0.3.0 | 2026-08-09 | failed |
| 32828505888 | v0.3.0 (re-run) | 2026-08-25 | failed |

Cause, from the last run's log:

```
400 Filename 'TranslonScorer-0.3.0-py3-none-any.whl' should
contain the normalized project name 'translonscorer', not 'TranslonScorer'.
```

`pyproject.toml` had `name = "TranslonScorer"`. The **sdist** was already
normalised (`translonscorer-0.3.0.tar.gz`), so only the wheel tripped it — and
`twine check` passes a non-normalised wheel, so the build looked green right up
to the upload. Fixed on main in `88adeb8` (`name = "translonscorer"`).

Two lasting consequences:

- **Never re-run `release.yml` on the `v0.3.0` tag.** That tag's tree still
  carries the capitalised name; it will fail identically forever. The fix ships
  in the next tag cut from `main`.
- The `build` job and `make preflight` now assert the filename prefix, so this
  class of failure is caught in CI rather than at upload.
