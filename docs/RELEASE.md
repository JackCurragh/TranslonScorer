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

  The `make` targets call the `bump2version` **executable**, so they fail with
  "command not found" unless it is on `PATH`. `pip install --user bump2version`
  puts it in the user scripts dir, which often is not. Either add that dir to
  `PATH` or call the module directly — equivalent, and PATH-independent:
  ```bash
  python3 -m bumpversion --dry-run --verbose minor   # check
  python3 -m bumpversion minor                       # commit + tag
  ```

## 2. CI (`.github/workflows/ci.yml`)

Runs on push/PR to `main` and `dev` (matches RiboMetric):

- **test** job — matrix Python 3.11 + 3.12; installs samtools + htslib system
  libs (for pysam / pyBigWig); `pip install -e '.[full]'`; `python -m build`;
  `pytest -n auto`; coverage + Codecov upload on 3.12.
- **type-check** job — `mypy` (`mypy.ini`). Report-only while alpha; drop the
  `|| true` to make it fail-on-error before 1.0 (RiboMetric fails on warnings).

## 3. Release (`.github/workflows/release.yml`)

Triggered by pushing a `v*` tag (which `bump2version` creates):

1. **build** — `python -m build` → sdist + wheel; `twine check`; upload artifact.
2. **publish** — downloads the artifact and publishes with
   `pypa/gh-action-pypi-publish` from the `pypi` GitHub environment using OIDC
   Trusted Publishing (no API token in the repo).

## 4. One-time setup (before the first PyPI release)

Status as of 2026-08-09: **steps 1 and 2 are done**; the first publish has not
run yet.

1. ✅ **PyPI Trusted Publisher (PENDING flow — first release).** Registered
   2026-08-09. Until the first publish creates the project, there is no project
   page to add a publisher on, so this is a *pending* publisher at the account
   level: https://pypi.org/manage/account/publishing/ → "Add a new pending
   publisher", with PyPI project name `TranslonScorer`, owner `JackCurragh`,
   repo `TranslonScorer`, workflow `release.yml`, environment `pypi`. On the
   first successful publish PyPI creates the project and converts this into a
   normal publisher, managed at
   https://pypi.org/manage/project/TranslonScorer/.

   Two fields are easy to get wrong and both fail as an opaque "not authorised"
   at publish time, long after a green build: **"Workflow name" wants the
   filename** (`release.yml`, not the `name:` inside the YAML), and the
   environment name must match `environment: pypi` in `release.yml` exactly.
2. ✅ **GitHub Environment `pypi`** (Settings → Environments). Created. The
   `publish` job declares `environment: pypi`; if it does not exist the job
   fails before it ever contacts PyPI, with an error that does not mention PyPI.
   Optionally add a required reviewer to make each publish pause for approval.
3. **Codecov:** enable the repo (no token needed for public repos).
4. Confirm the repo's default branches are `main`/`dev` (CI triggers).

**TestPyPI rehearsal is *not* a useful pre-flight for this path.** `make
release-test` uploads with a stored TestPyPI token, which exercises neither
OIDC nor the pending-publisher binding — TestPyPI needs its own separate
pending publisher. To catch metadata problems before tagging, run the build
locally instead; that is exactly what the `build` job does:
```bash
python3 -m build && python3 -m twine check dist/*
```

## 5. Release runbook

```bash
# 0. Ensure working tree is clean and on the release branch (main)
git switch main && git pull && git status   # clean

# 1. Update the changelog: move [Unreleased] items under the new version
$EDITOR CHANGELOG.md

# 2. Pre-flight locally
make clean && make test && make lint
make dist && twine check dist/*

# 3. Bump (commits + creates the vX.Y.Z tag) and push
make bump-minor                 # or bump-patch / bump-major
git push origin main --follow-tags

# 4. CI runs on the push; the tag triggers release.yml -> PyPI.
#    Verify: https://pypi.org/project/TranslonScorer/  and the GitHub Release.

# 5. Containers: the existing GHCR workflows
#    (publish-translonscorer.yml / docker-latest.yml) build images on push to
#    main. Tag-pin the image to the release if desired.
```

## 6. Gap closed vs. RiboMetric

| capability | RiboMetric | TranslonScorer (now) |
|---|---|---|
| Semantic version, single bump command | bump2version | ✅ `.bumpversion.cfg`, `make bump-*` |
| CI tests (py3.11/3.12) + coverage | ci.yml + Codecov | ✅ ci.yml |
| Type checking | mypy (fail) | ✅ mypy (report-only → tighten to fail) |
| Tag → PyPI Trusted Publishing | release.yml | ✅ release.yml |
| Changelog / citation | CHANGELOG.md / CITATION.cff | ✅ added |
| Make targets | Makefile | ✅ Makefile |
| Container image | GH Actions | ✅ pre-existing GHCR workflows |

## 7. Remaining to wire up (not code — config/secrets)

- ✅ PyPI Trusted Publisher + `pypi` GitHub environment registered 2026-08-09
  (§4). Not yet exercised — the first tag push is the real test.
- ✅ The workflows do live at the root of `JackCurragh/TranslonScorer` (`origin`
  points there and `.github/workflows/` is at the repo root), so the subtree
  concern is resolved for CI/release purposes.
- **Decide which branch releases are cut from.** §5 below says `main`, but
  `v0.3.0` was tagged on `dev` (which was 14 commits ahead of `main` at the
  time). `release.yml` triggers on any `v*` tag regardless of branch, so both
  work — but the runbook and practice should agree before this becomes a habit.
- Before 1.0: make `mypy` fatal, ensure `pytest` is green in CI, and pin a
  `requires-python` upper bound if needed.
