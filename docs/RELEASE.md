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

## 2. CI (`.github/workflows/ci.yml`)

Runs on push/PR to `main` and `dev` (matches RiboMetric):

- **test** job — matrix Python 3.10 + 3.12; installs samtools + htslib system
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

1. **Reserve the name / first upload via TestPyPI** to validate metadata:
   `make release-test` (uploads to TestPyPI; needs a TestPyPI token locally).
2. **PyPI Trusted Publisher (PENDING flow — first release).** The project is
   NOT yet on PyPI, so its project page does not exist and you cannot add a
   publisher there. Instead register a *pending* publisher at the account level:
   https://pypi.org/manage/account/publishing/ → "Add a new pending publisher":
   PyPI project name `TranslonScorer`, owner `JackCurragh`, repo `TranslonScorer`,
   workflow `release.yml`, environment `pypi`. On the first successful publish
   PyPI creates the project and binds it to this publisher. (After that, manage
   it at https://pypi.org/manage/project/TranslonScorer/ as usual.)
3. **GitHub:** create an Environment named `pypi` (Settings → Environments);
   optionally require a reviewer for publish.
4. **Codecov:** enable the repo (no token needed for public repos).
5. Confirm the repo's default branches are `main`/`dev` (CI triggers).

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
| CI tests (py3.10/3.12) + coverage | ci.yml + Codecov | ✅ ci.yml |
| Type checking | mypy (fail) | ✅ mypy (report-only → tighten to fail) |
| Tag → PyPI Trusted Publishing | release.yml | ✅ release.yml |
| Changelog / citation | CHANGELOG.md / CITATION.cff | ✅ added |
| Make targets | Makefile | ✅ Makefile |
| Container image | GH Actions | ✅ pre-existing GHCR workflows |

## 7. Remaining to wire up (not code — config/secrets)

- Register the PyPI Trusted Publisher + `pypi` GitHub environment (§4).
- Decide the public Git history: TranslonScorer currently lives as a subtree of
  `all-RiboSeq`. The CI/release workflows assume it is its own repo at
  `JackCurragh/TranslonScorer` (per `pyproject` URLs); confirm that repo exists
  and these `.github/workflows/` live at its root.
- Before 1.0: make `mypy` fatal, ensure `pytest` is green in CI, and pin a
  `requires-python` upper bound if needed.
