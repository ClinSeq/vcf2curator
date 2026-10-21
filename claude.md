## VCF2CURATOR 

A Python CLI utility that transforms raw bioinformatics pipeline results (VCFs) into clean, structured tables, ready for immediate ingestion into the Curator web app.

## Environment & dependency management — uv only

### Do:

- uv init to scaffold a new project.
- uv add <package> to add a runtime dependency.
- uv add --dev <package> to add a dev-only dependency (test/lint tools).
- uv remove <package> to remove a dependency.
- uv run <command> to run anything inside the project's environment (scripts, pytest, ruff, the CLI itself) — this guarantees the lockfile's resolved versions are used.
- uv sync to install/update the environment from uv.lock.
- uv lock after manually editing pyproject.toml dependency constraints.
- Commit both pyproject.toml and uv.lock to version control.
- Pin the Python version in .python-version and requires-python in pyproject.toml.
- Use uv tool install . (or uv pip install -e .) to install the CLI itself as an executable entry point during development.

### Don't:

- Don't create or edit a requirements.txt as the source of truth for dependencies — pyproject.toml + uv.lock is the source of truth.
- Don't manually create a venv with python -m venv — let uv manage it (.venv/ is created and managed automatically).
- Don't install packages globally or with sudo pip install.
- Don't mix package managers in one repo (no poetry.lock alongside uv.lock, etc.).
- Don't hand-edit uv.lock.

## Project layout

The package lives at the repo root (not under src/); pyproject.toml sets
`[tool.uv.build-backend] module-root = ""` accordingly:

```
vcf2curator/              # repo root
├── pyproject.toml
├── uv.lock
├── .python-version
├── README.md
├── CLAUDE.md
├── .gitignore
├── vcf2curator/          # the package
│   ├── __init__.py
│   └── ...
└── tests/
    └── ...
```

Keep the CLI simple: plain subcommands with explicit arguments, one module per
command under vcf2curator/commands/, registered in cli.py. No plugin/format
registry or `list-formats` style indirection.

## Code style & static analysis
- Formatting & linting: ruff, added as a dev dependency (uv add --dev ruff). Run uv run ruff format . and uv run ruff check . before committing. Fix lint issues rather than suppressing them; use # noqa sparingly and with a reason comment.
- Type hints are required on all function signatures (parameters and return types) in library code. Use from __future__ import annotations where useful for forward references, or built-in generics (list[str], dict[str, int]) directly on Python 3.9+.
- Consider mypy or pyright as a dev dependency for a project of any significant size; run it in CI once configured.
- Follow PEP 8 naming: snake_case for functions/variables, PascalCase for classes, UPPER_SNAKE_CASE for constants.
- Keep functions small and single-purpose. If a function needs a comment explaining "and then it also does X", it should probably be split.
- Prefer explicit over clever. Avoid deep nesting, dense one-liners, and metaprogramming unless there's a clear, documented reason.

## Docstrings & comments
- Use docstrings (Google or NumPy style — pick one per project and be consistent) on all public modules, classes, and functions. Document parameters, return values, and exceptions raised.
- Comments should explain why, not what — the code itself should make the "what" clear through good naming.
- Keep a project README.md current: install instructions (via uv), usage examples, and how to run tests/lint.

## Error handling & logging
- Never use bare except:. Catch specific exceptions.
- Define custom exception classes for domain-specific error conditions rather than raising generic Exception or ValueError everywhere.
- Fail loudly and clearly on bad input — a clear error message with context (what file, what field, what was expected) beats a stack trace or a silent None.
- Use the standard logging module, not print(), for anything beyond a quick throwaway script. Configure log level via a CLI flag (--verbose/--quiet) rather than hardcoding.
- Don't swallow exceptions silently ("except: pass"). If you must ignore an error, log why.

## Testing
- Use pytest (uv add --dev pytest).
- Mirror the package structure under tests/ (e.g. vcf2curator/commands/snvs.py → tests/test_commands/test_snvs.py).
- Use fixtures (conftest.py) for shared setup, and small realistic sample data files under tests/fixtures/ rather than large real datasets.
- Test both the happy path and realistic failure modes (malformed input, missing fields, empty files).
- Aim for tests to run fast and deterministically — no network calls, no reliance on external services, no reliance on wall-clock time without mocking it.
- Run tests with uv run pytest (never python -m pytest with a manually-activated venv).

## CLI tools specifically
- Prefer Typer or Click over hand-rolled argparse wrappers for anything with subcommands, multiple flags, or that needs --help to be genuinely useful.
- Every subcommand should have a clear one-line help string and, ideally, an example in --help or the README.
- Validate inputs early (file exists, expected format) and fail with a clear message before doing expensive work.
- Support --verbose/-v and --quiet/-q where output volume matters. Support --dry-run for anything destructive or that writes files, where feasible.
- Exit codes matter: 0 for success, non-zero (and distinct codes where it's useful to distinguish failure modes) for failure — don't let the process exit 0 after a caught, logged error.

## Version control hygiene
- .gitignore should exclude .venv/, __pycache__/, *.pyc, .pytest_cache/, .ruff_cache/, .mypy_cache/, and any local data or output files — but not uv.lock.
- Write commit messages that explain intent ("why"), not just a restatement of the diff.
- Keep commits focused — one logical change per commit, so history is reviewable and git bisect-friendly.
- Don't commit secrets, credentials, or real patient/sample data into the repo, even in test fixtures — synthesize or anonymize sample data.

## General do's and don'ts summary

### Do:

- Use uv for everything environment/dependency related.
- Use type hints and docstrings everywhere in library code.
- Write tests alongside new features, not after.
- Keep functions small, named clearly, and single-purpose.
- Log meaningfully; fail with clear, actionable error messages.
- Keep the README and CLAUDE.md up to date as the project evolves.

### Don't:

- Don't use pip/conda/poetry in this project.
- Don't commit requirements.txt as a dependency source of truth.
- Don't use bare except: or silently swallow errors.
- Don't use print() for anything beyond a quick debug session.
- Don't hardcode file paths, credentials, or environment-specific values — use config files, environment variables, or CLI arguments.
- Don't silently change established output formats/schemas without flagging the change — ask first if behavior seems like it should differ from what existing code does.
