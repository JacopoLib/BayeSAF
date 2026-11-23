# Development

## Repository layout

```
bayesaf/
├── python/            package, examples, pyproject.toml
├── database/          species database
├── exp/               experimental targets
├── docs/
├── LICENSE
└── README.md
```

## Workflow

- Work on a branch, open a merge request into `main`, and have one person review it.
- Commit messages: short imperative summary, then a paragraph on the reason when it is not
  obvious.
- Do not commit run outputs (`runs/`, `*.npz`, `*.log`, `Figures/`) or `.DS_Store`; they are in
  `.gitignore`.
- Tag the commit used for each paper or report (`git tag paper1-v1`).

## Tests

The repository does not contain a `tests/` directory yet. `pytest` is declared in the `dev`
extra. A useful first step is a regression test that runs a short CHRJ configuration with a
fixed seed and compares the MAP and the chains with a stored reference.

## Reproducibility

- The sampler takes a `seed`. Record it with the run.
- Record the commit hash and the `DATASET`/`DATA_DIR` used. A one-line `run_info.txt` next to
  each run is enough.

## License

Prosperity Public License 3.0.0 (non-commercial). See `LICENSE`.

## Maintaining the documentation

The documentation lives in `docs/` (Markdown, rendered by GitLab in the repository browser) and
is versioned together with the code. A merge request that changes behaviour or parameters
should update the matching page.
