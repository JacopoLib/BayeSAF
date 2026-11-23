# Development

## Repository layout

```
BayeSAF/
├── matlab/            MATLAB implementation, with examples
├── python/            Python package, examples, pyproject.toml
├── database/          species database
├── exp/               experimental targets
├── docs/              documentation
├── img/
├── LICENSE
└── README.md
```

## Workflow

- Commit messages: a short imperative summary, then a paragraph on the reason when it is not
  obvious.
- Do not commit run outputs (`Figures/`, `MAP.txt`, log files, `.mat` and `.npz` files) or `.DS_Store`; they are in
  `.gitignore`.

## Checking a change

After modifying the sampler, the likelihood or the post-processing, run a short configuration of the example (for instance
1500 iterations with a burn-in of 500) and check that the run completes, and that the MAP, the effective sample size and the
figures are what you expect. With the same seed, a change that does not alter the method must give the same chains.

## Reproducibility

- The sampler takes a seed (`seed` in Python, `'Seed'` in MATLAB). Record it with the run.
- Record the commit hash and the data set (`DATASET`/`DATA_DIR` in Python, `Dataset` in MATLAB) used. A one-line note next to
  each run is enough.

## License

Prosperity Public License 3.0.0 (non-commercial). See `LICENSE`.

## Maintaining the documentation

The documentation lives in `docs/` (Markdown, rendered in the repository browser) and
is versioned together with the code. A change that alters the behaviour or the parameters
should update the matching page.
