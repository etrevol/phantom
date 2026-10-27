# Projects on this phantom fork

Our work, kept apart from the upstream tree. Fortran stays where phantom's Makefile looks
for it (`src/setup/`, `src/utils/`, flat), named after the project; everything else lives
here, one directory per project.

| project | status | what | SETUP |
|---|---|---|---|
| [`ppi/`](ppi/) | frozen 2026-10-26 (tag `ppi/frozen-2026-10-26`) | Papaloizou-Pringle torus, unseeded and seeded | `ppitorus` |
| [`m1torus/`](m1torus/) | active | the slow m = 1 mode of Bannikova et al. (2026) in a self-gravitating torus | `m1torus` |
| [`common/`](common/) | - | shared Python: readers of phantom dumps, .ev and .setup files | - |

## Layout

```
src/setup/setup_<project>.f90, src/utils/{analysis,moddump}_<project>.f90
build/Makefile_setups                  one SETUP block per project
projects/<project>/README.md           what it is and how to run it
projects/<project>/PLAN.md             questions, checks, campaigns (active projects)
projects/<project>/formulary/          LaTeX derivations
projects/<project>/scripts/            theory/, run/, analysis/, vis/
results/<project>/<campaign>/          gitignored run output: tests/, runs/
materials/papers/                      gitignored, shared papers
materials/<project>/                   gitignored, project data
```

## Commits

`[<project>] area: what changed` - lowercase, no period, at most 72 characters, a body of
one or two sentences. Tags: `[ppi]`, `[m1torus]`, `[common]` for `projects/common/`,
`[repo]` for layout and ignore rules. Areas: `setup`, `analysis`, `build`, `theory`,
`run`, `vis`, `docs`. The PPI commits before the tag `ppi/frozen-2026-10-26` predate this
convention and carry no project tag.
