# Project Guidelines

## Code Style
- Preserve the current academic writing style and mathematical notation. Keep matrix and vector symbols consistent with the existing document, including forms such as $\mathbf{Y}$, $\mathbf{L}$, and $D_{it}$.
- Prefer display math for formal definitions and derivations, and keep tables formatted with `booktabs` conventions.

## Architecture
- This workspace is a single-source LaTeX project. Edit `literature-review.tex` and treat `literature-review.aux`, `literature-review.log`, `literature-review.out`, and `literature-review.pdf` as generated artifacts.
- Preserve the existing sectioned paper structure unless the task explicitly requires reorganizing the document.

## Build and Test
- Compile from the workspace root with `pdflatex -interaction=nonstopmode literature-review.tex`.
- There is no automated test suite. After edits, check `literature-review.log` for new LaTeX errors. The current build succeeds with some existing underfull or overfull box warnings and float placement warnings, so only treat them as regressions if your changes make layout materially worse.

## Conventions
- Keep the current preamble lean. Only add LaTeX packages when they are necessary for the requested change.
- Preserve the existing inline numeric citation style such as `[1, 2]`. There is no `.bib` or BibTeX workflow in this workspace.
- Use `\texttt{}` for package, library, and algorithm names.
- The document relies on manual prose references rather than `\label{}` and `\ref{}`. If you rename or reorder sections, tables, or algorithms, update the surrounding prose references as part of the same change.