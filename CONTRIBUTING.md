# Contributing to Physics Lab Reports

Thank you for your interest in improving and contributing to this open-source physics laboratory collection!

## Ways to Contribute

1. **Bug Reports & Data Corrections:** If you spot a typographical error, sign error, or calculation discrepancy in any of the analysis scripts or LaTeX files, please open an Issue with clear details.
2. **Additional Datasets:** Submitting alternative measurement trials or extended datasets from your own lab sessions is warmly welcomed.
3. **Persian & English Documentation Improvements:** Refinements to theoretical derivations, additional mock oral exam questions, or enhanced visual illustrations.
4. **Code Quality:** Optimizing NumPy/SciPy regressions or enhancing Matplotlib styling and aesthetics.

## Guidelines

- All Python scripts must be standalone, execute headlessly (`matplotlib.use('Agg')`), and produce vector output (PDF format) in a local `plots/` subdirectory.
- All derived parameters must be accompanied by rigorous first-order partial-derivative quadrature uncertainty propagation:
  $$\delta q = \sqrt{\sum \left(\frac{\partial q}{\partial x_i} \delta x_i\right)^2}$$
- LaTeX files should follow standard typesetting conventions (`fouriernc`, `booktabs`, `siunitx`, `xepersian`).
- Commit messages should follow the [Conventional Commits](https://www.conventionalcommits.org/) specification (`feat:`, `fix:`, `docs:`, `style:`, `refactor:`).
