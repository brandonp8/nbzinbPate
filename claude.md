 STA 651 final project: NB and ZINB regression R package
## Rules (must follow)
- Fitting code in R/ must NOT call optim, optimize, nlm, nlminb, uniroot,
glm, glm.fit, lm, lm.fit, lm.wfit, MASS::glm.nb, MASS::theta.ml,
pscl::zeroinfl, or any other optimizer package.
- NB: IRLS for beta alternated with Newton-Raphson on log(theta).
- ZINB: EM algorithm; M-steps are IRLS (weighted logistic + weighted NB).
- MASS::glm.nb and pscl::zeroinfl may be used ONLY in tests/,
inst/benchmark/, and the report, for comparison.
- Follow the equations in project-docs/STA651_Final_Project.pdf, Section 3.
- Comment estimation code with the equation it implements.
- Explain your reasoning to me; I must be able to explain every line.
## Workflow
- After editing roxygen comments: Rscript -e 'devtools::document()'
- Run tests: Rscript -e 'devtools::test()'
- Full check (must be 0 errors, 0 warnings): Rscript -e 'devtools::check()'
- Report PDF: quarto render vignettes/articles/report.qmd --to pdf
- Website: Rscript -e 'pkgdown::build_site()'
- Make small git commits with clear messages after each working step.