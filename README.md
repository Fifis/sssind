# Step size selection in numerical differences

Replication codes for the article ‘[Step size selection in numerical differences using a regression kink](https://orbilu.uni.lu/handle/10993/64958)’.

Tested with R 4.5.0 and the [pnd](https://CRAN.R-project.org/package=pnd) package. The function `step.K()` in this package specifically implements the algorithm from the article with extra safety checks and reasonable fall-back returns. The `plot()` method produces truncation--round-off error plots that provide a reliable visual guidance on the location of the optimal step size.
