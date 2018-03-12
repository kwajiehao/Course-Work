# On Optimal and Objective-Specific Selection of Sequential Offers

A group project with Alessandro de Sanctis and Jan Rademacher. We look at two cases: optimal selection in the multivariate case, and top a% selection in the univariate case. The reason we do so is because there are already many known algorithms which solve the univariate case (for e.g. Bruss, 2000).

For the univariate case, we look at a nonparametric quantile estimator to solve a variant of the secretary problem: accepting an offer
within the top a% of the true distribution with high certainty. We obtain an upper bound using Chernoff Bounding.

For the multivariate case, we use the upper confidence bound algorithm from the family of multi-armed bandit techniques to make selections
on various (possibly correlated stocks). Our algorithm manages to outperform slightly our chosen benchmark of an equally-weighted portfolio.

# Contents

* Paper - Kwa, Alessandro, Jan, Guillem.pdf - Paper summarizing our results
* code - Code used in analysis
