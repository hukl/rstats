# Description

The ```rstats``` module implements random sampling of the Poisson (http://en.wikipedia.org/wiki/Poisson_distribution) and exponential distribution (http://en.wikipedia.org/wiki/Exponential_distribution). The algorithms are re-implemented from R's native C implementions which
themselves are based on the following papers:

Ahrens, J.H. and Dieter, U. (1982). Computer generation of Poisson deviates from modified normal distributions.
ACM Trans. Math. Software 8, 163-179.

and

Ahrens, J.H. and Dieter, U. (1972). Computer methods for sampling from the exponential and normal distributions.
Comm. ACM, 15, 873-882.

Besides that there are a few helper functions for generating CSV files to import them into R for
comparison and validation as well as common copies of floor, ceiling, fsign and fact.


### Usage

```erlang
rstats:rpois(Lamda) - returns random sample from Poisson distribution
rstats:rexp()       - returns random sample from Exponential Distribution
```

### Create Large Samples

```erlang
rstats:write_csv("/tmp/samples.csv", [rstats:rpois(32.0) || _ <- lists:seq(1,1000000)]).

rstats:write_float_csv("/tmp/samples.csv", [rstats:rexp()) || _ <- lists:seq(1,1000000)]).

```

#### Analyze Data in R with:

```R
# Native Implementation
rpois(1000000, 20.0)
rexp(1000000)

# Read CSV
results <- read.table("/tmp/samples.csv")

# Analyze Samples
mean(results(1:1000000, 1))
sd(results(1:1000000, 1))
plot(hist(results(1:1000000, 1)))

# Compare with native implementation
plot(hist(rpois(1000000, 20.0)))
plot(hist(rexp(1000000)))
```

# NOTE

I've included the original native C implementation of rpois in R which is licensed under the GNU General Public License
