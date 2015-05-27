# Poisson Experiments
Re-implementing R's native rpois function in Ruby and Erlang for Lamda &lt;= 10

Analyze Data in R with:

```R
# Native Implementation

rpois(1000000, 20.0)

# Read results from Ruby / Erlang
results <- read.table("/tmp/samples")
mean(results(1:1000000, 1))
sd(results(1:1000000, 1))
plot(hist((results(1:1000000, 1))))

```

# NOTE

I've included the original native C implementation of rpois in R which is licensed under the GNU General Public License
