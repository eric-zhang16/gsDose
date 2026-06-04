# Why Is the Correlation 0.5 Between Two Treatment-vs-Control Statistics?

## Setup

Consider a three-arm trial:

* Treatment A
* Treatment B
* Common Control C

Suppose we compute incremental log-rank statistics for:

* A vs C
* B vs C

For intuition, write the statistics as:

```text
U_AC = X_A - X_C
U_BC = X_B - X_C
```

where:

```text
Var(X_A) = Var(X_B) = Var(X_C) = σ²
```

and

```text
Cov(X_A, X_B) = 0
Cov(X_A, X_C) = 0
Cov(X_B, X_C) = 0
```

That is, the arm-specific contributions are independent and have equal variance.

## Step 1: Compute the Variances

For A vs C:

```text
Var(U_AC)
= Var(X_A - X_C)
= Var(X_A) + Var(X_C)
= 2σ²
```

Similarly:

```text
Var(U_BC)
= 2σ²
```

## Step 2: Compute the Covariance

```text
Cov(U_AC, U_BC)
= Cov(X_A - X_C, X_B - X_C)
```

Expanding:

```text
= Cov(X_A, X_B)
  - Cov(X_A, X_C)
  - Cov(X_C, X_B)
  + Var(X_C)
```

Since all cross-arm covariances are zero:

```text
Cov(U_AC, U_BC)
= σ²
```

## Step 3: Compute the Correlation

```text
Corr(U_AC, U_BC)
=
Cov(U_AC, U_BC)
/
sqrt(Var(U_AC) * Var(U_BC))
```

Substituting the values:

```text
Corr(U_AC, U_BC)
=
σ² / sqrt((2σ²)(2σ²))
```

Since:

```text
sqrt((2σ²)(2σ²))
= 2σ²
```

we obtain:

```text
Corr(U_AC, U_BC)
=
σ² / (2σ²)
=
0.5
```

## Intuition

Each treatment-vs-control comparison consists of:

```text
Treatment contribution
−
Control contribution
```

The two comparisons share the same control contribution.

Under equal allocation:

```text
Variance from treatment arm = σ²
Variance from control arm   = σ²
Total variance             = 2σ²
```

The shared variance is therefore:

```text
σ²
```

out of a total variance of:

```text
2σ²
```

giving:

```text
Correlation = σ² / (2σ²) = 0.5
```

## Application to Log-Rank Statistics

Under large-sample theory, incremental log-rank statistics from multiple treatment-vs-control comparisons are approximately multivariate normal.

With:

* equal allocation,
* a common control arm, and
* comparable information,

the pairwise correlation between treatment-vs-control statistics is approximately:

```text
0.5
```

because the comparisons share the same control-arm information.
