# Correlation of Log-Rank Statistics in a Multi-Arm Trial with a Common Control

## Setup

Consider a three-arm survival trial:

* Treatment A
* Treatment B
* Common Control C

Suppose we compare:

* A vs C
* B vs C

using log-rank tests.

We would like to determine the correlation between the corresponding log-rank statistics.

## Log-Rank Statistic as a Sum of Event Contributions

The log-rank score statistic can be written as

```text
U = Σ(O - E)
```

where:

* `O` is the observed number of events,
* `E` is the expected number of events under the null hypothesis.

Under large-sample theory, the log-rank statistic is asymptotically normal:

```text
U ~ N(0, V)
```

where `V` is the variance of the score statistic.

The standardized log-rank statistic is

```text
Z = U / √V
```

## Relationship Between Information and Variance

For the log-rank test, the quantity

```text
I = V
```

is referred to as the *information*.

Thus:

```text
Information = Variance of the score statistic
```

and

```text
Z = U / √I.
```

Under the null hypothesis,

```text
Z ~ N(0,1).
```

For a two-arm log-rank comparison, information is approximately proportional to the total number of events. Consequently, group sequential survival designs are often described in terms of information fractions rather than sample size fractions.

Throughout this note, we use the terms *information* and *variance of the score statistic* interchangeably.

## Information Contributions by Arm

For intuition, suppose each arm contributes the same amount of information:

```text
I_A = I_B = I_C = I.
```

Let:

```text
Var(U_A) = I_A
Var(U_B) = I_B
Var(U_C) = I_C
```

where:

* `U_A` is the score contribution from arm A,
* `U_B` is the score contribution from arm B,
* `U_C` is the score contribution from the common control.

The treatment-versus-control score statistics can be written as

```text
U_AC = U_A - U_C
U_BC = U_B - U_C.
```

Under the null hypothesis, assume

```text
Cov(U_A,U_B) = 0
Cov(U_A,U_C) = 0
Cov(U_B,U_C) = 0.
```

This assumption provides the standard intuition behind the correlation calculation.

> **Technical note.**
>
> Strictly speaking, efficient score contributions from different arms are not completely independent because all comparisons share the same risk sets. However, the covariance matrix of the efficient score vector yields the same asymptotic result. Under equal allocation, the correlation between the treatment-versus-control log-rank statistics is 0.5.

## Variance of Each Log-Rank Statistic

For A vs C:

```text
Var(U_AC)
=
Var(U_A - U_C)
=
Var(U_A) + Var(U_C)
=
I + I
=
2I.
```

Similarly,

```text
Var(U_BC)
=
2I.
```

## Covariance Between the Log-Rank Statistics

```text
Cov(U_AC,U_BC)
=
Cov(U_A-U_C, U_B-U_C).
```

Expanding,

```text
=
Cov(U_A,U_B)
- Cov(U_A,U_C)
- Cov(U_C,U_B)
+ Var(U_C).
```

The first three terms are zero, leaving

```text
Cov(U_AC,U_BC)
=
I.
```

## Correlation

Therefore,

```text
Corr(U_AC,U_BC)
=
I / √((2I)(2I))
=
1/2.
```

Thus,

```text
Corr(U_AC,U_BC)
=
0.5.
```

## Standardized Log-Rank Statistics

The standardized statistics are

```text
Z_AC = U_AC / √(2I)
Z_BC = U_BC / √(2I).
```

Therefore,

```text
Corr(Z_AC,Z_BC)
=
Corr(U_AC,U_BC)
=
0.5.
```

The standardization does not change the correlation because both statistics are scaled by constants.

## Extension to Incremental Log-Rank Statistics

Suppose the trial is monitored at multiple interim analyses.

Let

```text
ΔU_AC,k
```

denote the incremental log-rank contribution for A vs C from event interval `k`, and

```text
ΔU_BC,k
```

the corresponding contribution for B vs C.

### Same Event Interval

The statistics

```text
ΔU_AC,k
```

and

```text
ΔU_BC,k
```

share the same control-arm events.

Under equal allocation,

```text
Corr(ΔU_AC,k, ΔU_BC,k)
≈ 0.5.
```

### Different Event Intervals

For non-overlapping event intervals,

```text
ΔU_AC,k
```

and

```text
ΔU_BC,l
```

with `k ≠ l` are asymptotically independent because of the independent-increment property of the log-rank score process.

Therefore,

```text
Corr(ΔU_AC,k, ΔU_BC,l)
≈ 0.
```

## Summary

Under equal allocation and a common control arm:

```text
Corr(U_AC,U_BC)
=
Corr(Z_AC,Z_BC)
=
0.5.
```

The same intuition applies to incremental log-rank statistics within a common event interval. The correlation arises entirely from the shared control-arm information.
