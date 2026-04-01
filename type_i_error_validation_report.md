# Type I Error Rate Validation Report

Generated: 2026-04-01 15:33:13.798263

## Executive Summary

**Overall Weighted phyper() Type I Error Control: ✗ FAIL**
**Overall Standard phyper() Type I Error Control: ✗ FAIL**

## Test Results

### 1. Multiple Alpha Levels

| Alpha | Weighted Rate | Standard Rate | W_Controlled | S_Controlled |
|-------|---------------|---------------|--------------|--------------|
| 0.010 | 0.1360 | 0.0060 | FALSE | TRUE |
| 0.050 | 0.2300 | 0.0200 | FALSE | FALSE |
| 0.100 | 0.2740 | 0.0420 | FALSE | FALSE |

### 2. Copy Number Distributions

| Scenario | Weighted Rate | Standard Rate | Mean CN | W_Controlled | S_Controlled |
|----------|---------------|---------------|---------|--------------|--------------|
| All_CN_1 | 0.0150 | 0.0225 | 1.0 | FALSE | FALSE |
| All_CN_5 | 0.1725 | 0.0300 | 5.0 | FALSE | FALSE |
| Uniform_1_3 | 0.1350 | 0.0250 | 2.0 | FALSE | FALSE |
| Uniform_1_8 | 0.1975 | 0.0225 | 4.4 | FALSE | FALSE |
| Uniform_1_20 | 0.2675 | 0.0150 | 10.6 | FALSE | FALSE |
| High_Variance | 0.3375 | 0.0100 | 9.8 | FALSE | FALSE |
| Bimodal | 0.2150 | 0.0300 | 3.6 | FALSE | FALSE |

### 3. Dataset Sizes

| Scenario | Background | Pathway | Query | Weighted Rate | Standard Rate | W_Controlled | S_Controlled |
|----------|------------|---------|-------|---------------|---------------|--------------|--------------|
| small | 200 | 20 | 15 | 0.2000 | 0.0433 | FALSE | TRUE |
| medium | 1000 | 100 | 50 | 0.2200 | 0.0367 | FALSE | FALSE |
| large | 5000 | 500 | 200 | 0.2567 | 0.0267 | FALSE | FALSE |
| wide_pathway | 1000 | 300 | 50 | 0.2267 | 0.0500 | FALSE | TRUE |
| narrow_pathway | 1000 | 25 | 50 | 0.2000 | 0.0700 | FALSE | FALSE |
| large_query | 1000 | 100 | 150 | 0.1967 | 0.0500 | FALSE | TRUE |
| small_query | 1000 | 100 | 20 | 0.2067 | 0.0400 | FALSE | FALSE |

## Conclusions

- Type I error rates were tested at tolerance level of ±1%
- Weighted method overall pass rate: Failed
- Standard method overall pass rate: Failed
