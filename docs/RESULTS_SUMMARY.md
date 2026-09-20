# Model 3 results summary

Sample: 72 participants. Each model has residual df=62. `q(12)` is BH-FDR across the 12 prespecified interactions.

| Edge | B | Standardized B | 95% CI | t | Partial r | p | q(12) | HC3 p | HC3 q(12) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Hippocampus-right medial OFC | 0.428 | 0.130 | [-0.188, 1.044] | 1.39 | 0.174 | .1698 | .2264 | .3944 | .5258 |
| Hippocampus-right anterior OFC | 1.195 | 0.198 | [-0.070, 2.459] | 1.89 | 0.233 | .0637 | .0955 | .3332 | .5017 |
| Hippocampus-right posterior OFC | 0.691 | 0.072 | [-1.003, 2.386] | 0.82 | 0.103 | .4179 | .4559 | .6890 | .6890 |
| Hippocampus-right insula | 0.543 | 0.251 | [0.161, 0.925] | 2.84 | 0.339 | .0061 | .0146 | .2757 | .5017 |
| Amygdala-right medial OFC | 3.712 | 0.457 | [1.713, 5.710] | 3.71 | 0.426 | .000443 | .00177 | .00732 | .0439 |
| Amygdala-right anterior OFC | 16.478 | 0.801 | [9.219, 23.737] | 4.54 | 0.499 | .0000267 | .000160 | .00396 | .0439 |
| Amygdala-right posterior OFC | 63.271 | 0.342 | [27.470, 99.072] | 3.53 | 0.409 | .000783 | .00235 | .3344 | .5017 |
| Amygdala-right insula | 1.764 | 0.461 | [0.991, 2.537] | 4.56 | 0.501 | .0000244 | .000160 | .0132 | .0527 |
| dmPFC-right medial OFC | -0.106 | -0.057 | [-0.394, 0.182] | -0.73 | -0.093 | .4654 | .4654 | .6196 | .6760 |
| dmPFC-right anterior OFC | 0.439 | 0.204 | [0.073, 0.804] | 2.40 | 0.291 | .0195 | .0335 | .1028 | .2467 |
| dmPFC-right posterior OFC | 0.355 | 0.205 | [0.097, 0.612] | 2.75 | 0.330 | .00774 | .0155 | .0708 | .2123 |
| dmPFC-right insula | -0.462 | -0.122 | [-1.181, 0.258] | -1.28 | -0.161 | .2045 | .2454 | .5298 | .6358 |

## Conventional OLS inference

Seven interactions pass the 12-test BH-FDR:

1. Hippocampus-right insula.
2. Amygdala-right medial OFC.
3. Amygdala-right anterior OFC.
4. Amygdala-right posterior OFC.
5. Amygdala-right insula.
6. dmPFC-right anterior OFC.
7. dmPFC-right posterior OFC.

## HC3 robust sensitivity inference

Two interactions pass a separate HC3 12-test BH-FDR:

1. Amygdala-right medial OFC, HC3 q=.0439.
2. Amygdala-right anterior OFC, HC3 q=.0439.

Amygdala-right insula is just above the HC3 threshold, q=.0527.

## Interpretation

The significant conventional interactions are positive: higher baseline KLDS shifts the UCLA-LS-change association with follow-up MoCA toward a more positive slope. This is compatible with resilience-like or cross-over moderation, but it is not evidence that worsening loneliness improves cognition.

The strongest defensible findings are Amygdala-right medial OFC and Amygdala-right anterior OFC because they survive both conventional and HC3 FDR. The full coefficients, diagnostics, and all 36 simple slopes are available in the `output/` directory.

