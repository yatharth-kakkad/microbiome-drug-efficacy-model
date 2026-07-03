# Microbiome-modulated irinotecan efficacy: a drug-metabolism control channel

A minimal, reproducible model of how a patient's gut microbiome can decide the outcome of
chemotherapy through drug metabolism. This is a different channel from the immune-modulation one
used in recent hierarchical microbiome-control frameworks for cancer.

The drug is irinotecan (colorectal cancer). Its efficacy is set inside the Kuznetsov et al. (1994)
tumor-immune dynamics by letting the delivered drug kill-rate depend on the patient's Loop-1
beta-glucuronidase (GUS) load, using real published metagenomic data. At the same dose, a low-GUS
patient clears the tumor and a high-GUS patient does not. The model gives a closed-form threshold
and the minimal microbiome change that flips the outcome.

![Bifurcation diagram](fig1_bifurcation.png)

![Trajectories](fig2_trajectories.png)

## Mechanism

Irinotecan is activated to SN-38, the cytotoxin that kills tumor cells. The liver inactivates
SN-38 by attaching a glucuronide group, producing SN-38G, which is passed into the gut. Gut
bacterial beta-glucuronidase (GUS) removes the glucuronide and reactivates SN-38 in the gut,
causing dose-limiting toxicity. That toxicity forces dose reductions, so less drug reaches the
tumor. The Loop-1 (L1) GUS structural class are the efficient SN-38G reactivators (Pollet et al.
2017), so a patient's L1-GUS carrier abundance sets how much of the delivered dose is effectively
lost.

## Model

Fast tumor-immune dynamics (Kuznetsov 1994), for tumor cells C_T and effector immune cells C_I:

```
dC_T/dt = a*C_T - b*C_T^2 - n*C_T*C_I - kill(G)*C_T
dC_I/dt = s - d*C_I + p*C_T*C_I/(g + C_T) - m*C_T*C_I
```

Microbiome drug-efficacy channel. The delivered kill-rate is eroded by the gut Loop-1 GUS load G:

```
kill(G) = D / (1 + k*G)
```

This saturating form is the same as the drug-efficacy function used in the reference control
framework, but here it is driven by GUS load rather than tumor resistance.

Closed-form threshold. The tumor-free state is stable when kill(G) > a - n*s/d, which gives a
transcritical bifurcation at:

```
G_crit = ( D / (a - n*s/d) - 1 ) / k
```

Below G_crit the patient is curable; above it the tumor persists. The minimal intervention is to
reduce the L1-GUS carriers below G_crit.

## Parameters and data

| Quantity | Value / source |
|---|---|
| `a, b, n, s, d, p, g, m` (tumor-immune) | Kuznetsov et al. 1994; curated values from BioModels `BIOMD0000000762` |
| `C_T0, C_I0` (initial condition) | BioModels `BIOMD0000000762` |
| `G` = Loop-1 GUS carrier abundance | Sun et al. 2022 (60 human gut metagenomes): typical near 1.3%, maximum 24.3% |
| `D` (delivered dose), `k` (GUS sensitivity) | the only two tuned constants; illustrative, not fitted |

Patients A and B use the reported low (1.3%) and high (24.3%) L1-GUS endpoints from Sun et al.
2022. The raw metagenomes are public on NCBI SRA.

## Result

| | Patient A | Patient B | B after minimal flip |
|---|---|---|---|
| Loop-1 GUS abundance | 1.3% | 24.3% | 4.5% |
| Outcome (same dose) | tumor cleared | tumor persists | tumor cleared |

G_crit is about 5.0%. The high-GUS patient is flipped to cure by reducing the validated L1-GUS
carriers (E. coli, F. prausnitzii, B. ovatus/dorei/fragilis, R. gnavus) below the threshold.

## Run

```bash
pip install -r requirements.txt
python model.py
```

This prints the fates and the closed-form threshold (with three self-check assertions) and writes
`fig1_bifurcation.png` (final tumor burden versus Loop-1 GUS load) and `fig2_trajectories.png`
(tumor size over time for both patients and the intervention). It runs in a few seconds on a
laptop.

## Scope and limitations

This is a mechanistic illustration, not a validated clinical predictor.

- The two calibration constants (`D`, `k`) are chosen for illustration, not fitted to
  pharmacokinetic or toxicity data.
- `G` is the scalar Loop-1 GUS category load, not a full per-taxon community composition.
- The immune influx `s` is held at baseline to isolate the drug-metabolism channel.

Natural extensions: fit `D`, `k`, and `G` to real drug-metabolism and metagenomic data, and
replace the scalar `G` with a per-taxon composition and a learned residual on top of the
mechanistic prior.

## References

- V. Kuznetsov, I. Makalkin, M. Taylor, A. Perelson. Nonlinear dynamics of immunogenic tumors:
  parameter estimation and global bifurcation analysis. Bull. Math. Biol. 56(2):295-321, 1994.
  Model: BioModels `BIOMD0000000762`.
- Y. Sun et al. Beta-glucuronidase pattern predicted from gut metagenomes indicates potentially
  diversified pharmacomicrobiomics. Front. Microbiol. 13:826994, 2022 (CC BY).
- S. Pollet et al. An atlas of beta-glucuronidases in the human intestinal microbiome. Structure
  25(7):967-977, 2017.
- A. E. Gollwitzer, D. A. Subramanian, I. Tucker, G. Traverso. Steering the Evolutionary Game:
  Hierarchical Control of Therapeutic Resistance in Cancer Treatment. NeurIPS 2025 (AI4Science).
  This drug-metabolism channel complements the immune channel of that framework.
