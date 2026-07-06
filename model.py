"""
Microbiome control of chemotherapy outcome through TWO parallel input channels.

Channel 1 - drug metabolism (this repo's original idea):
    Loop-1 gut beta-glucuronidase carriers G reactivate SN-38 in the gut, force
    irinotecan dose reductions, and so ERODE the delivered tumor kill-rate.
        kill(G) = D / (1 + k*G)        # saturating efficacy form of Gollwitzer et al.'s kappa(p)=k0/(1+p)

Channel 2 - immune recruitment (straight from Gollwitzer et al. 2025):
    Immunogenic commensals B (Faecalibacterium/Ruminococcaceae, Akkermansia,
    Bifidobacterium; Gopalakrishnan 2018, Routy 2018) raise effector-immune
    recruitment. This is the report's own affine recruitment law.
        s(M) = s0 + s1*B               # report's s(M) = s0 + sum_i s_i * M_i

Both feed the same Kuznetsov (1994) fast tumor-immune system, Eq.(1) of the
report. Linearising Eq.(1) at the tumor-free state (C_T=0, C_I=s(M)/d) gives the
tumor-free stability margin

        margin = n*s(M)/d + kill(G) - a          # >0 => cleared, <0 => tumor escapes

i.e. combined immune + drug pressure must beat tumor growth. This one line
contains both channels: with the drug fully eroded (kill->0) it collapses to the
report's pure-immune threshold s(M) > (d/n)*a; with B held fixed it is the pure
drug threshold. So the microbiome is a genuine two-lever control input, and each
lever alone can flip a patient's fate.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.integrate import solve_ivp

plt.rcParams.update({'font.size': 11, 'axes.spines.top': False,
                     'axes.spines.right': False, 'figure.dpi': 130})

# --- Kuznetsov 1994 tumor-immune constants (BioModels BIOMD0000000762) ---
a = 0.18          # tumor growth rate
b = 2.0e-9        # tumor competition
n = 1.101e-7      # immune kill rate
d = 0.0412        # immune death rate
r = 0.1245        # immune stimulation rate    (report's r; was 'p')
h = 2.019e7       # stimulation half-saturation (report's h; was 'g')
m = 3.422e-10     # immune inactivation rate
C_T0, C_I0 = 5.0e6, 3.2e5

# --- microbiome channels (D, k, s0, s1 are the only tuned constants; illustrative, not fitted) ---
D  = 0.20         # delivered irinotecan kill strength
k  = 7.5          # GUS sensitivity of delivered dose
s0 = 13000.0      # baseline immune recruitment (few immunogenic commensals)
s1 = 3.0e5        # recruitment gain per unit immunogenic-commensal abundance

kill = lambda G: D / (1 + k * G)          # drug channel: GUS erodes kill-rate
s_of = lambda B: s0 + s1 * B              # immune channel: report's s(M)=s0+s1*B

def margin(B, G):
    """Tumor-free stability margin (/day). >0 => tumor cleared, <0 => tumor escapes.
    Jacobian of Eq.(1) at (C_T=0, C_I=s(M)/d): a - n*s(M)/d - kill(G), sign-flipped."""
    return n * s_of(B) / d + kill(G) - a

def final_tumor(B, G, T_end=1500):
    rhs = lambda t, y: [
        a*y[0]*(1 - b*y[0]) - n*y[0]*y[1] - kill(G)*y[0],
        s_of(B) - d*y[1] + r*y[0]*y[1]/(h+y[0]) - m*y[0]*y[1],
    ]
    sol = solve_ivp(rhs, (0, T_end), [C_T0, C_I0], method='LSODA',
                    rtol=1e-9, atol=1e-2, t_eval=np.linspace(0, T_end, 400))
    return sol.t, np.maximum(sol.y[0], 0.0)

cured = lambda B, G: final_tumor(B, G)[1][-1] < 1e3

# --- patients: (immunogenic commensals B, Loop-1 GUS carriers G), Sun 2022 GUS range ---
B_EUB, G_EUB = 0.20, 0.02     # eubiotic:  strong immunity + intact drug
B_DYS, G_DYS = 0.05, 0.243    # dysbiotic: weak immunity + GUS-eroded drug -> fails
G_DRUG = 0.06                 # rescue 1 (drug channel):  cut GUS carriers, restore irinotecan
B_IMM  = 0.16                 # rescue 2 (immune channel): boost commensals, raise s(M); drug UNCHANGED

s_crit = (d / n) * a                       # report's pure-immune threshold s(M) > (d/n)*a
B_selfsuff = (s_crit - s0) / s1            # commensal level at which immunity alone clears

print(f"pure-immune threshold  s_crit = (d/n)*a = {s_crit:,.0f}")
print(f"  immunity alone clears the tumor once B > {B_selfsuff*100:.1f}%  (recovers report regime)")
print(f"{'patient':24s} {'B':>6s} {'G':>7s} {'margin/day':>11s}  outcome")
patients = [("eubiotic", B_EUB, G_EUB),
            ("dysbiotic", B_DYS, G_DYS),
            ("dysbiotic + drug channel", B_DYS, G_DRUG),
            ("dysbiotic + immune channel", B_IMM, G_DYS)]
for name, B, G in patients:
    print(f"{name:24s} {B*100:5.1f}% {G*100:6.1f}% {margin(B, G):>+11.4f}  "
          f"{'CLEARED' if cured(B, G) else 'tumor escapes'}")

# --- rigor: the analytic stability margin must match the simulated fate, every time ---
for _, B, G in patients:
    assert (margin(B, G) > 0) == cured(B, G), "stability margin must match simulation"
assert not cured(B_DYS, G_DYS), "dysbiotic patient must fail on the same dose"
assert cured(B_DYS, G_DRUG), "drug channel (lower GUS) alone must rescue"
assert cured(B_IMM, G_DYS), "immune channel (raise s(M)) alone must rescue with the drug unchanged"
print("self-check passed: each channel independently flips the outcome")

# ------------------------------------------------------------------ figure 1: two-channel fate map
C_EUB, C_DYS, C_FIX = '#1565C0', '#C62828', '#2E7D32'
gg, bb = np.meshgrid(np.linspace(0, 0.30, 300), np.linspace(0, 0.30, 300))
M = margin(bb, gg)

fig1, ax = plt.subplots(figsize=(7.6, 6))
fig1.subplots_adjust(left=0.12, right=0.99, top=0.92, bottom=0.12)
pc = ax.pcolormesh(gg*100, bb*100, M, shading='auto', cmap='RdBu',
                   norm=TwoSlopeNorm(vcenter=0, vmin=M.min(), vmax=M.max()))
ax.contour(gg*100, bb*100, M, levels=[0], colors='k', linewidths=2)
cb = fig1.colorbar(pc, ax=ax, pad=0.02)
cb.set_label('tumor-free stability margin  n·s(M)/d + kill(G) − a   (/day)')

ax.axhline(B_selfsuff*100, color='#333', ls=':', lw=1.3)
ax.text(0.5, B_selfsuff*100 + 0.5, 'immunity alone clears tumor:  s(M) > (d/n)·a',
        fontsize=8.5, color='#222')
ax.text(21, 2.0, 'tumor escapes', color='#7d1a1a', fontsize=11, ha='center')
ax.text(6.5, 25, 'tumor cleared', color='#123f78', fontsize=11, ha='center')

ax.plot(G_EUB*100, B_EUB*100, 'o', color=C_EUB, ms=13, zorder=5, label='eubiotic (cleared)')
ax.plot(G_DYS*100, B_DYS*100, 'o', color=C_DYS, ms=13, zorder=5, label='dysbiotic (escapes)')
ax.plot(G_DRUG*100, B_DYS*100, 's', color=C_FIX, ms=11, zorder=5)
ax.plot(G_DYS*100, B_IMM*100, '^', color=C_FIX, ms=12, zorder=5)
ax.annotate('', xy=(G_DRUG*100, B_DYS*100), xytext=(G_DYS*100, B_DYS*100),
            arrowprops=dict(arrowstyle='->', color=C_FIX, lw=2))
ax.annotate('', xy=(G_DYS*100, B_IMM*100), xytext=(G_DYS*100, B_DYS*100),
            arrowprops=dict(arrowstyle='->', color=C_FIX, lw=2))
ax.text(13.5, 3.6, 'drug channel:\ncut GUS carriers', color=C_FIX, fontsize=8.5)
ax.text(20.5, 11, 'immune channel:\nboost commensals', color=C_FIX, fontsize=8.5)

ax.set_xlabel('Loop-1 GUS carriers  G  (%)   —  drug-metabolism channel')
ax.set_ylabel('immunogenic commensals  B  (%)   —  immune-recruitment channel')
ax.set_title('Microbiome sets tumor fate through two parallel channels')
ax.set_xlim(0, 30); ax.set_ylim(0, 30)
fig1.savefig('fig1_bifurcation.png', bbox_inches='tight', dpi=300)
plt.show(); print("Saved: fig1_bifurcation")

# ------------------------------------------------------------------ figure 2: trajectories
MC = 1e6
fig2, ax = plt.subplots(figsize=(9, 5))
fig2.subplots_adjust(left=0.11, right=0.97, top=0.9, bottom=0.14)
for B, G, c, ls, lab in [
    (B_EUB, G_EUB, C_EUB, '-',  f'eubiotic (B={B_EUB*100:.0f}%, G={G_EUB*100:.0f}%): cleared'),
    (B_DYS, G_DYS, C_DYS, '-',  f'dysbiotic (B={B_DYS*100:.0f}%, G={G_DYS*100:.0f}%): persists'),
    (B_DYS, G_DRUG, C_FIX, '-', f'+ drug channel  (G {G_DYS*100:.0f}%→{G_DRUG*100:.0f}%): cleared'),
    (B_IMM, G_DYS, C_FIX, '--', f'+ immune channel (B {B_DYS*100:.0f}%→{B_IMM*100:.0f}%, drug fixed): cleared'),
]:
    t, T = final_tumor(B, G)
    ax.plot(t, T/MC, color=c, lw=2.6, ls=ls, label=lab)
ax.axhline(0, color='gray', ls=':', lw=1)
ax.set_xlabel('time  (days-scale, Kuznetsov nondimensional)')
ax.set_ylabel('tumor burden  (million cells)')
ax.set_title('Same tumor, same irinotecan dose — microbiome decides fate', fontweight='bold')
ax.legend(fontsize=9, loc='upper right')
ax.set_xlim(0, 1500); ax.set_ylim(-0.5, 16)
fig2.savefig('fig2_trajectories.png', bbox_inches='tight', dpi=300)
plt.show(); print("Saved: fig2_trajectories")
