import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

plt.rcParams.update({
    'font.size': 11,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'figure.dpi': 130,
})

a = 0.18
b = 2.0e-9
n = 1.101e-7
s = 13000.0
d = 0.0412
p = 0.1245
g = 2.019e7
m = 3.422e-10
C_T0, C_I0 = 5.0e6, 3.2e5

D = 0.20
k = 7.5

def kill(G):
    return D / (1 + k * G)

G_LOW = 0.013
G_HIGH = 0.243
L1_TAXA = "E. coli, F. prausnitzii, B. ovatus/dorei/fragilis, R. gnavus"

net_pressure = a - n * s / d
G_crit = (D / net_pressure - 1) / k

def final_tumor(G, T_end=1500):
    rhs = lambda t, y: [
        a*y[0] - b*y[0]**2 - n*y[0]*y[1] - kill(G)*y[0],
        s - d*y[1] + p*y[0]*y[1]/(g+y[0]) - m*y[0]*y[1],
    ]
    sol = solve_ivp(rhs, (0, T_end), [C_T0, C_I0], method='LSODA',
                    rtol=1e-9, atol=1e-2, t_eval=np.linspace(0, T_end, 400))
    return sol.t, np.maximum(sol.y[0], 0.0)

cured = lambda G: final_tumor(G)[1][-1] < 1e3
G_flip = 0.9 * G_crit

print(f"net tumor pressure  a - n*s/d = {net_pressure:.4f} /day")
print(f"critical L1-GUS load  G_crit  = {G_crit*100:.1f}%")
print(f"Patient LOW   G={G_LOW*100:4.1f}%  -> {'CURED' if cured(G_LOW) else 'tumor escapes'}")
print(f"Patient HIGH  G={G_HIGH*100:4.1f}%  -> {'CURED' if cured(G_HIGH) else 'tumor escapes'}")
print(f"Minimal flip: reduce L1-GUS {G_HIGH*100:.1f}% -> {G_flip*100:.1f}%  "
      f"(dG={G_HIGH-G_flip:.3f}); targets: {L1_TAXA}")
print(f"HIGH after flip  -> {'CURED' if cured(G_flip) else 'tumor escapes'}")

assert G_LOW < G_crit < G_HIGH, "real data must straddle the threshold"
assert cured(G_LOW) and not cured(G_HIGH), "same drug must split fate on real data"
assert cured(G_flip), "minimal flip must cure the high-GUS patient"
print("self-check passed")

C_CURE, C_ESC, C_FLIP = '#1565C0', '#C62828', '#2E7D32'
MC = 1e6

G_scan = np.linspace(0, 0.30, 60)
tumor_M = np.array([final_tumor(G)[1][-1] for G in G_scan]) / MC
TB_LOW  = final_tumor(G_LOW)[1][-1]  / MC
TB_HIGH = final_tumor(G_HIGH)[1][-1] / MC
TB_FLIP = final_tumor(G_flip)[1][-1] / MC

fig1, ax = plt.subplots(figsize=(8, 5))
fig1.subplots_adjust(left=0.12, right=0.96, top=0.9, bottom=0.14)

ax.axvspan(0, G_crit, alpha=0.08, color=C_CURE)
ax.axvspan(G_crit, 0.30, alpha=0.08, color=C_ESC)
ax.axvspan(0, G_LOW, alpha=0.16, color='#999', zorder=0)

ax.plot(G_scan, tumor_M, 'k-', lw=2.5, zorder=3)
ax.axvline(G_crit, color='#666', ls='--', lw=1.5)
ax.text(G_crit+0.005, 6.6, f'G* = {G_crit*100:.1f}%\ntranscritical\nbifurcation', color='#555', fontsize=9)

ax.plot(G_LOW,  TB_LOW,  'o', color=C_CURE, ms=13, zorder=5)
ax.plot(G_HIGH, TB_HIGH, 'o', color=C_ESC,  ms=13, zorder=5)
ax.plot(G_flip, TB_FLIP, 's', color=C_FLIP, ms=11, zorder=5)

ax.text(0.0015, 7.7, 'typical range\n<1.3%', fontsize=8.5, color='#555', va='top')

ax.set_xlabel('Gut Loop-1 GUS carrier abundance  G', fontsize=11.5)
ax.set_ylabel('Final tumor burden  (million cells)', fontsize=11.5)
ax.set_title('Tumor fate versus gut Loop-1 GUS abundance', fontsize=12)
ax.set_xlim(0, 0.30); ax.set_ylim(-0.3, 8.0)
fig1.savefig('fig1_bifurcation.png', bbox_inches='tight', dpi=300)
plt.show()
print("Saved: fig1_bifurcation")

fig2, ax = plt.subplots(figsize=(9, 5))
fig2.subplots_adjust(left=0.12, right=0.97, top=0.9, bottom=0.14)

for G, c, ls, lab in [
    (G_LOW,  C_CURE, '-',  f'Patient A ({G_LOW*100:.1f}%): cleared'),
    (G_HIGH, C_ESC,  '-',  f'Patient B ({G_HIGH*100:.1f}%): persists'),
    (G_flip, C_FLIP, '--', f'Patient B reduced ({G_flip*100:.1f}%): cleared'),
]:
    t, T = final_tumor(G)
    ax.plot(t, T / MC, color=c, lw=2.6, ls=ls, label=lab)

ax.axhline(0, color='gray', ls=':', lw=1)
ax.set_xlabel('Time  (days-scale, Kuznetsov nondimensional)', fontsize=11.5)
ax.set_ylabel('Tumor burden  (million cells)', fontsize=11.5)
ax.set_title('Same irinotecan dose, different gut microbiome, different fate',
             fontsize=12, fontweight='bold')
ax.legend(fontsize=9.5, loc='upper right')
ax.set_xlim(0, 1500); ax.set_ylim(-0.5, 16)
fig2.savefig('fig2_trajectories.png', bbox_inches='tight', dpi=300)
plt.show()
print("Saved: fig2_trajectories")
