import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root_scalar

# Input parameters
st.title("Cracked Section Analysis")

b = st.number_input("Section width b (mm)", value=1000)
h = st.number_input("Section height h (mm)", value=500)
cover = st.number_input("Concrete cover (mm)", value=60)
dia = st.number_input("Rebar diameter ∅ (mm)", value=20)
spacing = st.number_input("Rebar spacing (mm)", value=180)

E_s = st.number_input("Steel modulus E_s (MPa)", value=200000)
E_concrete = st.number_input("Concrete modulus E_concrete (MPa)", value=35000)
E_cm = E_concrete / 4  # corrected as per long-term concrete modulus
st.write(f"Calculated E_cm (MPa) = {E_cm:.2f}")
alpha = E_s / E_cm

# Areas and positions
A_sc = (b / spacing) * (np.pi * dia**2 / 4)
A_s1 = 0
A_s2 = A_sc
d_sc = cover + 1.5 * dia
d1 = h - (cover + 1.5 * dia) - dia
d2 = h - (cover + 1.5 * dia)

# First moment function for cracked section
def S_T(x):
    return (
        b * x * (x / 2)
        + (alpha - 1) * A_sc * (x - d_sc)
        + alpha * A_s1 * (x - d1)
        + alpha * A_s2 * (x - d2)
    )

sol = root_scalar(S_T, bracket=[0.1 * h, h - 1], method='bisect')
x = sol.root

# Moment of inertia for cracked section
I_T = (
    (1/12) * b * x**3 + b * x * (x / 2)**2 +
    (alpha - 1) * A_sc * (x - d_sc)**2 +
    alpha * A_s1 * (x - d1)**2 +
    alpha * A_s2 * (x - d2)**2
)

A_T = b * x + (alpha - 1) * A_sc + alpha * A_s1 + alpha * A_s2
EA = A_T * E_cm
EI = I_T * E_cm

# Outputs
st.subheader("Results")
st.write(f"Neutral axis depth x = {x:.2f} mm")
st.write(f"Cracked section moment of inertia I_T = {I_T:.2e} mm⁴")
st.write(f"Axial stiffness EA = {EA:.2e} N")
st.write(f"Flexural stiffness EI = {EI:.2e} N·mm²")
st.write(f"Modular ratio alpha = {alpha:.2f}")

# Compute reinforcement positions for plot
n_bars = b // spacing
bar_x = np.linspace(cover + dia/2, b - cover - dia/2, int(n_bars))
top_y = cover + dia/2
bot_y = h - cover - dia/2

# Draw section
fig, ax = plt.subplots(figsize=(6, 3))

# Concrete outline
ax.plot([0, b, b, 0, 0], [0, 0, h, h, 0], 'k', lw=2)

# Reinforcement
ax.scatter(bar_x, [top_y]*len(bar_x), color='blue', label='Top rebar')
ax.scatter(bar_x, [bot_y]*len(bar_x), color='red', label='Bottom rebar')

# Annotations
ax.text(b, h+10, f"b = {b} mm", va='bottom', ha='right')
ax.text(b+10, top_y, f"cover = {cover} mm", va='center', fontsize=9)
ax.text(b+10, bot_y, f"∅ = {dia} mm", va='center', fontsize=9)

# Formatting
ax.set_xlim(-50, b + 100)
ax.set_ylim(-50, h + 100)
ax.set_aspect('equal')
ax.axis('off')
ax.legend(loc='upper right')

# Display in Streamlit
st.pyplot(fig)
