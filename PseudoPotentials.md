
# 1. Summary of Pseudopotential Types in DFT

In Density Functional Theory (DFT), pseudopotentials replace core electrons to simplify calculations.  
There are several main types, each with their own strengths and weaknesses.

---

## Main Types

| Type | Name | Main Feature | Pros | Cons |
|:---|:---|:---|:---|:---|
| **NCPP** | Norm-Conserving Pseudopotentials | Charge outside the core matches the real atom exactly | Simple, very transferable | Requires **high** kinetic energy cutoff (80–100 Ry) |
| **USPP** | Ultrasoft Pseudopotentials | Relaxes norm conservation to reduce cutoff | Lower cutoff (more efficient) | Complex, needs augmentation charge treatment |
| **PAW** | Projector Augmented-Wave | Reconstructs full all-electron wavefunction inside spheres | Very accurate, behaves like all-electron | Slightly heavier computations |
| **ONCV** | Optimized Norm-Conserving Vanderbilt | NCPP optimized for lower cutoff energy | Good balance between accuracy and efficiency | Newer, not available for every element |
| **AE** | All-Electron Methods (e.g., LAPW) | No pseudopotential at all, full real potential used | Ultimate accuracy | Extremely computationally expensive |

---

## Brief Descriptions

### NCPP (Norm-Conserving Pseudopotentials)
- Classical method (e.g., Troullier–Martins, Hamann–Schlüter–Chiang).
- Accurate but needs **higher plane-wave energy cutoff**.
- Still important for very precise applications.

### USPP (Ultrasoft Pseudopotentials)
- Developed to **reduce** the required cutoff.
- Introduces **augmentation charges** inside atomic spheres.
- Faster, especially for large and metallic systems.

### PAW (Projector Augmented-Wave)
- Combines the benefits of **pseudopotential** and **all-electron** calculations.
- Provides **very accurate forces**, **phonons**, and **magnetic properties**.
- Slightly larger memory and CPU usage.

### ONCV (Optimized Norm-Conserving Vanderbilt)
- A modern improvement of NCPP.
- Maintains norm-conservation but **optimizes** for lower cutoff energies (~50–60 Ry).
- Available in libraries like **PseudoDojo** and **SG15**.

### AE (All-Electron Methods)
- No approximation for core electrons.
- Solved using methods like **LAPW** (e.g., WIEN2k) or **FLAPW**.
- Needed for **extreme precision** studies, but extremely costly computationally.

---

## Practical Guide: Which to Choose?

| Material | Recommended Pseudopotential |
|:---|:---|
| Simple metals (Al, Na) | USPP or ONCV |
| Transition metals (Nb, Mo, Fe) | PAW |
| Magnetic materials | PAW preferred |
| Molecules and quantum chemistry | NCPP |
| High-pressure phases | NCPP or PAW |

---

> **Tip:**  
> In Quantum ESPRESSO, pseudopotential files have `.UPF` extension and include metadata about their type (NCPP, USPP, PAW, etc.).

---






# 2. Side by side comparison of pseudo potential of Nb
## Difference between Nb.pz-spn-kjpaw_psl.1.0.0.UPF and Nb.pbe-spn-rrkjus_psl.1.0.0.UPF

| Feature | **Nb.pz-spn-kjpaw_psl.1.0.0.UPF** | **Nb.pbe-spn-rrkjus_psl.1.0.0.UPF** |
|:---|:---|:---|
| Functional | **PZ** (Perdew–Zunger, LDA) | **PBE** (Perdew–Burke–Ernzerhof, GGA) |
| Pseudopotential type | **PAW** (Projected Augmented Wave) | **USPP** (Ultrasoft Pseudopotential) |
| Norm-conserving? | Effectively yes (through PAW augmentation) | No (true USPP, allows lower cutoff) |
| Kinetic energy cutoff needs | Moderate (~50-80 Ry) | Higher (~60 Ry for wavefunctions, ~600 Ry for charge density) |
| Transferability (accuracy across environments) | Very high | Good but slightly lower than PAW |
| Best for | **High-accuracy** DFT (especially forces and phonons) | **Efficiency**, larger systems |
| File name hints | `kjpaw_psl`: "Kleinman–Bylander Projector + PAW" | `rrkjus_psl`: "Rappe–Rabe–Kaxiras–Joannopoulos Ultrasoft" |

---

## In words

- **Nb.pz-spn-kjpaw_psl** is **PAW** and uses the **PZ (LDA)** exchange-correlation functional.
  - LDA generally **overbinds** (gives smaller lattice constants).
  - PAW is **very accurate** but a little heavier to run than USPP.

- **Nb.pbe-spn-rrkjus_psl** is **Ultrasoft**, with **PBE (GGA)** functional.
  - PBE tends to **underbind** slightly (larger lattice constants), but it's more realistic for many materials.
  - Ultrasoft pseudopotentials **allow lower plane-wave cutoffs** for faster calculations.

---

## Important practical points

- **You can't mix them**. If you're using PBE pseudopotentials, your calculation should also use `input_dft = 'PBE'`.
- **PAW** pseudopotentials behave like norm-conserving, but still benefit from slightly higher cutoffs than USPP.
- **USPP** are more efficient for large systems (saves memory and CPU time).

---

## Which one should you use?

✅ If you want **speed and reasonable accuracy** → use your **PBE USPP** (`Nb.pbe-spn-rrkjus_psl`).  
✅ If you want **very accurate phonons or forces**, and are okay with slower jobs → use the **PAW** one.

---

> **Extra tip:**  
> If you're studying *phonons* and care about **dynamic stability**, **PAW** can sometimes be more reliable.  
> But PBE-USPP is often *good enough*, especially for early-stage calculations!

---





# 3. Key Terms and Features in Pseudopotentials

When choosing a pseudopotential for your DFT calculations, it's important to understand several key features that can impact your results. Below is a summary of terms and concepts related to **relativistic effects**, **spin-orbit coupling (SOC)**, and more.

---

## Pseudopotential Features to Look For

| **Feature** | **Description** | **Impact on Calculations** | **Examples** |
|:---|:---|:---|:---|
| **Relativistic Treatment** | Whether the pseudopotential includes relativistic effects (spin-orbit, mass-velocity correction, etc.). | Relativistic treatments become more important for heavier elements. A **fully relativistic** pseudopotential accounts for all effects, while **scalar relativistic** only considers mass-velocity correction without spin-orbit. | - **Scalar Relativistic**: `Nb.pbe-dn-kjpaw_psl.0.2.2.UPF` <br> - **Fully Relativistic**: `Te.rel-pbe-dn-kjpaw_psl.0.2.2.UPF` |
| **Spin-Orbit Coupling (SOC)** | A relativistic effect where the electron's spin interacts with its orbital motion. | Necessary for accurate description of **magnetic materials**, **topological states**, or when **heavy elements** (like Te, Bi, Pb) are involved. If enabled in your calculations, use fully relativistic potentials for both elements. | - **SOC Included**: `Te.rel-pbe-dn-kjpaw_psl.0.2.2.UPF` |
| **Norm-Conserving** | Refers to pseudopotentials that maintain **total charge** in the core region exactly as in the all-electron case. | High accuracy but requires **high cutoffs** (50–100 Ry), making computations slower. | - `Nb.pbe-dn-kjpaw_psl.0.2.2.UPF` |
| **Ultrasoft** | **Ultrasoft Pseudopotentials (USPP)** relax the norm-conserving condition for faster calculations. | **Faster computations** but at the expense of slightly reduced accuracy. **Lower cutoffs** (typically ~40-60 Ry) are possible. | - `Nb.pbe-spn-rrkjus_psl.1.0.0.UPF` |
| **PAW** (Projector Augmented Wave) | Pseudopotential method combining **all-electron behavior** inside a region around the atom and a **pseudo-core outside**. | Very accurate for forces, phonons, and electronic structure, especially for **transition metals** or **materials with complex bonding**. | - `Nb.pbe-dn-kjpaw_psl.0.2.2.UPF` |
| **On-the-Fly** | A feature where pseudopotentials are generated during the calculation (on-the-fly generation). | Efficient for large systems where a custom pseudopotential may be necessary. Typically used for **non-standard elements** or when **no suitable database exists**. | - Quantum ESPRESSO supports generating pseudopotentials on the fly. |
| **Ultra-Soft Core** | Allows **fewer plane waves** to describe the core region, which is beneficial for larger systems. | **Faster calculations** for large systems at the expense of accuracy in core states. | - `Nb.pbe-spn-rrkjus_psl.1.0.0.UPF` |
| **All-Electron (AE)** | Refers to pseudopotentials that model **core electrons** explicitly, providing the **most accurate description** of an atom. | Very accurate but requires significantly more computational resources. Used for **high-accuracy applications**. | - Typically used in **LAPW** or **APW** methods. |
| **Relativistic Effects** | Pseudopotentials that include **relativity** (mass-velocity, spin-orbit effects) to simulate relativistic phenomena. | Necessary for heavy atoms where **relativistic effects** significantly alter electronic behavior (e.g., transition metals, lanthanides). | - `Te.rel-pbe-dn-kjpaw_psl.0.2.2.UPF` |
| **Pseudopotential Type (PAW, USPP)** | The form of the pseudopotential, i.e., **PAW** (Projector Augmented-Wave) or **USPP** (Ultrasoft). | **PAW** is more accurate, but **USPP** is faster, especially in large systems. | - **PAW**: `Nb.pbe-dn-kjpaw_psl.0.2.2.UPF` <br> - **USPP**: `Nb.pbe-spn-rrkjus_psl.1.0.0.UPF` |

---

## Choosing the Right Pseudopotential

| **Scenario** | **Recommended Pseudopotential** |
|:---|:---|
| **Fast, large system calculations (supercells, defects)** | USPP (e.g., `Nb.pbe-spn-rrkjus_psl.1.0.0.UPF`) |
| **Precise calculation of phonons or forces** | PAW (e.g., `Nb.pbe-dn-kjpaw_psl.0.2.2.UPF`) |
| **Magnetic materials or systems with strong SOC** | Fully relativistic (e.g., `Te.rel-pbe-dn-kjpaw_psl.0.2.2.UPF`) |
| **Transition metals** | PAW or NCPP (e.g., `Nb.pbe-dn-kjpaw_psl.0.2.2.UPF`) |
| **Light elements (C, N, O)** | USPP or NCPP for speed, PAW for accuracy |

---

## Key Considerations for Pseudopotential Choice

- **Relativistic Treatment**: For heavier elements (Te, Pb, Bi), always consider fully relativistic pseudopotentials.
- **SOC**: If you're working with **topological materials** or need accurate **spin-orbit coupling**, use fully relativistic pseudopotentials for both elements.
- **Accuracy vs Efficiency**: Choose PAW for **accurate forces/phonons** or USPP for **faster calculations** (especially for large systems).
- **Cutoff Energy**: For **USPP**, lower cutoffs (~40–60 Ry) are typical. For **PAW**, you will need higher cutoffs (~60–100 Ry).

---

> **Tip**: Always test your pseudopotential choice on a small, well-studied system (e.g., bulk Nb or Te) before using it for more complex materials like NbTe₂.

---




# Download links
https://pseudopotentials.quantum-espresso.org/upf_files/Te.rel-pbe-n-kjpaw_psl.1.0.0.UPF
https://pseudopotentials.quantum-espresso.org/upf_files/Nb.rel-pbe-spn-kjpaw_psl.1.0.0.UPF
