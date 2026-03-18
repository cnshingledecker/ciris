```
================================================================================

        ██████╗██╗██████╗ ██╗███████╗
       ██╔════╝██║██╔══██╗██║██╔════╝
       ██║     ██║██████╔╝██║███████╗
       ██║     ██║██╔══██╗██║╚════██║
       ╚██████╗██║██║  ██║██║███████║
        ╚═════╝╚═╝╚═╝  ╚═╝╚═╝╚══════╝

   Chemistry  of  Ionizing  Radiation  In  Solids
   -------------------------------------------------------
              >> KINETIC MONTE CARLO ENGINE <<
                      V E R S I O N   1 . 1

================================================================================
```

---

## ═══════════════════  W H A T   I S   C I R I S  ═══════════════════

**CIRIS** is a three-dimensional kinetic Monte Carlo (KMC) simulation of the
radiation chemistry of astrophysical ices under proton / cosmic-ray bombardment.
Given a crystal lattice of ice molecules, CIRIS follows an energetic proton as it
punches through the solid, generates a shower of secondary electrons, and watches
the resulting radical chemistry unfold — site by site, hop by hop — until a
prescribed fluence is reached.


> **Cite this code:**
> C. N. Shingledecker, R. Le Gal & E. Herbst,
> *"A new model of the chemistry of ionizing radiation in solids: CIRIS"*,
> **Phys. Chem. Chem. Phys.**, 2017, **19**, 11043–11056.
> DOI: [10.1039/C7CP01472D](https://doi.org/10.1039/C7CP01472D)

---

## ════════════════════  P H Y S I C A L   M O D E L  ════════════════════

```
  ┌─────────────────────────────────────────────────────────────────────────┐
  │                                                                         │
  │   INCOMING PROTON  ══════════════════════════════════════►              │
  │                         │         │         │                           │
  │                   [collision] [collision] [collision]  ...              │
  │                         │         │         │                           │
  │                        (e⁻)      (e⁻)      (e⁻)   << secondaries >>    │
  │                         │         │         │                           │
  │              [ionize] [excite] [ionize] [excite] ...                    │
  │                         │         │         │                           │
  │               O•  O₂⁺  O•   O₂*  O•   O₃•  O•  ...   << products >>   │
  │                    │              │                                      │
  │               [hop / react]  [hop / react]         << KMC loop >>       │
  │                    │              │                                      │
  │                   O₃            O₂                  << chemistry >>     │
  └─────────────────────────────────────────────────────────────────────────┘
```

### PRIMARY ION TRACK

The proton track is modelled using the **Biersack–Haggmark (1980) MAGIC formula**
for nuclear stopping and the **Green–McNeal (1971)** scaled cross-sections for
electronic (inelastic) stopping.  At each collision site the code decides between
a nuclear elastic event and an inelastic event using the respective cross-sections.

| Quantity | Formula / Reference |
|---|---|
| Nuclear stopping | MAGIC formula — Biersack & Haggmark 1980 |
| Electronic stopping | Green–McNeal 1971; Edgar, Porter & Green 1975 |
| Reduced energy ε | Lindhard–Scharff–Sigmund |
| Impact parameter | Uniform random sampling over collision cylinder |

### SECONDARY ELECTRONS

Each inelastic event spawns a secondary electron.  CIRIS transports each
secondary until its energy falls below the cutoff `ECUTOFF` (default **4.5 eV**).
At every sub-step the electron may:

- **Ionise** the target via the **Green–Sawada (1973)** cross-sections
- **Excite** an allowed transition via the **Porter–Jackman–Green (1976)** formula
- **Excite** a forbidden transition via the **Green–Dutta (1967)** formula

Secondary electron energies are drawn from a **Gamma distribution** (Marsaglia &
Tsang 2000) parameterised by `AVAL`.

### LATTICE & DIFFUSION

The ice is represented as a **3-D integer matrix** of dimension
`NTHICK × NEDGE × NEDGE`.  Each site is either empty (`0`), occupied by a
mobile species (positive index → wait-list entry), or occupied by a frozen
lattice molecule (negative species index).

Mobile species diffuse via a **waiting-time KMC** algorithm (Chang & Herbst 2014):

```
   Rate surface:   b_hop = ν₀ · exp(−E_b · α_surf / T)
   Rate desorb:    b_des = ν₀ · exp(−E_b         / T)
   Rate bulk:      b_bk  = ν₀ · exp(−E_b · α_bulk / T)

   Waiting time:   Δt = −ln(ξ) / (b_hop + b_des)        [surface]
                   Δt = −ln(ξ) / b_bk                    [bulk]
```

where `ξ` is a uniform random deviate, `ν₀ = TRL_NU`, and `α_surf`, `α_bulk`
are the diffusion energy fractions `E_SURF` and `E_BULK`.

---

## ═══════════════════  I N P U T   F I L E S  ═══════════════════════════

CIRIS reads three plain-text input files at run time.

### `species.dat`

One species per line.  Comment lines begin with `!`.

```
! name    Ed (eV)
O         0.56
O2        0.45
O3        0.70
...
```

| Column | Description |
|---|---|
| 1 | Species name (≤ 10 chars).  Anions end in `-` (e.g. `O2-`) |
| 2 | Binding / desorption energy `E_d` in **eV** |

### `reactions.dat`

One reaction per line.  Comment lines begin with `!`.

```
! R1   R2   P1   P2   P3
  O    O    O2   0    0
  O    O2   O3   0    0
  ...
```

| Column | Description |
|---|---|
| 1–2 | Reactant species names |
| 3–5 | Product species names (`0` = no product) |

### `gp.dat`

Key-value pairs used to override module-level defaults at run time.
Unrecognised keys produce a warning and are ignored.

```
TRL_NU          2.6E11
DISPROB         0.0
NSUBEX          3
STEPFAC         0.01
AVAL            13.0
O_O2_BRANCHING  0.5
```

| Parameter | Default | Description |
|---|---|---|
| `TRL_NU` | 2.6×10¹¹ s⁻¹ | Trial / attempt frequency |
| `DISPROB` | 0.0 | Excitative dissociation probability |
| `NSUBEX` | 3 | Sub-excitation interaction count |
| `STEPFAC` | 0.01 | Collision step size factor |
| `AVAL` | 13.0 | Gamma-distribution shape parameter |
| `O_O2_BRANCHING` | 0.5 | O + O₂ → O₃ branching ratio |
| `O2_ION_BRANCHING` | 0.0 | O₂⁻ + O₂⁺ → O₃ + O branching ratio |
| `FRAGILE` | 21 | Species index that dissociates easily |
| `FAST_REACTS` | 21 | Species index that reacts upon formation |

---

## ════════════════════  K E Y   P A R A M E T E R S  ════════════════════

Hard-coded physical conditions live in `parameters.f03` and are recompiled
to change.

```
  ┌──────────────────────────────────────────────────────────────────────┐
  │  PHYSICAL CONDITIONS                                                 │
  │  ─────────────────────────────────────────────────────────────────   │
  │  Ice thickness          THICK     =  2.0 × 10⁻⁶ cm                 │
  │  Lateral edge           EDGE      =  2.0 × 10⁻⁶ cm                 │
  │  Ice density            RHO       =  1.313 × 10²² cm⁻³             │
  │  Kinetic temperature    KIN_TEMP  =  5.0 K                          │
  │  CR / proton flux       CR_FLUX   =  1.0 × 10¹¹ cm⁻² s⁻¹          │
  │  Total fluence          FLUENCE   =  1.0 × 10¹⁶ cm⁻²              │
  │                                                                      │
  │  ION PARAMETERS                                                      │
  │  ─────────────────────────────────────────────────────────────────   │
  │  Initial ion energy     EINIT     =  100 keV                        │
  │  Ion atomic number      ZP        =  1  (proton)                    │
  │  Ion mass               MP        =  1 amu                          │
  │  Target mass            MO2       =  16 amu                         │
  │  Secondary cutoff       ECUTOFF   =  4.5 eV                         │
  │  Primary cutoff         PCUTOFF   =  5.0 eV                         │
  └──────────────────────────────────────────────────────────────────────┘
```

---

## ═══════════════════  O U T P U T   F I L E S  ═══════════════════════

### `abundance.csv`

Main simulation output.  Written every `TIME_FREQ` (default: 100) KMC steps.

```
TIME, FLUENCE, [O], [O3], O3_PROD, O3_DEST, WAIT_LENGTH, N_PROTONS
```

| Column | Units | Description |
|---|---|---|
| `TIME` | s | Elapsed simulation time |
| `FLUENCE` | cm⁻² | Integrated proton fluence |
| `[O]` | molecules cm⁻³ / 10²⁰ | Atomic oxygen column abundance |
| `[O3]` | molecules cm⁻³ / 10²⁰ | Ozone column abundance |
| `O3_PROD` | count | Cumulative ozone production events |
| `O3_DEST` | count | Cumulative ozone destruction events |
| `WAIT_LENGTH` | count | Current length of the KMC wait list |
| `N_PROTONS` | count | Total protons processed |

### `ozone_reactions.wsv`  *(when `O3_ANALYTICS = .TRUE.`)*

Whitespace-delimited log of individual ozone formation and destruction events
for post-processing and reaction-network analysis.

### `reaction_analytics.csv`  *(when `DEBUG = .TRUE.`)*

Verbose per-event log of every reaction, excitation, and ionisation event.
**Warning:** enabling `DEBUG` incurs significant I/O overhead.

### `wait_list.txt`

Snapshot of the full KMC wait list at the end of the simulation run.
Useful for post-run inspection and restart development.

---

## ═════════════════  C O M P I L A T I O N  &  U S A G E  ══════════════

### Requirements

- **gfortran ≥ 9** (or Intel `ifort`)
- GNU `make`

### Build

```bash
# Standard optimised build
make

# Debug build (bounds checking, all warnings)
# Edit Makefile: uncomment FCFLAGS debug line, comment out -O3 line
make clean && make

# Static binary
make static
```

### Run

```bash
./ciris
```

Place `species.dat`, `reactions.dat`, and `gp.dat` in the working directory
before running.  Output files are written to the same directory.

### Compiler switch

To use Intel Fortran, edit the top of `Makefile`:

```makefile
#FC = gfortran
FC = ifort
```

---

## ════════════════════  S O U R C E   F I L E S  ════════════════════════

```
  ciris/
  ├── main.f03          Main program — simulation loop, I/O, initialisation
  ├── parameters.f03    All physical constants and compile-time parameters
  ├── typedefs.f03      Derived types: wait_info, se_info, ionstate, …
  ├── subroutines.f03   Core KMC engine: hopping, reactions, transport
  ├── functiondefs.f03  Physics functions: MAGIC, Green-Sawada, Green-McNeal
  ├── mc_toolbox.f03    Random-number utilities: Gamma, normal distributions
  ├── specdata.f03      Spectroscopic data for O₂ (ionisation, excitation)
  ├── qbert.f03         Reaction-network parser → 3-D reaction cube
  ├── gp.f03            Fitness / objective function (for parameter optimisation)
  ├── bresenham.f03     Bresenham line algorithm (track geometry)
  ├── species.dat       Species list with binding energies
  ├── reactions.dat     Reaction network
  └── gp.dat            Run-time parameter overrides
```

---

## ═══════════════════  R E F E R E N C E S  ══════════════════════════════

```
  ╔══════════════════════════════════════════════════════════════════════╗
  ║  REFERENCES                                                          ║
  ╠══════════════════════════════════════════════════════════════════════╣
  ║                                                                      ║
  ║  Biersack & Haggmark (1980)                                          ║
  ║    A Monte Carlo computer program for the transport of energetic     ║
  ║    ions in amorphous targets.  Nucl. Instrum. Methods 174, 257.      ║
  ║    >> MAGIC formula for nuclear stopping                             ║
  ║                                                                      ║
  ║  Green & Sawada (1973)                                               ║
  ║    Ionization cross sections and secondary electron distributions.   ║
  ║    J. Atmos. Terr. Phys. 34, 1719.                                   ║
  ║    >> Electron-impact ionisation cross-sections                      ║
  ║                                                                      ║
  ║  Porter, Jackman & Green (1976)                                      ║
  ║    Efficiencies for production of atomic nitrogen and oxygen by      ║
  ║    relativistic proton impact.  J. Chem. Phys. 65, 154.             ║
  ║    >> Allowed excitation cross-sections                              ║
  ║                                                                      ║
  ║  Jackman, Garvey & Green (1977)                                      ║
  ║    Electron impacts on atmospheric gases I.                          ║
  ║    J. Geophys. Res. 82, 5081.                                        ║
  ║    >> Ionisation & excitation cross-section parameterisations        ║
  ║                                                                      ║
  ║  Green & Dutta (1967)                                                ║
  ║    Electron impact cross sections for optically forbidden            ║
  ║    transitions.  J. Geophys. Res. 72, 3933.                         ║
  ║    >> Forbidden excitation cross-sections                            ║
  ║                                                                      ║
  ║  Edgar, Porter & Green (1975)                                        ║
  ║    Proton impact cross sections for O₂.                             ║
  ║    Planet. Space Sci. 23, 787.                                       ║
  ║    >> Proton inelastic cross-sections                                ║
  ║                                                                      ║
  ║  Miller & Green (1971)                                               ║
  ║    Proton ionisation cross sections.                                 ║
  ║    >> Green-McNeal scaled proton cross-sections                      ║
  ║                                                                      ║
  ║  Chang & Herbst (2014)                                               ║
  ║    A unified microscopic–macroscopic Monte Carlo simulation of       ║
  ║    gas-grain chemistry in cold dense interstellar clouds.            ║
  ║    ApJ 787, 135.                                                     ║
  ║    >> KMC waiting-time formalism                                     ║
  ║                                                                      ║
  ║  Marsaglia & Tsang (2000)                                            ║
  ║    A simple method for generating gamma variables.                   ║
  ║    ACM TOMS 26, 363.                                                 ║
  ║    >> Gamma-distribution random sampling                             ║
  ║                                                                      ║
  ╚══════════════════════════════════════════════════════════════════════╝
```

---

## ══════════════════════  V E R S I O N   L O G  ═══════════════════════

```
  v1.1  ──  Secondary electron production enabled (SECELEC = .TRUE.)
            Ozone analytics output (O3_ANALYTICS = .TRUE.)
            Performance: precomputed rate constants (~2.5× speedup)
            Bug fixes: degrees→radians in crystal geometry, file-unit
            conflict resolved, memory leaks patched

  v1.0  ──  Initial release
            Kinetic Monte Carlo ice-chemistry engine
            Proton + secondary-electron track simulation
            Oxygen reaction network (O, O₂, O₃, ions)
```

---

```
================================================================================
  Developed by  C. N. Shingledecker
  Fortran 2003  //  GNU General Public License
================================================================================
```
