
# Spatial pharmacokinetic and pharmacodynamic modeling in airway mucus

These scripts were used to study the spatial distribution of drugs in
airway mucus and the spatial antimicrobial pharmacodynamics using
spatial modeling.

### Scripts

The folder scripts contains the following files:

- **01 diffusion_Radius.R**, contains the calculation of diffusion
  coefficients using modified Stokes-Einstein equation.
- **02 tutorial.R**, uses small molecules as an example to showcase the
  spatial pharmacokinetics in the airway mucus (Fig.2)
- **03 impact factor.R**, identifies the influential factors of the
  spatial PK from molecule/particle size (radius, r), mucin binding
  affinity, and half-lives, and under systems with different mucin
  concentration and muco-ciliary clearance (Fig.3, Fig.S2)
- **04 impact dosing.R**, investigates the impact of infusion duration
  and frequency on drug spatial PK in mucus (Fig.4)
- **05 spatialPD.R**, investigates the impact of drug distribution on
  drug efficacy using imipenem as a case study (Fig.5, Fig.S3).
- **06 thin source.R**, compares the spatial pharmacokinetics with thin
  source and infinite source (Fig.S5).

### Functions

Above scripts depend on certain functions written in separate files.

- **pde.R**, contains the discretized partial differential equations for
  describing the spatial pharmacokinetics in airway mucus.
- **pde_pkpd_emax.R**, contains the discretized partial differential
  equations for describing the spatial pharmacokinetics and
  pharmacodynamics in airway mucus.
