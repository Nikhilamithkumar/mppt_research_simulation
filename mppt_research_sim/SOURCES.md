# Sources and provenance

1. pvlib Python documentation, “Single diode models”: https://pvlib-python.readthedocs.io/en/stable/user_guide/modeling_topics/singlediode.html
   - Explains the single-diode equation, temperature/irradiance parameter models, and numerical solutions for I–V curves.
2. Ndegwa et al., “A Fast and Accurate Analytical Method for Parameter Determination of a Photovoltaic System Based on Manufacturer’s Data” (2020): https://onlinelibrary.wiley.com/doi/10.1155/2020/7580279
   - Reports model parameter sets and KC200GT performance checks.
3. “A New Simplified Five-Parameter Estimation Method for Single-Diode Model of Photovoltaic Panels,” Energies (2019): https://www.mdpi.com/1996-1073/12/22/4271
   - Reports KC200GT STC values: Isc 8.21 A, Voc 32.9 V, Vmp 26.3 V, Imp 7.61 A and 54 cells, plus temperature coefficients.
4. PV Performance Modeling Collaborative, Single-Diode Equivalent Circuit: https://pvpmc.sandia.gov/modeling-guide/2-dc-module-iv/single-diode-equivalent-circuit-models/
   - Background on the standard single-diode module equivalent circuit.

The exact parameter set is intended to reproduce the nameplate characteristics closely. It is not a manufacturer-supplied certificate of every off-STC behavior.

5. “Variation of the performance of a PV panel with the number of bypass diodes and partial shading patterns”: https://scispace.com/pdf/variation-of-the-performance-of-a-pv-panel-with-the-number-33rrsxfc3x.pdf
   - Describes the KC200GT 54-cell module as three strings of 18 cells with a bypass diode across each string. The code models each substring separately.
