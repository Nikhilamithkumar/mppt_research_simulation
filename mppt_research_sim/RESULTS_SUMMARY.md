# Results summary (current run)

## PV module model check

The single-diode KC200GT model gives: Isc 8.2096 A, Voc 32.8834 V, Vmp 26.3488 V, Imp 7.5956 A, and Pmp 200.1357 W at 1000 W/m² and 25 °C. These are close to the published STC points (8.21 A, 32.9 V, 26.3 V, 7.61 A, 200 W).

## Theoretical MPPs for the three-module string

- Uniform 1000/1000/1000 W/m²: about 600.41 W.
- Uniform 600/600/600 W/m²: about 356.93 W.
- Whole-module irradiance 1000/1000/300 W/m²: global peak about 384.34 W at 50.74 V; local peak about 207.63 W at 86.64 V.
- One of nine 18-cell substrings at 300 W/m² and the other eight at 1000 W/m²: global peak about 528.38 W; the P–V curve is multi-peaked.

## Ideal voltage-actuation comparison

| Scenario | Algorithm | Efficiency | Final power | Global peak found |
|---|---|---:|---:|---|
| 1000→600 W/m² step | P&O | 99.48% | 356.80 W | Yes |
| 1000→600 W/m² step | Incremental Conductance | 99.50% | 356.83 W | Yes |
| 1000→600 W/m² step | PSO | 96.46% | 356.91 W | Yes |
| Whole-module partial shading | P&O | 53.97% | 207.47 W | No |
| Whole-module partial shading | Incremental Conductance | 53.96% | 207.42 W | No |
| Whole-module partial shading | PSO | 97.05% | 384.34 W | Yes |
| One substring shaded | P&O | 42.05% | 222.25 W | No |
| One substring shaded | Incremental Conductance | 42.05% | 222.25 W | No |
| One substring shaded | PSO | 98.22% | 528.37 W | Yes |

Thirty PSO random-seed repetitions in the whole-module shading case found the global peak in 30/30 runs. Mean efficiency was 97.30%, with a sample standard deviation of 0.72%. This empirical success rate is limited to this selected configuration and should not be generalized to all PSO designs.

## Averaged boost converter

The separately simulated averaged CCM boost-converter run gave approximately:

| Scenario | Algorithm | Efficiency | Final PV power | Final MPP reference |
|---|---|---:|---:|---:|
| Irradiance step | P&O | 95.81% | 356.88 W | 356.93 W |
| Irradiance step | Incremental Conductance | 95.83% | 356.88 W | 356.93 W |
| Irradiance step | PSO | 89.24% | 356.88 W | 356.93 W |
| Whole-module partial shading | P&O | 53.94% | 207.63 W | 384.34 W |
| Whole-module partial shading | Incremental Conductance | 53.98% | 207.63 W | 384.34 W |
| Whole-module partial shading | PSO | 94.80% | 384.12 W | 384.34 W |

The averaged converter results are a simulated controller/model combination, not hardware validation. They are separated from the ideal-actuation comparison because converter dynamics change the evaluation conditions.
