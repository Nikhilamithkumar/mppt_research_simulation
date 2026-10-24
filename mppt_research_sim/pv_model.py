"""Physics-based single-diode PV module and series-string model.

Reference module: Kyocera KC200GT (54 series cells), with commonly used
single-diode parameter values calibrated to its STC datasheet ratings.
The array model supports different irradiance per module or per 18-cell
substring, with an idealized bypass diode across each substring.

This is an engineering model, not a claim of manufacturer-certified accuracy.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
from typing import Iterable
import numpy as np
from scipy.optimize import brentq, minimize_scalar

Q = 1.602176634e-19       # C
K_B = 1.380649e-23        # J/K
K_B_EV = 8.617333262e-5   # eV/K
T_REF_K = 298.15
G_REF = 1000.0
# KC200GT open-circuit-voltage temperature coefficient reported in published
# module model references. The saturation-current temperature law below is
# calibrated to this coefficient, rather than assuming an unverified band-gap law.
BETA_VOC_V_PER_K = -0.123

@dataclass(frozen=True)
class PVModule:
    name: str = "Kyocera KC200GT"
    cells_in_series: int = 54
    p_mp_w: float = 200.0
    v_mp_v: float = 26.3
    i_mp_a: float = 7.61
    v_oc_v: float = 32.9
    i_sc_a: float = 8.21
    alpha_sc_a_per_k: float = 0.00318
    beta_voc_v_per_k: float = BETA_VOC_V_PER_K
    # Single-diode parameter set at 1000 W/m², 25 °C.
    il_ref_a: float = 8.214
    i0_ref_a: float = 9.825e-8
    rs_ohm: float = 0.221
    rsh_ref_ohm: float = 415.405
    diode_ideality: float = 1.3

    @property
    def a_ref_v(self) -> float:
        """n*N_series*k*T/q, the modified ideality/thermal voltage (V)."""
        return self.diode_ideality * self.cells_in_series * K_B * T_REF_K / Q

    def parameters_at(self, irradiance_w_m2: float, temp_c: float) -> dict[str, float]:
        """Return single-diode parameters with datasheet-calibrated temperature dependence."""
        # The PV model uses the common single-diode form: 
        # I = IL - I0*(exp((V + I*Rs)/a)-1) - (V + I*Rs)/Rsh.
        # This method scales the photocurrent, diode saturation current, and shunt
        # resistance to the current irradiance and cell temperature so the same
        # module equation can be solved across a wide operating range.
        G = float(max(0.0, irradiance_w_m2))
        T = float(temp_c + 273.15)
        if T <= 0:
            raise ValueError("Cell temperature must be above absolute zero.")
        t_ratio = T / T_REF_K
        il_at_ref_irradiance = self.il_ref_a + self.alpha_sc_a_per_k * (T - T_REF_K)
        photocurrent = (G / G_REF) * il_at_ref_irradiance
        a = self.a_ref_v * t_ratio
        rsh = self.rsh_ref_ohm * G_REF / max(G, 1.0)
        # Empirical one-diode temperature correction calibrated to the published
        # open-circuit-voltage coefficient. At the reference temperature this
        # ratio is 1; at other temperatures it adjusts I0 so Voc follows beta_Voc.
        voc_target = max(1.0, self.v_oc_v + self.beta_voc_v_per_k * (T - T_REF_K))
        i0_cal = max(il_at_ref_irradiance - voc_target / self.rsh_ref_ohm, 1e-12) / np.expm1(voc_target / a)
        i0_cal_ref = max(self.il_ref_a - self.v_oc_v / self.rsh_ref_ohm, 1e-12) / np.expm1(self.v_oc_v / self.a_ref_v)
        saturation = self.i0_ref_a * (i0_cal / i0_cal_ref)
        return {"IL": float(photocurrent), "I0": float(saturation),
                "Rs": self.rs_ohm, "Rsh": float(rsh), "a": float(a),
                "G": G, "T_K": T}

    def current_at_voltage(self, voltage_v, irradiance_w_m2: float, temp_c: float):
        """Solve the implicit single-diode equation for I(V) with Newton iteration.

        `voltage_v` may be a scalar or NumPy array. Supports modest negative
        module voltages so an external bypass-diode clamp can be represented.
        """
        x = np.asarray(voltage_v, dtype=float)
        p = self.parameters_at(irradiance_w_m2, temp_c)
        if p["G"] <= 0:
            out = np.zeros_like(x)
            return float(out) if out.ndim == 0 else out
        IL, I0, Rs, Rsh, a = p["IL"], p["I0"], p["Rs"], p["Rsh"], p["a"]
        # Good initial value on the PV forward quadrant; Newton then corrects it.
        cur = np.maximum(-0.02, IL - x / Rsh)
        for _ in range(30):
            z = np.clip((x + cur * Rs) / a, -100.0, 80.0)
            ez = np.exp(z)
            f = cur - IL + I0 * np.expm1(z) + (x + cur * Rs) / Rsh
            df = 1.0 + I0 * ez * Rs / a + Rs / Rsh
            step = f / df
            cur_next = cur - step
            if np.all(np.abs(cur_next - cur) < 1e-11):
                cur = cur_next
                break
            cur = cur_next
        out = np.asarray(cur, dtype=float)
        return float(out) if out.ndim == 0 else out

    def voc(self, irradiance_w_m2: float, temp_c: float) -> float:
        if irradiance_w_m2 <= 0:
            return 0.0
        return float(brentq(lambda v: self.current_at_voltage(v, irradiance_w_m2, temp_c),
                            0.0, 1.5 * self.v_oc_v, xtol=1e-10))

    def voltage_at_current(self, current_a, irradiance_w_m2: float, temp_c: float,
                           bypass_v: float | None = None):
        """Invert I(V) by vectorized bisection; optionally clamp at -bypass_v."""
        # For a given operating current, the model finds the corresponding PV
        # terminal voltage. This is used to build the full I-V curve by sweeping
        # current rather than voltage.
        target = np.asarray(current_a, dtype=float)
        if irradiance_w_m2 <= 0:
            out = np.where(target <= 0, 0.0, -(bypass_v or 0.0))
            return float(out) if out.ndim == 0 else out
        v_high = self.voc(irradiance_w_m2, temp_c)
        low = np.full_like(target, -(bypass_v if bypass_v is not None else 0.0), dtype=float)
        high = np.full_like(target, v_high, dtype=float)
        i_at_low = self.current_at_voltage(low, irradiance_w_m2, temp_c)
        if bypass_v is None:
            # No bypass diode: let the inverse curve continue into negative voltage.
            bypass_mask = np.zeros_like(target, dtype=bool)
        else:
            bypass_mask = target > i_at_low
        # Root exists for targets between I(Voc)=0 and I(-Vbypass).
        solve_mask = ~bypass_mask
        for _ in range(58):
            mid = (low + high) / 2.0
            imid = self.current_at_voltage(mid, irradiance_w_m2, temp_c)
            # I(V) decreases with voltage: if I(mid) > requested I, raise voltage.
            move_low = (imid > target) & solve_mask
            move_high = (~move_low) & solve_mask
            low = np.where(move_low, mid, low)
            high = np.where(move_high, mid, high)
        out = (low + high) / 2.0
        if bypass_v is not None:
            out = np.where(bypass_mask, -float(bypass_v), out)
        return float(out) if out.ndim == 0 else out

@dataclass
class IVCurve:
    voltage_v: np.ndarray  # strictly increasing, normally [0,Voc]
    current_a: np.ndarray
    irradiances_w_m2: tuple[float, ...]
    temp_c: float
    bypass_v: float

    @property
    def power_w(self) -> np.ndarray:
        return self.voltage_v * self.current_a

    @property
    def voc_v(self) -> float:
        return float(self.voltage_v[-1])

    def current(self, voltage_v):
        """Interpolate the array I-V curve; current is clipped to its valid range."""
        v = np.asarray(voltage_v, dtype=float)
        out = np.interp(v, self.voltage_v, self.current_a,
                        left=self.current_a[0], right=0.0)
        out = np.where(v < 0, self.current_a[0], out)
        return float(out) if out.ndim == 0 else out

    def power(self, voltage_v):
        v = np.asarray(voltage_v, dtype=float)
        out = v * self.current(v)
        return float(out) if out.ndim == 0 else out

    def mpps(self) -> list[dict[str, float]]:
        """Find local peaks with bounded refinement; returns peaks sorted by power."""
        p = self.power_w
        peaks: list[dict[str, float]] = []
        for j in range(1, len(p) - 1):
            if p[j] >= p[j-1] and p[j] >= p[j+1] and p[j] > 1e-8:
                lo, hi = self.voltage_v[j-1], self.voltage_v[j+1]
                opt = minimize_scalar(lambda v: -self.power(v), bounds=(lo, hi), method="bounded",
                                      options={"xatol": 1e-9})
                vmp = float(opt.x)
                peaks.append({"voltage_v": vmp, "current_a": self.current(vmp),
                              "power_w": self.power(vmp)})
        if not peaks:
            j = int(np.argmax(p))
            peaks.append({"voltage_v": float(self.voltage_v[j]),
                          "current_a": float(self.current_a[j]), "power_w": float(p[j])})
        # De-duplicate neighbors/refinement duplicates, and order by power descending.
        uniq: list[dict[str, float]] = []
        for item in sorted(peaks, key=lambda d: d["power_w"], reverse=True):
            if not any(abs(item["voltage_v"] - q["voltage_v"]) < 0.1 for q in uniq):
                uniq.append(item)
        return uniq

    def mpp(self) -> dict[str, float]:
        return self.mpps()[0]


def make_substring_string_iv_curve(substring_irradiances_w_m2: Iterable[float],
                                   temp_c: float = 25.0,
                                   module: PVModule | None = None,
                                   bypass_v: float = 0.7,
                                   substrings_per_module: int = 3,
                                   points: int = 8000) -> IVCurve:
    """Build a series-string curve using one bypass diode per cell substring.

    Each PVModule is decomposed into equal series-cell substrings. Series
    resistance, shunt resistance, open-circuit voltage, MPP voltage, and Voc
    temperature coefficient are scaled with substring length; photocurrent and
    diode saturation current stay the same for series-connected cells. Each
    substring has an idealized bypass clamp at -bypass_v volts.

    `substring_irradiances_w_m2` contains one irradiance value for every
    substring in the entire string. For three 54-cell modules with three
    substrings per module, provide nine values.
    """
    # The string is assembled by summing each substring's voltage contribution at
    # a common current. This creates multiple local maxima when one substring is
    # heavily shaded while the rest remain well illuminated.
    module = module or PVModule()
    gs = tuple(float(g) for g in substring_irradiances_w_m2)
    if not gs or any(g < 0 for g in gs):
        raise ValueError("Provide one or more non-negative substring irradiances.")
    if substrings_per_module < 1 or module.cells_in_series % substrings_per_module != 0:
        raise ValueError("substrings_per_module must divide module.cells_in_series.")
    fraction = 1.0 / float(substrings_per_module)
    substring_module = replace(
        module,
        cells_in_series=module.cells_in_series // substrings_per_module,
        p_mp_w=module.p_mp_w * fraction,
        v_mp_v=module.v_mp_v * fraction,
        v_oc_v=module.v_oc_v * fraction,
        rs_ohm=module.rs_ohm * fraction,
        rsh_ref_ohm=module.rsh_ref_ohm * fraction,
        beta_voc_v_per_k=module.beta_voc_v_per_k * fraction,
    )
    # Sweep the string current; all series substrings carry the same current.
    il_values = [substring_module.parameters_at(g, temp_c)["IL"] for g in gs]
    i_max = max(il_values) * 1.035
    i_grid = np.linspace(0.0, i_max, int(points))
    string_v = np.zeros_like(i_grid)
    for g in gs:
        string_v += substring_module.voltage_at_current(i_grid, g, temp_c,
                                                         bypass_v=bypass_v)
    keep = string_v >= 0.0
    if not np.any(keep):
        raise RuntimeError("PV string curve has no non-negative voltage points.")
    v = string_v[keep]
    i = i_grid[keep]
    # Estimate short-circuit current at V=0 from the last crossing of V(I).
    crossing = np.flatnonzero((string_v[:-1] >= 0.0) & (string_v[1:] < 0.0))
    if len(crossing):
        j0 = int(crossing[-1])
        v0, v1 = string_v[j0], string_v[j0 + 1]
        i0, i1 = i_grid[j0], i_grid[j0 + 1]
        isc_array = float(i0 + (0.0 - v0) * (i1 - i0) / (v1 - v0))
        v = np.append(v, 0.0)
        i = np.append(i, isc_array)
    if len(v) < 2:
        raise RuntimeError("PV curve has too few points.")
    order = np.argsort(v)
    v = v[order]
    i = i[order]
    v_unique, first = np.unique(v, return_index=True)
    i_unique = i[first]
    if v_unique[0] > 1e-9:
        raise RuntimeError("PV curve did not include the short-circuit endpoint.")
    return IVCurve(v_unique, i_unique, gs, float(temp_c), float(bypass_v))


def make_string_iv_curve(irradiances_w_m2: Iterable[float], temp_c: float = 25.0,
                         module: PVModule | None = None, bypass_v: float = 0.7,
                         points: int = 8000, substrings_per_module: int = 3) -> IVCurve:
    """Build a series-string I-V curve from one irradiance value per module.

    Each module is represented by `substrings_per_module` series substrings
    with an idealized bypass diode per substring. Module-level irradiances are
    repeated over that module's substrings; for within-module shading use
    `make_substring_string_iv_curve` and pass separate values for each substring.
    """
    module = module or PVModule()
    module_gs = tuple(float(g) for g in irradiances_w_m2)
    if not module_gs or any(g < 0 for g in module_gs):
        raise ValueError("Provide one or more non-negative module irradiances.")
    # Module-level irradiance values are repeated across each substring in the
    # module so a user can specify one number per module while the lower-level
    # model still sees one irradiance value per shaded substring.
    flat_gs = tuple(g for g in module_gs for _ in range(substrings_per_module))
    curve = make_substring_string_iv_curve(flat_gs, temp_c, module=module,
                                           bypass_v=bypass_v,
                                           substrings_per_module=substrings_per_module,
                                           points=points)
    # Preserve the user-facing module irradiance profile, not the flattened list.
    curve.irradiances_w_m2 = module_gs
    return curve


def validate_module_stc(module: PVModule | None = None) -> dict[str, float]:
    """Evaluate one module at STC and return datasheet-vs-model check values."""
    module = module or PVModule()
    voc = module.voc(1000.0, 25.0)
    v = np.linspace(0.0, voc, 10000)
    i = module.current_at_voltage(v, 1000.0, 25.0)
    p = v * i
    idx = int(np.argmax(p))
    return {"model_isc_a": float(module.current_at_voltage(0.0, 1000.0, 25.0)),
            "model_voc_v": voc, "model_vmp_v": float(v[idx]),
            "model_imp_a": float(i[idx]), "model_pmp_w": float(p[idx]),
            "datasheet_isc_a": module.i_sc_a, "datasheet_voc_v": module.v_oc_v,
            "datasheet_vmp_v": module.v_mp_v, "datasheet_imp_a": module.i_mp_a,
            "datasheet_pmp_w": module.p_mp_w}

if __name__ == "__main__":
    import json
    print(json.dumps(validate_module_stc(), indent=2))
    for name, gs in [("uniform_1000", [1000, 1000, 1000]),
                     ("uniform_600", [600, 600, 600]),
                     ("partial_shading", [1000, 1000, 300])]:
        c = make_string_iv_curve(gs, 25, points=12000)
        print(name, json.dumps(c.mpps(), indent=2))
