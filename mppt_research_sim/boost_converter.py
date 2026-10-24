"""Continuous-conduction averaged boost converter model with PV voltage control.

Averaged equations (ideal switch/diode, with a small inductor resistance):
  C_in * dV_pv/dt = I_pv(V_pv) - I_L
  L * dI_L/dt = V_pv - (1-D)*V_out - r_L*I_L
  C_out * dV_out/dt = (1-D)*I_L - V_out/R_load

The controller is an outer MPPT voltage reference plus a feed-forward and
PD/slow-integral duty controller. This is an averaged model: it does not resolve
switching ripple or semiconductor losses.
"""
from __future__ import annotations
from dataclasses import dataclass
import math
import numpy as np
from pv_model import IVCurve

@dataclass(frozen=True)
class BoostParameters:
    inductance_h: float = 1.5e-3
    input_capacitance_f: float = 2.2e-3
    output_capacitance_f: float = 1.0e-3
    load_resistance_ohm: float = 64.0
    inductor_resistance_ohm: float = 0.05
    duty_min: float = 0.03
    duty_max: float = 0.85
    kp_duty_per_v: float = 3.0e-4
    kd_duty_s_per_v: float = 5.0e-6
    ki_duty_per_vs: float = 0.01

class AveragedBoostConverter:
    def __init__(self, curve: IVCurve, initial_pv_voltage_v: float,
                 params: BoostParameters | None = None, dt_s: float = 25e-6):
        self.params = params or BoostParameters()
        self.dt_s = float(dt_s)
        if self.dt_s <= 0:
            raise ValueError("dt_s must be positive")
        self.state = self._steady_state(curve, initial_pv_voltage_v)
        self.integral_error = 0.0
        self.last_duty = None

    def _steady_state(self, curve: IVCurve, voltage_v: float):
        # Start the converter in a physically meaningful operating point by
        # selecting a PV voltage near a realistic operating region and computing
        # the corresponding load-side output voltage estimate.
        v = float(np.clip(voltage_v, max(0.5, 0.01*curve.voc_v), 0.98*curve.voc_v))
        i = max(0.0, float(curve.current(v)))
        p = max(1.0, v*i)
        r = self.params.load_resistance_ohm
        r_l = self.params.inductor_resistance_ohm
        vout_sq = max(25.0**2, (p-r_l*i*i)*r)
        vout = math.sqrt(vout_sq)
        return np.array([v, i, vout], dtype=float)

    def measured(self, curve: IVCurve) -> tuple[float, float, float, float]:
        v, il, vout = (float(x) for x in self.state)
        ipv = max(0.0, float(curve.current(v)))
        return v, ipv, v*ipv, vout

    def _rhs(self, state, duty: float, curve: IVCurve):
        # These equations are the continuous-time averaged model of the boost
        # converter: input capacitor dynamics, inductor current dynamics, and output
        # capacitor dynamics. The duty cycle is the controlled switching ratio.
        p = self.params
        vpv, il, vout = (float(x) for x in state)
        vpv_q = float(np.clip(vpv, 0.0, 1.15*curve.voc_v))
        il_q = max(0.0, il)
        vout_q = max(1.0, vout)
        ipv = max(0.0, float(curve.current(vpv_q)))
        dvpv = (ipv - il_q) / p.input_capacitance_f
        dil = (vpv_q - (1.0-duty)*vout_q - p.inductor_resistance_ohm*il_q) / p.inductance_h
        dvout = ((1.0-duty)*il_q - vout_q/p.load_resistance_ohm) / p.output_capacitance_f
        return np.array([dvpv, dil, dvout], dtype=float)

    def advance(self, pv_voltage_reference_v: float, curve: IVCurve, duration_s: float):
        # Advance the averaged converter over a finite time interval by solving the
        # plant model with a controller that compares the current PV voltage to the
        # MPPT setpoint. This step represents the dynamic response rather than the
        # instantaneous actuation used in the ideal benchmark.
        p = self.params
        n_steps = max(1, int(round(duration_s / self.dt_s)))
        dt = duration_s / n_steps
        vref = float(np.clip(pv_voltage_reference_v, 0.0, curve.voc_v * 1.02))
        for _ in range(n_steps):
            vpv, il, vout = (float(x) for x in self.state)
            ipv = max(0.0, float(curve.current(max(0.0, vpv))))
            dv_est = (ipv - max(0.0, il)) / p.input_capacitance_f
            err = vpv - vref
            self.integral_error = float(np.clip(self.integral_error + err*dt, -10.0, 10.0))
            # Feed-forward gives the nominal boost duty; feedback compensates for
            # PV input-capacitor dynamics. Duty increases when Vpv is too high.
            vout_safe = max(vout, 20.0)
            d_ff = 1.0 - (vref - p.inductor_resistance_ohm*max(0.0, il))/vout_safe
            duty = float(np.clip(d_ff + p.kp_duty_per_v*err
                                 + p.kd_duty_s_per_v*dv_est
                                 + p.ki_duty_per_vs*self.integral_error,
                                 p.duty_min, p.duty_max))
            x = self.state
            k1 = self._rhs(x, duty, curve)
            k2 = self._rhs(x + 0.5*dt*k1, duty, curve)
            k3 = self._rhs(x + 0.5*dt*k2, duty, curve)
            k4 = self._rhs(x + dt*k3, duty, curve)
            x = x + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4)
            # Enforce physical state bounds for the averaged model.
            x[0] = float(np.clip(x[0], 0.0, 1.15*curve.voc_v))
            x[1] = float(np.clip(x[1], 0.0, 1.5*max(curve.current_a[0], 1.0)))
            x[2] = float(np.clip(x[2], 1.0, 250.0))
            self.state = x
            self.last_duty = duty
        return self.state.copy()
