"""Reference implementations of three MPPT algorithms operating on PV voltage.

The algorithms request a PV terminal-voltage setpoint. The baseline plant model
assumes an ideal, sufficiently fast voltage-controlled DC-DC interface. This
isolates MPPT logic; a dynamic averaged boost-converter model is supplied
separately for later validation.
"""
from __future__ import annotations
import numpy as np

class PerturbAndObserve:
    """Classic hill-climbing MPPT: nudge the setpoint and reverse on power loss."""
    def __init__(self, v_min: float, v_max: float, start_v: float,
                 step_v: float = 0.7):
        self.v_min, self.v_max = float(v_min), float(v_max)
        self.start_v, self.step_v = float(start_v), float(step_v)
        self.direction = 1.0
        self.prev_power = None
        self.prev_voltage = None
        self.first = True

    def initial_voltage(self) -> float:
        return float(np.clip(self.start_v, self.v_min, self.v_max))

    def next_voltage(self, voltage_v: float, current_a: float, power_w: float) -> float:
        if self.first:
            # The first sample establishes the baseline power and then makes the
            # first deterministic perturbation away from the initial setpoint.
            self.first = False
            self.prev_power = float(power_w)
            self.prev_voltage = float(voltage_v)
            return float(np.clip(self.start_v + self.direction * self.step_v,
                                 self.v_min, self.v_max))
        # If the last perturbation reduced power, the previous step moved away
        # from the maximum-power point and the search direction must reverse.
        if power_w < self.prev_power:
            self.direction *= -1.0
        self.prev_power = float(power_w)
        self.prev_voltage = float(voltage_v)
        return float(np.clip(voltage_v + self.direction * self.step_v,
                             self.v_min, self.v_max))

    def reset(self, start_v: float | None = None):
        self.direction = 1.0
        self.prev_power = None
        self.prev_voltage = None
        self.first = True
        if start_v is not None:
            self.start_v = float(start_v)


class IncrementalConductance:
    """Uses dP/dV = I + V*dI/dV to know which side of the MPP the operating point is on."""
    def __init__(self, v_min: float, v_max: float, start_v: float,
                 step_v: float = 0.7, deadband_a_per_v: float = 0.005,
                 dv_epsilon: float = 1e-6):
        self.v_min, self.v_max = float(v_min), float(v_max)
        self.start_v, self.step_v = float(start_v), float(step_v)
        self.deadband = float(deadband_a_per_v)
        self.dv_epsilon = float(dv_epsilon)
        self.prev_v = None
        self.prev_i = None
        self.first = True

    def initial_voltage(self) -> float:
        return float(np.clip(self.start_v, self.v_min, self.v_max))

    def next_voltage(self, voltage_v: float, current_a: float, power_w: float) -> float:
        if self.first:
            # Seed the derivative estimate with the first current/voltage sample.
            self.first = False
            self.prev_v, self.prev_i = float(voltage_v), float(current_a)
            return float(np.clip(self.start_v + self.step_v, self.v_min, self.v_max))
        dv = float(voltage_v - self.prev_v)
        di = float(current_a - self.prev_i)
        self.prev_v, self.prev_i = float(voltage_v), float(current_a)
        if abs(dv) < self.dv_epsilon:
            # Startup/saturation: use current change as a conservative direction cue
            # because the derivative estimate is numerically unstable at zero voltage swing.
            if abs(di) < 1e-9:
                return float(np.clip(voltage_v, self.v_min, self.v_max))
            delta = self.step_v if di > 0 else -self.step_v
            return float(np.clip(voltage_v + delta, self.v_min, self.v_max))
        # The slope test compares the instantaneous conductance (dI/dV) to the
        # array's conductance (I/V). The sign of dP/dV tells us whether to raise
        # or lower the operating voltage to reach the maximum-power point.
        slope = di / dv
        condition = slope + current_a / max(abs(voltage_v), 1e-9)
        if abs(condition) <= self.deadband:
            nxt = voltage_v
        elif condition > 0:
            # dP/dV > 0: move to higher voltage.
            nxt = voltage_v + self.step_v
        else:
            nxt = voltage_v - self.step_v
        return float(np.clip(nxt, self.v_min, self.v_max))

    def reset(self, start_v: float | None = None):
        self.prev_v = self.prev_i = None
        self.first = True
        if start_v is not None:
            self.start_v = float(start_v)


class ParticleSwarmMPPT:
    """Sequential-evaluation PSO with optional reinitialization after power drops."""
    def __init__(self, v_min: float, v_max: float, start_v: float,
                 particles: int = 6, w: float = 0.55, c1: float = 1.4,
                 c2: float = 1.4, max_iterations: int = 12,
                 seed: int = 0, restart_drop_fraction: float = 0.18):
        if particles < 3:
            raise ValueError("Use at least three particles for this PSO implementation.")
        self.v_min, self.v_max = float(v_min), float(v_max)
        self.start_v = float(start_v)
        self.n = int(particles)
        self.w, self.c1, self.c2 = float(w), float(c1), float(c2)
        self.max_iterations = int(max_iterations)
        self.restart_drop_fraction = float(restart_drop_fraction)
        self.rng = np.random.default_rng(seed)
        self.seed = int(seed)
        self.max_velocity = 0.25 * (self.v_max - self.v_min)
        self._initialize_swarm()
        self.locked = False
        self.locked_power = None
        self.last_voltage = None

    def _initialize_swarm(self):
        # A PSO run starts with a random cloud of candidate setpoints across the
        # allowable voltage range. The swarm then explores the search space and
        # retains the best candidate found so far.
        self.positions = self.rng.uniform(self.v_min, self.v_max, self.n)
        self.velocities = np.zeros(self.n, dtype=float)
        self.personal_positions = self.positions.copy()
        self.personal_scores = np.full(self.n, -np.inf)
        self.global_position = float(self.start_v)
        self.global_score = -np.inf
        self.particle_index = 0
        self.iteration = 0
        self.locked = False
        self.locked_power = None

    def reset(self, start_v: float | None = None):
        if start_v is not None:
            self.start_v = float(start_v)
        self._initialize_swarm()

    def initial_voltage(self) -> float:
        return float(self.positions[0])

    def next_voltage(self, voltage_v: float, current_a: float, power_w: float) -> float:
        p = float(max(0.0, power_w))
        # If the swarm had already converged and the PV conditions suddenly change,
        # we reset the search to avoid the algorithm staying trapped at the old MPP.
        if self.locked and self.locked_power is not None and p < (1.0 - self.restart_drop_fraction) * self.locked_power:
            self._initialize_swarm()
            # This sample belonged to the old converged point; do not score it as
            # the first particle in the restarted swarm.
            return float(self.positions[0])
        if self.locked:
            self.locked_power = p
            return float(np.clip(self.global_position, self.v_min, self.v_max))
        j = self.particle_index
        # PSO fitness belongs to the candidate setpoint being evaluated, not
        # necessarily the measured terminal voltage. The distinction matters
        # when a dynamic converter cannot track the setpoint instantaneously.
        candidate_v = float(self.positions[j])
        if p > self.personal_scores[j]:
            self.personal_scores[j] = p
            self.personal_positions[j] = candidate_v
        if p > self.global_score:
            self.global_score = p
            self.global_position = candidate_v
        self.particle_index += 1
        if self.particle_index >= self.n:
            # One full swarm sweep is complete; the particles move toward their
            # personal bests and the best point found by the whole swarm.
            self.particle_index = 0
            self.iteration += 1
            r1 = self.rng.random(self.n)
            r2 = self.rng.random(self.n)
            self.velocities = (self.w * self.velocities
                               + self.c1 * r1 * (self.personal_positions - self.positions)
                               + self.c2 * r2 * (self.global_position - self.positions))
            self.velocities = np.clip(self.velocities, -self.max_velocity, self.max_velocity)
            self.positions = np.clip(self.positions + self.velocities, self.v_min, self.v_max)
            if self.iteration >= self.max_iterations:
                self.locked = True
                self.locked_power = self.global_score
                return float(np.clip(self.global_position, self.v_min, self.v_max))
        return float(self.positions[self.particle_index])


def make_algorithm(name: str, v_min: float, v_max: float, start_v: float,
                   *, po_step: float = 0.7, inc_deadband: float = 0.005,
                   pso_particles: int = 6, pso_iterations: int = 12,
                   seed: int = 0):
    key = "".join(ch for ch in name.lower() if ch.isalnum())
    if key in {"po", "perturbandobserve"}:
        return PerturbAndObserve(v_min, v_max, start_v, step_v=po_step)
    if key in {"inccond", "incrementalconductance"}:
        return IncrementalConductance(v_min, v_max, start_v, step_v=po_step,
                                      deadband_a_per_v=inc_deadband)
    if key in {"pso", "particleswarmoptimization"}:
        return ParticleSwarmMPPT(v_min, v_max, start_v, particles=pso_particles,
                                 max_iterations=pso_iterations, seed=seed)
    raise ValueError(f"Unknown algorithm: {name}")
