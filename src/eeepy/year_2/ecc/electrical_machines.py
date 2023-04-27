#!/usr/bin/env python3
"""Functions and Equations related to electrical machines."""
import math


class ElectricalLoadingEquation:
    """Equation to calculate Electrical Loading (A).

    A = N Z I / (2 π R)
    ---------------------------------------
    A is the elctrical loading.
    N is the total number of slots.
    Z is the number of conductors per slot.
    R is the airgap radius.
    """
    @staticmethod
    def electrical_loading(N: int, Z: int, I: float, R: float) -> float:
        """Calculates the electrical loading (A)."""
        return N * Z * I / (2 * math.pi * R)

    @staticmethod
    def slots(A: float, Z: int, I: float, R: float) -> float:
        """Calculates the number of slots (N)."""
        return A * 2 * math.pi * R / (Z * I)

    @staticmethod
    def cunductors(A: float, N: int, I: float, R: float) -> int:
        """Calculates the number of conductors per slot (Z).

        This value is rounded up.
        """
        return math.ceil(A * 2 * math.pi * R / (N * I))

    @staticmethod
    def current(A: float, N: int, Z: int, R: float) -> int:
        """Calculates the current (I).

        This value is rounded up.
        """
        return math.ceil(A * 2 * math.pi * R / (N * Z))

    @staticmethod
    def airgap_radius(A: float, N: int, Z: int, I: float) -> float:
        """Calculates the airgap radius (R)."""
        return N * Z * I / (A * 2 * math.pi)


class TorqueDeveloped:
    """Equation to calculate the torque developed by an ideal motor.

    T = 2 V B A
    ---------------------------------------------------------------------------
    T is the torque developed.
    V is the volume considered at the stator bore. (V must be calculated at the
    center of the airgap if both airgap and stacking factors are given.)
    B is the magnetic loading.
    A is the current loading.
    """
    @staticmethod
    def torque_developed(V: float, B: float, A: float) -> float:
        """Calculates the torque developed (T)."""
        return 2 * V * B * A

    @staticmethod
    def volume(T: float, B: float, A: float) -> float:
        """Calculates the volume considered at the stator bore (V)."""
        return T / (2 * B * A)

    @staticmethod
    def magnetic_loading(T: float, V: float, A: float) -> float:
        """Calculates the magnetic loading (B)."""
        return T / (2 * V * A)

    @staticmethod
    def current_loading(T: float, V: float, B: float) -> float:
        """Calculates the current loading (A)."""
        return T / (2 * V * B)
