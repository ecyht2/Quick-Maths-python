#!/usr/bin/env python3
"""Functions and equations related to DC machines."""
import math


class PoleFlux:
    """The total flux of the poles.

    ɸ_p = B l 2πR / (2P) λ_f / λ_p
    ----------------------------------------------
    ɸ_p is the flux of the poles.
    B is the flux density of the poles.
    l is the length of the
    R is the radius of the poles.
    P is the number of pole pairs.
    λ_f is the salient field winding tooth (pole).
    λ_p is the pole pitch.
    """
    @staticmethod
    def flux(B: float, R: float, P: int, l: float,
             lambda_f: float, lambda_p: float) -> float:
        """Calculates the total flux of the poles (ɸ_p).

        :param ɸ_p: The flux of the poles.
        :param B: is the flux density of the poles.
        :param l: is the length of the
        :param R: is the radius of the poles.
        :param P: is the number of pole pairs.
        :param λ_f: is the salient field winding tooth (pole).
        :param λ_p: is the pole pitch.
        """
        return B * l * math.pi * R / P * lambda_f / lambda_p

    @staticmethod
    def flux_poles(B: float, R: float, p: int, l: float,
                   lambda_f: float, lambda_p: float) -> float:
        """Calculates the total flux of the poles (ɸ_p).

        This method differs from flux as it uses the number of poles (p)
        instead of the number of pole pairs (P).

        :param ɸ_p: The flux of the poles.
        :param B: is the flux density of the poles.
        :param l: is the length of the
        :param R: is the radius of the poles.
        :param p: is the number of poles.
        :param λ_f: is the salient field winding tooth (pole).
        :param λ_p: is the pole pitch.
        """
        return B * l * 2 * math.pi * R / p * lambda_f / lambda_p


class InducedTorque:
    """Induced torque in an ideal DC machine.

    T = ɸ_p N Z I_a / (2 π)
    -------------------------------------
    T is the induced torque.
    N is the number of slots.
    Z is the number of windings per slot.
    I_a is the armature current.
    ɸ_p is the flux of the poles.

    From Electrical Loading
    T = 2 V B A λ_f / λ_p
    -------------------------------------
    T is the induced torque.
    V is the volume considered at the stator bore. (V must be calculated at the
    center of the airgap if both airgap and stacking factors are given.)
    B is the magnetic loading.
    A is the current loading.
    λ_f is the salient field winding tooth (pole).
    λ_p is the pole pitch.
    """
    @staticmethod
    def torque(N: int, Z: int, I: float, flux: float) -> float:
        """Calcultes the induced torque (T).

        :param N: The number of slots.
        :param Z: The number of windings per slot.
        :param I: The armature current (I_a).
        :param flux: is the flux of the poles (ɸ_p).
        :returns: The torque induced (T).
        """
        return flux * N * Z * I / 2 / math.pi

    @staticmethod
    def torque_electrical_loading(V: float, B: float, A: float,
                                  lambda_f: float, lambda_p: float) -> float:
        """Calcultes the induced using electrical loading torque (T).

        :param V: The volume considered at the stator bore. (V must be
        calculated at the center of the airgap if both airgap and stacking
        factors are given.).
        :param B: The magnetic loading.
        :param A: The the current loading.
        :param λ_f: The salient field winding tooth (pole).
        :param λ_p: The pole pitch.
        :returns: The torque induced.
        """
        return 2 * V * B * A * lambda_f / lambda_p

    @staticmethod
    def armature_current(T: float, N: int, Z: int, flux: float) -> float:
        """Calculates the armature current (I_a).

        :param T: The induced torque.
        :param N: The number of slots.
        :param Z: The number of windings per slot.
        :param flux: is the flux of the poles (ɸ_p).
        :returns: The armature current (I_a).
        """
        return T * 2 * math.pi / N / Z / flux


class EMFInduced:
    """Induced EMF on the armature.

    Ea = N Z ɸ_p f_m
    -------------------------------------
    N is the number of slots.
    Z is the number of windings per slot.
    ɸ_p is the flux of the poles.
    f_m is the mechanical frequency.
    """
    @staticmethod
    def emf_induced(N: int, Z: int, flux: float, f_m: float) -> float:
        """EMF induced on the armature (Ea).

        :param N: is the number of slots.
        :param Z: is the number of windings per slot.
        :param flux: is the flux of the poles (ɸ_p).
        :param f_m: is the mechanical frequency.
        :returns: The EMF induced.
        """
        return N * Z * flux * f_m

    @staticmethod
    def mechanical_frequency(N: int, Z: int, flux: float, f_m: float) -> float:
        """EMF induced on the armature (f_m).

        :param N: is the number of slots.
        :param Z: is the number of windings per slot.
        :param flux: is the flux of the poles (ɸ_p).
        :param f_m: is the mechanical frequency.
        :returns: The EMF induced.
        """
        return N * Z * flux * f_m
