#!/usr/bin/env python3
"""Functions and equations related to DC machines."""
import math


class PoleFlux:
    """The total flux of the poles.

    Formula:
    :math:`\\Phi_p = B l \\frac{2πR}{2P} \\frac{\\lambda_f}{\\lambda_p}`

    :math:`\\Phi_p` is the flux of the poles.

    B is the flux density of the poles.

    l is the length of the

    R is the radius of the poles.

    P is the number of pole pairs.

    :math:`\\lambda_f` is the salient field winding tooth (pole).

    :math:`\\lambda_p` is the pole pitch.
    """
    @staticmethod
    def flux(B: float, R: float, P: int, l: float,
             lambda_f: float, lambda_p: float) -> float:
        """Calculates the total flux of the poles (:math:`\\Phi_p`).

        :param B: The flux density of the poles.
        :param R: The radius of the poles.
        :param P: The number of pole pairs.
        :param l: The length of the
        :param lambda_f: The salient field winding tooth (pole).
        :param lambda_p: The pole pitch.
        :returns: The flux of the poles.
        """
        return B * l * math.pi * R / P * lambda_f / lambda_p

    @staticmethod
    def flux_poles(B: float, R: float, p: int, l: float,
                   lambda_f: float, lambda_p: float) -> float:
        """Calculates the total flux of the poles (:math:`\\Phi_p`).

        This method differs from flux as it uses the number of poles (p)
        instead of the number of pole pairs (P).

        :param B: The flux density of the poles.
        :param R: The radius of the poles.
        :param p: The number of poles.
        :param l: The length of the
        :param lambda_f: The salient field winding tooth (pole).
        :param lambda_p: The pole pitch.
        :returns: The flux of the poles.
        """
        return B * l * 2 * math.pi * R / p * lambda_f / lambda_p


class InducedTorque:
    """Induced torque in an ideal DC machine.

    Formula: :math:`T = \\Phi_p \\frac{N Z I_a}{2 \\pi}`

    T is the induced torque.

    N is the number of slots.

    Z is the number of windings per slot.

    :math:`I_a` is the armature current.

    :math:`\\Phi_p` is the flux of the poles.

    From Electrical Loading:

    Formula: :math:`T = 2 V B A \\frac{\\lambda_f}{\\lambda_p}`

    T is the induced torque.

    V is the volume considered at the stator bore. (V must be calculated at the
    center of the airgap if both airgap and stacking factors are given.)

    B is the magnetic loading.

    A is the current loading.

    :math:`\\lambda_f` is the salient field winding tooth (pole).

    :math:`\\lambda_p` is the pole pitch.
    """
    @staticmethod
    def torque(N: int, Z: int, I: float, flux: float) -> float:
        """Calcultes the induced torque (T).

        :param N: The number of slots.
        :param Z: The number of windings per slot.
        :param I: The armature current (:math:`I_a`).
        :param flux: The flux of the poles (:math:`\\Phi_p`).
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
        :param lambda_f: The salient field winding tooth (pole).
        :param lambda_p: The pole pitch.
        :returns: The torque induced.
        """
        return 2 * V * B * A * lambda_f / lambda_p

    @staticmethod
    def armature_current(T: float, N: int, Z: int, flux: float) -> float:
        """Calculates the armature current (:math:`I_a`).

        :param T: The induced torque.
        :param N: The number of slots.
        :param Z: The number of windings per slot.
        :param flux: The flux of the poles (:math:`ɸ_p`).
        :returns: The armature current (:math:`I_a`).
        """
        return T * 2 * math.pi / N / Z / flux


class EMFInduced:
    """Induced EMF on the armature.

    Formula: :math:`E_a = N Z \\Phi_p f_m`

    N is the number of slots.

    Z is the number of windings per slot.

    :math:`\\Phi_p` is the flux of the poles.

    :math:`f_m` is the mechanical frequency.
    """
    @staticmethod
    def emf_induced(N: int, Z: int, flux: float, f_m: float) -> float:
        """EMF induced on the armature (:math:`Ea`).

        :param N: The number of slots.
        :param Z: The number of windings per slot.
        :param flux: The flux of the poles (:math:`ɸ_p`).
        :param f_m: The mechanical frequency.
        :returns: The EMF induced.
        """
        return N * Z * flux * f_m

    @staticmethod
    def mechanical_frequency(Ea: float, N: int, Z: int, flux: float) -> float:
        """EMF induced on the armature (:math:`f_m`).

        :param Ea: The EMF induced.
        :param N: The number of slots.
        :param Z: The number of windings per slot.
        :param flux: The flux of the poles (:math:`\\Phi_p`).
        :returns: The mechanical frequency.
        """
        return Ea / (N * Z * flux)
