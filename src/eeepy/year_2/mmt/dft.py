#!/usr/bin/env python3
"""Module containing functions and equations relating to discrete fourier
 transform."""
import math


class SampleTime(float):
    """The sample time needed parameter of a DFT.

    :param f_max: The maximum frequency of the sampled signal.
    """
    def __new__(cls, f_max: float):
        return super().__new__(cls, 1 / (2 * f_max))

    def __init__(self, f_max: float):
        super().__init__()

        self._f_sample: float = self**-1

    def calculate_f_max(self) -> float:
        """Calculates the f_max of a given sample time.

        :return: The f_max.
        """
        return (2 * self)**-1

    @property
    def f_sample(self):
        """The sampling frequency of the sample time."""
        return self._f_sample

    @classmethod
    def from_t_min(cls, t_min: float):
        """Calculates the sample time using minimum period instead.

        :param t_min: Minimum period of the sampled signal.
        """
        return cls(t_min**-1)


class SampleNumber(float):
    """The number of samples parameter of a DFT.

    :param sample_time: The sample time of the sampled signal.
    :param f_min: The minimum frequency of the sampled signal.
    """
    def __new__(cls, dt: SampleTime, f_min: float):
        return super().__new__(cls, math.ceil(1 / (dt * f_min)))

    def __init__(self, dt: SampleTime, f_min: float):
        super().__init__()

    @classmethod
    def from_f_max(cls, f_min: float, f_max: float):
        """Calculate the SampleNumber using f_min and f_max.

        :param f_min: The minimum frequency of the sampled signal.
        :param f_max: The maximum frequency of the sampled signal.
        """
        dt: SampleTime = SampleTime(f_max)
        return cls(dt, f_min)

    def calculate_sample_time(self, f_min: float):
        """Calculates the sample time.

        :param f_min: The minimum frequency of the sampled signal.
        :return: The sample time Δt
        """
        return 1 / (self * f_min)

    def calculate_f_min(self, dt: SampleTime):
        """Calculates the f_min.

        :param sample_time: The sample time of the sampled signal.
        :return: The minimum frequency of the sampled signal.
        """
        return 1 / (self * dt)
