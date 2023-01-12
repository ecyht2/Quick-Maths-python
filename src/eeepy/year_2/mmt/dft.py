#!/usr/bin/env python3
"""Module containing functions and equations relating to discrete fourier
 transform."""


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

        :return: The f_max
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
