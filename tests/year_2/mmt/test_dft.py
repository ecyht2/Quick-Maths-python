#!/usr/bin/env python3
"""This file contains test for dft module of eeepy.year_2.mmt"""
from eeepy.year_2.mmt.dft import (SampleTime)


class TestSampleTime:
    """Tests for SampleTime equation."""
    f_max: float = 20e3

    def test_sample_time(self):
        """Testing sample time Δt calculation."""
        dt: SampleTime = SampleTime(self.f_max)
        assert dt == 25e-6

    def test_f_sample(self):
        """Testing sample frequency f_sample."""
        dt: SampleTime = SampleTime(self.f_max)
        assert dt.f_sample == 40e3

    def test_from_t_min(self):
        """Testing sample time calculation using t_min."""
        t_min: float = self.f_max**-1
        dt: SampleTime = SampleTime.from_t_min(t_min)
        assert dt == 25e-6

    def test_f_max(self):
        """Testing f_max calulation."""
        dt: SampleTime = SampleTime(self.f_max)
        assert dt.calculate_f_max() == self.f_max
