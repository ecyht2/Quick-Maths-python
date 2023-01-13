#!/usr/bin/env python3
"""This file contains test for dft module of eeepy.year_2.mmt"""
from eeepy.year_2.mmt.dft import (SampleTime, SampleNumber,
                                  FrequencySamplingInterval)


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
        assert dt.calculate_f_max() == 20e3


class TestSampleNumber:
    """Tests for SampleNumber equation."""
    f_min: float = 10
    f_max: float = 20e3

    def test_sample_number(self):
        """Testing sample number calculation."""
        dt: SampleTime = SampleTime(self.f_max)
        N: SampleNumber = SampleNumber(dt, self.f_min)
        assert N == 4000

    def test_from_f_max(self):
        """Testing sample number calculation using f_max."""
        N: SampleNumber = SampleNumber.from_f_max(self.f_min, self.f_max)
        assert N == 4000

    def test_sample_time(self):
        """Testing sample time calculation."""
        N: SampleNumber = SampleNumber.from_f_max(self.f_min, self.f_max)
        assert N.calculate_sample_time(self.f_min) == 25e-6

    def test_f_min(self):
        """Testing f_min calculation."""
        N: SampleNumber = SampleNumber.from_f_max(self.f_min, self.f_max)
        dt: SampleTime = SampleTime(self.f_max)
        assert N.calculate_f_min(dt) == 10


class TestFrequencySamplingInterval:
    """Tests for FrequencySamplingInterval equation."""
    f_min: float = 10

    def test_sample_interval(self):
        """Testing frequency sample interval Δɷ calculation."""
        dt: FrequencySamplingInterval = FrequencySamplingInterval(self.f_min)
        assert dt == 62.83185307179586

    def test_from_t_max(self):
        """Testing frequency sample interval calculation using t_max."""
        t_max: float = self.f_min**-1
        dt: FrequencySamplingInterval = FrequencySamplingInterval.from_t_max(
            t_max
        )
        assert dt == 62.83185307179586

    def test_f_min(self):
        """Testing f_min calulation."""
        dt: FrequencySamplingInterval = FrequencySamplingInterval(self.f_min)
        assert dt.calculate_f_min() == 10
