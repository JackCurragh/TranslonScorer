"""Test package marker.

Makes ``tests`` an importable package so cross-module fixtures resolve, e.g.
``from tests.test_golden import _gapdh_events`` in test_bam_provider.py.
"""
