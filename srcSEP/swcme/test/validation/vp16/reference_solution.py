#!/usr/bin/env python3
"""Outcome-blind skill-score utilities for VP16."""
import statistics
def skill_score(model_error: float, baseline_error: float) -> float:
    if baseline_error<=0: raise ValueError("baseline error must be positive")
    return 1.0-model_error/baseline_error
def median(values): return statistics.median(float(v) for v in values)

