#!/usr/bin/env python3
"""Frozen qualitative connection histories for the VP11 literature benchmark."""
HISTORIES={"STEREO-A":[1,1,1,1,1,1],"STEREO-B":[1,1,1,0,0,0],"SOHO":[0,0,1,1,0,0]}
def expected_history(observer: str):
    """Return a copy so a runner cannot mutate the frozen oracle."""
    return list(HISTORIES[observer])

