from .configs import LogOptions as opt

__all__ = ["SILENT", "MINIMAL", "INFO", "DEBUG", "PB_WIDTH"]

SILENT = 0
MINIMAL = 1
INFO = 2
DEBUG = 3

PB_START = -1
PB_STOP = -2
PB_WIDTH = 25

def Log(level, *args):
    if level <= opt.LEVEL:
        print(*args, flush=True)

def ProgressBar(level, val):
    """
    Starts, advances, or stops a progress bar.
    
    Arguments:
    
        - level: log level (one of SILENT, MINIMAL, INFO, DEBUG).
        - val: a float in [0,1] indicating progress, or one of PB_START, PB_STOP.
    """
    if level > opt.LEVEL:
        return
    if val == PB_START:
        print("Progress: ["+ (" "*PB_WIDTH) +"]", end="\r", flush=True)
        ProgressBar.last = 0
        return
    count = round(val * PB_WIDTH) if val != PB_STOP else PB_WIDTH
    if count != ProgressBar.last:
        print("Progress: ["+ "="*count + " "*(PB_WIDTH-count) +"]", end="\r", flush=True)
        ProgressBar.last = count
    if val == PB_STOP:
        print(flush=True)
