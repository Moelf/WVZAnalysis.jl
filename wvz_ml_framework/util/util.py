import numpy as np

def significance(s: float, b: float) -> float:
    '''
    Compute the asymptotic significance.

    Parameters
    ----------
    s: float
        Number of signal events.
    b: float
        Number of background events.
    '''
    if s == 0:
        return 0
    return np.sqrt(2 * ((s + b) * np.log(1 + s / b) - s))