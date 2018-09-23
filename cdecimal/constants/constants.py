import decimal as dec
from .. import CDecimal

getcontext = dec.getcontext

def PI():
    """Returns the ratio of a circle's circumference to its diameter."""
    getcontext().prec += 2
    agm0 = agm1 = pow2 = dec.Decimal(1)
    denm = dec.Decimal(0.25)
    agm2 = dec.Decimal(0.5).sqrt()
    # Uses the elliptic integral relation to the agm to calculate pi.
    while True:
        agm1 = (agm1 + agm2) / 2
        agm2 = (agm0 * agm2).sqrt()
        diff = agm1 - agm0
        agm0 = agm1
        denm -= pow2 * diff * diff
        if diff == 0:
            ave = (agm1 + agm2) / 2
            arc = ave * ave / denm
            getcontext().prec -= 2
            return CDecimal(+arc)
        pow2 *= 2

def E():
    """Returns the base of the natural logarithm."""
    return CDecimal(1).exp()

def I():
    """Returns the imaginary unit."""
    return CDecimal(0, 1)
