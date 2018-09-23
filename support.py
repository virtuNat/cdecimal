"""Support functions for cdecimal, these should not be directly importable nor interactable."""
import decimal as dec
import re

getcontext = dec.getcontext

imag_regex = re.compile(
    r"""
    (
        [-+]? # Sign is optional
        (?:   # Manitssa
            \d+          # Integer part: 345
            (?:\.\d*)?   # Decimal part after integer: 345.123
            |\.\d+       # Integer part optional: .123
        )
        (?:   # Exponent is optional
            [eE][-+]?\d+ # 1E-7
        )?
    )[jJ]     # Signifies this is an imaginary number
    """,
    re.VERBOSE
    )
cplx_regex = re.compile(
    r"""
    (   # Real part, sign is optional
        [-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?
    )
    (   # Imaginary part, sign is mandatory
        [-+](?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?
    )[jJ]
    """,
    re.VERBOSE
    )

def pi():
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
            return +arc
        pow2 *= 2

def sin(ang):
    """Returns the signed distance from the x axis of a point on the unit circle given
    the angle in radians that it makes with the positive x axis with the origin as its
    vertex.
    """
    # Note: Pretty fucking slow.
    if not isinstance(ang, dec.Decimal):
        ang = dec.Decimal(ang)
    pi_val = pi()
    pi_2 = pi_val / 2
    context = getcontext()
    context.prec += 2
    # Clamp all equivalent angles to the interval [-2*pi, 2*pi]
    # sin(x + 2*pi) = sin(x)
    ang = ang % (2 * pi_val)
    # Clamp all equivalent angles to the interval [0, 2*pi]
    # sin(x + 2*pi) = sin(x)
    if ang < 0:
        ang = pi_val.fma(2, ang)
    # Clamp all equivalent angles to the interval [-pi/2, pi/2]
    # sin(pi - x) = sin(x)
    if pi_2 < ang < 3 * pi_2:
        ang = pi_val - ang
    # sin(-x) = -sin(x)
    sgn = dec.Decimal(1).copy_sign(ang)
    ang = abs(ang)
    # Calculates sine using the taylor series expansion at 0
    ctr = 5
    val = -ang * ang * ang / 6
    total = ang + val
    while True:
        val = -val * ang * ang / ctr / (ctr - 1)
        if val == 0 or val.logb() < total.logb() - context.prec:
            context.prec -= 2
            if total == 0 or abs(total).logb() < -context.prec:
                return sgn * dec.Decimal(0)
            return sgn * total
        total += val
        ctr += 2

def cos(ang):
    """Returns the signed distance from the y axis of a point on the unit circle given
    the angle in radians that it makes with the positive x axis with the origin as its
    vertex.
    """
    # Note: Pretty fucking slow.
    if not isinstance(ang, dec.Decimal):
        ang = dec.Decimal(ang)
    pi_val = pi()
    pi_2 = pi_val / 2
    context = getcontext()
    context.prec += 2
    # Clamp all equivalent angles to the interval [-2*pi, 2*pi]
    # cos(x + 2*pi) = cos(x)
    ang = ang % (2 * pi_val)
    # Clamp all equivalent angles to the interval [0, 2*pi]
    # cos(x + 2*pi) = cos(x)
    if ang < 0:
        ang = pi_val.fma(2, ang)
    # Clamp all equivalent angles to the interval [0, pi]
    # cos(2*pi - x) = cos(x)
    if ang > pi_val:
        ang = pi_val.fma(2, -ang)
    # cos(pi - x) = -cos(x)
    if ang > pi_2:
        sgn = -1 
        ang = pi_val - ang
    elif ang < pi_2:
        sgn = 1
    else:
        return dec.Decimal(0)
    # Calculates cosine using the taylor series expansion at 0
    ctr = 4
    val = -ang * ang / 2
    total = 1 + val
    while True:
        val = -val * ang * ang / ctr / (ctr - 1)
        if val == 0 or val.logb() < total.logb() - context.prec:
            context.prec -= 2
            if total == 0 or abs(total).logb() < -context.prec:
                return sgn * dec.Decimal(0)
            return sgn * total
        total += val
        ctr += 2

def atan_half():
    """Returns atan(0.5)."""
    context = getcontext()
    context.prec += 2
    prec = dec.Decimal('1E-'+str(context.prec))
    val = dec.Decimal(0.5)
    num1 = 1
    num2 = val * val
    den1 = 3
    den2 = 1 + num2
    term = val / den2
    total = term
    while True:
        prev = term
        term = term * 4 * num1 * num1 * num2 / den1 / (den1 - 1) / den2
        if term == 0 or term.logb() < total.logb() - context.prec:
            context.prec -= 2
            return +total
        total += term
        num1 += 1
        den1 += 2

def atan(val):
    """Returns the principal angle in radians for which, on the unit circle,
    a point whose ray casted from the origin makes that angle with respect
    to the positive x axis forms a right triangle with respect to the coordinate
    axes such that the ratio of the length of the leg parallel to the y axis
    to the length of the leg along the x axis is the given value.
    """
    if not isinstance(val, dec.Decimal):
        val = dec.Decimal(val)
    # atan(-x) = -atan(x)
    sgn = dec.Decimal(1).copy_sign(val)
    val = abs(val)
    pi_val = pi()
    context = getcontext()
    context.prec += 2
    if val == dec.Decimal('Infinity'):
        ans = (pi_val / 2).copy_sign(sgn)
        context.prec -= 2
        return +ans
    # atan(x) = pi/2 - atan(1/x)
    if val > 1:
        off = pi_val / 2
        val = 1 / val
    else:
        off = 0
    # atan(x) = atan(y) + atan((x - y) / (1 + x*y))
    if val > 0.5:
        at_hlf = atan_half()
        val = (val - dec.Decimal(0.5)) / (1 + val/2)
    else:
        at_hlf = 0
    num1 = 1
    num2 = val * val
    den1 = 3
    den2 = 1 + num2
    term = val / den2
    total = term
    while True:
        term *= 4 * num1 * num1 * num2 / den1 / (den1 - 1) / den2
        if term == 0 or term.logb() < total.logb() - context.prec:
            if total == 0 or abs(total).logb() < -context.prec:
                context.prec -= 2
                return sgn * dec.Decimal(0)
            total += at_hlf
            if off != 0:
                total = off - total 
            context.prec -= 2
            return +(sgn * total)
        total += term
        num1 += 1
        den1 += 2

def atan2(y, x):
    """Returns the principal angle in radians that a given point on the cartesian
    plane makes with respect to the positive x axis with the origin as the vertex.
    """
    if x == 0:
        if y == 0:
            raise ValueError('2-arg inverse tangent not defined at the origin')
        return pi().copy_sign(dec.Decimal(y)) / 2
    elif y == 0:
        return pi() if x < 0 else dec.Decimal(0)
    getcontext().prec += 2
    arg = dec.Decimal(y) / dec.Decimal(x)
    if arg != 1:
        if x < 0:
            val = atan(arg) - pi().copy_sign(arg)
        else:
            val = atan(arg)
    else:
        if x < 0:
            val = -3 * pi().copy_sign(arg) / 4
        else:
            val = pi().copy_sign(arg) / 4
    getcontext().prec -= 2
    return +val

def sinh(ang):
    """Returns the signed distance from the x axis of a point on the unit hyperbola given
    the hyperbolic angle in radians that it makes with the positive x axis with the origin
    as its vertex.
    """
    if not isinstance(ang, dec.Decimal):
        ang = dec.Decimal(ang)
    context = getcontext()
    context.prec += 2
    num = ang * ang
    den = 5
    term = ang * num / 6
    total = ang + term
    while True:
        prev = dec.Decimal(term)
        term *= num / den / (den - 1)
        if term == 0 or term.logb() < total.logb() - context.prec:
            context.prec -= 2
            return +total
        total += term
        den += 2

def cosh(ang):
    """Returns the absolute distance from the y axis of a point on the unit hyperbola given
    the hyperbolic angle in radians that it makes with the positive x axis with the origin
    as its vertex.
    """
    if not isinstance(ang, dec.Decimal):
        ang = dec.Decimal(ang)
    context = getcontext()
    context.prec += 2
    num = ang * ang
    den = 4
    term = num / 2
    total = 1 + term
    while True:
        prev = dec.Decimal(term)
        term *= num / den / (den - 1)
        if term == 0 or term.logb() < total.logb() - context.prec:
            context.prec -= 2
            return +total
        total += term
        den += 2

def cbrt(num):
    """Returns the principal number for which, when multiplied to itself twice
    would have a value equal to the argument.
    """
    if not isinstance(num, dec.Decimal):
        num = dec.Decimal(num)
    if num == 0:
        return dec.Decimal(0)
    context = getcontext()
    context.prec += 2
    x = num / 3
    while True:
        p = x
        x = (2*x + num/x/x) / 3
        if p == x:
            context.prec -= 2
            return x
