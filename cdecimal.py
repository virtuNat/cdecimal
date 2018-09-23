import decimal as dec
import re

getcontext = dec.getcontext
setcontext = dec.setcontext
localcontext = dec.localcontext

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
            raise InvalidOperationError('2-arg inverse tangent not defined at the origin')
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

def is_close(num1, num2, prec=dec.Decimal('1E-9')):
    """Returns True if two numbers are sufficiently close to each other in value,
    which is to say the absolute value of the relative error is less than some
    referential benchmark value.
    """
    if not isinstance(num1, dec.Decimal):
        num1 = dec.Decimal(num1)
    if not isinstance(num2, dec.Decimal):
        num2 = dec.Decimal(num2)
    err = abs(num1 - num2)
    if num1 == 0:
        if num2 == 0:
            return True
        return err < prec
    if num2 == 0:
        return err < prec
    return 2 * err / (num1 + num2) < prec


class DivisionByZero(ZeroDivisionError):
    """Thrown when a division by zero occurs."""


class InvalidOperationError(ArithmeticError):
    """Thrown when an operation with an undefined value occurs."""


class CDecimal(object):
    """Constructs a new CDecimal object. 'real' and 'imag' can be any value
    type accepted by decimal.Decimal objects for construction. If no value
    is given, returns CDecimal('0', '0').
    """
    __slots__ = ('_real', '_imag')
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

    def __init__(self, real='0', imag='0', context=None):
        self.real = dec.Decimal(real, context)
        self.imag = dec.Decimal(imag, context)

    @classmethod
    def from_string(cls, value='0+0j', context=None):
        """Alternate constructor that constructs from a single string that
        represents a complex number. If no value is given, returns the same
        default value as the __init__ constructor.
        """
        value = value.strip()
        match = cls.imag_regex.match(value)
        if match:
            return cls(0, match[1], context)
        match = cls.cplx_regex.match(value)
        if match:
            return cls(match[1], match[2], context)
        raise ValueError('CDecimal.from_string argument is a malformed string')

    @classmethod
    def from_polar(cls, mod='0', arg='0', context=None):
        """Alternate constructor that constructs from a polar coordinate system
        instead of the default cartesian system. The unit of 'arg' is taken to be
        in radians. If no value is given, returns the same default value as the
        __init__ constructor.
        """
        mod = dec.Decimal(mod)
        arg = dec.Decimal(arg)
        return cls(mod * cos(arg), mod * sin(arg), context)

    def __getattr__(self, name):
        if name == 'real':
            return self._real
        elif name == 'imag':
            return self._imag

    def __setattr__(self, name, value):
        if name == 'real':
            self._real = +dec.Decimal(value)
        elif name == 'imag':
            self._imag = +dec.Decimal(value)
        else:
            super().__setattr__(name, value)

    def __complex__(self):
        """Returns the nearest complex float representation."""
        return eval(str(self))

    def __str__(self):
        return '({!s}{}{!s}j)'.format(
            self._real, '' if self._imag < 0 else '+', self._imag
            )

    def __repr__(self):
        return '{}(\'{!s}\', \'{!s}\')'.format(
            self.__class__.__name__,
            self._real,
            self._imag,
            )

    def __hash__(self):
        return super().__hash__()

    def __eq__(self, other):
        """Return self == other."""
        return self._real == other.real and self._imag == other.imag

    def __ne__(self, other):
        """Return self != other."""
        return self._real != other.real or self._imag != other.imag

    def __gt__(self, other):
        raise TypeError('complex numbers are not an ordered set')

    def __ge__(self, other):
        raise TypeError('complex numbers are not an ordered set')

    def __lt__(self, other):
        raise TypeError('complex numbers are not an ordered set')

    def __le__(self, other):
        raise TypeError('complex numbers are not an ordered set')

    def __pos__(self):
        """Return +self, useful when one wants to force rounding to current precision."""
        return self.__class__(+self._real, +self._imag)

    def __neg__(self):
        """Return -self"""
        return self.__class__(-self._real, -self._imag)

    def __abs__(self):
        """Returns the absolute distance of the complex number from zero."""
        return (self._real.fma(self._real, self._imag*self._imag)).sqrt()

    def __invert__(self):
        """Returns the complex conjugate of the number."""
        return self.__class__(self._real, -self._imag)

    def __add__(self, value):
        """Return self + value."""
        if isinstance(value, (int, dec.Decimal)):
            return self.__class__(self._real + value, self._imag)
        elif isinstance(value, self.__class__):
            return self.__class__(self._real + value._real, self._imag + value._imag)
        raise TypeError(
            'unsupported operand type(s) for +: {!r} and {!r}'.format(
                self.__class__.__name__, value.__class__.__name__
                )
            )

    def __radd__(self, value):
        """Return value + self."""
        if isinstance(value, (int, dec.Decimal)):
            return self.__class__(value + self._real, self._imag)
        elif isinstance(value, self.__class__):
            return self.__class__(value._real + self._real, value._imag + self._imag)
        raise TypeError(
            'unsupported operand type(s) for +: {!r} and {!r}'.format(
                value.__class__.__name__, self.__class__.__name__
                )
            )

    def __iadd__(self, value):
        """self += value"""
        if isinstance(value, (int, dec.Decimal)):
            self._real += value
        elif isinstance(value, self.__class__):
            self._real += value._real
            self._imag += value._imag
        raise TypeError(
            'unsupported operand type(s) for +: {!r} and {!r}'.format(
                self.__class__.__name__, value.__class__.__name__
                )
            )

    def __sub__(self, value):
        """Return self - value."""
        if isinstance(value, (int, dec.Decimal)):
            return self.__class__(self._real - value, self._imag)
        elif isinstance(value, self.__class__):
            return self.__class__(self._real - value._real, self._imag - value._imag)
        raise TypeError(
            'unsupported operand type(s) for -: {!r} and {!r}'.format(
                self.__class__.__name__, value.__class__.__name__
                )
            )

    def __rsub__(self, value):
        """Return value - self."""
        if isinstance(value, (int, dec.Decimal)):
            return self.__class__(value - self._real, -self._imag)
        elif isinstance(value, self.__class__):
            return self.__class__(value._real - self._real, value._imag - self._imag)
        raise TypeError(
            'unsupported operand type(s) for -: {!r} and {!r}'.format(
                value.__class__.__name__, self.__class__.__name__
                )
            )

    def __isub__(self, value):
        """self -= value"""
        if isinstance(value, (int, dec.Decimal)):
            self._real -= value
        elif isinstance(value, self.__class__):
            self._real -= value._real
            self._imag -= value._imag
        raise TypeError(
            'unsupported operand type(s) for -: {!r} and {!r}'.format(
                self.__class__.__name__, value.__class__.__name__
                )
            )

    def __mul__(self, value):
        """Return self * value."""
        if isinstance(value, (int, dec.Decimal)):
            return self.__class__(self._real * value, self._imag * value)
        elif isinstance(value, self.__class__):
            # (a + bj)*(c + dj) = (ac - bd) + (ad + bc)*j
            return self.__class__(
                self._real.fma(value._real, -self._imag*value._imag),
                self._real.fma(value._imag, self._imag*value._real)
                )
        raise TypeError(
            'unsupported operand type(s) for *: {!r} and {!r}'.format(
                self.__class__.__name__, value.__class__.__name__
                )
            )

    def __rmul__(self, value):
        """Return value * self."""
        if isinstance(value, (int, dec.Decimal)):
            return self.__class__(value * self._real , value * self._imag)
        elif isinstance(value, self.__class__):
            # (a + bj)*(c + dj) = (ac - bd) + (ad + bc)*j
            return self.__class__(
                value._real.fma(self._real, -value._imag*self._imag), 
                value._imag.fma(self._real, value._real*self._imag)
                )
        raise TypeError(
            'unsupported operand type(s) for *: {!r} and {!r}'.format(
                value.__class__.__name__, self.__class__.__name__
                )
            )

    def __imul__(self, value):
        """self *= value."""
        if isinstance(value, (int, dec.Decimal)):
            self._real *= value
            self._imag *= value
        elif isinstance(value, self.__class__):
            # (a + bj)*(c + dj) = (ac - bd) + (ad + bc)*j
            r = self._real.fma(value._real, -self._imag*value._imag)
            i = self._real.fma(value._imag, self._imag*value._real)
            self._real = r
            self._imag = i
        raise TypeError(
            'unsupported operand type(s) for *: {!r} and {!r}'.format(
                self.__class__.__name__, value.__class__.__name__
                )
            )

    def __truediv__(self, value):
        """Return self / value."""
        if isinstance(value, (int, dec.Decimal)):
            if value == 0:
                if self.is_zero():
                    raise InvalidOperationError(
                        'zero divided by zero is indeterminate'
                        )
                raise DivisionByZero('division by zero')
            return self.__class__(self._real / value, self._imag / value)
        elif isinstance(value, self.__class__):
            # (a + bj)/(c + dj) = ((ac + bd) + (bc - ad)*j) / (c*c + d*d)
            dr = value._real
            di = value._imag
            hy = dr.fma(dr, di*di)
            if hy == 0:
                if self.is_zero():
                    raise InvalidOperationError(
                        'zero divided by zero is indeterminate'
                        )
                raise DivisionByZero('division by zero')
            qr = self._real.fma(dr, self._imag*di) / hy
            qi = self._imag.fma(dr, -self._real*di) / hy
            return self.__class__(qr, qi)
        raise TypeError(
            'unsupported operand type(s) for /: {!r} and {!r}'.format(
                self.__class__.__name__, value.__class__.__name__
                )
            )

    def __rtruediv__(self, value):
        """Return value / self."""
        if isinstance(value, (int, dec.Decimal)):
            rr = dec.Decimal(value)
            ri = dec.Decimal('0')
        elif isinstance(value, self.__class__):
            rr = value._real
            ri = value._imag
        else:
            raise TypeError(
                'unsupported operand type(s) for /: {!r} and {!r}'.format(
                    value.__class__.__name__, self.__class__.__name__
                    )
                )
        # (a + bj)/(c + dj) = ((ac + bd) + (bc - ad)*j) / (c*c + d*d)
        dr = self._real
        di = self._imag
        hy = dr.fma(dr, di*di)
        if hy == 0:
            if rr == 0 and ri == 0:
                raise InvalidOperationError(
                    'zero divided by zero is indeterminate'
                    )
            raise DivisionByZero('division by zero')
        qr = rr.fma(dr, ri*di) / hy
        qi = ri.fma(dr, -rr*di) / hy
        return self.__class__(qr, qi)

    def __itruediv__(self, value):
        """self /= value."""
        if isinstance(value, (int, dec.Decimal)):
            if value == 0:
                if self.is_zero():
                    raise InvalidOperationError(
                        'zero divided by zero is indeterminate'
                        )
                raise DivisionByZero('division by zero')
            self._real /= value
            self._imag /= value
        elif isinstance(value, self.__class__):
            # (a + bj)/(c + dj) = ((ac + bd) + (bc - ad)*j) / (c*c + d*d)
            dr = value._real
            di = value._imag
            hy = dr.fma(dr, di*di)
            if hy == 0:
                if self.is_zero():
                    raise InvalidOperationError(
                        'zero divided by zero is indeterminate'
                        )
                raise DivisionByZero('division by zero')
            qr = self._real.fma(dr, self._imag*di) / hy
            qi = self._imag.fma(dr, -self._real*di) / hy
            self._real = qr
            self._imag = qi
        raise TypeError(
            'unsupported operand type(s) for /: {!r} and {!r}'.format(
                self.__class__.__name__, value.__class__.__name__
                )
            )

    def __pow__(self, value):
        """Return self ** value."""
        if isinstance(value, (int, dec.Decimal)):
            getcontext().prec += 2
            mod = abs(self)
            if mod.is_zero():
                if value == 0:
                    raise InvalidOperationError(
                        'zero raised to zero is indeterminate'
                        )
                arg = 0
            else:
                arg = atan2(self._imag, self._real) * value
            ans = self.__class__.from_polar(mod ** value, arg)
            getcontext().prec -= 2
            return +ans
        elif isinstance(value, self.__class__):
            if self.is_zero() and value.is_zero():
                raise InvalidOperationError(
                    'zero raised to zero is indeterminate'
                    )
            getcontext().prec += 2
            log_val = self.ln() * value
            mod = log_val.real.exp()
            getcontext().prec -= 2
            return self.__class__.from_polar(+mod, +log_val.imag)
        raise TypeError(
            'unsupported operand type(s) for **: {!r} and {!r}'.format(
                self.__class__.__name__, value.__class__.__name__
                )
            )

    def __rpow__(self, value):
        """Return value ** self."""
        if isinstance(value, (int, dec.Decimal)):
            if self.is_zero() and value == 0:
                raise InvalidOperationError(
                    'zero raised to zero is indeterminate'
                    )
            getcontext().prec += 2
            ans = self.__class__.from_polar(
                value ** self._real,
                value.ln() * self._imag,
                )
            getcontext().prec -= 2
            return +ans
        elif isinstance(value, self.__class__):
            if self.is_zero() and value.is_zero():
                raise InvalidOperationError(
                    'zero raised to zero is indeterminate'
                    )
            getcontext().prec += 2
            log_val = value.ln() * self
            mod = log_val.real.exp()
            getcontext().prec -= 2
            return self.__class__.from_polar(mod, log_val.imag)
        raise TypeError(
            'unsupported operand type(s) for **: {!r} and {!r}'.format(
                value.__class__.__name__, self.__class__.__name__
                )
            )

    def __ipow__(self, value):
        """self **= value."""
        if isinstance(value, (int, dec.Decimal)):
            getcontext().prec += 2
            mod = abs(self) ** value
            if mod.is_zero():
                if value == 0:
                    raise InvalidOperationError(
                        'zero raised to zero is indeterminate'
                        )
                arg = 0
            else:
                arg = atan2(self._imag, self._real) * value
            x = mod * cos(arg)
            y = mod * sin(arg)
            getcontext().prec -= 2
            self._real = +x
            self._imag = +y
        elif isinstance(value, self.__class__):
            if self.is_zero() and value.is_zero():
                raise InvalidOperationError(
                    'zero raised to zero is indeterminate'
                    )
            getcontext().prec += 2
            log_val = self.ln() * value
            mod = log_val.real.exp()
            arg = log_val.imag
            x = mod * cos(arg)
            y = mod * sin(arg)
            getcontext().prec -= 2
            self._real = +x
            self._imag = +y
        raise TypeError(
            'unsupported operand type(s) for **: {!r} and {!r}'.format(
                self.__class__.__name__, value.__class__.__name__
                )
            )

    def __floordiv__(self, value):
        raise TypeError('can\'t take the floor of a complex number')

    def __rfloordiv__(self, value):
        raise TypeError('can\'t take the floor of a complex number')

    def __ifloordiv__(self, value):
        raise TypeError('can\'t take the floor of a complex number')

    def __floor__(self):
        raise TypeError('can\'t take the floor of a complex number')

    def __ceil__(self):
        raise TypeError('can\'t take the ceiling of a complex number')

    def __trunc__(self):
        raise TypeError('can\'t truncate complex number to an integer')

    # Data representation functions
    def re(self):
        """Returns the real component of self."""
        return self._real

    def im(self):
        """Returns the imaginary component of self."""
        return self._imag

    def arg(self):
        """Returns the argument t such that r*cis(t) returns the complex number
        for some positive real value r, given that t is real.
        """
        if self.is_zero():
            raise InvalidOperationError(
                'argument of a complex number is not defined at 0'
                )
        return atan2(self._imag, self._real)

    def copy_abs(self):
        """Returns the absolute distance of the complex number from zero."""
        return (self._real.fma(self._real, self._imag*self._imag)).sqrt()

    def copy_negate(self):
        """Return -self"""
        return self.__class__(-self._real, -self._imag)

    # Format representation functions
    def as_tuple(self):
        """Returns the real and imaginary components as a tuple."""
        return (self._real, self._imag)

    def conjugate(self):
        """Returns the complex conjugate of the number."""
        return self.__class__(self._real, -self._imag)

    def to_float_complex(self):
        """Returns the nearest complex float representation."""
        return eval(str(self))

    def to_polar_tuple(self):
        """Returns the modulus and argument as a tuple."""
        return (abs(self), self.arg())

    def copy_sign(self, other, context=None):
        """Returns a copy of the first operand with the argument set to be
        the same as that of the second.
        """
        if isinstance(other, (int, dec.Decimal)):
            return self.__class__(
                self._real.copy_sign(other),
                self._imag.copy_sign(other),
                context
                )
        elif isinstance(other, self.__class__):
            try:
                arg = other.arg()
            except InvalidOperationError:
                arg = 0
            return self.__class__.from_polar(abs(self), arg)

    def is_zero(self):
        """Returns True if the argument is a zero, and False otherwise."""
        return self._real.is_zero() and self._imag.is_zero()
    
    def is_finite(self):
        """Returns True if the argument is a finite number, and False if either
        the real or imaginary components are infinity or NaN.
        """
        return self._real.is_finite() and self._imag.is_finite()

    def is_infinite(self):
        """Returns True if the argument is infinite, and False otherwise."""
        return self._real.is_infinite() or self._imag.is_infinite()

    def is_nan(self):
        """Returns True if the argument contains a NaN component."""
        return self._real.is_nan() or self._imag.is_nan()

    def is_real(self):
        """Returns True if the argument is equivalent to a real number and
        False otherwise.
        """
        return self._imag.is_zero()

    def is_imaginary(self):
        """Returns True if the argument is equivalent to a pure imaginary number
        and False otherwise.
        """
        return self._real.is_zero()

    def is_close(self, other, context=None):
        """Returns True if both arguments have a relative error that is less than
        the the last place unit in a given context.
        """
        if context is None:
            context = getcontext()
        prec = dec.Decimal('1E-'+str(context.prec))
        if isinstance(other, (int, dec.Decimal)):
            return (
                is_close(self._real, other, prec)
                and is_close(self._imag, dec.Decimal('0'), prec)
                )
        elif not isinstance(other, self.__class__):
            raise TypeError(
                'cannot compare {} and {} objects'.format(
                    self.__class__.__name__, other.__class__.__name__
                    )
                )
        return (
            is_close(self._real, other._real, prec)
            and is_close(self._imag, other._imag, prec)
            )

    # Power functions
    def sqrt(self):
        """Returns the two square roots of the operand."""
        getcontext().prec += 2
        mod = abs(self).sqrt()
        try:
            arg = atan2(self._imag, self._real) / 2
        except InvalidOperationError:
            arg = 0
        val = self.__class__.from_polar(mod, arg)
        getcontext().prec -= 2
        return (+val, -val)

    def sqrt1(self):
        """Returns the principal square root of the operand."""
        getcontext().prec += 2
        mod = abs(self).sqrt()
        try:
            arg = atan2(self._imag, self._real) / 2
        except InvalidOperationError:
            arg = 0
        val = self.__class__.from_polar(mod, arg)
        getcontext().prec -= 2
        return +val

    def cbrt(self):
        """Returns the three cube roots of the operand."""
        getcontext().prec += 2
        off = self.__class__(-0.5, dec.Decimal(0.75).sqrt()) # (-0.5+0.866j)
        mod = cbrt(abs(self))
        try:
            arg = atan2(self._imag, self._real) / 3
        except InvalidOperationError:
            arg = 0
        rt1 = self.__class__.from_polar(mod, arg)
        rt2 = rt1 * off
        rt3 = rt2 * off
        getcontext().prec -= 2
        return (+rt1, +rt2, +rt3)

    def cbrt1(self):
        """Returns the principal cube root of the operand."""
        getcontext().prec += 2
        mod = cbrt(abs(self))
        try:
            arg = atan2(self._imag, self._real) / 3
        except InvalidOperationError:
            arg = 0
        val = self.__class__.from_polar(mod, arg)
        getcontext().prec -= 2
        return +val

    def n_roots(self, index):
        """Returns the n nth roots of the complex number as a tuple.
        The index of the radicand must be a natural number.
        """
        if isinstance(index, int):
            pass
        elif isinstance(index, dec.Decimal):
            if index % 1 != 0:
                raise TypeError('radical index must be an integer')
            index = int(index)
        else:
            raise TypeError('radical index must be an integer')
        if index < 1:
            raise ValueError('radical index must be greater than zero')
        getcontext().prec += 2
        if index == 1:
            getcontext().prec -= 2
            return (self,)
        elif index == 2:
            mod = abs(self).sqrt()
        elif index == 3:
            mod = abs(self).cbrt()
        else:
            mod = abs(self) ** (dec.Decimal(1) / index)
        try:
            arg = atan2(self._imag, self._real) / index
        except InvalidOperationError:
            arg = 0
        val = self.__class__.from_polar(mod, arg)
        off = self.__class__.from_polar(1, 2 * pi() / index)
        tmp = val
        lst = [val]
        for _ in range(index - 1):
            tmp = tmp * off
            lst.append(tmp)
        return tuple(lst)

    def exp(self):
        """Returns the value of e raised to the power of the operand."""
        mod = abs(self).exp()
        try:
            arg = atan2(self._imag, self._real)
        except InvalidOperationError:
            arg = 0
        return self.__class__.from_polar(mod, arg)

    def ln(self):
        """Returns the natural logarithm of the operand."""
        mod = abs(self)
        if mod == 0:
            raise InvalidOperationError('logarithm is not defined at 0')
        rad = 0 if mod == 1 else mod.ln()
        return self.__class__(
            rad, atan2(self._imag, self._real)
            )

    def log10(self):
        """Returns the decimal logarithm of the operand."""
        base = dec.Decimal(10).ln()
        mod = abs(self)
        if mod == 0:
            raise InvalidOperationError('logarithm is not defined at 0')
        rad = 0 if mod == 1 else mod.log10()
        return self.__class__(
            rad, atan2(self._imag, self._real) / base
            )

    def log_base(self, base):
        """Returns the logarithm of self to the base 'base'."""
        if isinstance(base, int):
            base = dec.Decimal(base)
        elif not isinstance(base, (dec.Decimal, self.__class__)):
            raise TypeError((
                'base must be integer, decimal.Decimal instance,'
                'or cdecimal.CDecimal instance.'
                ))
        mod = abs(self)
        if mod == 0:
            raise InvalidOperationError('logarithm is not defined at 0')
        elif mod == 1:
            return self.__class__(0, 0)
        return self.__class__(mod.ln(), atan2(self._imag, self._real)) / base.ln()

    # Trigonometric functions
    def sin(self):
        """Returns the complex sine of self."""
        getcontext().prec += 2
        re = sin(self._real) * cosh(self._imag)
        im = cos(self._real) * sinh(self._imag)
        ans = self.__class__(re, im)
        getcontext().prec -= 2
        return +ans

    def cos(self):
        """Returns the complex cosine of self."""
        getcontext().prec += 2
        re = cos(self._real) * cosh(self._imag)
        im = sin(self._real) * sinh(self._imag)
        ans = self.__class__(re, -im)
        getcontext().prec -= 2
        return +ans

    def tan(self):
        """Returns the complex tangent of self."""
        getcontext().prec += 2
        re2 = 2 * self._real
        im2 = 2 * self._imag
        den = cos(re2) + cosh(im2)
        ans = self.__class__(sin(re2) / den, sinh(im2) / den)
        getcontext().prec -= 2
        return +ans

    def asin(self):
        """Returns the complex inverse sine of self."""
        getcontext().prec += 2
        im1 = self.__class__(0, 1)
        arg = im1*self + (1 - self*self).sqrt1()
        ans = -im1 * arg.ln()
        getcontext().prec -= 2
        return +ans

    def acos(self):
        """Returns the complex inverse cosine of self."""
        getcontext().prec += 2
        arg = self + (self*self - 1).sqrt1()
        ans = self.__class__(0, -1) * arg.ln()
        getcontext().prec -= 2
        return +ans

    def atan(self):
        """Returns the complex inverse tangent of self."""
        getcontext().prec += 2
        im1 = self.__class__(0, 1) * self
        arg = (1 - im1) / (1 + im1)
        ans = self.__class__(0, 0.5) * arg.ln()
        getcontext().prec -= 2
        return +ans

    # Hyperbolic trigonometric functions
    def sinh(self):
        """Returns the complex hyperbolic sine of self."""
        getcontext().prec += 2
        re = sinh(self._real) * cos(self._imag)
        im = cosh(self._real) * sin(self._imag)
        ans = self.__class__(re, im)
        getcontext().prec -= 2
        return +ans

    def cosh(self):
        """Returns the complex hyperbolic cosine of self."""
        getcontext().prec += 2
        re = cosh(self._real) * cos(self._imag)
        im = sinh(self._real) * sin(self._imag)
        ans = self.__class__(re, im)
        getcontext().prec -= 2
        return +ans

    def tanh(self):
        """Returns the complex hyperbolic tangent of self."""
        getcontext().prec += 2
        re2 = 2 * self._real
        im2 = 2 * self._imag
        den = cosh(re2) + cos(im2)
        ans = self.__class__(sinh(re2) / den, sin(im2) / den)
        getcontext().prec -= 2
        return +ans

    def asinh(self):
        """Returns the complex inverse hyperbolic sine of self."""
        getcontext().prec += 2
        arg = self + (self*self + 1).sqrt1()
        ans = arg.ln()
        getcontext().prec -= 2
        return +ans

    def acosh(self):
        """Returns the complex inverse hyperbolic cosine of self."""
        getcontext().prec += 2
        arg = self + (self*self - 1).sqrt1()
        ans = arg.ln()
        getcontext().prec -= 2
        return +ans

    def atanh(self):
        """Returns the complex inverse hyperbolic tangent of self."""
        getcontext().prec += 2
        arg = (1 - self) / (1 + self)
        ans = arg.ln() / 2
        getcontext().prec -= 2
        return +ans

    # Special functions
    def gamma(self):
        """Returns the result of the gamma function on self, which is the
        analytic continuation of the factorial on all complex numbers.
        """
        raise NotImplementedError('not implemented yet, will use spouge approximation')

    def zeta(self):
        """Returns the result of the riemann zeta function on self."""
        raise NotImplementedError('not implemented yet')

    # decimal.Decimal accuracy effecting functions
    def normalize(self, context=None):
        """Normalizes both components of the argument."""
        self._real.normalize(context)
        self._imag.normalize(context)

    def quantize(self, exp, rounding=None, context=None):
        """Quantizes both components of the argument."""
        self._real.quantize(exp, rounding, context)
        self._imag.quantize(exp, rounding, context)

    def same_quantum(self, context=None):
        """Returns True if both components have the same exponent or if
        both components are NaN.
        """
        return self._real.same_quantum(self._imag, context)

    def radix(self):
        """Returns the radix of mathematical operation."""
        return self.__class__(10, 0)

IMAG_UNIT = CDecimal(0, 1)
