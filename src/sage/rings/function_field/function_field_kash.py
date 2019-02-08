# -*- coding: utf-8 -*-
"""
Function Fields implemented with an expect interface to the optional
package kash.

This implementation is incomplete.  Many standard methods have not
yet been implemented.

A major feature of this implemention is its support for QQbar as
a base field.  Since kash only supports number fields, this is
implemented by keeping a working_constant_field, which is a number
field, and extending it as needed to include new algebraic numbers.
Most dependent objects keep a local variable kash_constant_field
and check it against the base function field's kash_constant_field
any time an operation is requested, recomputing all of the kash
objects if the kash_constant_field has changed.

EXAMPLES::

    sage: F.<x> = FunctionField(QQ, implementation='kash')
    sage: R.<Y> = F[]
    sage: L.<y> = F.extension(Y^2 - x^8 - 1)
    sage: O = L.maximal_order()
    sage: I = O.ideal(x, y-1)
    sage: P = I.place()
    sage: D = P.divisor()
    sage: D.basis_function_space()
    [1]
    sage: (2*D).basis_function_space()
    [1]
    sage: (3*D).basis_function_space()
    [1]
    sage: (4*D).basis_function_space()
    [1, 1/x^4*y + 1/x^4]


    sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<t> = PolynomialRing(K)
    sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
    sage: O = F.maximal_order()
    sage: I = O.ideal(y)
    sage: I.divisor()
    2*Place (x, (1/(x^3 + x^2 + x))*y^2)
     + 2*Place (x^2 + x + 1, (1/(x^3 + x^2 + x))*y^2)

    sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
    sage: L.<y> = K.extension(Y^2+Y+x+1/x)
    sage: O = L.maximal_order()
    sage: I = O.ideal(y)
    sage: I.divisor()
    -1*Place (x, x*y)
     + Place (x^2 + 1, x*y)

    sage: R.<x> = FunctionField(QQbar, implementation='kash')
    sage: L.<y> = R[]
    sage: F.<y> = R.extension(y^2 - (x^2+1))
    sage: D = (x/y).divisor()
    sage: D
    -1*Place (x - I, y)
     + Place (x, y + x - 1)
     + Place (x, y + x + 1)
     - Place (x + I, y)
    sage: D2 = (x/y*x.differential()).divisor()
    sage: pl1 = D.support()[0]
    sage: m = F.completion(pl1)
    sage: m(x/y)
    1/2*s^-1 + 1/2*s + O(s^20)
    sage: m(x/y*x.differential())
    2*I + 6*I*s^2 + 10*I*s^4 + 14*I*s^6 + 18*I*s^8 + 22*I*s^10 + 26*I*s^12 + 30*I*s^14 + 34*I*s^16 + 38*I*s^18 + O(s^20)
    sage: pl2 = D.support()[1]
    sage: m2 = F.completion(pl2)
    sage: m2(x)
    s + 1/4*s^3 + 1/16*s^5 + 1/64*s^7 + 1/256*s^9 + 1/1024*s^11 + 1/4096*s^13 + 1/16384*s^15 + 1/65536*s^17 + 1/262144*s^19 + O(s^20)
    sage: m2(x/y)
    s - 1/4*s^3 + 1/16*s^5 - 1/64*s^7 + 1/256*s^9 - 1/1024*s^11 + 1/4096*s^13 - 1/16384*s^15 + 1/65536*s^17 - 1/262144*s^19 + O(s^20)
    sage: m2(x.differential())
    1 + 3/4*s^2 + 5/16*s^4 + 7/64*s^6 + 9/256*s^8 + 11/1024*s^10 + 13/4096*s^12 + 15/16384*s^14 + 17/65536*s^16 + 19/262144*s^18 + O(s^20)
    sage: m2(x * x.differential())
    s + s^3 + 9/16*s^5 + 1/4*s^7 + 25/256*s^9 + 9/256*s^11 + 49/4096*s^13 + 1/256*s^15 + 81/65536*s^17 + 25/65536*s^19 + O(s^20)
"""

import itertools

from sage.misc.cachefunc import cached_method
from sage.misc.flatten import flatten

from sage.interfaces.kash import kash, KashElement

from sage.structure.richcmp import richcmp
from sage.structure.factorization import Factorization

from sage.rings.rational_field import QQ
from sage.rings.fraction_field import FractionField
from sage.rings.finite_rings.finite_field_base import FiniteField
from sage.rings.number_field.number_field_base import NumberField
from sage.rings.qqbar import QQbar, number_field_elements_from_algebraics

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

from .function_field import RationalFunctionField
from .function_field import FunctionFieldElement_rational
from .function_field import FunctionField_polymod
from .function_field import FunctionFieldElement_polymod
from .differential import DifferentialsSpace, FunctionFieldDifferential_global
from .order import FunctionFieldMaximalOrder
from .order import FunctionFieldMaximalOrderInfinite
from .ideal import FunctionFieldIdeal
from .place import FunctionFieldPlace
from .divisor import FunctionFieldDivisor
from .constructor import FunctionField
from .maps import FunctionFieldCompletion

class RationalFunctionField_kash(RationalFunctionField):
    """
    Rational function field `K(t)` in one variable, over the rationals,
    implemented using kash.

    EXAMPLES::

        sage: K.<t> = FunctionField(QQ, implementation='kash'); K
        Rational function field in t over Rational Field
        sage: K.gen()
        t
        sage: 1/t + t^3 + 5
        (t^4 + 5*t + 1)/t

    There are various ways to get at the underlying fields and rings
    associated to a rational function field::

        sage: K.<t> = FunctionField(QQ, implementation='kash')
        sage: K.base_field()
        Rational function field in t over Rational Field
        sage: K.field()
        Fraction Field of Univariate Polynomial Ring in t over Rational Field
        sage: K.constant_field()
        Rational Field
        sage: K.maximal_order()
        Maximal order of Rational function field in t over Rational Field

    We define a morphism::

        sage: K.<t> = FunctionField(QQ, implementation='kash')
        sage: L = FunctionField(QQ, 'tbar') # give variable name as second input
        sage: K.hom(L.gen())
        Function Field morphism:
          From: Rational function field in t over Rational Field
          To:   Rational function field in tbar over Rational Field
          Defn: t |--> tbar

    Here's a calculation over a number field::

        sage: R.<x> = FunctionField(QQ, implementation='kash')
        sage: L.<y> = R[]
        sage: F.<y> = R.extension(y^2 - (x^2+1))
        sage: (y/x).divisor()
        -1*Place (x, y + x - 1)
         - Place (x, y + x + 1)
         + Place (x^2 + 1, y)

        sage: A.<z> = QQ[]
        sage: NF.<i> = NumberField(z^2+1)
        sage: R.<x> = FunctionField(NF, implementation='kash')
        sage: L.<y> = R[]
        sage: F.<y> = R.extension(y^2 - (x^2+1))

        sage: (x/y*x.differential()).divisor()
        -2*Place (1/x, 1/x*y + (-x + 1)/x)
         - 2*Place (1/x, 1/x*y + (x + 1)/x)
         + Place (x, y + x - 1)
         + Place (x, y + x + 1)

        sage: (x/y).divisor()
        -1*Place (x - i, y)
         + Place (x, y + x - 1)
         + Place (x, y + x + 1)
         - Place (x + i, y)

    We try the same calculation over QQbar::

        sage: R.<x> = FunctionField(QQbar, implementation='kash')
        sage: L.<y> = R[]
        sage: F.<y> = R.extension(y^2 - (x^2+1))
        sage: (y/x).divisor()
        Place (x - I, y)
         - Place (x, y + x - 1)
         - Place (x, y + x + 1)
         + Place (x + I, y)

    """

    def __init__(self, constant_field, names, category=None):
        """
        Create a rational function field in one variable over the rational
        field.

        INPUT:

        - ``constant_field`` -- QQ

        - ``names`` -- string or tuple of length 1

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ, implementation='kash'); K
            Rational function field in t over Rational Field
            sage: K.category()
            Category of function fields
            sage: FunctionField(QQ[I], 'alpha', implementation='kash')   # not implemented
            Rational function field in alpha over Number Field in I with defining polynomial x^2 + 1

        Must be over a field::

            sage: FunctionField(ZZ, 't')
            Traceback (most recent call last):
            ...
            TypeError: constant_field must be a field
        """

        RationalFunctionField.__init__(self, constant_field, names, category)

        # the field, without the function field structure (referenced by divisor.py)
        self._field = constant_field[names[0]].fraction_field()

        if constant_field is QQbar:
            self._set_working_constant_field(QQ)
        else:
            self._set_working_constant_field(constant_field)

    def _set_working_constant_field(self, constant_field):

        self.working_constant_field = constant_field

        if constant_field is QQ:
            self.kash_constant_field = kash.RationalField()
            self.reverse_map = {}
        elif isinstance(constant_field, FiniteField):
            assert constant_field.order().is_prime()
            self.kash_constant_field = kash.GaloisField(constant_field.order())
            self.reverse_map = {}
        elif isinstance(constant_field, NumberField):
            kZa = kash.IntegerRing().PolynomialAlgebra()
            ka = kZa.Element(list(constant_field.defining_polynomial()))
            self.kash_constant_field = ka.NumberField()
            # If constant_field is equipped with an embedding (almost surely to QQbar),
            # we map the kash generator to that target.  Otherwise, we map the kash
            # generator to the number field generator.
            if constant_field.gen_embedding():
                self.reverse_map = {self.kash_constant_field.gen(1) : constant_field.gen_embedding()}
            else:
                self.reverse_map = {self.kash_constant_field.gen(1) : constant_field.gen(0)}
        else:
            raise ValueError("The constant field must be either QQ, QQbar, a number field, or a finite field.")

        # we seem to need this to avoid getting variable names like '$.1' Kash's output
        self.kash_constant_field.PolynomialAlgebra().AssignNames_(['"x"'])

        self._kash_ = self.kash_constant_field.RationalFunctionField()

    def _extend_constant_field(self, arg):
        """
        Extend the constant field to express `arg`.
        (only for function fields over QQbar)

        Currently, `arg` can only be a FunctionFieldDivisor.
        """

        # XXX - should also be able to factor polynomials in the extension function field C(x,y)

        algebraics = []
        if self.working_constant_field != QQ:
            algebraics.append(self.working_constant_field.gen())

        if isinstance(arg, FunctionFieldDivisor):
            arg = itertools.chain.from_iterable([pls.prime_ideal().gens() for pls in arg.support()])

        for g in arg:
            if g.parent() is self:
                for r,m in g.element().numerator().change_ring(QQbar).roots():
                    algebraics.append(r)

        (constant_field, new_algebraics, nftoQQbar) = number_field_elements_from_algebraics(algebraics)

        # If we decided to expand our number field, then redo the divisor
        # computation in a new function field

        if constant_field.degree() != self.working_constant_field.degree():
            import sage.rings.number_field.number_field as number_field
            constant_field = number_field.NumberField(constant_field.polynomial(), constant_field.gen(),
                                         embedding = nftoQQbar(constant_field.gen()))
            self._set_working_constant_field(constant_field)

    def kash(self):
        return self._kash_

    def to_kash(self, c):
        x = self.kash().gen(1)
        # c will be an element of a fraction field with polynomial coefficients
        c = c.element()
        # This isn't just an optimization.  Kash won't build multivariate polynomials
        # if a constant coefficient is in a polynomial ring / fraction field.  The
        # main line code below does that; this code makes sure that constants
        # are presented as constants, and not as zero-degree polynomials.
        # XXX probably still has problems with algebraic constants like sqrt(2),
        # which fall through into the code below
        try:
            if self.constant_field() is QQbar:
                return QQ(c)
            else:
                return self.constant_field()(c)
        except:
            pass
        # Sage considers QQ to be a number field, which is why I check
        # to see if there's a reverse_map instead of just checking
        # to see if self.constant_field() is a NumberField.
        if len(self.reverse_map) > 0:
            # this is the number field case
            if c.numerator().parent().base_ring() is QQbar:
                # XXX assumes that the polynomial coefficients are in working_constant_field
                # If they're not, we should extend working_constant_field, but instead
                # we just throw an exception.  This isn't a problem if we work with
                # divisors and ideals generated by the code itself, but if the user tries to
                # create an ideal himself with a previously unseen algebraic number,
                # we'll get an exception here.
                n = c.numerator().change_ring(self.working_constant_field)
                d = c.denominator().change_ring(self.working_constant_field)
            else:
                n = c.numerator()
                d = c.denominator()
            # P will be a polynomial ring with coeffs in a number field
            P = c.numerator().parent()
            # create a bivariate ring, one variable for the polynomial var, and one for the number field generator
            R2 = QQ[P.gens() + ('a',)]
            # the number field generator will map to 'a'
            a = R2.gen(1)
            # for each coefficient, convert it to a polynomial in the number field generator, then map it into R2
            n = n.map_coefficients(lambda v: v.polynomial().subs({v.polynomial().parent().gen(0) : a}), new_base_ring = R2)
            # now map n's polynomial variable to R2's polynomial variable (n is now in R2)
            n = n.subs({n.parent().gen(0): R2.gen(0)})
            # now map n into kash
            n = n.subs({R2.gen(0) : x, R2.gen(1) : self.kash_constant_field.gen(1)})
            # ditto for the denominator
            d = d.map_coefficients(lambda v: v.polynomial().subs({v.polynomial().parent().gen(0) : a}), new_base_ring = R2)
            d = d.subs({d.parent().gen(0): R2.gen(0)})
            d = d.subs({R2.gen(0) : x, R2.gen(1) : self.kash_constant_field.gen(1)})
            return n/d
        else:
            # P will be a polynomial ring with coeffs in either QQ or QQbar
            # XXX this code assumes that even if we're in QQbar, actual coeffs are only in QQ
            # (same issue as the comment above)
            P = c.numerator().parent()
            R2 = QQ[P.gens()]
            n = c.numerator().change_ring(QQ)(x)
            d = c.denominator().change_ring(QQ)(x)
            return n/d

    def extension(self, f, names=None):
        """
        Create an extension `L = K[y]/(f(y))` of the rational function field.

        INPUT:

        - ``f`` -- univariate polynomial over self

        - ``names`` -- string or length-1 tuple

        OUTPUT:

        - a function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); R.<y> = K[]
            sage: K.extension(y^5 - x^3 - 3*x + x*y)
            Function field in y defined by y^5 + x*y - x^3 - 3*x

        A nonintegral defining polynomial::

            sage: K.<t> = FunctionField(QQ, implementation='kash'); R.<y> = K[]
            sage: K.extension(y^3 + (1/t)*y + t^3/(t+1))
            Function field in y defined by y^3 + 1/t*y + t^3/(t + 1)

        The defining polynomial need not be monic or integral::

            sage: K.extension(t*y^3 + (1/t)*y + t^3/(t+1))
            Function field in y defined by t*y^3 + 1/t*y + t^3/(t + 1)
        """
        from . import constructor
        return constructor.FunctionFieldExtension(f, names, implementation='kash')

    @cached_method
    def place_set(self):
        """
        Return the set of all places of the function field.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ, implementation='kash')
            sage: K.place_set()
            Set of places of Rational function field in t over Rational Field
        """
        from .place import PlaceSet
        return PlaceSet(self)

    @cached_method
    def space_of_differentials(self):
        """
        Return the space of differentials attached to the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(5), implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 - (x^3 - 1)/(x^3 - 2))
            sage: L.space_of_differentials()
            Space of differentials of Function field in y defined by y^3 + (4*x^3 + 1)/(x^3 + 3)

            sage: K.<t> = FunctionField(QQ, implementation='kash')
            sage: K.space_of_differentials()
            Space of differentials of Rational function field in t over Rational Field
        """
        return DifferentialsSpace_kash(self)

def divisor_decorator(f):

    def wrapper(self):

        if self.is_zero():
            raise ValueError("divisor not defined for zero")

        # if we're working over QQbar, this divisor is over the current constant field
        D = getattr(super(FunctionFieldElement_polymod_kash, self), f.__name__)()

        if self.parent().constant_base_field() is not QQbar:
            return D

        # QQbar case - compute a number field that can factor all of our support
        # polynomials in the base function field C(x).

        # If we decided to expand our number field, then redo the divisor
        # computation in a new function field

        self.parent().base_field()._extend_constant_field(D)

        D = getattr(super(FunctionFieldElement_polymod_kash, self), f.__name__)()

        return D

    return wrapper

class FunctionFieldDifferential_kash(FunctionFieldDifferential_global):

    def divisor_of_poles(self):
        """
        Return the divisor of poles of the differential.

        NOTE::
            To avoid extending the constant field unnecessarily, we first
            compute the full divisor in the current constant field, then
            discard any positive places to get the divisor of poles,
            then extend the constant field so that all of of those places
            factor, then repeat the calculation for the final result.
        """
        from .divisor import FunctionFieldDivisor
        F = self._field
        x = F.base_field().gen()
        D = self._f._divisor() + (-2) * F(x)._divisor_of_poles() + F.different()
        D = FunctionFieldDivisor(D.parent().function_field(), {p:-m for p,m in D.list() if m < 0})
        F.base_field()._extend_constant_field(D)
        D = self._f._divisor() + (-2) * F(x)._divisor_of_poles() + F.different()
        D = FunctionFieldDivisor(D.parent().function_field(), {p:-m for p,m in D.list() if m < 0})
        return D

    def residue(self, place):
        return self._field.completion(place, prec=0)(self)[-1]

class DifferentialsSpace_kash(DifferentialsSpace):
    """
    Space of differentials of a function field.
    """
    Element = FunctionFieldDifferential_kash

class FunctionFieldElement_polymod_kash(FunctionFieldElement_polymod):

    def valuation(self, place):
        """
        Return the valuation of the element at the place.

        INPUT:

        - ``place`` -- a place of the function field

        EXAMPLES::

            sage: R.<x> = FunctionField(QQbar, implementation='kash')
            sage: L.<y> = R[]
            sage: F.<y> = R.extension(y^2 - (x^2+1))
            sage: pl = F.maximal_order().ideal(x-QQbar(sqrt(-1)),y).place()
            sage: pl
            Place (x - I, y)
            sage: x.valuation(pl)
            0
            sage: y.valuation(pl)
            1
            sage: (x-QQbar(sqrt(-1))).valuation(pl)
            2
        """
        prime = place.prime_ideal()
        ideal = prime.ring().ideal(self)
        return prime.valuation(ideal)

    @divisor_decorator
    def divisor(self):
        pass

    @divisor_decorator
    def divisor_of_zeros(self):
        pass

    @divisor_decorator
    def divisor_of_poles(self):
        pass

    def _divisor(self):
        """
        Return the divisor of the element, without extending the constant subfield.
        """
        return super(FunctionFieldElement_polymod_kash, self).divisor()

    def _divisor_of_poles(self):
        """
        Return the divisor of poles of the element, without extending the constant
        subfield.
        """
        return super(FunctionFieldElement_polymod_kash, self).divisor_of_poles()


class FunctionField_polymod_kash(FunctionField_polymod):
    """
    Function fields defined by a univariate polynomial, as an extension of the
    base field, implemented using kash.

    EXAMPLES:

    We make a function field defined by a degree 5 polynomial over the
    rational function field over the rational numbers::

        sage: K.<x> = FunctionField(QQ, implementation='kash')
        sage: R.<y> = K[]
        sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x)); L
        Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x

    We next make a function field over the above nontrivial function
    field L::

        sage: S.<z> = L[]
        sage: M.<z> = L.extension(z^2 + y*z + y); M     # not tested
        Function field in z defined by z^2 + y*z + y
        sage: 1/z                                       # not tested
        ((x/(-x^4 - 1))*y^4 - 2*x^2/(-x^4 - 1))*z - 1
        sage: z * (1/z)                                 # not tested
        1

    We drill down the tower of function fields::

        sage: M.base_field()                            # not tested
        Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x
        sage: M.base_field().base_field()               # not tested
        Rational function field in x over Rational Field
        sage: M.base_field().base_field().constant_field() # not tested
        Rational Field
        sage: M.constant_base_field()                   # not tested
        Rational Field

    The polynomial must be irreducible::

        sage: K.<x>=FunctionField(QQ, implementation='kash')
        sage: R.<y> = K[]
        sage: L.<y>=K.extension(x^2-y^2)
        Traceback (most recent call last):
        ...
        ValueError: The polynomial must be irreducible
    """
    Element = FunctionFieldElement_polymod_kash

    def __init__(self, polynomial, names, category=None):
        """
        Create a function field defined as an extension of another function
        field by adjoining a root of a univariate polynomial.

        INPUT:

        - ``polynomial`` -- univariate polynomial over a function field

        - ``names`` -- tuple of length 1 or string; variable names

        - ``category`` -- category (default: category of function fields)

        EXAMPLES:

        We create an extension of a function field::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); R.<y> = K[]
            sage: L = K.extension(y^5 - x^3 - 3*x + x*y); L
            Function field in y defined by y^5 + x*y - x^3 - 3*x

        Note the type::

            sage: type(L)
            <class 'sage.rings.function_field.function_field_kash.FunctionField_polymod_kash_with_category'>

        We can set the variable name, which doesn't have to be y::

            sage: L.<w> = K.extension(y^5 - x^3 - 3*x + x*y); L
            Function field in w defined by w^5 + x*w - x^3 - 3*x

        TESTS:

        Test that :trac:`17033` is fixed::

            sage: K.<t> = FunctionField(QQ)
            sage: R.<x> = QQ[]
            sage: M.<z> = K.extension(x^7-x-t)  # not tested
            sage: M(x)                          # not tested
            z
            sage: M('z')                        # not tested
            z
            sage: M('x')                        # not tested
            Traceback (most recent call last):
            ...
            TypeError: unable to evaluate 'x' in Fraction Field of Univariate
            Polynomial Ring in t over Rational Field
        """

        FunctionField_polymod.__init__(self, polynomial, names, category)

        assert isinstance(polynomial.base_ring(), RationalFunctionField_kash)

        if not polynomial.is_irreducible():
            raise ValueError("The polynomial must be irreducible")

        # the field, without the function field structure (referenced by divisor.py)

        # don't set these here; they were set in the superclass constructor
        # self._base_field = polynomial.base_ring()
        # self._polynomial = polynomial

        # self._field = self._base_field[self._gen].extension(polynomial)
        # self._field = self._base_field[polynomial.parent().gen(0)].extension(polynomial)
        self._field = self

        self._place_class = FunctionFieldPlace_kash

        self.kash_constant_field = None

    def kash(self):

        if self.kash_constant_field != self.base_field().kash_constant_field:

            self.kash_constant_field = self.base_field().kash_constant_field

            polynomial = self.polynomial()

            kash_base_field = self.base_field().kash()
            kTy = kash_base_field.PolynomialAlgebra()

            x = kash_base_field.gen(1)
            y = kTy.gen(1)

            # Goofiness because 'subs' doesn't recurse into coefficient
            # rings, so we convert our polynomial to a list of
            # coefficients, map them to Kash, and pass the list to Kash's
            # Element constructor, which builds a polynomial from a list.

            # Also, FunctionFieldElement's aren't callable, and thus the
            # extra call to element() to convert the FunctionFieldElement
            # to a FractionFieldElement

            #kash_polynomial = kTy.Element(map(lambda c: c.element()(x), list(polynomial)))
            kash_polynomial = kTy.Element([self.base_field().to_kash(c) for c in polynomial])

            self._kash_ = kash_polynomial.FunctionField()

            ybar = self._kash_.gen(1)

            # Set up reverse map to convert elements back from kash

            self.reverse_map = {x : polynomial.base_ring().gen(0), ybar : self.gen(0)}
            self.reverse_map.update(polynomial.base_ring().reverse_map)

        return self._kash_

    def to_kash(self, polynomial):

        # Same kind of goofiness as above.  When building an element
        # of an algebraic extension, Kash requires the list to be the
        # same length as the degree of the extension.  Also, Kash's
        # Element function doesn't accept constant polynomials - they
        # have to actually be constants.

        x = self.base_ring().kash().gen(1)
        coeffs = list(self(polynomial).element())
        coeffs += [self.base_ring().zero()] * (self.degree() - len(coeffs))
        coeffs = [self.base_ring().to_kash(c) for c in coeffs]
        return self.kash().Element(coeffs)

    @cached_method
    def place_set(self):
        """
        Return the set of all places of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: L.place_set()
            Set of places of Function field in y defined by y^2 + y + (x^2 + 1)/x
        """
        from .place import PlaceSet
        return PlaceSet(self)

    @cached_method
    def space_of_differentials(self):
        """
        Return the space of differentials attached to the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(5), implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 - (x^3 - 1)/(x^3 - 2))
            sage: L.space_of_differentials()
            Space of differentials of Function field in y defined by y^3 + (4*x^3 + 1)/(x^3 + 3)

            sage: K.<t> = FunctionField(QQ, implementation='kash')
            sage: K.space_of_differentials()
            Space of differentials of Rational function field in t over Rational Field
        """
        return DifferentialsSpace_kash(self)

    @cached_method
    def divisor_group(self):
        """
        Return the group of divisors attached to the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 - (x^3 - 1)/(x^3 - 2))
            sage: L.divisor_group()
            Divisor group of Function field in y defined by y^3 + (-x^3 + 1)/(x^3 - 2)
        """
        from .divisor import DivisorGroup
        return DivisorGroup(self)

    @cached_method
    def maximal_order(self):
        """
        Return the maximal order of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash');
            sage: R.<t> = PolynomialRing(K);
            sage: F.<y> = K.extension(t^4 + x^12*t^2 + x^18*t + x^21 + x^18);
            sage: O = F.maximal_order()
            sage: O.basis()
            (1, 1/x^4*y, 1/x^9*y^2, 1/x^13*y^3)

            sage: K.<x> = FunctionField(GF(2), implementation='kash');
            sage: R.<t> = PolynomialRing(K);
            sage: F.<y> = K.extension(t^4 + x^12*t^2 + x^18*t + x^21 + x^18);
            sage: O = F.maximal_order()
            sage: O.basis()
            (1, 1/x^4*y, 1/x^11*y^2 + 1/x^2, 1/x^15*y^3 + 1/x^6*y)

        The basis of the maximal order *always* starts with 1. This is assumed
        in some algorithms.
        """
        return FunctionFieldMaximalOrder_kash(self)

    @cached_method
    def maximal_order_infinite(self):
        """
        Return the maximal infinite order of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<t> = K[]
            sage: F.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)
            sage: F.maximal_order_infinite()
            Maximal infinite order of Function field in y defined by y^3 - x^6 - 2*x^5 - 3*x^4 - 2*x^3 - x^2

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: L.maximal_order_infinite()
            Maximal infinite order of Function field in y defined by y^2 + y + (x^2 + 1)/x
        """
        return FunctionFieldMaximalOrderInfinite_kash(self)

    def different(self):
        """
        Return the different divisor of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2), implementation='kash'); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
            sage: F.different()
            2*Place (x, (1/(x^3 + x^2 + x))*y^2)
             + 2*Place (x^2 + x + 1, (1/(x^3 + x^2 + x))*y^2)
        """
        D = self.kash().DifferentDivisor()
        support = D.Support()
        data = {FunctionFieldPlace_kash(self, place) : D.Valuation(place) for place in support}
        return FunctionFieldDivisor(self, data)

    def completion(self, place, name=None, prec=None):
        """
        Return the completion of the function field at the place.

        INPUT:

        - ``place`` -- place

        - ``name`` -- string; name of the series variable

        - ``prec`` -- integer; default precision

        EXAMPLES::

            sage: R.<x> = FunctionField(QQbar, implementation='kash')
            sage: L.<y> = R[]
            sage: F.<y> = R.extension(y^2 - (x^2+1))
            sage: D = (y/x).divisor()
            sage: p = D.support()[0]
            sage: m = F.completion(p)
            sage: m
            Completion map:
              From: Function field in y defined by y^2 - x^2 - 1
              To:   Laurent Series Ring in s over Algebraic Field
            sage: m(x, 10)
            I + 2*I*s^2 + 2*I*s^4 + 2*I*s^6 + 2*I*s^8 + O(s^10)
            sage: m(y, 10)
            2*I*s + 2*I*s^3 + 2*I*s^5 + 2*I*s^7 + 2*I*s^9 + O(s^10)

        """
        return FunctionFieldCompletion_kash(self, place, name=name, prec=prec)

    def places_infinite(self, degree=1):
        """
        Return a list of the infinite places of degree ``degree``.

        INPUT:

        - ``degree`` -- positive integer (default: `1`)

        EXAMPLES::

            sage: F.<a>=GF(2)
            sage: K.<x>=FunctionField(F)
            sage: R.<t>=PolynomialRing(K)
            sage: L.<y>=K.extension(t^4+t-x^5)
            sage: L.places_infinite(1)
            [Place (1/x, 1/x^4*y^3)]
        """
        return [place for place in self._places_infinite(degree)]

    def _places_infinite(self, degree):
        """
        Return a generator of *infinite* places of the function field of the degree.

        INPUT:

        - ``degree`` -- positive integer

        EXAMPLES::

            sage: F.<a>=GF(2)
            sage: K.<x>=FunctionField(F)
            sage: R.<t>=PolynomialRing(K)
            sage: L.<y>=K.extension(t^4+t-x^5)
            sage: L._places_infinite(1)
            <generator object ...>
        """
        Oinf = self.maximal_order_infinite()
        for prime,_,_ in Oinf.decomposition():
            place = prime.place()
            if place.degree() == degree:
                yield place

class FunctionFieldCompletion_kash(FunctionFieldCompletion):
    """
    Completions on kash function fields.  Currently only supports
    QQbar as the field of constants.

    EXAMPLES::

        sage: R.<x> = FunctionField(QQbar, implementation='kash')
        sage: L.<y> = R[]
        sage: F.<y> = R.extension(y^2 - (x^2+1))
        sage: D = (y/x).divisor()
        sage: p = D.support()[0]
        sage: m = F.completion(p)
        sage: m
        Completion map:
          From: Function field in y defined by y^2 - x^2 - 1
          To:   Laurent Series Ring in s over Algebraic Field
        sage: m(x)
        I + 2*I*s^2 + 2*I*s^4 + 2*I*s^6 + 2*I*s^8 + 2*I*s^10 + 2*I*s^12 + 2*I*s^14 + 2*I*s^16 + 2*I*s^18 + O(s^20)
        sage: m(y)
        2*I*s + 2*I*s^3 + 2*I*s^5 + 2*I*s^7 + 2*I*s^9 + 2*I*s^11 + 2*I*s^13 + 2*I*s^15 + 2*I*s^17 + 2*I*s^19 + O(s^20)
        sage: m(x*y) == m(x) * m(y)
        True
        sage: m(x+y) == m(x) + m(y)
        True
        sage: m(y)^2 == m(x)^2 + 1
        True

    The variable name of the series can be supplied, as can the default precision.

        sage: p2 = D.support()[1]
        sage: p2
        Place (x, y + x - 1)
        sage: m2 = F.completion(p2, 't', prec=10)
        sage: m2(x)
        t + 1/4*t^3 + 1/16*t^5 + 1/64*t^7 + 1/256*t^9 + O(t^10)
        sage: m2(y)
        1 + 1/2*t^2 + 1/8*t^4 + 1/32*t^6 + 1/128*t^8 + O(t^10)
    """
    def __init__(self, field, place, name=None, prec=None):
        """
        Initialize.

        INPUT:

        - ``field`` -- function field

        - ``place`` -- place of the function field

        - ``name`` -- string for the name of the series variable

        - ``prec`` -- positive integer; default precision

        EXAMPLES::

            sage: R.<x> = FunctionField(QQbar, implementation='kash')
            sage: L.<y> = R[]
            sage: F.<y> = R.extension(y^2 - (x^2+1))
            sage: D = (y/x).divisor()
            sage: p = D.support()[0]
            sage: m = F.completion(p)
            sage: m
            Completion map:
              From: Function field in y defined by y^2 - x^2 - 1
              To:   Laurent Series Ring in s over Algebraic Field

        """
        from sage.rings.laurent_series_ring import LaurentSeriesRing

        if name is None:
            name = 's' # default

        # Currently we only work on places of degree one, where our residue field
        # is simply the constant base field.  In particular, this always works over QQbar.

        if place.degree() != 1:
            raise NotImplementedError("series expansions not implemented at places of degree > 1")

        # if prec is None, the Laurent series ring provides default
        # precision
        codomain = LaurentSeriesRing(field.constant_base_field(), name=name, default_prec=prec)

        FunctionFieldCompletion.__init__(self, field, codomain)

        self._place = place
        self._precision = codomain.default_prec()

    def _call_(self, f):
        """
        Call the completion for f

        EXAMPLES::

            sage: R.<x> = FunctionField(QQbar, implementation='kash')
            sage: L.<y> = R[]
            sage: F.<y> = R.extension(y^2 + y + x + 1/x)
            sage: D = (x*y).divisor()
            sage: p = D.support()[2]
            sage: m = F.completion(p, prec=10)
            sage: m
            Completion map:
              From: Function field in y defined by y^2 + y + (x^2 + 1)/x
              To:   Laurent Series Ring in s over Algebraic Field
            sage: m(x)
            -s^2 - s^3 - s^4 - s^5 - 2*s^6 - 4*s^7 - 7*s^8 - 11*s^9 + O(s^10)
        """
        return self._expand(f, prec=None)

    def _call_with_args(self, f, args=(), kwds={}):
        """
        Call the completion with ``args`` and ``kwds``.

        EXAMPLES::

            sage: R.<x> = FunctionField(QQbar, implementation='kash')
            sage: L.<y> = R[]
            sage: F.<y> = R.extension(y^2 + y + x + 1/x)
            sage: D = (x*y).divisor()
            sage: p = D.support()[2]
            sage: m = F.completion(p)
            sage: m(x+y, 10)  # indirect doctest
            -s^-1 - s^2 - s^3 - s^4 - s^5 - 2*s^6 - 4*s^7 - 7*s^8 - 11*s^9 + O(s^10)
        """
        return self._expand(f, *args, **kwds)

    def _expand(self, f, prec=None):
        """
        Return the power series representation of f of precision prec.

        INPUT:

        - ``f`` -- element of the function field

        - ``prec`` -- integer; absolute precision of the series

        OUTPUT:

        - a series of precision prec

        EXAMPLES::

            sage: R.<x> = FunctionField(QQbar, implementation='kash')
            sage: L.<y> = R[]
            sage: F.<y> = R.extension(y^2 + y + x + 1/x)
            sage: D = (x*y).divisor()
            sage: p = D.support()[2]
            sage: m = F.completion(p)
            sage: m(x, 10)  # indirect doctest
            -s^2 - s^3 - s^4 - s^5 - 2*s^6 - 4*s^7 - 7*s^8 - 11*s^9 + O(s^10)
        """
        if prec is None:
            prec = self._precision

        place = self._place
        F = place.function_field()

        kash_series = F.to_kash(f).Expand(place.kash(), AbsPrec=prec)

        val = kash_series.Valuation()
        coeffs = [c.sage(F.reverse_map) for c in list(kash_series.Coefficients())]

        return self.codomain()(coeffs, val).add_bigoh(prec)

    def pushforward(self, f, prec=None):
        """
        Allows the map to work on differentials in addition to function field elements.

        EXAMPLES::

            sage: R.<x> = FunctionField(QQbar, implementation='kash')
            sage: L.<y> = R[]
            sage: F.<y> = R.extension(y^2 - (x^2+1))
            sage: D = (x/y).divisor()
            sage: D
            -1*Place (x - I, y)
             + Place (x, y + x - 1)
             + Place (x, y + x + 1)
             - Place (x + I, y)
            sage: pl = D.support()[0]
            sage: m = F.completion(pl)
            sage: m(x.differential())
            4*I*s + 8*I*s^3 + 12*I*s^5 + 16*I*s^7 + 20*I*s^9 + 24*I*s^11 + 28*I*s^13 + 32*I*s^15 + 36*I*s^17 + 40*I*s^19 + O(s^20)
            sage: m(x/y)
            1/2*s^-1 + 1/2*s + O(s^20)
            sage: m(x/y*x.differential())
            2*I + 6*I*s^2 + 10*I*s^4 + 14*I*s^6 + 18*I*s^8 + 22*I*s^10 + 26*I*s^12 + 30*I*s^14 + 34*I*s^16 + 38*I*s^18 + O(s^20)

            sage: K.<x> = FunctionField(QQbar, implementation='kash')
            sage: L.<y> = K[]
            sage: F.<y> = K.extension(y^2-x+1)
            sage: pl = y.divisor().support()[1]
            sage: m = F.completion(pl, prec=1)
            sage: m(y*x.differential())
            O(s^1)

        """
        if prec is None:
            prec = self._precision

        t = f._field.base_field().gen()       # all differentials construction w.r.t this differential
        vt = t.valuation(self._place)         # t's valuation
        vf = f._f.valuation(self._place)      # f's valuation
        st = self._expand(t, max(prec-vf,vt)+1) # t's series
        sdt = st.derivative()                 # dt's series
        sf = self._expand(f._f, prec-(vt-1))  # f's series
        return sf * sdt

    def default_precision(self):
        """
        Return the default precision.

        EXAMPLES::

            sage: R.<x> = FunctionField(QQbar, implementation='kash')
            sage: L.<y> = R[]
            sage: F.<y> = R.extension(y^2 + y + x + 1/x)
            sage: D = (x*y).divisor()
            sage: p = D.support()[2]
            sage: m = F.completion(p)
            sage: m.default_precision()
            20
        """
        return self._precision


class FunctionFieldMaximalOrder_kash(FunctionFieldMaximalOrder):
    """
    Base class of kash-implemented maximal orders of function fields.
    """

    def __init__(self, field, category=None):
        """
        Initialize.

        INPUT:

        - ``field`` -- function field

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ, implementation='kash'); K
            Rational function field in t over Rational Field
            sage: R = K.maximal_order(); R
            Maximal order of Rational function field in t over Rational Field
        """

        FunctionFieldMaximalOrder.__init__(self, field, category)
        self.kash_constant_field = None

    def kash(self):
        if self.kash_constant_field != self.function_field().base_field().kash_constant_field:
            self.kash_constant_field = self.function_field().base_field().kash_constant_field
            self._kash_ = self.function_field().kash().MaximalOrderFinite()
        return self._kash_

    def _element_constructor_(self, f, check=True):
        """
        Construct an element of this order from ``f``.

        INPUT:

        - ``f`` -- element convertible to the function field

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2), implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2-x*Y+x^2+1)
            sage: O = L.maximal_order()
            sage: y in O
            True
            sage: 1/y in O
            False
            sage: x in O
            True
            sage: 1/x in O
            False
            sage: L.<y>=K.extension(Y^2+Y+x+1/x)
            sage: O = L.maximal_order()
            sage: 1 in O
            True
            sage: y in O
            False
            sage: x*y in O
            True
            sage: x^2*y in O
            True
        """
        field = self.function_field()

        #if f.parent() is field:
        #    f = f.element()
        #f = self._field._ring(f)
        if check:
            if not kash._contains(self._field.to_kash(f).name(), self.kash().name()):
                raise TypeError("%r is not an element of %r"%(f,self))
        # return f and not field._element_class(self, f) because
        # 1. that's what FunctionFieldMaximalOrder_global does
        # 2. that's what makes "d in O" work (it tests if O(d) == d)
        # 3. but it's probably not right (O(d).parent() should report the order, not the field)
        # return field._element_class(self, f)
        return f

    def polynomial(self):
        """
        Return the defining polynomial of the function field of which this is an order.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.equation_order()
            sage: O.polynomial()
            y^4 + x*y + 4*x + 1
        """
        return self._field.polynomial()

    def basis(self):
        """
        Return the basis of the order as a module over the polynomial ring.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.equation_order()
            sage: O.basis()
            (1, y, y^2, y^3)
        """

        return tuple(self.kash().Basis().sage(self._field.reverse_map))

    def ideal(self, *gens):
        """
        Return the fractional ideal generated by ``gens``.

        INPUT:

        - ``gens`` -- list of generators or an ideal in a ring which
          coerces to this order

        EXAMPLES::

            sage: K.<y> = FunctionField(QQ, implementation='kash')
            sage: O = K.maximal_order()
            sage: O.ideal(y)
            Ideal (y) of Maximal order of Rational function field in y over Rational Field
            sage: O.ideal([y,1/y]) == O.ideal(y,1/y) # multiple generators may be given as a list
            True

        A fractional ideal of a nontrivial extension::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); R.<y> = K[]
            sage: O = K.maximal_order()
            sage: I = O.ideal(x^2-4)
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: S = L.equation_order()
            sage: S.ideal(1/y)
            Ideal (1, (-1/(x^3 + 1))*y) of Order in Function field in y defined by y^2 - x^3 - 1
            sage: I2 = S.ideal(x^2-4); I2
            Ideal (x^2 - 4, (x^2 - 4)*y) of Order in Function field in y defined by y^2 - x^3 - 1
            sage: I2 == S.ideal(I)
            True

            sage: K.<x> = FunctionField(GF(7), implementation='kash'); R.<y> = K[]
            sage: O = K.maximal_order()
            sage: I = O.ideal(x^2-4)
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: S = L.equation_order()
            sage: S.ideal(1/y)
            Ideal (1, (6/(x^3 + 1))*y) of Order in Function field in y defined by y^2 + 6*x^3 + 6
            sage: I2 = S.ideal(x^2-4); I2
            Ideal (x^2 + 3, (x^2 + 3)*y) of Order in Function field in y defined by y^2 + 6*x^3 + 6
            sage: I2 == S.ideal(I)
            True
        """

        return FunctionFieldIdeal_kash(self, gens)

class FunctionFieldMaximalOrderInfinite_kash(FunctionFieldMaximalOrderInfinite):
    """
    Base class of kash-implemented maximal infinite orders of function fields.
    """

    def __init__(self, field, category=None):
        """
        Initialize.

        INPUT:

        - ``field`` -- function field

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ, implementation='kash'); K
            Rational function field in t over Rational Field
            sage: R = K.maximal_order_infinite(); R
            Maximal infinite order of Rational function field in t over Rational Field
        """

        FunctionFieldMaximalOrderInfinite.__init__(self, field, category)
        self.kash_constant_field = None

    def kash(self):
        if self.kash_constant_field != self.function_field().base_field().kash_constant_field:
            self.kash_constant_field = self.function_field().base_field().kash_constant_field
            self._kash_ = self.function_field().kash().MaximalOrderInfinite()
        return self._kash_

    def _element_constructor_(self, f):
        """
        Make ``f`` an element of this order.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2), implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: Oinf.basis()
            (1, 1/x*y)
            sage: 1 in Oinf
            True
            sage: 1/x*y in Oinf
            True
            sage: x*y in Oinf
            False
            sage: 1/x in Oinf
            True
        """
        if not f.parent() is self.function_field():
            f = self.function_field()(f)

        if not kash._contains(self._field.to_kash(f).name(), self.kash().name()):
            raise TypeError("%r is not an element of %r"%(f,self))

        return f

    def basis(self):
        """
        Return a basis of this order as a module over the maximal order
        of the base function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<t> = K[]
            sage: L.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)
            sage: Oinf = L.maximal_order_infinite()
            sage: Oinf.basis()
            (1, 1/x^2*y, 1/x^4*y^2)

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: Oinf.basis()
            (1, 1/x*y)

            sage: K.<x> = FunctionField(GF(2), implementation='kash'); _.<t> = K[]
            sage: L.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)
            sage: Oinf = L.maximal_order_infinite()
            sage: Oinf.basis()    # not tested - kash returns a different (but equivalent) basis
            (1, 1/x^2*y, 1/x^4*y^2)

            sage: K.<x> = FunctionField(GF(2), implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: Oinf.basis()
            (1, 1/x*y)
        """

        return tuple(self.kash().Basis().sage(self._field.reverse_map))

    def ideal(self, *gens):
        """
        Return the fractional ideal generated by ``gens``.

        INPUT:

        - ``gens`` -- tuple of elements of the function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<t> = K[]
            sage: F.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)
            sage: Oinf = F.maximal_order_infinite()
            sage: I = Oinf.ideal(x,y); I
            Ideal (x, y) of Maximal infinite order of Function field
            in y defined by y^3 - x^6 - 2*x^5 - 3*x^4 - 2*x^3 - x^2

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: I = Oinf.ideal(x,y); I
            Ideal (x, y) of Maximal infinite order of Function field
            in y defined by y^2 + y + (x^2 + 1)/x
        """

        return FunctionFieldIdeal_kash(self, gens)

    def decomposition(self):
        """
        Return prime ideal decomposition of `pO_\infty` where `p` is the unique
        prime ideal of the maximal infinite order.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<t> = K[]
            sage: F.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)
            sage: Oinf = F.maximal_order_infinite()
            sage: Oinf.decomposition()
            [(Ideal (1/x, 1/x^2*y - 1) of Maximal infinite order
             of Function field in y defined by y^3 - x^6 - 2*x^5 - 3*x^4 - 2*x^3 - x^2, 1, 1),
             (Ideal (1/x, 1/x^4*y^2 + 1/x^2*y + 1) of Maximal infinite order
             of Function field in y defined by y^3 - x^6 - 2*x^5 - 3*x^4 - 2*x^3 - x^2, 2, 1)]

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: Oinf.decomposition()
            [(Ideal (1/x, 1/x*y) of Maximal infinite order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x, 1, 2)]
        """

        return [ ( FunctionFieldIdeal_kash(self, i), i.Degree().sage(), i.RamificationIndex().sage() )
                 for i in self.kash().Decomposition()]

class FunctionFieldIdeal_kash(FunctionFieldIdeal):
    """
    Fractional ideal of the maximal order of a kash-implemented function field.
    """

    def __init__(self, ring, gens):
        """
        Initialize.

        INPUT:

        - ``ring`` -- maximal order

        - ``gens``-- a list of generators, or a kash ideal
        """

        FunctionFieldIdeal.__init__(self, ring)

        if isinstance(gens, KashElement):
            self._kash_ = gens
            self._gens = tuple(gens.Generators().sage(ring._field.reverse_map))
        else:
            if ring.function_field().constant_base_field() is QQbar:
                ring.function_field().base_field()._extend_constant_field(gens)
            self._gens = tuple(flatten(gens))
            self._kash_ = ring.kash().Ideal(map(ring._field.to_kash, self._gens))

        self.kash_constant_field = ring.function_field().kash_constant_field

    def kash(self):
        if self.kash_constant_field != self.ring().function_field().base_field().kash_constant_field:
            self.kash_constant_field = self.ring().function_field().base_field().kash_constant_field
            self._kash_ = self.ring().kash().Ideal(map(self.ring().function_field().to_kash, self._gens))
        return self._kash_

    def __hash__(self):
        """
        Return hash computed from the data.
        """

        return hash( (self._ring, self._gens) )

    def _richcmp_(self, other, op):
        """
        Compare the ideal with the other ideal with respect to ``op``.

        INPUT:

        - ``other`` -- ideal

        - ``op`` -- comparison operator

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(1/y)
            sage: I == I + I
            True
            sage: I == I * I
            False

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(1/y)
            sage: I == I + I
            True
            sage: I == I * I
            False
            sage: I < I * I
            True
            sage: I > I * I
            False
        """

        return richcmp((self.denominator(), self.gens_over_base()), (other.denominator(), other.gens_over_base()), op)

    def __repr__(self):
        """
        Return a string representation of the ideal.
        """
        return "Ideal (%s) of %s"%(', '.join([repr(g) for g in self.gens()]), self.ring())

    def __contains__(self, x):
        """
        Return ``True`` if ``x`` is in the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3 - 1)
            sage: O = L.maximal_order()
            sage: I = O.ideal([y]); I
            Ideal (y) of Maximal order of Function field in y
            defined by y^2 - x^3 - 1
            sage: x * y in I
            True
            sage: y / x in I
            False
            sage: y^2 - 2 in I
            False

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal([y]); I
            Ideal (y) of Maximal order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x
            sage: x * y in I
            True
            sage: y / x in I
            False
            sage: y^2 - 2 in I
            False
        """
        kashx = self._ring._field.to_kash(x)
        return kash.eval("%s in %s" % (kashx.name(), self.kash().name())) == 'TRUE'

    def __invert__(self):
        """
        Return the inverse fractional ideal of the ideal.

        EXAMPLES::
            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3 - 1)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: ~I
            Ideal (1, (1/(x^3 + 1))*y) of Maximal order of Function field in y defined by y^2 - x^3 - 1
            sage: I^(-1)
            Ideal (1, (1/(x^3 + 1))*y) of Maximal order of Function field in y defined by y^2 - x^3 - 1
            sage: ~I * I
            Ideal (1) of Maximal order of Function field in y defined by y^2 - x^3 - 1

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y>=K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: ~I
            Ideal (x, (x/(x^2 + 1))*y + x/(x^2 + 1)) of Maximal order
            of Function field in y defined by y^2 + y + (x^2 + 1)/x
            sage: I^(-1)
            Ideal (x, (x/(x^2 + 1))*y + x/(x^2 + 1)) of Maximal order
            of Function field in y defined by y^2 + y + (x^2 + 1)/x
            sage: ~I * I
            Ideal (1) of Maximal order of Function field in y defined by y^2 + y + (x^2 + 1)/x
        """

        return FunctionFieldIdeal_kash(self._ring, 1 / self.kash())

    def _add_(self, other):
        """
        Add with other ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: J = O.ideal(x+y)
            sage: I + J
            Ideal (x, y) of Maximal order of Function field in y defined by y^2 - x^3*y - x

            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: J = O.ideal(x+y)
            sage: I + J
            Ideal (1, y) of Maximal order of Function field in y defined by y^2 + y + (x^2 + 1)/x
        """

        return FunctionFieldIdeal_kash(self._ring, self.kash() + other.kash())

    def _mul_(self, other):
        """
        Multiply with other ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: J = O.ideal(x+y)
            sage: I * J
            Ideal (x^4 + x^2 - x, x*y + x^2) of Maximal order
            of Function field in y defined by y^2 - x^3*y - x

            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: J = O.ideal(x+y)
            sage: I * J
            Ideal ((x^5 + x^3 + x^2 + 1)/x,
                   y + (1/2*x^4 + 1/2*x^3 + x^2 + 1/2*x + 1/2)/x)
                of Maximal order of Function field in y
                defined by y^2 + y + (x^2 + 1)/x

        TESTS:

        Verify the examples using Sage's standard ideal operations::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: J = O.ideal(x+y)
            sage: IJ = I * J

            sage: R = PolynomialRing(QQ, K.gens() + L.gens())
            sage: Q = R.quo(L.polynomial())
            sage: I1 = ideal(Q(g*b) for g in I.gens() for b in O.basis())
            sage: J1 = ideal(Q(g*b) for g in J.gens() for b in O.basis())
            sage: IJ = ideal(Q(g*b) for g in (IJ).gens() for b in O.basis())
            sage: I1 * J1 == IJ
            True

            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: J = O.ideal(x+y)
            sage: IJ = I * J

            sage: Q = R.quo(L.polynomial().numerator())
            sage: I1 = ideal(Q(g*b*x) for g in I.gens() for b in O.basis())
            sage: J1 = ideal(Q(g*b) for g in J.gens() for b in O.basis())
            sage: IJ = ideal(Q(g*b*x) for g in (IJ).gens() for b in O.basis())
            sage: I1 * J1 == IJ
            True
        """

        return FunctionFieldIdeal_kash(self._ring, self.kash() * other.kash())

    def denominator(self):
        """
        Return the denominator of the fractional ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y/(y+1))
            sage: d = I.denominator(); d
            x^3
            sage: d in O
            True
        """
        return self.kash().Denominator().sage(self._ring._field.reverse_map)

    def factor(self):
        """
        Return the factorization of the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(y)
            sage: I == I.factor().prod()
            True

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: I == I.factor().prod()
            True
        """
        factors = self._factor()
        return Factorization(factors, cr=True)

    def _factor(self):
        """
        Return the list of prime and multiplicity pairs of the
        factorization of the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<t> = K[]
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(y)
            sage: I == I.factor().prod()  # indirect doctest
            True
        """
        factors = []
        for f,m in self.kash().Factorization():
            factors.append( (FunctionFieldIdeal_kash(self._ring, f), m) )
        return factors

    def gens(self):
        """
        Return the generators of the ideal.

        This provides whatever set of generators as quickly
        as possible.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(x+y)
            sage: I.gens()
            (y + x,)

            sage: L.<y> = K.extension(Y^2 +Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(x+y)
            sage: I.gens()
            (y + x,)
        """
        return self._gens

    def gens_over_base(self):
        """
        Return the generators of the ideal as a module over the
        maximal order of the base rational function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(x+y)
            sage: I.gens_over_base()
            (x^4 + x^2 - x, y + x)

            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(x+y)
            sage: I.gens_over_base()
            (x^3 + 1, y + x)
        """
        return tuple(self.kash().Basis().sage(self._ring._field.reverse_map))

    def is_prime(self):
        """
        Return ``True`` if the ideal is a prime ideal.

        If checked to be a prime ideal, then the ideal can be used
        as a prime ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(y)
            sage: I.is_prime()
            False
            sage: [f.is_prime() for f,_ in I.factor()]
            [True, True]

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: I.is_prime()
            False
            sage: [f.is_prime() for f,_ in I.factor()]
            [True, True]
        """

        return bool(self.kash().IsPrime())

    def place(self):
        """
        Return the place corresponding to the prime ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(y)
            sage: [f.place() for f,_ in I.factor()]
            [Place (x, (1/(x^3 + x^2 + x))*y^2),
             Place (x^2 + x + 1, (1/(x^3 + x^2 + x))*y^2)]

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: [f.place() for f,_ in I.factor()]
            [Place (x, x*y), Place (x^2 + 1, x*y)]
        """
        if not self.is_prime():
            raise TypeError("not a prime ideal")
        return FunctionFieldPlace_kash(self._ring._field, self)

    def valuation(self, ideal):
        """
        Return the valuation of the ideal at this prime ideal.

        INPUT:

        - ``ideal`` -- fractional ideal

        EXAMPLES::

            sage: F.<x> = FunctionField(QQ, implementation='kash')
            sage: O = F.maximal_order()
            sage: I = O.ideal(x^2*(x^2+x+1)^3)
            sage: [f.valuation(I) for f,_ in I.factor()]
            [2, 3]
        """
        if not self.is_prime():
            raise TypeError("not a prime ideal")

        return ideal.kash().Valuation(self.kash())


class FunctionFieldPlace_kash(FunctionFieldPlace):
    """
    Places of kash-implemented function field.
    """

    def __init__(self, field, arg):
        """
        Initialize the place.

        INPUT:

        - ``field`` -- function field

        - ``arg`` -- prime ideal associated with the place, or a kash place
        """

        if isinstance(arg, KashElement):
            ideal = arg.Ideal()
            order = field.maximal_order()
            prime = FunctionFieldIdeal_kash(order, order.kash().CoerceIdeal(ideal))
            self._kash_ = arg
        else:
            prime = arg
            self._kash_ = arg.kash().Place()

        FunctionFieldPlace.__init__(self, field, prime)
        self.kash_constant_field = field.base_field().kash_constant_field

    def kash(self):
        if self.kash_constant_field != self.function_field().base_field().kash_constant_field:
            self.kash_constant_field = self.function_field().base_field().kash_constant_field
            self._kash_ = self.prime_ideal().kash().Place()
        return self._kash_

    def is_infinite_place(self):
        """
        Return ``True`` if the place is at infinity.
        """
        F = self.function_field()
        return self.prime_ideal().ring() == F.maximal_order_infinite()

    def degree(self):
        return self.kash().Degree()
