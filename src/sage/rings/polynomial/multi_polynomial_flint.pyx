"""
Sparse multivariate polynomials over `\ZZ`, implemented using FLINT
"""

import sys

import itertools

from collections import Counter

from sage.rings.integer_ring import ZZ

from sage.rings.polynomial.polydict cimport ETuple

from sage.rings.fraction_field_element import is_FractionFieldElement

from sage.structure.element import coerce_binop
from sage.structure.element cimport Element, CommutativeRingElement
from sage.structure.richcmp cimport rich_to_bool
from sage.structure.category_object cimport CategoryObject

from sage.structure.factorization import Factorization

from sage.libs.gmp.mpz cimport mpz_get_ui, mpz_get_si

from sage.libs.flint.mpoly cimport *
from sage.libs.flint.fmpz cimport *
from sage.libs.flint.fmpz_mpoly cimport *
from sage.libs.flint.fmpz_mpoly_factor cimport *

from libc.stdlib cimport malloc, free

from sage.misc.sage_eval import sage_eval
from sage.misc.misc_c import prod

from sage.cpython.string cimport char_to_str, str_to_bytes

from functools import reduce

cdef class MPolynomialRing_flint(MPolynomialRing_base):

    def __init__(self, base_ring, n, names, order='degrevlex'):
        """
        Construct a multivariate polynomial ring subject to the
        following conditions:

        INPUT:

        - ``base_ring`` - base ring (must be ZZ)

        - ``n`` - number of variables (must be at least 1)

        - ``names`` - names of ring variables, may be string of list/tuple

        - ``order`` - term order (must be ``degrevlex``; FLINT also supports ``lex`` and ``deglex``)

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: P
            Multivariate Polynomial Ring in x, y, z over Rational Field

            sage: f = 27/113 * x^2 + y*z + 1/2; f
            27/113*x^2 + y*z + 1/2

            sage: P.term_order()
            Degree reverse lexicographic term order

            sage: P = PolynomialRing(GF(127),3,names='abc', order='lex')
            sage: P
            Multivariate Polynomial Ring in a, b, c over Finite Field of size 127

            sage: a,b,c = P.gens()
            sage: f = 57 * a^2*b + 43 * c + 1; f
            57*a^2*b + 43*c + 1

            sage: P.term_order()
            Lexicographic term order

            sage: z = QQ['z'].0
            sage: K.<s> = NumberField(z^2 - 2)
            sage: P.<x,y> = PolynomialRing(K, 2)
            sage: 1/2*s*x^2 + 3/4*s
            (1/2*s)*x^2 + (3/4*s)

            sage: P.<x,y,z> = ZZ[]; P
            Multivariate Polynomial Ring in x, y, z over Integer Ring

            sage: P.<x,y,z> = Zmod(2^10)[]; P
            Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 1024

            sage: P.<x,y,z> = Zmod(3^10)[]; P
            Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 59049

            sage: P.<x,y,z> = Zmod(2^100)[]; P
            Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 1267650600228229401496703205376

            sage: P.<x,y,z> = Zmod(2521352)[]; P
            Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 2521352
            sage: type(P)
            <type 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomialRing_libsingular'>

            sage: P.<x,y,z> = Zmod(25213521351515232)[]; P
            Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 25213521351515232
            sage: type(P)
            <class 'sage.rings.polynomial.multi_polynomial_ring.MPolynomialRing_polydict_with_category'>

            sage: P.<x,y,z> = PolynomialRing(Integers(2^32),order='lex')
            sage: P(2^32-1)
            4294967295

        TESTS:

        Make sure that a faster coercion map from the base ring is used;
        see :trac:`9944`::

            sage: R.<x,y> = PolynomialRing(ZZ)
            sage: R.coerce_map_from(R.base_ring())
            Polynomial base injection morphism:
              From: Integer Ring
              To:   Multivariate Polynomial Ring in x, y over Integer Ring

        Check some invalid input::

            sage: from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular
            sage: MPolynomialRing_libsingular(Zmod(1), 1, ["x"], "lex")
            Traceback (most recent call last):
            ...
            NotImplementedError: polynomials over Ring of integers modulo 1 are not supported in Singular
            sage: MPolynomialRing_libsingular(SR, 1, ["x"], "lex")
            Traceback (most recent call last):
            ...
            NotImplementedError: polynomials over Symbolic Ring are not supported in Singular
            sage: MPolynomialRing_libsingular(QQ, 0, [], "lex")
            Traceback (most recent call last):
            ...
            NotImplementedError: polynomials in 0 variables are not supported in Singular
            sage: MPolynomialRing_libsingular(QQ, -1, [], "lex")
            Traceback (most recent call last):
            ...
            NotImplementedError: polynomials in -1 variables are not supported in Singular
        """

        if base_ring != ZZ:
            raise NotImplementedError("FLINT base_ring must be ZZ")

        if order != 'degrevlex':
            raise NotImplementedError("FLINT orderings must be degrevlex")

        fmpz_mpoly_ctx_init(self._ctx, n, ORD_DEGREVLEX);

        self._cnames = <const char **>malloc(sizeof(char *) * n)
        self._bnames = []
        for i in range(n):
            self._bnames.append(str_to_bytes(names[i]))
            self._cnames[i] = self._bnames[i]

        MPolynomialRing_base.__init__(self, base_ring, n, names, order)

    def __dealloc__(self):
        free(<void *>self._cnames)

    cpdef _coerce_map_from_(self, other):
        """
        Return True if and only if there exists a coercion map from
        ``other`` to ``self``.

        TESTS::

            sage: R.<x,y> = QQ[]
            sage: type(R)
            <type 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomialRing_libsingular'>
            sage: R.has_coerce_map_from(ZZ['t'])
            False
            sage: R.coerce_map_from(ZZ['x'])
            Coercion map:
              From: Univariate Polynomial Ring in x over Integer Ring
              To:   Multivariate Polynomial Ring in x, y over Rational Field

        """
        if other == ZZ:
            return True
        else:
            return False
#        base_ring = self.base_ring()
#        if other is base_ring:
#            # Because this parent class is a Cython class, the method
#            # UnitalAlgebras.ParentMethods.__init_extra__(), which normally
#            # registers the coercion map from the base ring, is called only
#            # when inheriting from this class in Python (cf. Trac #26958).
#            return self._coerce_map_from_base_ring()
#        f = self._coerce_map_via([base_ring], other)
#        if f is not None:
#            return f
#
#        if isinstance(other, MPolynomialRing_libsingular):
#            if self is other:
#                return True
#            n = other.ngens()
#            if(other.base_ring is base_ring and self.ngens() >= n and
#               self.variable_names()[:n] == other.variable_names()):
#                return True
#            elif base_ring.has_coerce_map_from(other._mpoly_base_ring(self.variable_names())):
#                return True
#        elif isinstance(other, MPolynomialRing_polydict):
#            if self == other:
#                return True
#            elif other.ngens() == 0:
#                return True
#            elif base_ring.has_coerce_map_from(other._mpoly_base_ring(self.variable_names())):
#                return True
#        elif is_PolynomialRing(other):
#            if base_ring.has_coerce_map_from(other._mpoly_base_ring(self.variable_names())):
#                return True

    Element = MPolynomial_flint

    def _element_constructor_(self, element, check=True):
        """
        Construct a new element in this polynomial ring by converting
        ``element`` into ``self`` if possible.

        INPUT:

        - ``element`` -- several types are supported, see below

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]

        We can coerce elements of self to self::

            sage: P._coerce_(x*y + 1/2)
            x*y + 1/2

        We can coerce elements for a ring with the same algebraic properties::

            sage: from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular
            sage: R.<x,y,z> = MPolynomialRing_libsingular(QQ,3)
            sage: P == R
            True

            sage: P is R
            False

            sage: P._coerce_(x*y + 1)
            x*y + 1

        We can coerce base ring elements::

            sage: P._coerce_(3/2)
            3/2

        and all kinds of integers::

            sage: P._coerce_(ZZ(1))
            1

            sage: P._coerce_(int(1))
            1

            sage: k.<a> = GF(2^8)
            sage: P.<x,y> = PolynomialRing(k,2)
            sage: P._coerce_(a)
            (a)

            sage: z = QQ['z'].0
            sage: K.<s> = NumberField(z^2 - 2)
            sage: P.<x,y> = PolynomialRing(K, 2)
            sage: P._coerce_(1/2*s)
            (1/2*s)

        TESTS::

            sage: P.<x,y> = PolynomialRing(GF(127))
            sage: P("111111111111111111111111111111111111111111111111111111111")
            21
            sage: P.<x,y> = PolynomialRing(QQ)
            sage: P("111111111111111111111111111111111111111111111111111111111")
            111111111111111111111111111111111111111111111111111111111
            sage: P("31367566080")
            31367566080

        Check if :trac:`7582` is fixed::

            sage: R.<x,y,z> = PolynomialRing(CyclotomicField(2),3)
            sage: R.coerce(1)
            1

        Check if :trac:`6160` is fixed::

            sage: x=var('x')
            sage: K.<j> = NumberField(x-1728)
            sage: R.<b,c> = K[]
            sage: R.coerce(1)
            1

        Check if coercion from zero variable polynomial rings work
        (:trac:`7951`)::

            sage: P = PolynomialRing(QQ,0,'')
            sage: R.<x,y> = QQ[]
            sage: P(5)*x
            5*x
            sage: P = PolynomialRing(ZZ,0,'')
            sage: R.<x,y> = GF(127)[]
            sage: R.coerce(P(5))
            5

        Conversion from strings::

            sage: P.<x,y,z> = QQ[]
            sage: P('x+y + 1/4')
            x + y + 1/4

        Coercion from SINGULAR elements::

            sage: P._singular_()
            polynomial ring, over a field, global ordering
            //   coefficients: QQ
            //   number of vars : 3
            //        block   1 : ordering dp
            //                  : names    x y z
            //        block   2 : ordering C

            sage: P._singular_().set_ring()
            sage: P(singular('x + 3/4'))
            x + 3/4

        Coercion from symbolic variables::

            sage: R = QQ['x,y,z']
            sage: var('x')
            x
            sage: R(x)
            x

        Coercion from 'similar' rings, which maps by index::

            sage: P.<x,y,z> = QQ[]
            sage: R.<a,b,c> = ZZ[]
            sage: P(a)
            x

        ::

            sage: P.<x,y> = QQ[]
            sage: R.<a,b,c> = QQ[]
            sage: R(x)
            a

        Coercion from PARI objects::

            sage: P.<x,y,z> = QQ[]
            sage: P(pari('x^2 + y'))
            x^2 + y
            sage: P(pari('x*y'))
            x*y

        Coercion from boolean polynomials, also by index::

            sage: B.<x,y,z> = BooleanPolynomialRing(3)
            sage: P.<x,y,z> = QQ[]
            sage: P(B.gen(0))
            x

        If everything else fails, we try to coerce to the base ring::

            sage: R.<x,y,z> = GF(3)[]
            sage: R(1/2)
            -1

        Finally, conversions from other polynomial rings which are not
        coercions are provided. Variables are mapped as follows. Say,
        we are mapping an element from `P` to `Q` (this ring). If the
        variables of `P` are a subset of `Q`, we perform a name
        preserving conversion::

            sage: P.<y_2, y_1, z_3, z_2, z_1> = GF(3)[]
            sage: Q = GF(3)['y_4', 'y_3', 'y_2', 'y_1', 'z_5', 'z_4', 'z_3', 'z_2', 'z_1']
            sage: Q(y_1*z_2^2*z_1)
            y_1*z_2^2*z_1

        Otherwise, if `P` has less than or equal the number of
        variables as `Q`, we perform a conversion by index::

            sage: P.<a,b,c> = GF(2)[]
            sage: Q = GF(2)['c','b','d','e']
            sage: f = Q.convert_map_from(P)
            sage: f(a), f(b), f(c)
            (c, b, d)

        ::

            sage: P.<a,b,c> = GF(2)[]
            sage: Q = GF(2)['c','b','d']
            sage: f = Q.convert_map_from(P)
            sage: f(a),f(b),f(c)
            (c, b, d)

        In all other cases, we fail::

            sage: P.<a,b,c,f> = GF(2)[]
            sage: Q = GF(2)['c','d','e']
            sage: f = Q.convert_map_from(P)
            sage: f(a)
            Traceback (most recent call last):
            ...
            TypeError: Could not find a mapping of the passed element to this ring.

        Coerce in a polydict where a coefficient reduces to 0 but isn't 0. ::

            sage: R.<x,y> = QQ[]; S.<xx,yy> = GF(5)[]; S( (5*x*y + x + 17*y)._mpoly_dict_recursive() )
            xx + 2*yy

        Coerce in a polynomial one of whose coefficients reduces to 0. ::

            sage: R.<x,y> = QQ[]; S.<xx,yy> = GF(5)[]; S(5*x*y + x + 17*y)
            xx + 2*yy

        Some other examples that illustrate the same coercion idea::

            sage: R.<x,y> = ZZ[]
            sage: S.<xx,yy> = GF(25,'a')[]
            sage: S(5*x*y + x + 17*y)
            xx + 2*yy

            sage: S.<xx,yy> = Integers(5)[]
            sage: S(5*x*y + x + 17*y)
            xx + 2*yy

        See :trac:`5292`::

            sage: R.<x> = QQ[]; S.<q,t> = QQ[]; F = FractionField(S)
            sage: x in S
            False
            sage: x in F
            False

        Check if :trac:`8228` is fixed::

            sage: P.<x,y> = Zmod(10)[]; P(0)
            0
            sage: P.<x,y> = Zmod(2^10)[]; P(0)
            0

        And :trac:`7597` is fixed if this does not segfault::

            sage: F2 = GF(2)
            sage: F.<x> = GF(2^8)
            sage: R4.<a,b> = PolynomialRing(F)
            sage: R.<u,v> = PolynomialRing(F2)
            sage: P = a
            sage: (P(0,0).polynomial()[0])*u
            0
            sage: P(a,b)
            a

        Check that :trac:`15746` is fixed::

            sage: R.<x,y> = GF(7)[]
            sage: R(2^31)
            2

        Check that :trac:`17964` is fixed::

            sage: K.<a> = QuadraticField(17)
            sage: Q.<x,y> = K[]
            sage: f = (-3*a)*y + (5*a)
            sage: p = K.primes_above(5)[0]
            sage: R = K.residue_field(p)
            sage: S = R['x','y']
            sage: S(f)
            (2*abar)*y

        """

        cdef MPolynomial_flint p

        base_ring = self.base_ring()

        if isinstance(element, MPolynomial_flint):
            n = (<MPolynomial_flint>element)._parent.ngens()
            if element.parent() is self:
                return element
#            elif(base_ring is element.base_ring() and
#                 self.ngens() >= n and
#                 self.variable_names()[:n] == (<MPolynomial_libsingular>element)._parent.variable_names()):
#                if self.term_order() == (<MPolynomial_libsingular>element)._parent.term_order():
#                    _p = prCopyR_NoSort((<MPolynomial_libsingular>element)._poly,
#                                        (<MPolynomial_libsingular>element)._parent_ring,
#                                        _ring)
#                else:
#                    _p = prCopyR((<MPolynomial_libsingular>element)._poly,
#                                 (<MPolynomial_libsingular>element)._parent_ring, _ring)
#                return new_MP(self, _p)
#            elif base_ring.has_coerce_map_from(element.parent()._mpoly_base_ring(self.variable_names())):
#                return self(element._mpoly_dict_recursive(self.variable_names(), base_ring))
#
#        elif isinstance(element, MPolynomial_polydict):
#            if element.parent() == self:
#                bucket = sBucketCreate(_ring)
#                try:
#                    for (m,c) in element.element().dict().iteritems():
#                        mon = p_Init(_ring)
#                        p_SetCoeff(mon, sa2si(c, _ring), _ring)
#                        for pos in m.nonzero_positions():
#                            overflow_check(m[pos], _ring)
#                            p_SetExp(mon, pos+1, m[pos], _ring)
#                        p_Setm(mon, _ring)
#                        #we can use "_m" because we're merging a monomial and
#                        #"Merge" because this monomial is different from the rest
#                        sBucket_Merge_m(bucket, mon)
#                    e=0
#                    #we can use "Merge" because the monomials are distinct
#                    sBucketClearMerge(bucket, &_p, &e)
#                    sBucketDestroy(&bucket)
#                except:
#                     sBucketDeleteAndDestroy(&bucket)
#                     raise
#                return new_MP(self, _p)
#            elif element.parent().ngens() == 0:
#                # zero variable polynomials
#                _p = p_NSet(sa2si(base_ring(element[tuple()]), _ring),
#                        _ring)
#                return new_MP(self, _p)
#            elif base_ring.has_coerce_map_from(element.parent()._mpoly_base_ring(self.variable_names())):
#                return self(element._mpoly_dict_recursive(self.variable_names(), base_ring))
#
#        elif isinstance(element, polynomial_element.Polynomial):
#            if base_ring.has_coerce_map_from(element.parent()._mpoly_base_ring(self.variable_names())):
#                return self(element._mpoly_dict_recursive(self.variable_names(), base_ring))
#
#        if isinstance(element, (SingularElement, cypari2.gen.Gen)):
#            element = str(element)
#
#        if isinstance(element, MPolynomial_libsingular) and element.parent() is not self and element.parent() != self:
#            variable_names_s = element.parent().variable_names()
#            variable_names_t = self.variable_names()
#
#            if set(variable_names_s).issubset(variable_names_t):
#                for v in variable_names_s:
#                    ind_map.append(variable_names_t.index(v)+1)
#            else:
#                ind_map = [i+1 for i in range(_ring.N)]
#
#            if element.parent().ngens() <= self.ngens():
#                # Map the variables by indices
#                _p = p_ISet(0, _ring)
#                Element = <MPolynomial_libsingular>element
#                El_poly = Element._poly
#                El_parent = Element._parent
#                El_ring = Element._parent_ring
#                El_base = El_parent._base
#
#                #this loop needs improvement
#                while El_poly:
#                    c = si2sa(p_GetCoeff(El_poly, El_ring), El_ring, El_base)
#                    if check:
#                        try:
#                            c = base_ring(c)
#                        except TypeError:
#                            p_Delete(&_p, _ring)
#                            raise
#                    if c:
#                        mon = p_Init(_ring)
#                        p_SetCoeff(mon, sa2si(c, _ring), _ring)
#                        for j from 1 <= j <= El_ring.N:
#                            e = p_GetExp(El_poly, j, El_ring)
#                            if e:
#                                p_SetExp(mon, ind_map[j-1], e, _ring)
#                        p_Setm(mon, _ring)
#                        _p = p_Add_q(_p, mon, _ring)
#                    El_poly = pNext(El_poly)
#                return new_MP(self, _p)
#
#        if isinstance(element, MPolynomial_polydict):
#            variable_names_s = element.parent().variable_names()
#            variable_names_t = self.variable_names()
#
#            if set(variable_names_s).issubset(variable_names_t):
#                for v in variable_names_s:
#                    ind_map.append(variable_names_t.index(v)+1)
#            else:
#                ind_map = [i+1 for i in range(_ring.N)]
#
#            if element.parent().ngens() <= self.ngens():
#                bucket = sBucketCreate(_ring)
#                try:
#                    for (m,c) in element.element().dict().iteritems():
#                        if check:
#                            c = base_ring(c)
#                        if not c:
#                            continue
#                        mon = p_Init(_ring)
#                        p_SetCoeff(mon, sa2si(c , _ring), _ring)
#                        for pos in m.nonzero_positions():
#                            overflow_check(m[pos], _ring)
#                            p_SetExp(mon, ind_map[pos], m[pos], _ring)
#                        p_Setm(mon, _ring)
#                        sBucket_Merge_m(bucket, mon)
#                    e=0
#                    sBucketClearMerge(bucket, &_p, &e)
#                    sBucketDestroy(&bucket)
#                except TypeError:
#                    sBucketDeleteAndDestroy(&bucket)
#                    raise
#                return new_MP(self, _p)
#
#        from sage.rings.polynomial.pbori import BooleanPolynomial
#        if isinstance(element, BooleanPolynomial):
#            if element.constant():
#                if element:
#                    return self._one_element
#                else:
#                    return self._zero_element
#
#            variable_names_s = set(element.parent().variable_names())
#            variable_names_t = self.variable_names()
#
#            if variable_names_s.issubset(variable_names_t):
#                return eval(str(element),self.gens_dict(copy=False))
#
#            elif element.parent().ngens() <= self.ngens():
#                Q = element.parent()
#                gens_map = dict(zip(Q.variable_names(),self.gens()[:Q.ngens()]))
#                return eval(str(element),gens_map)

        if isinstance(element, str):
            # let python do the parsing
            d = self.gens_dict()
            if self.base_ring().gen() != 1:
                d[str(self.base_ring().gen())]=self.base_ring().gen()
            try:
                if '/' in element:
                    element = sage_eval(element,d)
                else:
                    element = element.replace("^","**")
                    element = eval(element, d, {})
            except (SyntaxError, NameError):
                raise TypeError("Could not find a mapping of the passed element to this ring.")

            # we need to do this, to make sure that we actually get an
            # element in self.
            return self._coerce_c(element)

#        if isinstance(element, dict):
#            if len(element)==0:
#                _p = p_ISet(0, _ring)
#            else:
#                bucket = sBucketCreate(_ring)
#                try:
#                    for (m,c) in element.iteritems():
#                        if check:
#                            c = base_ring(c)
#                        if not c:
#                            continue
#                        mon = p_Init(_ring)
#                        p_SetCoeff(mon, sa2si(c , _ring), _ring)
#                        if len(m) != self.ngens():
#                            raise TypeError("tuple key must have same length as ngens")
#                        for pos from 0 <= pos < len(m):
#                            if m[pos]:
#                                overflow_check(m[pos], _ring)
#                                p_SetExp(mon, pos+1, m[pos], _ring)
#                        p_Setm(mon, _ring)
#                        sBucket_Merge_m(bucket, mon)
#                    e=0
#                    sBucketClearMerge(bucket, &_p, &e)
#                    sBucketDestroy(&bucket)
#                except TypeError:
#                    sBucketDeleteAndDestroy(&bucket)
#                    raise
#            return new_MP(self, _p)

        try: #if hasattr(element,'_polynomial_'):
            # SymbolicVariable
            return element._polynomial_(self)
        except AttributeError:
            pass

#        if is_Macaulay2Element(element):
#            return self(element.external_string())
#        try:
#            return self(str(element))
#        except TypeError:
#            pass

        try:
            # now try calling the base ring's __call__ methods
            element = self.base_ring()(element)
            elementint = Integer(element)

            p = MPolynomial_flint.__new__(MPolynomial_flint)
            fmpz_mpoly_init(p._poly, self._ctx)
            p._parent = self

            fmpz_mpoly_set_ui(p._poly, mpz_get_ui(elementint.value), self._ctx)

            return p
        except (TypeError, ValueError):
            raise TypeError("Could not find a mapping of the passed element to this ring.")

    def gen(self, int n=0):
        """
        Returns the ``n``-th generator of this multivariate polynomial
        ring.

        INPUT:

        - ``n`` -- an integer ``>= 0``

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: P.gen(),P.gen(1)
            (x, y)

            sage: P = PolynomialRing(GF(127),1000,'x')
            sage: P.gen(500)
            x500

            sage: P.<SAGE,SINGULAR> = QQ[] # weird names
            sage: P.gen(1)
            SINGULAR
        """

        if n < 0 or n >= self.__ngens:
            raise ValueError("Generator not defined.")

        cdef MPolynomial_flint p = MPolynomial_flint.__new__(MPolynomial_flint)
        p._parent = self

        fmpz_mpoly_gen(p._poly, n, self._ctx)

        return p

    def addmul_multi(self, terms, verbose=False):
        """
        Accepts a list of terms, each of which is a list of factors, and constructs the sum
        of the products.  The factors can be either polynomials or rational functions in
        the polynomial ring's fraction field.
        """

        cdef MPolynomial_flint p = MPolynomial_flint.__new__(MPolynomial_flint)
        p._parent = self

        if verbose: print("fmpz_mpoly_addmul_multi entry", len(terms), file=sys.stderr)

        # Construct a list of polynomials that we will refer to by their index number.
        # Factoring each polynomial seems a bit too time consuming, so I just ensure
        # that they're relatively prime.

        # Each list element is either a polynomial or a list of indices indicating that
        # that polynomial has been factored.

        polys = list()

        # Add a polynomial to the polys list (if needed), ensuring that polys is maintained
        # relatively prime and return a list of indices into polys that when multiplied
        # together form the new polynomial.  Returns a list of indices whose polynomials,
        # when multiplied together, form the input polynomial.

        def add_to_polys(newpoly, startpoly=0):
            assert newpoly.parent() is self
            if verbose: print("add_to_polys: adding polynomial of length", len(newpoly))
            result = []
            factorization = newpoly.factor()
            for poly, power in factorization:
                if poly not in polys:
                    polys.append(poly)
                result += [polys.index(poly)] * power
            return (result, Integer(factorization.unit()))

        # Lists of Counters, one pair for each term, mapping polynomial indices
        # to counts of polynomials in each numerator and denominator
        n = list()
        d = list()

        # Lists of integer coefficients, one pair for each term
        nc = list()
        dc = list()

        # build a list of lists of indices into polys for numerators and denominators
        next_len = 10
        for term in terms:
            n.append(Counter())
            d.append(Counter())
            nc.append(Integer(1))
            dc.append(Integer(1))
            for factor in term:
               if is_FractionFieldElement(factor) or factor.parent() is self.base_ring().fraction_field():
                   assert(factor.numerator().parent() is self or factor.numerator().parent() is self.base_ring())
                   assert(factor.denominator().parent() is self or factor.denominator().parent() is self.base_ring())
                   l,u = add_to_polys(self(factor.numerator()))
                   n[-1].update(l)
                   nc[-1] *= u
                   if verbose: print("added polynomial to numerator", len(n), n[-1], nc[-1])
                   l,u = add_to_polys(self(factor.denominator()))
                   d[-1].update(l)
                   dc[-1] *= u
                   if verbose: print("added polynomial to denominator", len(d), d[-1], dc[-1])
               else:
                   if not (factor.parent() is self or factor.parent() is self.base_ring()):
                       print("factor", factor, factor.parent())
                   assert(factor.parent() is self or factor.parent() is self.base_ring())
                   l,u = add_to_polys(self(factor))
                   n[-1].update(l)
                   nc[-1] *= u
                   if verbose: print("added polynomial to numerator", len(n), n[-1], nc[-1])
               if verbose and len(polys) >= next_len:
                   next_len += 10
                   print("len(polys) =", len(polys))

        def flatten(t):
            return [item for sublist in t for item in sublist]

        # some of the earlier elements in n and d might have "split" into multiple polynomials.  Fix them.

        def expand_splits(p):
            if type(polys[p]) == list:
                return flatten(map(expand_splits, polys[p]))
            else:
                return [p]

        for i in range(len(n)):
            n[i] = Counter(flatten(map(expand_splits, n[i].elements())))
            d[i] = Counter(flatten(map(expand_splits, d[i].elements())))

        # construct the LCM of all the denominators, in the form of a list of indices

        def common_elements(lists):
            common_set = set(itertools.chain.from_iterable(lists))
            return flatten([[e]*max(map(lambda l: l[e], lists)) for e in common_set])

        #print("fmpz_mpoly_addmul_multi term_lengths =", term_lengths, file=sys.stderr)

        lcm = common_elements(d)
        lcmc = reduce(lambda a,b: a.lcm(b), dc)

        if verbose: print("fmpz_mpoly_addmul_multi lcm =", lcm, "lcmc =", lcmc, file=sys.stderr)

        #print("fmpz_mpoly_addmul_multi n =", n, file=sys.stderr)
        #print("fmpz_mpoly_addmul_multi d =", d, file=sys.stderr)

        if verbose: print("fmpz_mpoly_addmul_multi array building", file=sys.stderr)

        cdef slong * iptr = <slong *>malloc(sizeof(slong) * len(terms))

        num_polys = 0
        for t in range(len(terms)):
           term_len = sum(n[t].values())
           lcm2 = Counter(lcm)
           lcm2.subtract(d[t])
           term_len += sum(lcm2.values())
           if nc[t] * lcmc != dc[t]:
               term_len += 1
           iptr[t] = term_len
           num_polys += term_len

        if verbose: print("fmpz_mpoly_addmul_multi num_polys =", num_polys, file=sys.stderr)

        cdef const fmpz_mpoly_struct ** fptr = <const fmpz_mpoly_struct **>malloc(sizeof(fmpz_mpoly_struct *) * num_polys)
        cdef MPolynomial_flint ni
        constants = []
        k = 0
        try:
            for t in range(len(terms)):
                for i in n[t].elements():
                    ni = polys[i]
                    fptr[k] = <const fmpz_mpoly_struct *>ni._poly
                    k += 1
                lcm2 = Counter(lcm)
                lcm2.subtract(d[t])
                for i in lcm2.elements():
                    ni = polys[i]
                    fptr[k] = <const fmpz_mpoly_struct *>ni._poly
                    k += 1
                if nc[t] * lcmc != dc[t]:
                    constant = Integer(nc[t] * lcmc / dc[t])
                    #if verbose: print("constant", constant)
                    if constant not in constants:
                        constants.append(self(constant))
                    ni = constants[constants.index(constant)]
                    #if verbose: print("constant ni", ni)
                    fptr[k] = <const fmpz_mpoly_struct *>ni._poly
                    k += 1
        except Exception as ex:
            print(ex, file=sys.stdout)
            raise

        if verbose: print("fmpz_mpoly_addmul_multi", len(terms), file=sys.stderr)
        fmpz_mpoly_addmul_multi(p._poly, fptr, iptr, len(terms), self._ctx)

        if len(lcm) == 0:
            return p
        else:
            denom = prod(polys[i] for i in lcm)
            return p/denom

    # It is required in cython to redefine __hash__ when __richcmp__ is
    # overloaded. Also just writing
    #         __hash__ = CategoryObject.__hash__
    # doesn't work.
    def __hash__(self):
        """
        Return a hash for this ring, that is, a hash of the string
        representation of this polynomial ring.

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: hash(P)      # somewhat random output
            967902441410893180 # 64-bit
            -1767675994        # 32-bit
        """
        return CategoryObject.__hash__(self)

    def __reduce__(self):
        """
        Serializes self.

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQ, order='degrevlex')
            sage: P == loads(dumps(P))
            True

            sage: P.<x,y,z> = PolynomialRing(ZZ, order='degrevlex')
            sage: P == loads(dumps(P))
            True

            sage: P = PolynomialRing(GF(127), names='abc')
            sage: P == loads(dumps(P))
            True

            sage: P = PolynomialRing(GF(2^8,'F'), names='abc')
            sage: P == loads(dumps(P))
            True

            sage: P = PolynomialRing(GF(2^16,'B'), names='abc')
            sage: P == loads(dumps(P))
            True
            sage: z = QQ['z'].0
            sage: P = PolynomialRing(NumberField(z^2 + 3,'B'), names='abc')
            sage: P == loads(dumps(P))
            True
        """
        return unpickle_MPolynomialRing_flint, \
            (self.base_ring(), self.variable_names(), self.term_order())

def unpickle_MPolynomialRing_flint(base_ring, names, term_order):
    """
    inverse function for ``MPolynomialRing_flint.__reduce__``

    EXAMPLES::

        sage: P.<x,y> = PolynomialRing(QQ)
        sage: loads(dumps(P)) is P # indirect doctest
        True
    """
    from sage.rings.polynomial.polynomial_ring_constructor import _multi_variate
    # If libsingular would be replaced by a different implementation in future
    # sage version, the unpickled ring will belong the new implementation.
    return _multi_variate(base_ring, tuple(names), None, term_order, "FLINT")


cdef class MPolynomial_flint(MPolynomial):
    """
    A multivariate polynomial implemented using FLINT.
    """

    def __init__(self, MPolynomialRing_flint parent):
        """
        Construct a zero element in parent.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomial_libsingular
            sage: P = PolynomialRing(GF(32003),3,'x')
            sage: MPolynomial_libsingular(P)
            0
        """
        fmpz_mpoly_init(self._poly, parent._ctx)
        
        self._parent = parent
        self._ctx = & parent._ctx[0]

    def __dealloc__(self):
        # WARNING: the Cython class self._parent is now no longer accessible!
        fmpz_mpoly_clear(self._poly, self._ctx)

    cpdef _richcmp_(left, right, int op):
        """
        Compare left and right.

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQ,order='degrevlex')
            sage: x == x
            True

            sage: x > y
            True
            sage: y^2 > x
            True

            sage: (2/3*x^2 + 1/2*y + 3) > (2/3*x^2 + 1/4*y + 10)
            True

        TESTS::

            sage: P.<x,y,z> = PolynomialRing(QQ, order='degrevlex')
            sage: x > P(0)
            True

            sage: P(0) == P(0)
            True

            sage: P(0) < P(1)
            True

            sage: x > P(1)
            True

            sage: 1/2*x < 3/4*x
            True

            sage: (x+1) > x
            True

            sage: f = 3/4*x^2*y + 1/2*x + 2/7
            sage: f > f
            False
            sage: f < f
            False
            sage: f == f
            True

            sage: P.<x,y,z> = PolynomialRing(GF(127), order='degrevlex')
            sage: (66*x^2 + 23) > (66*x^2 + 2)
            True
        """
        assert (<MPolynomial_flint>left)._parent == (<MPolynomial_flint>right)._parent

        equal = fmpz_mpoly_equal((<MPolynomial_flint>left)._poly, (<MPolynomial_flint>right)._poly, (<MPolynomialRing_flint>left._parent)._ctx)

        if equal:
            return rich_to_bool(op, 0)
        else:
            # XXX wrong
            return rich_to_bool(op, 1)

    cpdef _add_(left, right):
        """
        Add left and right.

        EXAMPLES::

            sage: P.<x,y,z>=PolynomialRing(QQ,3)
            sage: 3/2*x + 1/2*y + 1 #indirect doctest
            3/2*x + 1/2*y + 1
        """

        assert (<MPolynomial_flint>left)._parent == (<MPolynomial_flint>right)._parent

        cdef MPolynomial_flint p = MPolynomial_flint.__new__(MPolynomial_flint)
        p._parent = (<MPolynomial_flint>left)._parent

        fmpz_mpoly_add(p._poly, (<MPolynomial_flint>left)._poly, (<MPolynomial_flint>right)._poly, (<MPolynomialRing_flint>left._parent)._ctx)

        return p

    cpdef _sub_(left, right):
        """
        Subtract left and right.

        EXAMPLES::

            sage: P.<x,y,z>=PolynomialRing(QQ,3)
            sage: 3/2*x - 1/2*y - 1 #indirect doctest
            3/2*x - 1/2*y - 1
        """

        assert (<MPolynomial_flint>left)._parent == (<MPolynomial_flint>right)._parent

        cdef MPolynomial_flint p = MPolynomial_flint.__new__(MPolynomial_flint)
        p._parent = (<MPolynomial_flint>left)._parent

        fmpz_mpoly_sub(p._poly, (<MPolynomial_flint>left)._poly, (<MPolynomial_flint>right)._poly, (<MPolynomialRing_flint>left._parent)._ctx)

        return p

    cpdef _lmul_(self, Element left):
        """
        Multiply self with a base ring element.

        EXAMPLES::

            sage: P.<x,y,z>=PolynomialRing(QQ,3)
            sage: 3/2*x # indirect doctest
            3/2*x

        ::

            sage: P.<x,y,z>=PolynomialRing(QQ,3)
            sage: (3/2*x - 1/2*y - 1) * (3/2) # indirect doctest
            9/4*x - 3/4*y - 3/2
        """

        cdef MPolynomial_flint p = MPolynomial_flint.__new__(MPolynomial_flint)
        p._parent = (<MPolynomial_flint>self)._parent

        leftint = Integer(left)

        fmpz_mpoly_scalar_mul_ui(p._poly, (<MPolynomial_flint>self)._poly, mpz_get_ui(leftint.value), (<MPolynomialRing_flint>self._parent)._ctx)

        return p


    cpdef _mul_(left, right):
        """
        Multiply left and right.

        EXAMPLES::

            sage: P.<x,y,z>=PolynomialRing(QQ,3)
            sage: (3/2*x - 1/2*y - 1) * (3/2*x + 1/2*y + 1) # indirect doctest
            9/4*x^2 - 1/4*y^2 - y - 1

            sage: P.<x,y> = PolynomialRing(QQ,order='lex')
            sage: (x^2^15) * x^2^15
            Traceback (most recent call last):
            ...
            OverflowError: exponent overflow (...)
        """
        # all currently implemented rings are commutative
        assert (<MPolynomial_flint>left)._parent == (<MPolynomial_flint>right)._parent

        cdef MPolynomial_flint p = MPolynomial_flint.__new__(MPolynomial_flint)
        p._parent = (<MPolynomial_flint>left)._parent

        fmpz_mpoly_mul(p._poly, (<MPolynomial_flint>left)._poly, (<MPolynomial_flint>right)._poly, (<MPolynomialRing_flint>left._parent)._ctx)

        return p

    cpdef _div_(left, right_ringelement):
        """
        Divide left by right

        EXAMPLES::

            sage: R.<x,y>=PolynomialRing(QQ,2)
            sage: f = (x + y)/3 # indirect doctest
            sage: f.parent()
            Multivariate Polynomial Ring in x, y over Rational Field

        Note that / is still a constructor for elements of the
        fraction field in all cases as long as both arguments have the
        same parent and right is not constant. ::

            sage: R.<x,y>=PolynomialRing(QQ,2)
            sage: f = x^3 + y
            sage: g = x
            sage: h = f/g; h
            (x^3 + y)/x
            sage: h.parent()
            Fraction Field of Multivariate Polynomial Ring in x, y over Rational Field

        If we divide over `\ZZ` the result is the same as multiplying
        by 1/3 (i.e. base extension). ::

            sage: R.<x,y> = ZZ[]
            sage: f = (x + y)/3
            sage: f.parent()
            Multivariate Polynomial Ring in x, y over Rational Field
            sage: f = (x + y) * 1/3
            sage: f.parent()
            Multivariate Polynomial Ring in x, y over Rational Field

        But we get a true fraction field if the denominator is not in
        the fraction field of the base ring.""

            sage: f = x/y
            sage: f.parent()
            Fraction Field of Multivariate Polynomial Ring in x, y over Integer Ring

        Division will fail for non-integral domains::

            sage: P.<x,y> = Zmod(1024)[]
            sage: x/P(3)
            Traceback (most recent call last):
            ...
            TypeError: self must be an integral domain.

            sage: x/3
            Traceback (most recent call last):
            ...
            TypeError: self must be an integral domain.

        TESTS::

            sage: R.<x,y>=PolynomialRing(QQ,2)
            sage: x/0
            Traceback (most recent call last):
            ...
            ZeroDivisionError: rational division by zero
        """
        assert (<MPolynomial_flint>left)._parent == (<MPolynomial_flint>right_ringelement)._parent

        return (left._parent).fraction_field()(left,right_ringelement)

    def _repr_(self):
        """
        EXAMPLES::

            sage: R.<x,y>=PolynomialRing(QQ,2)
            sage: f = x^3 + y
            sage: f # indirect doctest
            x^3 + y
        """
        cdef const char * cstr
        cstr = fmpz_mpoly_get_str_pretty(self._poly, (<MPolynomialRing_flint>self._parent)._cnames, (<MPolynomialRing_flint>self._parent)._ctx)
        pstr = char_to_str(cstr)
        free(<void *>cstr)
        return pstr

    def __len__(self):
        cdef fmpz_mpoly_struct A = self._poly[0]
        return Integer(A.length)

    def degree(self, x=None):
        """
        Return the degree of this polynomial.

        INPUT:

        - ``x`` -- (default: ``None``) a generator of the parent ring

        OUTPUT:

        If ``x`` is not given, return the maximum degree of the monomials of
        the polynomial. Note that the degree of a monomial is affected by the
        gradings given to the generators of the parent ring. If ``x`` is given,
        it is (or coercible to) a generator of the parent ring and the output
        is the maximum degree in ``x``. This is not affected by the gradings of
        the generators.

        EXAMPLES::

            sage: R.<x, y> = QQ[]
            sage: f = y^2 - x^9 - x
            sage: f.degree(x)
            9
            sage: f.degree(y)
            2
            sage: (y^10*x - 7*x^2*y^5 + 5*x^3).degree(x)
            3
            sage: (y^10*x - 7*x^2*y^5 + 5*x^3).degree(y)
            10

        The term ordering of the parent ring determines the grading of the
        generators. ::

            sage: T = TermOrder('wdegrevlex', (1,2,3,4))
            sage: R = PolynomialRing(QQ, 'x', 12, order=T+T+T)
            sage: [x.degree() for x in R.gens()]
            [1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4]

        A matrix term ordering determines the grading of the generators by the
        first row of the matrix. ::

            sage: m = matrix(3, [3,2,1,1,1,0,1,0,0])
            sage: m
            [3 2 1]
            [1 1 0]
            [1 0 0]
            sage: R.<x,y,z> = PolynomialRing(QQ, order=TermOrder(m))
            sage: x.degree(), y.degree(), z.degree()
            (3, 2, 1)
            sage: f = x^3*y + x*z^4
            sage: f.degree()
            11

        If the first row contains zero, the grading becomes the standard one. ::

            sage: m = matrix(3, [3,0,1,1,1,0,1,0,0])
            sage: m
            [3 0 1]
            [1 1 0]
            [1 0 0]
            sage: R.<x,y,z> = PolynomialRing(QQ, order=TermOrder(m))
            sage: x.degree(), y.degree(), z.degree()
            (1, 1, 1)
            sage: f = x^3*y + x*z^4
            sage: f.degree()
            5

        To get the degree with the standard grading regardless of the term
        ordering of the parent ring, use ``std_grading=True``. ::

            sage: f.degree(std_grading=True)
            5

        TESTS::

            sage: P.<x, y> = QQ[]
            sage: P(0).degree(x)
            -1
            sage: P(1).degree(x)
            0

        The following example is inspired by :trac:`11652`::

            sage: R.<p,q,t> = ZZ[]
            sage: poly = p + q^2 + t^3
            sage: poly = poly.polynomial(t)[0]
            sage: poly
            q^2 + p

        There is no canonical coercion from ``R`` to the parent of ``poly``, so
        this doesn't work::

            sage: poly.degree(q)
            Traceback (most recent call last):
            ...
            TypeError: argument is not coercible to the parent

        Using a non-canonical coercion does work, but we require this
        to be done explicitly, since it can lead to confusing results
        if done automatically::

            sage: poly.degree(poly.parent()(q))
            2
            sage: poly.degree(poly.parent()(p))
            1
            sage: T.<x,y> = ZZ[]
            sage: poly.degree(poly.parent()(x))  # noncanonical coercions can be confusing
            1

        The argument to degree has to be a generator::

            sage: pp = poly.parent().gen(0)
            sage: poly.degree(pp)
            1
            sage: poly.degree(pp+1)
            Traceback (most recent call last):
            ...
            TypeError: argument is not a generator

        Canonical coercions are used::

            sage: S = ZZ['p,q']
            sage: poly.degree(S.0)
            1
            sage: poly.degree(S.1)
            2
        """
        cdef slong deg
        if x == None:
            deg = fmpz_mpoly_total_degree_si(self._poly, (<MPolynomialRing_flint>self._parent)._ctx)
            return Integer(deg)
        else:
            for i,gen in enumerate(self.parent().gens()):
                if gen == x:
                    deg = fmpz_mpoly_degree_si(self._poly, i, (<MPolynomialRing_flint>self._parent)._ctx)
                    return Integer(deg)
            raise TypeError('argument is not a generator')

    def dict(self):
        """
        Return a dictionary representing self. This dictionary is in
        the same format as the generic MPolynomial: The dictionary
        consists of ``ETuple:coefficient`` pairs.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: f=2*x*y^3*z^2 + 1/7*x^2 + 2/3
            sage: f.dict()
            {(0, 0, 0): 2/3, (1, 3, 2): 2, (2, 0, 0): 1/7}
        """
        cdef fmpz_mpoly_struct A = self._poly[0]

        n = (<MPolynomialRing_flint>self._parent).ngens()
        cdef ulong *exp = <ulong *>malloc(sizeof(ulong) * n)

        result = {}

        for i in range(A.length):
            fmpz_mpoly_get_term_exp_ui(exp, self._poly, i, (<MPolynomialRing_flint>self._parent)._ctx)
            explist = []
            for j in range(n):
                explist.append(Integer(exp[j]))
            result[ETuple(explist)] = Integer(fmpz_mpoly_get_coeff_ui_ui(self._poly, exp, (<MPolynomialRing_flint>self._parent)._ctx))

        free(exp)

        return result

    def exponents(self, as_ETuples=True):
        """
        Return the exponents of the monomials appearing in this
        polynomial.

        INPUT:

        - ``as_ETuples`` -- (default: ``True``) if ``True`` returns the
          result as an list of ETuples, otherwise returns a list of tuples

        EXAMPLES::

            sage: R.<a,b,c> = QQ[]
            sage: f = a^3 + b + 2*b^2
            sage: f.exponents()
            [(3, 0, 0), (0, 2, 0), (0, 1, 0)]
            sage: f.exponents(as_ETuples=False)
            [(3, 0, 0), (0, 2, 0), (0, 1, 0)]
        """
        cdef fmpz_mpoly_struct A = self._poly[0]

        n = (<MPolynomialRing_flint>self._parent).ngens()
        cdef ulong *exp = <ulong *>malloc(sizeof(ulong) * n)

        result = []

        for i in range(A.length):
            fmpz_mpoly_get_term_exp_ui(exp, self._poly, i, (<MPolynomialRing_flint>self._parent)._ctx)
            explist = []
            for j in range(n):
                explist.append(Integer(exp[j]))
            if as_ETuples:
                result.append(ETuple(explist))
            else:
                result.append(tuple(explist))

        free(exp)

        return result

    def inverse_of_unit(self):
        """
        Return the inverse of this polynomial if it is a unit.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: x.inverse_of_unit()
            Traceback (most recent call last):
            ...
            ArithmeticError: Element is not a unit.

            sage: R(1/2).inverse_of_unit()
            2
        """


        if fmpz_mpoly_equal_si((<MPolynomial_flint>self)._poly, 1, (<MPolynomialRing_flint>self._parent)._ctx) or \
           fmpz_mpoly_equal_si((<MPolynomial_flint>self)._poly, -1, (<MPolynomialRing_flint>self._parent)._ctx):
            return self
        else:
            raise ArithmeticError("Element is not a unit.")

    def constant_coefficient(self):
        """
        Return the constant coefficient of this multivariate
        polynomial.

        EXAMPLES::

            sage: P.<x, y> = QQ[]
            sage: f = 3*x^2 - 2*y + 7*x^2*y^2 + 5
            sage: f.constant_coefficient()
            5
            sage: f = 3*x^2
            sage: f.constant_coefficient()
            0
        """
        n = (<MPolynomialRing_flint>self._parent).ngens()
        cdef ulong *exp = <ulong *>malloc(sizeof(ulong) * n)

        for i in range(n):
            exp[i] = 0

        cdef ulong coeff = fmpz_mpoly_get_coeff_ui_ui(self._poly, exp, (<MPolynomialRing_flint>self._parent)._ctx)

        free(exp)

        return Integer(coeff)

    cpdef _floordiv_(self, right):
        """
        Perform division with remainder and return the quotient.

        INPUT:

        - ``right`` - something coercible to an MPolynomial_flint
          in ``self.parent()``

        EXAMPLES::

            sage: R.<x,y,z> = GF(32003)[]
            sage: f = y*x^2 + x + 1
            sage: f//x
            x*y + 1
            sage: f//y
            x^2

            sage: P.<x,y> = ZZ[]
            sage: x//y
            0
            sage: (x+y)//y
            1

            sage: P.<x,y> = QQ[]
            sage: (x+y)//y
            1
            sage: (x)//y
            0

            sage: P.<x,y> = Zmod(1024)[]
            sage: (x+y)//x
            1
            sage: (x+y)//(2*x)
            Traceback (most recent call last):
            ...
            NotImplementedError: Division of multivariate polynomials over non fields by non-monomials not implemented.

        TESTS::

            sage: P.<x,y> = ZZ[]
            sage: p = 3*(-x^8*y^2 - x*y^9 + 6*x^8*y + 17*x^2*y^6 - x^3*y^2)
            sage: q = 7*(x^2 + x*y + y^2 + 1)
            sage: p*q//q == p
            True
            sage: p*q//p == q
            True
        """

        assert (<MPolynomial_flint>self)._parent == (<MPolynomial_flint>right)._parent

        cdef MPolynomial_flint p = MPolynomial_flint.__new__(MPolynomial_flint)
        p._parent = (<MPolynomial_flint>self)._parent

        fmpz_mpoly_div(p._poly, (<MPolynomial_flint>self)._poly, (<MPolynomial_flint>right)._poly, (<MPolynomialRing_flint>self._parent)._ctx)

        return p

#        if _right._poly == NULL:
#            raise ZeroDivisionError
#        elif p_IsOne(_right._poly, r):
#            return self

#        if r.cf.type != n_unknown:
#            if (singular_polynomial_length_bounded(_right._poly, 2) == 1
#                    and r.cf.cfIsOne(p_GetCoeff(_right._poly, r), r.cf)):
#                p = self._poly
#                quo = p_ISet(0,r)
#                while p:
#                    if p_DivisibleBy(_right._poly, p, r):
#                        temp = p_MDivide(p, _right._poly, r)
#                        p_SetCoeff0(temp, n_Copy(p_GetCoeff(p, r), r), r)
#                        quo = p_Add_q(quo, temp, r)
#                    p = pNext(p)
#                return new_MP(parent, quo)
#            if r.cf.type == n_Znm or r.cf.type == n_Zn or r.cf.type == n_Z2m :
#                raise NotImplementedError("Division of multivariate polynomials over non fields by non-monomials not implemented.")

#        count = singular_polynomial_length_bounded(self._poly, 15)

#        # fast in the most common case where the division is exact; returns zero otherwise
#        if count >= 15:  # note that _right._poly must be of shorter length than self._poly for us to care about this call
#            sig_on()
#        quo = p_Divide(p_Copy(self._poly, r), p_Copy(_right._poly, r), r)
#        if count >= 15:
#            sig_off()

#        if quo == NULL:
#            if r.cf.type == n_Z:
#                P = parent.change_ring(QQ)
#                f = (<MPolynomial_libsingular>P(self))._floordiv_(P(right))
#                return parent(sum([c.floor() * m for c, m in f]))
#            else:
#                sig_on()
#                quo = singclap_pdivide(self._poly, _right._poly, r)
#                sig_off()

#        return new_MP(parent, quo)


    @coerce_binop
    def gcd(left, right, algorithm=None, **kwds):
        """
        Return the greatest common divisor of left and right.

        INPUT:

        - ``right`` - polynomial
        - ``algorithm``
          - ``ezgcd`` - EZGCD algorithm
          - ``modular`` - multi-modular algorithm (default)
        - ``**kwds`` - ignored

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: f = (x*y*z)^6 - 1
            sage: g = (x*y*z)^4 - 1
            sage: f.gcd(g)
            x^2*y^2*z^2 - 1
            sage: GCD([x^3 - 3*x + 2, x^4 - 1, x^6 -1])
            x - 1

            sage: R.<x,y> = QQ[]
            sage: f = (x^3 + 2*y^2*x)^2
            sage: g = x^2*y^2
            sage: f.gcd(g)
            x^2

        We compute a gcd over a finite field::

            sage: F.<u> = GF(31^2)
            sage: R.<x,y,z> = F[]
            sage: p = x^3 + (1+u)*y^3 + z^3
            sage: q = p^3 * (x - y + z*u)
            sage: gcd(p,q)
            x^3 + (u + 1)*y^3 + z^3
            sage: gcd(p,q)  # yes, twice -- tests that singular ring is properly set.
            x^3 + (u + 1)*y^3 + z^3

        We compute a gcd over a number field::

            sage: x = polygen(QQ)
            sage: F.<u> = NumberField(x^3 - 2)
            sage: R.<x,y,z> = F[]
            sage: p = x^3 + (1+u)*y^3 + z^3
            sage: q = p^3 * (x - y + z*u)
            sage: gcd(p,q)
            x^3 + (u + 1)*y^3 + z^3

        TESTS::

            sage: Q.<x,y,z> = QQ[]
            sage: P.<x,y,z> = QQ[]
            sage: P(0).gcd(Q(0))
            0
            sage: x.gcd(1)
            1

            sage: k.<a> = GF(9)
            sage: R.<x,y> = PolynomialRing(k)
            sage: f = R.change_ring(GF(3)).gen()
            sage: g = x+y
            sage: g.gcd(f)
            1
            sage: x.gcd(R.change_ring(GF(3)).gen())
            x

            sage: Pol.<x,y,z> = ZZ[]
            sage: p = x*y - 5*y^2 + x*z - z^2 + z
            sage: q = -3*x^2*y^7*z + 2*x*y^6*z^3 + 2*x^2*y^3*z^4 + x^2*y^5 - 7*x*y^5*z
            sage: (21^3*p^2*q).gcd(35^2*p*q^2) == -49*p*q
            True
        """
#        cdef poly *_res
#        cdef ring *_ring = self._parent_ring
#        cdef MPolynomial_libsingular _right = <MPolynomial_libsingular>right
#
#        if _right._poly == NULL:
#            return self
#        elif self._poly == NULL:
#            return right
#        elif p_IsOne(self._poly, _ring):
#            return self
#        elif p_IsOne(_right._poly, _ring):
#            return right
#
#        res = new_MP(self._parent, _res)
#        return res

        assert (<MPolynomial_flint>left)._parent == (<MPolynomial_flint>right)._parent

        cdef MPolynomial_flint p = MPolynomial_flint.__new__(MPolynomial_flint)
        p._parent = (<MPolynomial_flint>left)._parent

        if not fmpz_mpoly_gcd(p._poly, (<MPolynomial_flint>left)._poly, (<MPolynomial_flint>right)._poly, (<MPolynomialRing_flint>left._parent)._ctx):
            raise RuntimeError("GCD failed")

        return p

    def factor(self, proof=None):
        cdef fmpz_mpoly_factor_t f
        fmpz_mpoly_factor_init(f, (<MPolynomialRing_flint>self._parent)._ctx)
        if not fmpz_mpoly_factor(f, self._poly, (<MPolynomialRing_flint>self._parent)._ctx):
            raise RuntimeError("multivariate factorization failed")
        factors = []
        cdef MPolynomial_flint p
        for i in range(f.num):
            p = MPolynomial_flint.__new__(MPolynomial_flint)
            p._parent = self._parent
            fmpz_mpoly_init(p._poly, (<MPolynomialRing_flint>self._parent)._ctx)
            # would like to fmpz_mpoly_swap, but this is an inline function in FLINT
            # and I'm not sure how to import that into Cython
            fmpz_mpoly_set(p._poly, &f.poly[i], (<MPolynomialRing_flint>self._parent)._ctx)
            factors.append((p, f.exp[i]))
        return Factorization(factors, unit=fmpz_get_ui(f.constant))

    @coerce_binop
    def lcm(left, right, algorithm=None, **kwds):
        """
        Return the least common multiple of left and right.
        """
#        cdef poly *_res
#        cdef ring *_ring = self._parent_ring
#        cdef MPolynomial_libsingular _right = <MPolynomial_libsingular>right
#
#        if _right._poly == NULL:
#            return self
#        elif self._poly == NULL:
#            return right
#        elif p_IsOne(self._poly, _ring):
#            return self
#        elif p_IsOne(_right._poly, _ring):
#            return right
#
#        res = new_MP(self._parent, _res)
#        return res

        return left*right/left.gcd(right)

    def __reduce__(self):
        """
        Serialize this polynomial.

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQ,3, order='degrevlex')
            sage: f = 27/113 * x^2 + y*z + 1/2
            sage: f == loads(dumps(f))
            True

            sage: P = PolynomialRing(GF(127),3,names='abc')
            sage: a,b,c = P.gens()
            sage: f = 57 * a^2*b + 43 * c + 1
            sage: f == loads(dumps(f))
            True

        TESTS:

        Verify that :trac:`9220` is fixed.

            sage: R=QQ['x']
            sage: S=QQ['x','y']
            sage: h=S.0^2
            sage: parent(h(R.0,0))
            Univariate Polynomial Ring in x over Rational Field
        """
        return unpickle_MPolynomial_flint, (self._parent, self.dict())

def unpickle_MPolynomial_flint(MPolynomialRing_flint R, d):
    """
    Deserialize an ``MPolynomial_flint`` object

    INPUT:

    - ``R`` - the base ring
    - ``d`` - a Python dictionary as returned by :meth:`MPolynomial_flint.dict`

    EXAMPLES::

        sage: P.<x,y> = PolynomialRing(QQ)
        sage: loads(dumps(x)) == x # indirect doctest
        True
    """

    cdef MPolynomial_flint p = MPolynomial_flint.__new__(MPolynomial_flint)
    p._parent = R
    cdef int N = R.ngens()
    cdef ulong *exp = <ulong *>malloc(sizeof(ulong) * N)
    cdef int _i, _e

    try:
        for mon,c in d.iteritems():
            for i in range(N):
                exp[i] = 0
            for i,e in mon.sparse_iter():
                _i = i
                if _i >= N:
                    raise TypeError("variable index too big")
                _e = e
                if _e <= 0:
                    raise TypeError("exponent too small")
                exp[_i] = _e

            fmpz_mpoly_push_term_si_ui(p._poly, mpz_get_si(Integer(c).value), exp, R._ctx)

    except:
        free(exp)
        raise

    free(exp)

    fmpz_mpoly_sort_terms(p._poly, R._ctx)
    fmpz_mpoly_combine_like_terms(p._poly, R._ctx)
    fmpz_mpoly_assert_canonical(p._poly, R._ctx)

    return p
