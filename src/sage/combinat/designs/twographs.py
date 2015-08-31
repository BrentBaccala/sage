r"""
Two-graphs

A two-graph on `n` points is a family `T \subset \binom {[n]}{3}`
of `3`-sets, such that any `4`-set `S\subset [n]` of size four
contains an even number of elements of `T`. Any graph `([n],E)`
gives rise to a two-graph
`T(E)=\{t \in \binom {[n]}{3} : \left| \binom {t}{2} \cap E \right|\ odd \}`,
and any two graphs with the same two-graph can be obtained one
from the other by :meth:`Seidel switching <Graph.seidel_switching>`.
This defines an equivalence relation on the graphs on `[n]`,
called Seidel switching equivalence.
Conversely, given a two-graph `T`, one can construct a graph
`\Gamma` in the corresponding Seidel switching class with an
isolated vertex `w`. The graph `\Gamma \setminus w` is called
the :meth:`descendant <TwoGraph.descendant>` of `T` w.r.t. `v`.

`T` is called regular if each two-subset of `[n]` is contained
in the same number alpha of triples of `T`.

This module implements a direct construction of a two-graph from a list of
triples, constrution of descendant graphs, regularity checking, and other
things such as constructing the complement two-graph, cf. [BH12]_.

AUTHORS:

- Dima Pasechnik (Aug 2015)

Index
-----

This module's methods are the following :

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~TwoGraph.is_regular_twograph` | tests if ``self`` is a regular two-graph, i.e. a 2-design
    :meth:`~TwoGraph.complement` | returns the complement of ``self``
    :meth:`~TwoGraph.descendant` | returns the descendant graph at `w`

This module's functions are the following :

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :func:`~is_twograph`         | checks that the incidence system is a two-graph
    :func:`~twograph_descendant`  | returns the descendant graph w.r.t. a given vertex of the two-graph of a given graph

Methods
---------
"""
from sage.combinat.designs.incidence_structures import IncidenceStructure
from itertools import combinations

class TwoGraph(IncidenceStructure):
    r"""
    Two-graphs class.

    A two-graph on `n` points is a 3-uniform hypergraph, i.e.  a family `T
    \subset \binom {[n]}{3}` of `3`-sets, such that any `4`-set `S\subset [n]`
    of size four contains an even number of elements of `T`. For more
    information, see the documentation of the
    :mod:`~sage.combinat.designs.twographs` module.

    """
    def __init__(self, points=None, blocks=None, incidence_matrix=None,
            name=None, check=False, copy=True):
        r"""
        Constructor of the class

        TESTS::

            sage: from sage.combinat.designs.twographs import TwoGraph
            sage: TwoGraph([[1,2]])
            Incidence structure with 2 points and 1 blocks
            sage: TwoGraph([[1,2]], check=True)
            Traceback (most recent call last):
            ...
            AssertionError: the structure is not a 2-graph!
            sage: p=graphs.PetersenGraph().twograph()
            sage: TwoGraph(p, check=True)
            Incidence structure with 10 points and 60 blocks
        """
        IncidenceStructure.__init__(self, points=points, blocks=blocks,
                                    incidence_matrix=incidence_matrix,
                                    name=name, check=False, copy=copy)
        if check:  # it is a very slow, O(|points|^4), test...
           from sage.combinat.designs.twographs import is_twograph
           assert is_twograph(self), "the structure is not a 2-graph!"

    def is_regular_twograph(self, alpha=False):
        r"""
        Tests if the :class:`TwoGraph` is regular, i.e. is a 2-design.

        Namely, each pair of elements of :meth:`ground_set` is contained in
        exactly ``alpha`` triples.

        INPUT:

        - ``alpha`` -- (optional, default is ``False``) return the value of
          ``alpha``, if possible.

        EXAMPLES::

            sage: p=graphs.PetersenGraph().twograph()
            sage: p.is_regular_twograph(alpha=True)
            4
            sage: p.is_regular_twograph()
            True
            sage: p=graphs.PathGraph(5).twograph()
            sage: p.is_regular_twograph(alpha=True)
            False
            sage: p.is_regular_twograph()
            False
        """
        r, (_,_,_,a) = self.is_t_design(t=2, k=3, return_parameters=True)
        if r and alpha:
            return a
        return r

    def descendant(self, v):
        """
        The descendant :class:`graph <sage.graphs.graph.Graph>` at ``v``

        The :mod:`switching class of graphs <sage.combinat.designs.twographs>`
        corresponding to ``self`` contains a graph ``D`` with ``v`` its own connected
        component; removing ``v`` from ``D``, one obtains the descendant graph of
        ``self`` at ``v``, which is constructed by this method.

        INPUT:

        - ``v`` -- an element of :meth:`ground_set`

        EXAMPLES::

            sage: p=graphs.PetersenGraph().twograph().descendant(0)
            sage: p.is_strongly_regular(parameters=True)
            (9, 4, 1, 2)
        """
        from sage.graphs.graph import Graph
        return Graph(map(lambda y: filter(lambda z: z != v, y),
                            filter(lambda x: v in x, self.blocks())))

    def complement(self):
        """
        The two-graph which is the complement of ``self``

        That is, the two-graph constisting exactly of triples not in ``self``.
        Note that this is different from :meth:`complement
        <sage.combinat.designs.incidence_structures.IncidenceStructure.complement>`
        of the :class:`parent class
        <sage.combinat.designs.incidence_structures.IncidenceStructure>`.

        EXAMPLES::

            sage: p=graphs.CompleteGraph(8).line_graph().twograph()
            sage: pc = p.complement(); pc
            Incidence structure with 28 points and 1260 blocks

        TESTS::

            sage: from sage.combinat.designs.twographs import is_twograph
            sage: is_twograph(pc)
            True
        """
        return super(TwoGraph, self).complement(uniform=True)

def taylor_twograph(q):
    r"""
    constructing Taylor's two-graph for U_3(q), q odd

    """
    from sage.rings.arith import is_prime_power
    p, k = is_prime_power(q,get_data=True)
    if k==0 or p==2:
       raise ValueError('q must be a an odd prime power')
    from sage.schemes.projective.projective_space import ProjectiveSpace
    from sage.rings.finite_rings.constructor import FiniteField
    from sage.modules.free_module_element import free_module_element as vector
    from sage.rings.finite_rings.integer_mod import mod
    from __builtin__ import sum
    Fq = FiniteField(q**2, 'a')
    PG = ProjectiveSpace(2, Fq)
    def S(xx,yy):
        x = vector(xx)
        y = vector(yy)
        return sum(map(lambda j: x[j]*y[2-j]**q, xrange(3)))

    V = filter(lambda x: S(x,x)==0, PG)
    def make_tester():
        if mod(q,4)==1:
            return lambda (x,y,z): not (S(x,y)*S(y,z)*S(z,x)).is_square()
        else:
            return lambda (x,y,z): (S(x,y)*S(y,z)*S(z,x)).is_square()
    f = make_tester()
    T = filter(f, combinations(V,3))
    return T

def is_twograph(T):
    r"""
    Checks that the incidence system `T` is a two-graph

    INPUT:

    - ``T`` -- an :class:`incidence structure <sage.combinat.designs.IncidenceStructure>`

    EXAMPLES:

    a two-graph from a graph::

        sage: from sage.combinat.designs.twographs import (is_twograph, TwoGraph)
        sage: p=graphs.PetersenGraph().twograph()
        sage: is_twograph(p)
        True

    a non-regular 2-uniform hypergraph which is a two-graph::

        sage: is_twograph(TwoGraph([[1,2,3],[1,2,4]]))
        True

    TESTS:

    wrong size of blocks::

        sage: is_twograph(designs.projective_plane(3))
        False

    a triple system which is not a two-graph::

        sage: is_twograph(designs.projective_plane(2))
        False
    """
    if not T.is_uniform(3):
        return False

    # A structure for a fast triple existence check
    v_to_blocks = {v:set() for v in range(T.num_points())}
    for B in T._blocks:
        B = frozenset(B)
        for x in B:
            v_to_blocks[x].add(B)

    has_triple = lambda (x,y,z) : bool(v_to_blocks[x]&v_to_blocks[y]&v_to_blocks[z])

    # Check that every quadruple contains an even number of triples
    from __builtin__ import sum
    for quad in combinations(range(T.num_points()),4):
        if sum(map(has_triple,combinations(quad,3))) % 2 == 1:
            return False

    return True

def twograph_descendant(G, v, name=None):
    r"""
    Returns the descendant graph w.r.t. vertex `v` of the two-graph of `G`

    In the :mod:`switching class <sage.combinat.designs.twographs>` of `G`,
    construct a graph `\Delta` with `v` an isolated vertex, and return the subgraph
    `\Delta \setminus v`. It is equivalent to, although much faster than, computing the
    :meth:`TwoGraph.descendant` of :meth:`two-graph of G <sage.graphs.graph.Graph.twograph>`, as the
    intermediate two-graph is not constructed.

    INPUT:

    - ``G`` -- a :class:`graph <sage.graphs.graph.Graph>`

    - ``v`` -- a vertex of ``G``

    - ``name`` -- (optional) ``None`` - no name, otherwise derive from the construction

    EXAMPLES:

    one of s.r.g.'s from the :mod:`database <sage.graphs.strongly_regular_db>`::

        sage: from sage.combinat.designs.twographs import twograph_descendant
        sage: A=graphs.strongly_regular_graph(280,135,70)                    # optional - gap_packages internet
        sage: twograph_descendant(A, 0).is_strongly_regular(parameters=True) # optional - gap_packages internet
        (279, 150, 85, 75)

    TESTS::

        sage: T8 = graphs.CompleteGraph(8).line_graph()
        sage: v = T8.vertices()[0]
        sage: twograph_descendant(T8, v)==T8.twograph().descendant(v)
        True
        sage: twograph_descendant(T8, v).is_strongly_regular(parameters=True)
        (27, 16, 10, 8)
        sage: p = graphs.PetersenGraph()
        sage: twograph_descendant(p,5)
        Graph on 9 vertices
        sage: twograph_descendant(p,5,name=True)
        descendant of Petersen graph at 5: Graph on 9 vertices
    """
    G = G.seidel_switching(G.neighbors(v),inplace=False)
    G.delete_vertex(v)
    if name:
        G.name('descendant of '+G.name()+' at '+str(v))
    else:
        G.name('')
    return G
