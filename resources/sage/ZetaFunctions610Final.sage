
################################################################################
#                                                                              #
# CALCULATION OF THE (LOCAL) IGUSA AND TOPOLOGICAL ZETA FUNCTIONS OF A         #
# NON-DEGENERATED POLYNOMIAL WITH RESPECT TO HIS NEWTON'S POLYHEDRON.          #
# For Sage v6.10.                                                              #
#                                                                              #
# Last update: 18-01-2016                                                      #
#                                                                              #
# These functions are based on the work of K. Hoornaert and D. Loots:          #
# "Computer program written  in Maple for the calculation of Igusa local zeta  #
# function".                                                                   #
# http://www.wis.kuleuven.ac.be/algebra/kathleen.htm, 2000.                    #
#                                                                              #
# For any bug or commentary, please contact me.                                #
#                                                                              #
# Juan Viu-Sos                                                                 #
# Universite de Pau et des Pays de l'Adour                                     #
# http://jviusos.perso.univ-pau.fr                                             #
#                                                                              #
#                                                                              #
# This program is free software; you can redistribute it and/or modify         #
# it under the terms of the GNU General Public License as published by         #
# the Free Software Foundation; either version 2 of the License, or            #
# (at your option) any later version.                                          #
#                                                                              #
# This program is distributed in the hope that it will be useful,              #
# but WITHOUT ANY WARRANTY; without even the implied warranty of               #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                #
# GNU General Public License for more details.                                 #
#                                                                              #
# You should have received a copy of the GNU General Public License            #
# along with this program; if not, write to the Free Software                  #
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,                   #
# MA 02110-1301, USA.                                                          #
#                                                                              #
################################################################################



class ZetaFunctions(object):
    """
    Class ``Zetafunctions`` takes a multivariate polynomial as argument for calculate
    their associated (local) Igusa and Topological zeta functions. This class allows us to get 
    information about: the associated Newton's polyhedron, their faces, the associated cones,...\n
    This class is composed by a multivariate polynomial `f` of degree `n` with non-constant 
    term and his associated Newton's polyhedron `\Gamma(f)`.\n
    Methods in ZetaFunctions:\n
    - give_info_facets(self)
    - give_info_newton(self, faces = False, cones = False)
    - newton_plot(self)
    - cones_plot(self)
    - give_expected_pole_info(self,d = 1, local = False, weights = None)
    - igusa_zeta(self, p = None, dict_Ntau = {}, local = False, weights = None, info = False, \ check = 'ideals')
    - topological_zeta(self, d = 1, local = False, weights = None, info = False, check = 'ideals')
    - monodromy_zeta(self, weights = None, char = False, info = False)

    .. WARNING::
        These formulas for the Igusa and Topological zeta functions only work when the \
        given polynomial is NOT DEGENERATED with respect to his associated Newton \
        Polyhedron (see [DenHoo]_, [DenLoe]_ and [Var]_).
    
    EXAMPLES::

        sage: R.<x,y,z> = QQ[]
        sage: zex = ZetaFunctions(x^2 + y*z)
        sage: zex.give_info_newton()
        Newton's polyhedron of x^2 + y*z:
            support points = [(2, 0, 0), (0, 1, 1)]
            vertices = [(0, 1, 1), (2, 0, 0)]
            number of proper faces = 13
            Facet 1: x >= 0
            Facet 2: y >= 0
            Facet 3: z >= 0
            Facet 4: x + 2*z - 2 >= 0
            Facet 5: x + 2*y - 2 >= 0
        sage: zex.topological_zeta()
        (s + 3)/((s + 1)*(2*s + 3))
        sage: zex.give_expected_pole_info()
        The candidate poles of the (local) topological zeta function (with d =
        1) of x^2 + y*z in function of s are:
        
        -3/2 with expected order: 2
        The responsible face of maximal dimension is ``tau_0`` = minimal face
        who intersects with the diagonal of ambient space:
             tau6: dim 1,  vertices = [(0, 1, 1), (2, 0, 0)],  rays = []
            generators of cone = [(1, 0, 2), (1, 2, 0)], partition into simplicial
            cones = [[(1, 0, 2), (1, 2, 0)]]
            
        -1 with expected order: 1
        (If all Vol(tau) are 0, where tau runs through the selected faces that
        are no vertices, then the expected order of -1 is 0).


    REFERENCES:

    .. [DenHoo] J . Denef and K . Hoornaert, "Newton Polyhedra and Igusa's Local Zeta Function.", 2001.

    .. [DenLoe] J . Denef and F . Loeser, "Caracteristiques d'Euler-Poincare;, fonctions zeta locales et modifications analytiques.", 1992.

    .. [HooLoo] K . Hoornaert and D . Loots, "Computer program written in Maple for the calculation of Igusa's local zeta function.", http://www.wis.kuleuven.ac.be/algebra/kathleen.htm, 2000.

    .. [Var] A . N . Varchenko, "Zeta-function of monodromy and Newton's diagram.", 1976.

    .. [Viu] J . Viu-Sos, "Funciones zeta y poliedros de Newton: Aspectos teoricos y computacionales.", 2012.


    AUTHORS:

    - Kathleen Hoornaert (2000): initial version for Maple

    - Juan Viu-Sos (2012): initial version for Sage
    """
    def __init__(self, poly):
        #Polynomial
        self._f = poly
        #Newton's polyhedron
        self._Gammaf = newton_polyhedron(poly)

    def give_info_facets(self):
        """
        Prints a relation of facets in Newton's polyhedron and the inequalities which define them.
        """
        give_all_facets_info(self._f,self._Gammaf)

    def give_info_newton(self, faces = False, cones = False):
        """
        Prints information about the the Newton's polyhedron of ``f``:\n
            - Support points of f.
            - Vertices of Newton's polyhedron.
            - Number of proper faces.
            - Inequalities defining facets.
        Options:
            - ``faces = True`` prints information about each face in polyhedron.
            - ``cones = True`` prints information about each cone associated to faces in \
                polyhedron.
        """
        print "Newton's polyhedron of " + str(self._f) + ":"
        print "\t" + "support points = " + str(support_points(self._f))
        print "\t" + "vertices = " + str(map(tuple,self._Gammaf.vertices()))
        print "\t" + "number of proper faces = " + str(len(proper_faces(self._Gammaf)))
        give_all_facets_info(self._f,self._Gammaf)
        if faces or cones:
            print "Information about faces:"
            faces_set = proper_faces(self._Gammaf)
            for tau in faces_set:
                face_info, cone_info = str(), str()
                i = faces_set.index(tau)
                if faces: face_info = face_info_output(tau) + "\n"
                if cones: cone_info = cone_info_output(cone_from_face(tau)) + "\n"
                print "tau" + str(i) + ": " + face_info + cone_info

    def newton_plot(self):
        """
        Plots the associated Newton's polyhedron (for n = 2 , 3).
        """
        show(self._Gammaf.plot())

    def cones_plot(self):
        """
        Plots the Fan of cones associated to Newton's polyhedron (for n = 2 , 3).
        """
        F = fan_all_cones(self._Gammaf)
        show(F.plot())

    def give_expected_pole_info(self,d = 1, local = False, weights = None):
        """
        Prints information about the candidate real poles of the topological zeta function 
        `Z_{top, f}^{(d)}(s)` like order of poles and responsible faces of highest dimension.\n
            - ``local = True`` calculates the local (at the origin) topological Zeta function.
            - ``weights`` -- a list `[k_1,\ldots,k_n]` of weights for the volume form.

        REFERENCES:

        .. [DenHoo] J . Denef and K . Hoornaert, "Newton Polyhedra and Igusa's Local Zeta Function.", 2001.
        """
        give_expected_pole_info(self._f,d, local, weights)

    def igusa_zeta(self, p = None, dict_Ntau = {}, local = False, weights = None, info = False, 
                   check = 'ideals'):
        """
        Returns the expression of the Igusa zeta function for ``p`` a prime number (explicit or \
        abstract), in terms of a symbolic variable ``s``.\n
        - ``local = True`` calculates the local Igusa zeta function (at the origin).\n
        - ``weights`` -- a list `[k_1,\ldots,k_n]` of weights for the volume form.\n
        - ``info = True`` gives information of each face `\\tau`, the associated cone of \
            `\\tau`, and the values ``L_tau`` and ``S_tau`` in the process.\n
        - ``check`` -- choose the method to check the non-degeneracy condition ('default' or \
            'ideals'). If ``check = 'no_check'``, degeneracy checking is omitted.

        .. WARNING::
            This formula is only valid when the the given polynomial is NOT DEGENERATED for \
            ``p`` with respect to his associated Newton Polyhedron (see [DenHoo]_).

        In the abstract case ``p = None``, you can give a dictionary ``dict_Ntau`` where:\n
        - The keys are the polynomials `f_{\\tau}` associated of each face `\\tau` of the \
          Newton Polyhedron.\n
        - The items are the associated abstract values \
          `N_{\\tau}= \#\{a\in(\mathbb{F}_p - 0)^d \mid f^*_{\\tau}(a)=0\}` with \
          `f^*_{\\tau}=\mathbb{F}_p(f_{\\tau})`, depending of a symbolic variable ``p``.
        If the value associated to a face `\\tau_k` is not in the dictionary, function \
        introduces a new symbolic variable ``N_tauk`` to represent `N_{\\tau_k}`.\n

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]; p = var('p')
            sage: zex1 = ZetaFunctions(x^2 - y^2 + z^3)
            sage: #For p=3 given
            sage: zex1.igusa_zeta(p = 3)
            2*(3^(2*s + 4) - 3^(s + 1) + 2)*3^(2*s)/((3^(s + 1) - 1)*(3^(3*s + 4) - 1))
            sage: #For p arbitrary, we can give the number of solutions over the faces
            sage: dNtau1 = { x^2-y^2+z^3 : (p-1)*(p-3), -y^2+z^3 : (p-1)^2, x^2+z^3 : (p-1)^2, x^2-y^2 : 2*(p-1)^2 }
            sage: zex1.igusa_zeta(p = None, dict_Ntau = dNtau1)
            (p - 1)*(p + p^(2*s + 4) - p^(s + 1) - 1)*p^(2*s)/((p^(s + 1) - 1)*(p^(3*s + 4) - 1))
            sage: #
            sage: zex2 = ZetaFunctions(x^2 + y*z + z^2)
            sage: #For p=3 mod 4, we can give the number of solutions over the faces
            sage: dNtau2 = { x^2+y*z+z^2 : (p-1)^2, y*z+z^2 : (p-1)^2,  x^2+y*z : (p-1)^2, x^2+z^2 : 0 }
            sage: zex2.igusa_zeta(p = None, dict_Ntau = dNtau2)
            (p^(s + 3) - 1)*(p - 1)*p^(2*s)/((p^(s + 1) - 1)*(p^(2*s + 3) - 1))
            sage: #For p=1 mod 4
            sage: dNtau2bis = { x^2+y*z+z^2 : (p-1)*(p-3), y*z+z^2 : (p-1)^2,  x^2+y*z : (p-1)^2, x^2+z^2 : 2*(p-1)^2 }
            sage: zex2.igusa_zeta(p = None, dict_Ntau = dNtau2bis)
            (p^(s + 3) - 1)*(p - 1)*p^(2*s)/((p^(s + 1) - 1)*(p^(2*s + 3) - 1))

        REFERENCES:

        .. [DenHoo] J . Denef and K . Hoornaert, "Newton Polyhedra and Igusa's Local Zeta Function.", 2001.
        """
        return igusa_zeta(self._f, p, dict_Ntau, local, weights, info, check)

    def topological_zeta(self, d = 1, local = False, weights = None, info = False, 
                         check = 'ideals'):
        """
        Returns the expression of the Topological zeta function `Z_{top, f}^{(d)}` for `d\geq 1`, 
        in terms of the symbolic variable ``s``:\n
            - ``local = True`` calculates the local Topological zeta function (at the origin).\n
            - ``weights`` -- a list `[k_1,\ldots,k_n]` of weights for the volume form.\n
            - ``d`` -- (default:1) an integer. We consider only the divisor whose multiplicity \
            is a multiple of ``d`` (see [DenLoe]_). 
            - ``info = True`` gives information of each face `\\tau`, the associated cone of \
            `\\tau`, and the values `J_\\tau` and `\mathop{dim}(\\tau)!\cdot\mathop{Vol}(\\tau)` \
            in the process (see [DenLoe]_).
            - ``check`` -- choose the method to check the non-degeneracy condition ('default' or \
            'ideals'). If ``check = 'no_check'``, degeneracy checking is omitted.

        .. WARNING::
            This formula is only valid when the the given polynomial is NOT DEGENERATED \
            with respect to his associated Newton Polyhedron (see [DenLoe]_).
            
        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: zex1 = ZetaFunctions(x^2 + y*z)
            sage: zex1.topological_zeta()
            (s + 3)/((s + 1)*(2*s + 3))
            sage: #For d = 2
            sage: zex1.topological_zeta(d = 2)
            1/(2*s + 3)

        REFERENCES:

        .. [DenLoe] J . Denef and F . Loeser, "Caracteristiques d'Euler-Poincare;, fonctions zeta locales et modifications analytiques.", 1992.
        """
        return topological_zeta(self._f, d, local, weights, info, check)

    def monodromy_zeta(self, weights = None, char = False, info = False, check = 'ideals'):
        """
        Returns the expression of the Monodromy zeta function at the origin, in terms of the 
        symbolic variable ``s``.\n
            - ``weights`` -- a list `[k_1,\ldots,k_n]` of weights for the volume form.\n
            - ``char = True`` prints the characteristic polynomial of the monodromy (only if \
            ``f`` has an isolated singularity at the origin).\n
            - ``info = True`` gives information of each face `\\tau`, the associated cone of \
            `\\tau`, and the values `J_\\tau` and `\mathop{dim}(\\tau)!*\mathop{Vol}(\\tau)` \
            in the process (see [Var]_).
            - ``check`` -- choose the method to check the non-degeneracy condition ('default' or \
            'ideals'). If ``check = 'no_check'``, degeneracy checking is omitted.

        .. WARNING::
            This formula is only valid when the the given polynomial is NOT DEGENERATED \
            with respect to his associated Newton Polyhedron (see [Var]_).

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: zex1 = ZetaFunctions(y^7+x^2*y^5+x^5*y^3)
            sage: zex1.monodromy_zeta(char = True)
            The characteristic polynomial of the monodromy is (T - 1)^3*(T^6 + T^5 + T^4 + T^3 +
            T^2 + T + 1)*(T^18 + T^17 + T^16 + T^15 + T^14 + T^13 + T^12 + T^11 + T^10 + T^9 +
            T^8 + T^7 + T^6 + T^5 + T^4 + T^3 + T^2 + T + 1)

            1/((t^7 - 1)*(t^19 - 1))
            sage: S.<x,y,z> = QQ[]
            sage: zex2 = ZetaFunctions(x*y + z^3)
            sage: zex2.monodromy_zeta(char = True)
            The characteristic polynomial of the monodromy is T^2 + T + 1

            -t^3 + 1

        REFERENCES:

        .. [Var] A . N . Varchenko, "Zeta-function of monodromy and Newton's diagram.", 1976.
        """
        return monodromy_zeta(self._f, weights, char, info, check)




###------------------------AUXILIARY FUNCTIONS------------------------###
##---NEWTON'S POLYHEDRON
def support_points(f):
    """
    Returns a list containing the support of ``f``: points of `\ZZ^n` corresponding to \
    exponents of each monomial into ``f``.
    """
    points = f.exponents()
    return points

def newton_polyhedron(f):
    """
    Constructs an object ``Polyhedra`` who represents the Newton's Polyhedron `\Gamma(f)` of \
    the polynomial ``f``.
    """
    P = Polyhedron(vertices = support_points(f), rays=VectorSpace(QQ,f.parent().ngens()).basis())
    return P


##--- FACES 
def faces(P):
    """
    Returns a ``Poset`` of the faces in the polyhedron ``P`` with a partial relation of order 
    defined by the content relation between faces.
    """
    P_poset = Poset(LatticePoset(P.face_lattice()), facade = False)  
    return P_poset

def proper_faces(P):
    """
    Returns a list with the proper faces of the polyhedron ``P`` sorted in increasing order 
    dimension.
    """
    L = faces(P).list()[1:-1]
    return L

#-- Informations about faces
def face_Vinfo(tau):
    """
    Returns a list containing strings with descriptions of the face ``tau`` in terms of vertices \
    and rays.
    """
    return tau.element.ambient_Vrepresentation()

def face_Hinfo(tau):
    """
    Returns a list containing strings with descriptions of the face ``tau`` in terms of the \
    inequalities of the facets which define ``tau`` by their intersection.
    """
    return tau.element.ambient_Hrepresentation()

def contains_a_ray(tau):
    """
    Checks if the face ``tau`` contains some ray (i.e, if is not bounded).
    """
    bool = False
    Vrep = face_Vinfo(tau)
    for e in Vrep:
        if e.is_ray() == True:
            bool = True
            break
    return bool

def compact_faces(P):
    """
    Returns a list with the compact faces of the polyhedron ``P`` sorted in increasing dimension \
    order.
    """
    pfaces = proper_faces(P)
    return filter(lambda i: not contains_a_ray(i), pfaces)

def vertices(tau):
    """
    Returns a list with the vertices of the face ``tau``.
    """
    L = map(lambda i:i.vector(),filter(lambda j:j.is_vertex(),face_Vinfo(tau)))
    return L

def rays(tau):
    """
    Returns a list with the rays of the face ``tau``.
    """
    L = map(lambda i:i.vector(),filter(lambda j:j.is_ray(),face_Vinfo(tau)))
    return L

def translate_points(points_list):
    """
    Returns a list of re-parametrized points by affine translation consider the first point in \
    points_list as the origin.
    """
    origin = points_list[0]
    L = map(lambda v: v - origin, points_list)
    return L

def dim_face(tau):
    """
    Gives the dimension of the face ``tau``.
    """
    vertices_tau = vertices(tau)
    rays_tau = rays(tau)
    if len(vertices_tau) == 0 and len(vertices_tau) == 0: dim_tau = -1
    else:
        v_list = translate_points(vertices_tau)
        dim_tau = matrix(v_list + rays_tau).rank()
    return dim_tau

def facets(P):
    """
    Returns a list of facets ((n-1)-dimensional faces) of the polyhedron ``P``.
    """
    dim = P.dim() 
    L = filter(lambda i: dim_face(i) == dim-1, proper_faces(P))
    return L

def facet_info(f, facet):
    """
    Returns a string with the inequality which define the facet, writed in form 
    `a_1 x_1 + a_2 x_2 + ... + a_n x_n + b \geq 0`.
    """
    rep = face_Hinfo(facet)[0]
    message = str(vector(rep.A()).dot_product(vector(f.parent().gens())) + rep.b())
    message = message  + " >= 0"
    return message

def give_all_facets_info(f,P):
    """
    Prints a relation of facets in ``P`` and their associated inequalities.
    """
    i = 1
    for facet in facets(P):
        print "\tFacet " + str(i) + ": " + facet_info(f, facet)
        i = i + 1

def face_info_output(tau):
    """
    Returns a string containing a description of vertices and rays in face ``tau``.
    """
    info =  "dim " + str(dim_face(tau)) + ",  vertices = " + str(vertices(tau)) + \
            ",  rays = " + str(rays(tau))
    return info

#-- Relations Polyhedron-points
def point_in_face(point,tau):
    """
    Checks if the tuple ``point`` belongs to the face ``tau``.
    """
    return Polyhedron(vertices = vertices(tau), rays = rays(tau)).contains(point)

def support_points_in_face(f, tau):
    """
    Returns a list of support points of ``f`` contained in the face ``tau``.
    """
    L = filter(lambda i: point_in_face(i,tau),support_points(f))
    return L


##---CONES, FANS AND SIMPLE CONES
def primitivize(v):
    """
    Returns a linearly dependent vector with ``v`` whose components are integers whose greatest \
    common divisor is 1.
    """
    return v/gcd(v)

def primitive_vectors(tau):
    """
    Returns a list of primitive vectors of a face: primitive normal vectors to the hyperplanes \
    defining ``tau``.
    """
    L = map(lambda i:primitivize(i.A()),face_Hinfo(tau))
    return L

def cone_from_face(tau):
    """
    Constructs an object ``Cone`` which represents the dual cone of the face ``tau``.\n
    Note: for the total face, it gives a cone generated by the zero vector.
    """
    gens = primitive_vectors(tau)
    if len(gens) == 0: cone = Cone([vertices(tau)[0].parent()(0)])
    else: cone = Cone(gens)
    return cone

def primitive_vectors_cone(cone):
    """
    Returns a list of primitive rays generators of ``cone``.
    """
    L = map(lambda i:primitivize(i.sparse_vector()),cone.rays())
    return L

def all_cones(P):
    """
    Returns a list with all the cones generated by the faces of ``P``.
    """
    L = map(cone_from_face, faces(P)[1:])
    return L

def fan_all_cones(P):
    """
    Construct an object ``Fan`` representing the fan of all associated cones of ``P``.
    """
    F = Fan(all_cones(P),discard_faces=True)
    return F

def same_facet(lcone_gens, p, bcone_gens):
    """
    Checks if ``lcone_gens`` (a cone represented by their generators) and the fixed point 
    ``p`` belongs to the same facet of ``bcone_gens``.
    """
    bool = False
    for face_cone in Cone(bcone_gens).facets():
        rays_lcone = set(map(tuple,lcone_gens))
        rays_face = set(map(tuple,primitive_vectors_cone(face_cone)))
        if ({tuple(p)}.union(rays_lcone)).issubset(rays_face):
            bool = True
            break
    return bool

def simplicial_partition(cone):
    """
    Returns a list with the subcones which forms the simplicial partition of ``cone``.
    """
    L = [cone]
    if not cone.is_simplicial():
        dict = {}
        F = Fan(L)
        list_subcones = []
        #Ordered list of subcones by increasing dimension
        cone_poset = Poset(F.cone_lattice(), facade = False).level_sets()
        for ls in cone_poset[1:-1]: list_subcones = list_subcones + ls 
        for subcone in list_subcones:
            if subcone.element.is_simplicial(): dict[subcone] = [set(subcone.element.rays())]
            else:
                partition = []
                fixpoint = subcone.element.rays()[0]
                for subsubcone in filter(lambda x: x<subcone, list_subcones):
                    if not same_facet(subsubcone.element.rays(), fixpoint, \
                                        subcone.element.rays()):
                        for part in dict[subsubcone]: partition = partition + \
                                                        [part.union({fixpoint})]
                dict[subcone] = partition
        total_cone = list_subcones[-1]
        L = map(Cone,dict[total_cone])
    return L

def cone_info_output(cone, fan_simplicial = None):
    """
    Returns a string containing information about the generators of the cone and his simplicial \
    partition.\n
    - ``fan_simplicial`` -- a fan with the simplicial partition of ``cone`` (if there is already \
    calculated).
    """
    F = fan_simplicial
    if  F == None: F = simplicial_partition(cone)
    info = "generators of cone = " + str(primitive_vectors_cone(cone)) + ", partition into "\
           "simplicial cones = " + str(map(primitive_vectors_cone,F))
    return info

def integral_vectors(scone):
    """
    Returns a list of integral vectors contained in \
    `\{\sum \lambda_j a_j \mid 0\leq \lambda_j <1\}` where `\{a_j\}` is a basis of the simple \
    cone ``scone``.
    """
    origin = VectorSpace(QQ,scone.lattice_dim()).zero_vector()
    if scone.dim() == 0: integrals = [origin]
    else:
        cone_gens = primitive_vectors_cone(scone)
        ngens = len(cone_gens)
        A = transpose(matrix(ZZ, cone_gens))
        D, U, V = A.smith_form()
        diag = D.diagonal()
        coords = mrange(diag, vector)
        #Aux function to scale the vectors on the list
        def escale(v):
            v = vector(QQ, v)
            for i in range(ngens):
                if diag[i] != 0: v[i] = v[i]/diag[i]
                else: v[i] = 0
            return v
        #Aux function 'floor' for vectors component by component
        def floor(v):
            for i in range(ngens):
                v[i] = v[i] - v[i].floor()
            return v
        #Now, we scale and we return to the canonical basis
        L = map(lambda v: V*v, map(escale, coords))
        #Finally, we find the integral vectors of own region
        integrals = map(lambda v: matrix(QQ,A)*v, map(floor,L))
    return integrals

def multiplicity(scone):
    """
    Returns the multiplicity of a simple cone.
    """
    L = primitive_vectors_cone(scone)
    A = matrix(ZZ, L)
    S = A.smith_form()[0]
    result = 1
    for i in range(len(L)):
        result = result*S[i,i]
    return result

#-- Sigma and m functions defined in the paper
def sigma(v, weights = None):
    """
    Returns the weighted sum of the components of ``v`` by ``weights``.
    """
    if weights == None: result = sum(v)
    else: result = vector(v).dot_product(vector(weights))
    return result

def m(v,P):
    """
    Returns `m(v):=\mathop{min}\{v\cdot x \mid x\in P\}` where ``v`` is a vector and ``P`` is a\
    ``Polyhedra`` in the affine space.
    """
    W = QQ^len(v); vrat=W(v)
    L = map(lambda x: vrat.dot_product(x.vector()), P.vertices())
    return min(L)


##--- MONOMIALS ASSOCIATED TO A FACE
def ftau(f,tau):
    """
    Returns the polynomial ``f_tau`` associated to the face ``tau`` of the Newton's Polyhedron \
    of ``f``.
    """
    from sage.rings.polynomial.polydict import ETuple
    # We take a dictionary between the exponents and the coefficients of monomials.
    D=f.dict()
    vars = f.parent().gens()
    nvars = len(vars)
    g = 0
    for v in support_points_in_face(f,tau):
        mon = 1
        for i in range(nvars): mon = mon*vars[i]^v[i]
        g = g + D[ETuple(v)]*mon
    return g

def solve_in_Fp_x(f,p):
    """
    For ``f`` be an integral polynomial, returns the list `[ \{a\in (\mathbb{F}_p-0)^d \mid 
    f^*(a)=0\}, \{` vars of `f^*\}]` where `f^*` is consider `f` with coefficients in \
    `\mathbb{F}_p`, for a given prime number ``p``.
    """
    g = f.change_ring(GF(p))    #We can lose variables in GF(p)
    vars = g.variables() 
    nvars = g.nvariables()
    h = (GF(p)[vars])(g)
    if len(h.exponents()) == 1: sols = []    #If f_tau is a monomial
    else:
        Fp_x_nvars = list(Tuples(range(1,p), nvars))    #(Fp-0)^nvars
        if h == 0: return [Fp_x_nvars, vars]
        sols = filter(lambda a:h(tuple(a))==0,Fp_x_nvars)
    return [sols,vars]

def is_degenerated(f_tau, p = None, method = 'default'):
    """
    Checks if the polynomial ``f_tau`` is degenerated over `\mathbb{F}_p`, for ``p`` a given 
    prime number ``p`` (see [DenHoo]_).\n
    If ``p = None``, checks degeneration  over `\CC` (which is equivalent to be degenerated \
    over `\mathbb{F}_p` with `p>>0`).\n
    For finite fields (``p`` is a given prime):\n
        - ``method = 'default'`` checks the condition using evaluation over the \
    `(\mathbb{F}_p-0)^n` in the system of equations.\n
        - ``method = 'ideals'`` checks the condition using ideals over the finite field.

    REFERENCES:

    .. [DenHoo] J . Denef and K . Hoornaert, "Newton Polyhedra and Igusa's Local Zeta Function.", 2001.
    """
    bool = False
    vars = f_tau.parent().gens()
    if type(p) != Integer:
        S = QQ[vars]
        I = S(f_tau).jacobian_ideal() + S*(S(f_tau))
        bool = prod(S.gens()) not in I.radical()
    else:
        if method == 'ideals':
            S = GF(p)[vars]
            I = S(f_tau).jacobian_ideal() + S*(f_tau)
            for xi in vars:
                I = I + S*(xi^(p-1)-1)    #xi unity in Fp iff xi^{(p-1)-1}=0
            bool = 1 not in I    #True if I in NOT the ring (ie, sist. has a solution)
        else:
            [candidates, vars] = solve_in_Fp_x(f_tau,p)
            if vars == []: bool = True
            else:
                S = GF(p)[vars]
                g = f_tau.change_ring(GF(p))
                for xi in S.gens():
                    df_tau = S(g).derivative(xi)
                    candidates = filter(lambda a: df_tau(tuple(a)) == 0, candidates)
                if len(candidates) != 0: bool = True
    return bool

def is_newton_degenerated(f,P, p = None, local = False, method = 'default'):
    """
    Checks if the polynomial ``f`` is degenerated over `\mathbb{F}_p` (``p`` prime) with respect \
    the faces of the polyhedron ``P`` (see [DenHoo]_).\n
    If ``p = None``, checks degeneration over `\CC` (which is equivalent to be 
    degenerated over `\mathbb{F}_p` with `p>>0`).\n
    ``local = True`` checks degeneration for local case (only with respect the compact faces).\n
    For finite fields (``p`` is a given prime):\n
        - ``method = 'default'`` checks the condition using evaluation over `(\mathbb{F}_p-0)^n` \
    in the system of equations.\n
        - ``method = 'ideals'`` checks the condition using ideals over the finite field.

    REFERENCES:

    .. [DenHoo] J . Denef and K . Hoornaert, "Newton Polyhedra and Igusa's Local Zeta Function.", 2001.
    """
    bool = False
    if local == True: faces_set = compact_faces(P)
    else: faces_set = faces(P)[1:]
    for tau in faces_set:
        f_tau = ftau(f,tau)
        if is_degenerated(f_tau, p, method) == True:
            bool = True
            print "The formula for Igusa Zeta function is not valid:"
            if type(p) != Integer: print "The polynomial is degenerated at least with respect "\
                                          "to the face tau = {" + face_info_output(tau) + "} "\
                                          "over the complex numbers!"
            else: print "The polynomial is degenerated at least with respect to the face tau = "\
                        "{" + face_info_output(tau) + "} over GF(" + str(p) + ")!"
            break
    return bool


##--- IGUSA ZETA FUNCTION
#-- Values Ntau, Ltau and Stau defined in paper [DenHoo]
def Ntau(f,tau,p):
    """
    Returns the number `N_{\\tau} = \#\{a\in(\mathbb{F}_p - 0)^d \mid f^*_{\\tau}(a)=0\}` with 
    `f^*_{\\tau}=\mathbb{F}_p(f_{\\tau}) for a given face `\\tau` and `p` a given prime number \
    (see [DenHoo]_).

    REFERENCES:

    .. [DenHoo] J . Denef and K . Hoornaert, "Newton Polyhedra and Igusa's Local Zeta Function.", 2001.
    """
    n = f.parent().ngens()
    f_tau = ftau(f,tau)
    if type(p) != Integer:
        print "You must to give a 'Dictionary' with the number of solutions in GF(" + str(p) + \
              ")^" + str(n) + " associated to each face."
    else:
        [sols,vars] = solve_in_Fp_x(f_tau,p)
        nsols = len(sols)*(p-1)^(n - len(vars))
    return nsols

def Ltau(f,tau,p,abs_Ntau,s):
    """
    Returns a list `[L_{\\tau}, N_{\\tau}]`` in terms of a symbolic variable ``s`.\n
    If `p` is an abstract prime number, ``abs_Ntau`` is the given symbolic expression of 
    `N_{\\tau}` in this case (see [DenHoo]_).

    REFERENCES:

    .. [DenHoo] J . Denef and K . Hoornaert, "Newton Polyhedra and Igusa's Local Zeta Function.", 2001.
    """
    n = f.parent().ngens()
    if type(p) != Integer: N_tau = abs_Ntau
    else: N_tau = Ntau(f,tau,p)
    result = p^(-n)*((p-1)^n - p*N_tau*((p^s-1)/(p^(s+1)-1)))
    result = factor(result)
    return [result, N_tau]

def Lgamma(f,p,abs_Ngamma,s):
    """
    Returns the value `N_{\\tau}` for the total polyhedron `\\Gamma` in terms of a symbolic \
    variable ``s`.\n
    ``abs_Ngamma`` is the corresponding ``Ngamma`` value for abstract prime ``p``  \
    (see [DenHoo]_).

    REFERENCES:

    .. [DenHoo] J . Denef and K . Hoornaert, "Newton Polyhedra and Igusa's Local Zeta Function.", 2001.
    """
    n = f.parent().ngens()
    if type(p) != Integer: N_gamma = abs_Ngamma
    else:
        [sols,vars] = solve_in_Fp_x(f,p)
        N_gamma = len(sols)*(p-1)^(n - len(vars))
    result = p^(-n)*((p-1)^n - p*N_gamma*((p^s-1)/(p^(s+1)-1)))
    return result

def Stau(f,P,tau,p, weights,s):
    """
    Returns a list ``[S_tau, cone_info]`` with ``cone_info`` containing a string of information 
    about the cones, simplicial partition, multiplicity and integral points (see [DenHoo]_).\n
    Value ``S_tau`` is expressed in terms of a symbolic variable ``s``.

    REFERENCES:

    .. [DenHoo] J . Denef and K . Hoornaert, "Newton Polyhedra and Igusa's Local Zeta Function.", 2001.
    """
    c = cone_from_face(tau)
    dim_cone = c.dim()
    F = simplicial_partition(c)
    result = 0
    for scone in F:
        num = 0
        den = 1
        for h in integral_vectors(scone): num = num + p^(sigma(h, weights) + m(h,P)*s)
        for a in primitive_vectors_cone(scone): den = den*(p^(sigma(a, weights) + m(a,P)*s) - 1)
        result = result + num/den
#        result = factor(simplify(expand(result + num/den)))
    info = cone_info_output(c,F)+ "\n" + "multiplicities = " + str(map(multiplicity,F)) + \
            ", integral points = " + str(map(integral_vectors,F))
    return [result, info]


def igusa_zeta(f, p = None, dict_Ntau = {}, local = False, weights = None, info = False, 
                check = 'ideals'):
    """
    Returns the expression of the Igusa zeta function for ``p`` a prime number (explicit or \
    abstract), in terms of a symbolic variable ``s``.\n
        - ``local = True`` calculates the local Igusa zeta function (at the origin).\n
        - ``weights`` -- a list `[k_1,\ldots,k_n]` of weights for the volume form.\n
        - ``info = True`` gives information of each face `\\tau`, the associated cone of \
            `\\tau`, and the values ``L_tau`` and ``S_tau`` in the process.\n
        - ``check`` -- choose the method to check the non-degeneracy condition ('default' or \
            'ideals'). If ``check = 'no_check'``, degeneracy checking is omitted.

    .. WARNING::
        This formula is only valid when the the given polynomial is NOT DEGENERATED for \
        ``p`` with respect to his associated Newton Polyhedron (see [DenHoo]_).

    In the abstract case ``p = None``, you can give a dictionary ``dict_Ntau`` where:\n
        - The keys are the polynomials `f_{\\tau}` associated of each face `\\tau` of the \
          Newton Polyhedron.\n
        - The items are the associated abstract values \
          `N_{\\tau}= \#\{a\in(\mathbb{F}_p - 0)^d \mid f^*_{\\tau}(a)=0\}` with \
          `f^*_{\\tau}=\mathbb{F}_p(f_{\\tau})`, depending of a symbolic variable ``p``.
    If the value associated to a face `\\tau_k` is not in the dictionary, function \
    introduces a new symbolic variable ``N_tauk`` to represent `N_{\\tau_k}`.\n

    EXAMPLES::

        sage: R.<x,y,z> = QQ[]; p = var('p')
        sage: f1 = x^2 - y^2 + z^3
        sage: #For p=3 given
        sage: igusa_zeta(f, p = 3)
        2*(3^(2*s + 4) - 3^(s + 1) + 2)*3^(2*s)/((3^(s + 1) - 1)*(3^(3*s + 4) - 1))
        sage: #For p arbitrary, we can give the number of solutions over the faces
        sage: dNtau1 = { x^2-y^2+z^3 : (p-1)*(p-3), -y^2+z^3 : (p-1)^2, x^2+z^3 : (p-1)^2, x^2-y^2 : 2*(p-1)^2 }
        sage: igusa_zeta(f1, p = None, dict_Ntau = dNtau1)
        (p - 1)*(p + p^(2*s + 4) - p^(s + 1) - 1)*p^(2*s)/((p^(s + 1) - 1)*(p^(3*s + 4) - 1))
        sage: 
        sage: f2 = x^2 + y*z + z^2
        sage: #For p=3 mod 4, we can give the number of solutions over the faces
        sage: dNtau2 = { x^2+y*z+z^2 : (p-1)^2, y*z+z^2 : (p-1)^2,  x^2+y*z : (p-1)^2, x^2+z^2 : 0 }
        sage: igusa_zeta(f2, p = None, dict_Ntau = dNtau2)
        (p^(s + 3) - 1)*(p - 1)*p^(2*s)/((p^(s + 1) - 1)*(p^(2*s + 3) - 1))
        sage: #For p=1 mod 4
        sage: dNtau2bis = { x^2+y*z+z^2 : (p-1)*(p-3), y*z+z^2 : (p-1)^2,  x^2+y*z : (p-1)^2, x^2+z^2 : 2*(p-1)^2 }
        sage: igusa_zeta(f2, p = None, dict_Ntau = dNtau2bis)
        (p^(s + 3) - 1)*(p - 1)*p^(2*s)/((p^(s + 1) - 1)*(p^(2*s + 3) - 1))

    REFERENCES:

    .. [DenHoo] J . Denef and K . Hoornaert, "Newton Polyhedra and Igusa's Local Zeta Function.", 2001.
    """
    s = var('s')
    if type(p) != Integer: p = var('p')
    P = newton_polyhedron(f)
    abs_Ngamma = None; abs_Ntau = None
    if check!='no_check':
        if is_newton_degenerated(f,P,p,local, method = check) == True: return NaN
    else: print "Warning: you are not checking the non-degeneracy condition for " + str(f) + \
                ", the obtained expression can be false!.\n"
    if local == True:
        faces_set = compact_faces(P)
        result = 0
    else:
        faces_set = proper_faces(P)
        if type(p) != Integer:
            abs_Ngamma = dict_Ntau.get(f)
            if abs_Ngamma == None: abs_Ngamma = var('N_Gamma')
        result = Lgamma(f,p,abs_Ngamma,s)
        if info == True: print "Gamma: total polyhedron\n" + "L_gamma = " + str(result) + "\n\n"
    for tau in faces_set:
        i = proper_faces(P).index(tau)
        if type(p) != Integer:
            f_tau = ftau(f, tau)
            if len(f_tau.exponents()) == 1: abs_Ntau = 0    #If f_tau is a monomial
            else:
                abs_Ntau = dict_Ntau.get(f_tau)
                if abs_Ntau == None: abs_Ntau = var('N_tau' + str(i))
        [L_tau, N_tau] = Ltau(f,tau,p,abs_Ntau,s)
        [S_tau, cone_info] = Stau(f,P,tau,p,weights,s)
        if info == True:
            print "tau" + str(i) + ": " + face_info_output(tau) + "\n" + cone_info + "\n" +\
                  "N_tau = " + str(N_tau) + ", L_tau = " + str(L_tau) + " , S_tau = " +\
                  str(S_tau) + "\n\n"
        result = result + L_tau*S_tau
    result = factor(result)
    return result


##--- TOPOLOGICAL ZETA FUNCTION
#-- Calculation of the numbers Jtau and Mtau defined in paper [DenLoe]
def Jtau(tau,P,weights,s):
    """
    Returns a list [J_tau, cone_info] with ``cone_info`` containing a string of information 
    about the cones, simplicial partition, multiplicity and integral points. (see [DenLoe]_)\n
    Value J_tau is defined in terms of a symbolic variable ``s``.

    REFERENCES:

    .. [DenLoe] J . Denef and F . Loeser, "Caracteristiques d'Euler-Poincare;, fonctions zeta locales et modifications analytiques.", 1992.
    """
    c = cone_from_face(tau)
    dim_cone = c.dim()
    F = simplicial_partition(c)
    if dim_cone == 0: result = 1
    else:
        result = 0
        for scone in filter(lambda i: i.dim()==dim_cone, F):
            num = multiplicity(scone)
            den = 1
            for a in primitive_vectors_cone(scone): den = den*(m(a,P)*s + sigma(a,weights))
            result = factor(simplify(expand(result + num/den)))
    cone_info = cone_info_output(c,F)+ "\n" + "multiplicities = " + str(map(multiplicity,F)) + \
                ", integral points = " + str(map(integral_vectors,F))
    return [result, cone_info]

def Mtau(tau,P,t):
    """
    Returns the value Mtau for ``tau`` face in ``P`` in terms of a symbolic variable ``s`` \
    (see [DenLoe]_).

    REFERENCES:

    .. [DenLoe] J . Denef and F . Loeser, "Caracteristiques d'Euler-Poincare;, fonctions zeta locales et modifications analytiques.", 1992.
    """
    s = var('s')
    c = cone_from_face(tau)
    dim_cone = c.dim()
    F = simplicial_partition(c)
    result = 1
    for scone in filter(lambda i: i.dim()==dim_cone, F):
        mult = multiplicity(scone)
        den = SR(1)
        for a in primitive_vectors_cone(scone): den = den*(m(a,P)*s + sigma(a,weights = None))
        if den.degree(s) == 1:
            M = factor(simplify(s*den.subs(s=1/s)))
            result = result*(1-t^(M.subs(s=0)/mult))
    return result

def face_volume(f,tau):
    """
    Returns the value `\mathop{Vol}(\\tau)\cdot(\mathop{dim}\\tau)!`, for a given face ``tau`` .\n
    `\mathop{Vol}(\\tau)` is defined as follows:
    Let `\omega_\\tau` be the volume form over `\mathop{Aff}(\\tau)`, the affine space generated \
    by `\\tau` such that the parallelepiped spanned by a lattice basis of \
    `\mathop{Aff}(\\tau)\cap\ZZ^n`` has volume 1. \
    Then `\mathop{Vol}(\\tau)` is the volume of `\\tau` intersection the Global Newton \
    Polyhedron of ``f`` with respect to `\omega_\\tau` (see [DenLoe]_).

    REFERENCES:

    .. [DenLoe] J . Denef and F . Loeser, "Caracteristiques d'Euler-Poincare;, fonctions zeta locales et modifications analytiques.", 1992.
    """
    n = f.parent().ngens()
    dim_tau = dim_face(tau)
    result = 0
    if dim_tau != 0:
        tau_in_global = Polyhedron(vertices = support_points_in_face(f,tau))
        vertices_in_global = map(vector,tau_in_global.vertices())
        trans_vertices = translate_points(vertices_in_global)
        if matrix(ZZ,trans_vertices).rank() == dim_tau:
            V = QQ^n
            basis_aff = V.submodule(trans_vertices).intersection(ZZ^n).basis()
            W = V.submodule_with_basis(basis_aff)
            coords_list = map(W.coordinate_vector, trans_vertices)
            p = PointConfiguration(coords_list)
            result = p.volume()    # Returns dimtau!*n-volumen de tau
    else:
        result = 1
    return result

def face_divisors(d,faces_set,P):
    """
    Returns a list of faces `\\tau` in ``faces_set`` such that ``d`` divides \
    `m(\Delta_\\tau) = \mathop{gcd}\{m(a) \mid a\in\Delta_\\tau\cap\ZZ^n\}` where `\Delta_\\tau` \
    is the associated cone of `\\tau`.
    """
    if d ==1: return faces_set
    L_faces = list()
    dim_total = P.dim()
    for tau in faces_set:
        c = cone_from_face(tau)
        F = simplicial_partition(c)
        L_vectors = list()
        #We need to evaluate m over the basis of the cone and the integral points views above.
        for scone in F:
            L_vectors = L_vectors + integral_vectors(scone) + primitive_vectors_cone(scone)
        l = gcd(map(lambda i: m(i,P), L_vectors))
        if d.divides(l): L_faces.append(tau)
    return L_faces

def is_global_degenerated(f, p = None, method = 'default'):
    """
    Checks if the polynomial ``f`` is degenerated over `\mathbb{F}_p` (``p`` prime) with respect \
    the faces of the of the Global Newton Polyhedron of ``f`` (see [DenLoe]_).\n
    If ``p = None``, checks degeneration over `\CC` (which is equivalent to be 
    degenerated over `\mathbb{F}_p` with `p>>0`).\n
    ``local = True`` checks degeneration for local case (only with respect the compact faces).\n
    For finite fields (``p`` is a given prime):\n
        - ``method = 'default'`` checks the condition using evaluation over `(\mathbb{F}_p-0)^n` \
    in the system of equations.\n
        - ``method = 'ideals'`` checks the condition using ideals over the finite field.

    REFERENCES:

    .. [DenLoe] J . Denef and F . Loeser, "Caracteristiques d'Euler-Poincare;, fonctions zeta locales et modifications analytiques.", 1992.
    """
    Q = f.newton_polytope()    #Global Newton Polyhedron of f
    bool = False
    for tau in faces(Q)[1:]:
        f_tau = ftau(f,tau)
        if is_degenerated(f_tau, p, method) == True:
            bool = True
            print "The formula for Topological Zeta function is not valid:"
            if type(p) != Integer: print "The polynomial is degenerated at least with respect "\
                                         "to the face tau = {" + face_info_output(tau) + \
                                         "} over the complex numbers!"
            else: print "The polynomial is degenerated at least with respect to the face tau "\
                        "= {" + face_info_output(tau) + "} over GF(" + str(p) + ")!"
            break
    return bool

#-- Topological Zeta Function of f,  Z_{top, f}^(d) for d>=1 and Monodromy Zeta:
def topological_zeta(f, d = 1, local = False, weights = None, info = False, check = 'ideals'):
    """
    Returns the expression of the Topological zeta function `Z_{top, f}^{(d)}` for `d\geq 1`, 
    in terms of the symbolic variable ``s``:\n
        - ``local = True`` calculates the local Topological zeta function (at the origin).\n
        - ``weights`` -- a list `[k_1,\ldots,k_n]` of weights for the volume form.\n
        - ``d`` -- (default:1) an integer. We consider only the divisor whose multiplicity \
        is a multiple of ``d`` (see [DenLoe]_). 
        - ``info = True`` gives information of each face `\\tau`, the associated cone of \
        `\\tau`, and the values `J_\\tau` and `\mathop{dim}(\\tau)!\cdot\mathop{Vol}(\\tau)` \
        in the process (see [DenLoe]_).
        - ``check`` -- choose the method to check the non-degeneracy condition ('default' or \
        'ideals'). If ``check = 'no_check'``, degeneracy checking is omitted.

    .. WARNING::
        This formula is only valid when the the given polynomial is NOT DEGENERATED \
        with respect to his associated Newton Polyhedron (see [DenLoe]_).

    REFERENCES:

    .. [DenLoe] J . Denef and F . Loeser, "Caracteristiques d'Euler-Poincare;, fonctions zeta locales et modifications analytiques.", 1992.
    """
    s = var('s')
    P = newton_polyhedron(f)
    if check!='no_check':
        if local == True:
            if is_newton_degenerated(f,P,local = True, method = check) == True: return NaN
        else:
            if is_global_degenerated(f, method = check) == True: return NaN
    else: print "Warning: you are not checking the non-degeneracy condition for " + str(f) + \
                ", the obtained expression can be false!.\n"
    result = 0
    if local == True: faces_set = compact_faces(P)
    else:
        faces_set = proper_faces(P)
        if d == 1:
            total_face = faces(P)[-1]
            dim_gamma = dim_face(total_face)
            vol_gamma = face_volume(f,total_face)
            result = (s/(s+1))*((-1)^dim_gamma)*vol_gamma
            if info == True: print "Gamma: total polyhedron\n" + "J_gamma = 1 , "\
                                   "dim_Gamma!*Vol(Gamma) = " + str(vol_gamma) +"\n\n"
    faces_set = face_divisors(d,faces_set,P)
    for tau in faces_set:
        [J_tau, cone_info] = Jtau(tau,P,weights,s)
        dim_tau = dim_face(tau)
        vol_tau = face_volume(f,tau)
        if info == True:
            i = proper_faces(P).index(tau)
            print "tau" + str(i) + ": " + face_info_output(tau) + "\n" + cone_info + "\n" +\
                  "J_tau = " + str(J_tau) + " , dim_tau!*Vol(tau) = " + str(vol_tau) + "\n\n"
        if d == 1:
            if dim_tau == 0: term = J_tau
            else: term = (s/(s+1))*((-1)^dim_tau)*vol_tau*J_tau
        else:
            term = ((-1)^dim_tau)*vol_tau*J_tau
        result = simplify(expand(result + term))
        if result != 0: result = factor(result)
    return result

def monodromy_zeta(f, weights = None, char = False, info = False, check = 'ideals'):
    """
    Returns the expression of the Monodromy zeta function at the origin, in terms of the 
    symbolic variable ``s``.\n
        - ``weights`` -- a list `[k_1,\ldots,k_n]` of weights for the volume form.\n
        - ``char = True`` prints the characteristic polynomial of the monodromy (only if \
        ``f`` has an isolated singularity at the origin).\n
        - ``info = True`` gives information of each face `\\tau`, the associated cone of \
        `\\tau`, and the values `J_\\tau` and `\mathop{dim}(\\tau)!*\mathop{Vol}(\\tau)` \
        in the process (see [Var]_).
        - ``check`` -- choose the method to check the non-degeneracy condition ('default' or \
        'ideals'). If ``check = 'no_check'``, degeneracy checking is omitted.

    .. WARNING::
        This formula is only valid when the the given polynomial is NOT DEGENERATED \
        with respect to his associated Newton Polyhedron (see [Var]_).

    REFERENCES:

    .. [Var] A . N . Varchenko, "Zeta-function of monodromy and Newton's diagram.", 1976.
    """
    n = f.parent().ngens()
    t = var('t')
    P = newton_polyhedron(f)
    if check!='no_check':
        if is_global_degenerated(f, method = check) == True: return NaN
    else: print "Warning: you are not checking the non-degeneracy condition for " + str(f) + \
                ", the obtained expression can be false!.\n"
    result = 1
    i = 0
    for tau in compact_faces(P):
        zeta_tau = Mtau(tau,P,t)
        dim_tau = dim_face(tau)
        vol_tau = face_volume(f,tau)
        if info == True:
            print "tau" + str(i) + ": " + str(face_Vinfo(tau)) + "\n" + "M_tau = " + \
                  str(zeta_tau) + " , dim_tau!*Vol(tau) = " + str(vol_tau) + "\n\n"
            i = i + 1
        result = result*zeta_tau^((-1)^(dim_tau)*vol_tau)
    if char == True:
        result_aux = copy(result)
        T = var('T')
        mu = (-1)^(n - 1)*(result_aux.numerator().degree(t) - \
            result_aux.denominator().degree(t) - 1)
        factor(result_aux)
        aux = result_aux.subs(t = 1/T)
        char = factor(simplify(expand(T^mu*(T/(T - 1)*aux)^((-1)^(n - 1)))))
        print "The characteristic polynomial of the monodromy is " + str(char) +\
              "\n"
    return result



# -- INFORMATION ABOUT POLES OF THE TOPOLOGICAL ZETA FUNCTION
def dict_info_poles(f,d = 1, weights = None, local = False):
    """
    Returns a dictionary where the keys are the candidate real poles of the chosen zeta function. 
    Items are lists containing, by order:\n
    1.- The list of perpendicular vectors to the facets that are responsible for the candidate
    real pole.\n
    2.- A list of the faces of maximal dimension that are responsible for the expected order.\n
    3.- The expected order.\n
    4.- Boolean: for the candidate pole -1 the factor L_tau or s/(s+1) can contribute to the 
    order. If this is is the case, we increase the expected order due to the S_Delta_tau by 1 
    and this record gets the value ``True``. In all other cases we don't increase this expected 
    order and this record gets the value ``False``.
    """
    P = newton_polyhedron(f)
    if local == True: faces_set = compact_faces(P)
    else: faces_set = proper_faces(P)
    faces_set = face_divisors(d,faces_set,P)
    dict_poles = {}
    all_prim_vect = set([])
    for tau in faces_set: all_prim_vect = all_prim_vect.union(set(map(tuple, \
                                            primitive_vectors(tau))))
    valid_prim_vect = filter(lambda v: m(v,P)!=0, all_prim_vect)
    for v in valid_prim_vect:
        realpole = -sigma(v,weights)/m(v,P)
        #We initialize a list of attributes if the pole is not detected yet
        if dict_poles.get(realpole) == None: dict_poles[realpole] = [set([v]), [], 0, False]
        else: dict_poles[realpole][0].add(v)
    #We calculate the maximal expected of each pole and the faces of higher dimension
    #responsibles of this order
    poles_set = dict_poles.keys()
    #If d=1, we have the face tau_0
    #(tau_0 is the smallest face who contains the intersecion between diagonal and polyhedron)
    if d == 1:
        max_pole = max(poles_set)
        for tau in faces_set:
            gens_cone = map(tuple,primitive_vectors(tau))
            if set(gens_cone) ==  dict_poles[max_pole][0]:
                dict_poles[max_pole][1] = [tau]
                dict_poles[max_pole][2] = rank(matrix(gens_cone))
                dict_poles[max_pole][3] = False
                break
        poles_set.remove(max_pole)
    for pole in poles_set:
        maxorder = 0
        responsible_faces = set([])
        for tau in faces_set:
            prim_vect = primitive_vectors(tau)
            interscone = set(map(tuple, prim_vect)).intersection(dict_poles[pole][0])
            if len(interscone) != 0:
                diminters = matrix(QQ, list(interscone)).rank()
                if diminters > maxorder:
                    maxorder = diminters
                    responsible_faces.add(tau)
                elif diminters == maxorder:
                    responsible_faces.add(tau)
        #We find the maximal elements in set of responsible faces
        max_faces = set(responsible_faces)
        for face in responsible_faces:
            if face in max_faces:
                if len(filter(lambda i: face<i, max_faces)): max_faces.remove(face)
        dict_poles[pole][1] = list(max_faces)
        #Max order of pole is max dim of the associated cones
        dict_poles[pole][2] = maxorder
        #We convert the set of vectors into a list
        dict_poles[pole][0] = map(vector, dict_poles[pole][0])
    #Special pole -1 has sometimes a larger order
    if -1 in dict_poles.keys():
        faces_minus_one = dict_poles[-1][1]
        if max(map(lambda tau: len(support_points_in_face(f,tau)), faces_minus_one)) > 1:
            dict_poles[-1][2] = dict_poles[-1][2] + 1
            dict_poles[-1][3] = True
    return dict_poles


def give_expected_pole_info(f,d = 1, local = False, weights = None):
    """
    Prints information about the candidate real poles of the topological zeta function 
    `Z_{top, f}^{(d)}(s)` like order of poles and responsible faces of highest dimension.\n
        - ``local = True`` calculates the local (at the origin) topological Zeta function.
        - ``weights`` -- a list `[k_1,\ldots,k_n]` of weights for the volume form.

    REFERENCES:

    .. [DenHoo] J . Denef and K . Hoornaert, "Newton Polyhedra and Igusa's Local Zeta Function.", 2001.
    """
    dict_poles = dict_info_poles(f,d, weights, local)
    P = newton_polyhedron(f)
    if local == True: faces_set = compact_faces(P)
    else: faces_set = proper_faces(P)
    faces_set = face_divisors(d,faces_set,P)
    n_supp_by_face = map(lambda tau: len(support_points_in_face(f,tau)), faces_set)
    if dict_poles == {}:
        if ( d == 1 and max(n_supp_by_face) == 1 ) or d!=1:
            print "There will be no poles for the (local) topological zeta function " + \
                    "(with d = " + str(d) + ") of " + str(f) + ".\n"
        else:
            print "The candidate poles of the (local) topological zeta function (with d = " + \
                    str(d) + ") of " + str(f) + " in function of s are:\n"
            print "-1 with expected order: 1"
            print "(If all Vol(tau) are 0, where tau runs through the selected faces " + \
                    "that are no vertices, then the expected order of -1 is 0)\n"
    else:
        poles_set = dict_poles.keys()
        #We reconstruct the list of all faces accessing to element
        some_face = dict_poles[poles_set[0]][1][0]
        list_all_faces = some_face.parent().list()[1:-1]
        if d == 1:
            print "The candidate poles of the (local) topological zeta function (with d = " + \
                    str(d) + ") of " + str(f) + " in function of s are:\n"
            max_pole = max(poles_set)
            print str(max_pole) + " with expected order: " + str(dict_poles[max_pole][2])
            if max_pole == -1:
                if dict_poles[-1][3] == True:
                    print "(If all the Vol(tau) of the faces tau that are no vertices and " + \
                            "contained in Gamma are 0, then the expected order of -1 is " + \
                            str(dict_poles[-1][3]-1) + ")."
            tau_0 = dict_poles[max_pole][1][0]
            print "The responsible face of maximal dimension is ``tau_0`` = minimal face " + \
                    "who intersects with the diagonal of ambient space:"
            i = list_all_faces.index(tau_0)
            print "\t tau" + str(i) + ": " +  face_info_output(tau_0) + "\n\t" + \
                    cone_info_output(cone_from_face(tau_0)) + "\n"
            poles_set.remove(max_pole)
            if -1 not in poles_set:
                print "-1 with expected order: 1"
                print "(If all Vol(tau) are 0, where tau runs through the selected faces that " +\
                        "are no vertices, then the expected order of -1 is 0).\n"
        elif local == True:
            print "The candidate poles of the local topological zeta function (with d = " +\
                    str(d) + ") of " + str(f) + " in function of s are:\n"
            if max(n_supp_by_face) > 1 and -1 not in poles_set:
                print "-1 with expected order: 1"
                print "(If all Vol(tau) are 0, where tau runs through the selected faces that " +\
                        "are no vertices, then the expected order of -1 is 0).\n"
        for pole in poles_set:
            print str(pole) + " with expected order: " + str(dict_poles[pole][2])
            if dict_poles[pole][3] == True:
                print "(If all the Vol(tau) of the faces that are no vertices and contained"
                print "one or more of the faces below are 0, then the expected order of -1 is " +\
                     str(dict_poles[pole][3]-1) + ")."
            print "The responsible face(s) of maximal dimension is/are:"
            for tau in dict_poles[pole][1]:
                i = list_all_faces.index(tau)
                print "\t tau" + str(i) + ": " + face_info_output(tau) + "\n\t" + \
                        cone_info_output(cone_from_face(tau)) + "\n"


