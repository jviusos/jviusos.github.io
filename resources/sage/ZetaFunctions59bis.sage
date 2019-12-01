
#################################################################################################
#                                                                                               #
# CALCULATION OF THE (LOCAL) IGUSA AND TOPOLOGICAL ZETA FUNCTIONS OF A POLYNOMIAL WITH          #
# NON-DEGENERATED NEWTON'S POLYHEDRON. For Sage v5.9.                                           #
#                                                                                               #
# These functions are based on the work of K. Hoornaert and D. Loots: "Computer program written #
# in Maple for the calculation of Igusa local zeta function".                                   #
# http://www.wis.kuleuven.ac.be/algebra/kathleen.htm, 2000.                                     #
#                                                                                               #
# For any bug or commentary, please contact me.                                                 #
#                                                                                               #
# Juan Viu-Sos                                                                                  #
# Universite de Pau et des Pays de l'Adour                                                      #
# http://jviusos.perso.univ-pau.fr                                                              #
#                                                                                               #
# Last update: 26-06-2013                                                                       #
#                                                                                               #
#################################################################################################

# (Docstring update in progress)

class ZetaFunctions(object):
	"""
	Class ``Zetafunctions`` gets an integral multivariate polynomial as argument for calculate
	their associated (local) Igusa and Topological zeta functions. This class allows to get 
	information about the associated Newton's polyhedron, their faces, the associated cones,...\n
	This class is composed by a integer multivariate polynomial with non-constant term and his 
	associated Newton's polyhedron.\n
	Methods in ZetaFunctions:\n
	- give_info_facets(self)
	- give_info_newton(self, faces = False, cones = False)
	- newton_plot(self)
	- cones_plot(self)
	- give_expected_pole_info(self,d = 1, local = False, weights = None)
	- igusa_zeta(self, p = None, dict_Ntau = {}, local = False, weights = None, info = False, \ check = 'ideals')
	- topological_zeta(self, d = 1, local = False, weights = None, info = False, check = 'ideals')
	- monodromy_zeta(self, weights = None, char = False, info = False)
	
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
	"""
	def __init__(self, poly):
		#Polynomial
		self._f = poly
		#Newton's polyhedron
		self._Gammaf = newton_polyhedron(poly)												   
			
	def give_info_facets(self):
		"""
		Prints a relation of facets in Newton's polyhedron and their inequalities.
		"""
		give_all_facets_info(self._f,self._Gammaf)

	def give_info_newton(self, faces = False, cones = False):
		"""
		Prints information about the the Newton's polyhedron of ``f``:\n
			- Support points of f.
			- Vertices of Newton's polyhedron.
			- Number of proper faces.
			- Inequalities defining facets.
		``faces = True`` prints information about each face of polyhedron.
		``cones = True`` prints information about each cone associated to faces of polyhedron.\n
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
		Prints information about the candidate real poles for the topological zeta function of 
		``f`` relative to orders and responsible faces of highest dimension.\n
		 - ``local = True`` calculates the local (at the origin) topological Zeta function.
		 - ``weights`` -- a list of weights for the volume form.
		"""
		give_expected_pole_info(self._f,d, local, weights)

	def igusa_zeta(self, p = None, dict_Ntau = {}, local = False, weights = None, info = False, check = 'ideals'):
		"""
		Returns the expression of the Igusa zeta function for ``p`` prime (given or abstract), 
		in terms of symbolic variable ``s``.\n
		For the abstract case (``p = None``), you must to give a dictionary ``dist_Ntau`` where 
		the	polynomials ftau for the faces of the Newton Polyhedron are the keys and the 
		abstract value ``N_tau`` (depending of var ``p``) as associated item. If ftau for face 
		``tauk`` is not in the	dictionary, program introduces a new variable ``N_tauk``.\n
			- ``local = True`` calculates the local Igusa zeta function (at the origin).\n
			- ``weights``-- a list of weights for the volume form.\n
			- ``info = True`` gives information of face ``tau``, cone of ``tau`` (all), \
			``L_tau`` and ``S_tau`` in the process.\n
			- ``check`` -- choose the method to check the non-degeneracy condition ('default' or \
			'ideals'). If ``check = 'no_check'``, degeneracy checking is omitted.
		"""
		return igusa_zeta(self._f, p, dict_Ntau, local, weights, info, check)

	def topological_zeta(self, d = 1, local = False, weights = None, info = False, check = 'ideals'):
		"""
		Returns the expression of the Topological zeta function `Z_{top, f}^{(d)}` for `d\geq 1`, 
		in terms of symbolic variable ``s``:\n
			- ``local = True`` calculates the local Topological zeta function (at the origin).\n
			- ``weights`` -- a list of weights for the volume form.\n
			- ``info = True`` gives information of face ``tau``, cone of ``tau``, ``J_tau`` and
			``dim_tau!*Vol(tau)`` in the process.
			- ``check`` -- choose the method to check the non-degeneracy condition ('default' or \
			'ideals'). If ``check = 'no_check'``, degeneracy checking is omitted.
	"""
		return topological_zeta(self._f, d, local, weights, info, check)

	def monodromy_zeta(self, weights = None, char = False, info = False):
		"""
		The Monodromy Zeta Function in the origin, in terms of variable ``s``.\n
		``weights`` is a list of weights if you want to considerate some weighting.\n
		``char = True`` prints the characteristic polynomial of the monodromy (only if ``f``
		has an isolated singularity in the origin).\n
		``info = True`` gives information of face ``tau``, cone of ``tau`` (all), ``J_tau`` and
		``dim_tau!*Vol(tau)`` in the process.
		"""
		return monodromy_zeta(self._f, weights, char, info)

###------------------------AUXILIARY FUNCTIONS------------------------###
##---NEWTON'S POLYHEDRON
def support_points(f):
	"""
	Support of f: points of ZZ^n corresponding to the exponents of the monomial into ``f``.
	"""
	points = f.exponents()
	return points

def newton_polyhedron(f):
	"""
	Construction of Newton's Polyhedron Gamma(f) for the polynomial ``f``.
	"""
	P = Polyhedron(vertices = support_points(f), rays=VectorSpace(QQ,f.parent().ngens()).basis())
	return P


##--- FACES 
def faces(P):
	"""
	Returns a Poset of the faces in the polyhedron ``P`` with a relation of order 
	(content between the faces).
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
	Returns a list containing the descriptions of the face in terms of vertices and rays.
	"""
	return tau.element.ambient_Vrepresentation()

def face_Hinfo(tau):
	"""
	Returns a list containing the descriptions of the face in terms of the inequalities of the 
	facets who intersects into the face.
	"""
	return tau.element.ambient_Hrepresentation()

def contains_a_ray(tau):
	"""
	Checks if the face contains some ray.
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
	Returns a list with the compact faces of the polyhedron ``P`` sorted in increasing order 
	dimension.
	"""
	pfaces = proper_faces(P)
	return filter(lambda i: not contains_a_ray(i), pfaces)

def vertices(tau):
	"""
	Returns a list with the vertices of the face.
	"""
	L = map(lambda i:i.vector(),filter(lambda j:j.is_vertex(),face_Vinfo(tau)))
	return L

def rays(tau):
	"""
	Returns a list with the rays of the face.
	"""
	L = map(lambda i:i.vector(),filter(lambda j:j.is_ray(),face_Vinfo(tau)))
	return L

def translate_points(points_list):
	"""
	Returns a list of points taking the first point in the original list how the origin and 
	rewriting the other points in terms of new origin.
	"""
	origin = points_list[0]
	L = map(lambda v: v - origin, points_list)
	return L

def dim_face(tau):
	"""
	Gives the dimension of the face.
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
	Returns a string with the inequation of the facet write in form 
	a_1*x_1 + a_2*x_2 + ... + a_n*x_n + b >= 0.
	"""
	rep = face_Hinfo(facet)[0]
	message = str(vector(rep.A()).dot_product(vector(f.parent().gens())) + rep.b())
	message = message  + " >= 0"
	return message

def give_all_facets_info(f,P):
	"""
	Prints a relation of facets in ``P`` and their inequalities.
	"""
	i = 1
	for facet in facets(P):
		print "\tFacet " + str(i) + ": " + facet_info(f, facet)
		i = i + 1

def face_info_output(tau):
	"""
	Returns a string containing a description of vertices and rays in face.
	"""
	info =  "dim " + str(dim_face(tau)) + ",  vertices = " + str(vertices(tau)) + \
			",  rays = " + str(rays(tau))
	return info

#-- Relations Polyhedron-points
def point_in_face(point,tau):
	"""
	Checks if point belongs to the face.
	"""
	return Polyhedron(vertices = vertices(tau), rays = rays(tau)).contains(point)

def support_points_in_face(f, tau):
	"""
	Returns a list of support points of ``f`` contained in the face.
	"""
	L = filter(lambda i: point_in_face(i,tau),support_points(f))
	return L


##---CONES, FANES AND SIMPLE CONES
def prim(v):
	"""
	Returns the primitivitation of an integral vector.
	"""
	return v/gcd(v)

def primitive_vectors(tau):
	"""
	Returns a list of primitive vectors of a face (normal vectors of the hyperplanes who defines 
	the face, components are relatively primes).
	"""
	L = map(lambda i:prim(i.A()),face_Hinfo(tau))
	return L

def cone_from_face(tau):
	"""
	Construction of the dual cone of the face. In particular, for the total face it gives a cone 
	generated by the zero vector.
	"""
	gens = primitive_vectors(tau)
	if len(gens) == 0: cone = Cone([vertices(tau)[0].parent()(0)])
	else: cone = Cone(gens)
	return cone

def primitive_vectors_cone(cone):
	"""
	Returns a list of primitive vectors (rays generators) of cone.
	"""
	L = map(lambda i:prim(i.sparse_vector()),cone.rays())
	return L

def all_cones(P):
	"""
	Returns a list with all the cones generated by the faces of ``P``.
	"""
	L = map(cone_from_face, faces(P)[1:])
	return L

def fan_all_cones(P):
	"""
	Fan of all cones of a Polyhedron.
	"""
	F = Fan(all_cones(P),discard_faces=True)
	return F

def same_facet(lcone_gens, p, bcone_gens):
	"""
	Checks if ``lcone_gens`` (a cone represented by their generators) and the fix point 
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
	Returns a list with the subcones who forms the simplicial partition of ``cone``.
	"""
	L = [cone]
	if not cone.is_simplicial():
		dict = {}
		F = Fan(L)
		list_subcones = []
		#Ordered list of subcones by ascending dimension
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

def cone_info_output(cone, F = None):
	"""
	Returns a string containing information about the generators of the cone and his simplicial 
	partition.
	"""
	if F == None: F = simplicial_partition(cone)
	info = "generators of cone = " + str(primitive_vectors_cone(cone)) + ", partition into "\
		   "simplicial cones = " + str(map(primitive_vectors_cone,F))
	return info

def integral_vectors(scone):
	"""
	Returns a list of integral vectors contained in ``{Sum lambda_j*a_j | 0<= lambda_j <1, a_j`` 
	basis of the simple cone ``}``.
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
	Returns the pondered sum of the components in vector.
	"""
	if weights == None: result = sum(v)
	else: result = vector(v).dot_product(vector(weights))
	return result

def m(v,P):
	"""
	Returns min{v.x | x in polyhedron P}.
	"""		
	L = [vector(v).dot_product(vector(x)) for x in P.vertices()]
	return min(L)


##--- MONOMIALS ASSOCIATED TO A FACE
def ftau(f,tau):
	"""
	Returns the polynomial ``f_tau`` associated to the face ``tau`` of the Newton's Polyhedron.
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
	f^*(a)=0\}, \{` vars of `f^*\}]` with	`f^*` being `f` with coefficients in `\mathbb{F}_p`, 
	for a given prime number ``p``.
	"""
	g = f.change_ring(GF(p))	#We can lose variables in GF(p)
	vars = g.variables()		 
	nvars = g.nvariables()
	h = (GF(p)[vars])(g)
	if len(h.exponents()) == 1: sols = []	#If f_tau is a monomial
	else:
		Fp_x_nvars = list(Tuples(range(1,p), nvars))	#(Fp-0)^nvars
		if h == 0: return [Fp_x_nvars, vars]
		sols = filter(lambda a:h(tuple(a))==0,Fp_x_nvars)
	return [sols,vars]

def is_degenerated(f_tau, p = None, method = 'default'):
	"""
	Checks if the polynomial ``f_tau`` is degenerated over `\mathbb{F}_p`, for ``p`` a given 
	prime number ``p``.\n
	If ``p = None``, checks degeneration  over `\mathbb{C}` (which is equivalent to be 
	degenerated	over `\mathbb{F}_p` with `p>>0`).\n
	For finite fields (``p`` is a given prime):\n
		- ``method = 'default'`` checks the condition using evaluation over the \
	`(\mathbb{F}_p-0)^n` in the system of equations.\n
		- ``method = 'ideals'`` checks the condition using ideals over the finite field.
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
				I = I + S*(xi^(p-1)-1)	#xi unity in Fp iff xi^{(p-1)-1}=0
			bool = 1 not in I	#True if I in NOT the ring (ie, sist. has a solution)
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

def is_all_degenerated(f,P, p = None, local = False, method = 'default'):
	"""
	Checks if own polynomial ``f`` is degenerated over F_p (``p`` prime) with respect the 
	faces of the polyhedron ``P``.\n
	If ``p = None``, checks degeneration over `\mathbb{CC}` (which is equivalent to be 
	degenerated over `F_p` with `p>>0`).\n
	``local = True`` checks degeneration for local case (only with respect the compact faces).\n
	For finite fields (``p`` is a given prime):\n
		- ``method = 'default'`` checks the condition using evaluation over `(F_p-0)^n` in the 
	system of equations.\n
		- ``method = 'ideals'`` checks the condition using ideals over the finite field.
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


##--- IGUSA'S ZETA FUNCTION
#-- Values Ntau, Ltau and Stau defined in the paper
def Ntau(f,tau,p):
	"""
	Returns the number ``Ntau = #{a in (F_p^x)^d | f*_tau(a)=0}`` with ``f*_tau`` being f_tau with 
	coefficients in F_p(f_tau) for tau face.
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
	Returns a list ``[L_tau, N_tau]`` in terms of variable ``s`.\n
	``abs_Ntau`` is the corresponding ``Ntau``'s values for abstract prime ``p``.
	"""
	n = f.parent().ngens()
	if type(p) != Integer: N_tau = abs_Ntau
	else: N_tau = Ntau(f,tau,p)
	result = p^(-n)*((p-1)^n - p*N_tau*((p^s-1)/(p^(s+1)-1)))
	result = factor(result)
	return [result, N_tau]

def Lgamma(f,p,abs_Ngamma,s):
	"""
	Returns the value ``Ntau`` for the total polyhedron in terms of variable ``s`.\n
	``abs_Ngamma`` is the corresponding ``Ngamma`` value for abstract prime ``p``.
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
	about the cones, simplicial partition, multiplicity and integral points.\n
	Value ``S_tau`` is in terms of variable ``s``.
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
#		result = factor(simplify(expand(result + num/den)))		
	info = cone_info_output(c,F)+ "\n" + "multiplicities = " + str(map(multiplicity,F)) + \
			", integral points = " + str(map(integral_vectors,F))
	return [result, info]


def igusa_zeta(f, p = None, dict_Ntau = {}, local = False, weights = None, info = False, check = 'ideals'):
	"""
	Returns the expression of the Igusa zeta function for ``p`` prime (given or abstract), 
	in terms of symbolic variable ``s``.\n
	For the abstract case (``p = None``), you must to give a dictionary ``dist_Ntau`` where the 
	polynomials ftau for the faces of the Newton Polyhedron are the keys and the abstract value 
	``N_tau`` (depending of var ``p``) as associated item. If ftau for face ``tauk`` is not in 
	the	dictionary, program introduces a new variable ``N_tauk``.\n
		- ``local = True`` calculates the local Igusa zeta function (at the origin).\n
		- ``weights``-- a list of weights for the volume form.\n
		- ``info = True`` gives information of face ``tau``, cone of ``tau`` (all), ``L_tau`` \
			and ``S_tau`` in the process.\n
		- ``check`` -- choose the method to check the non-degeneracy condition ('default' or \
			'ideals'). If ``check = 'no_check'``, degeneracy checking is omitted.
	"""
	s = var('s')
	if type(p) != Integer: p = var('p')
	P = newton_polyhedron(f)
	abs_Ngamma = None; abs_Ntau = None
	if check!='no_check':
		if is_all_degenerated(f,P,p,local, method = check) == True: return NaN
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
			if len(f_tau.exponents()) == 1: abs_Ntau = 0	#If f_tau is a monomial
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
#-- Calculation of the numbers Jtau and Mtau
def Jtau(tau,P,weights,s):
	"""
	Returns a list [J_tau, cone_info] with ``cone_info`` containing a string of information 
	about the cones, simplicial partition, multiplicity and integral points.\n
	Value J_tau is in terms of variable ``t``.
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
	Returns the value Mtau for ``tau`` face in ``P`` in terms of variable ``t``.
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
	Returns the value ``Vol(tau)*(dim tau)!``, for a given face ``tau``.\n
	The points of the face tau are contained in RR^dim and ``Vol(tau)`` is defined as follows:
	Let ``omega[tau]`` be the volume form on Aff(tau) such that the parallelepiped spanned by a 
	lattice  basis of ZZ^n intersect ``(aff tau)_0`` has volume 1. Then Vol(tau) is the volume of 
	tau intersection the Global Newton Polyhedron with respect to omega[tau].
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
			result = p.volume()	   # Returns dimtau!*n-volumen de tau
	else:
		result = 1
	return result

def face_divisors(d,faces_set,P):
	"""
	Returns a list of faces in ``faces_set`` such that d divides m(Delta_tau) = gcd{m(a)| a in 
	Delta_tau and ZZ^n}.
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
	Checks if own polynomial ``f`` over F_p with respect all the faces of the Global Newton's
	Polyhedron.
	If p = None, checks degeneration  over CC and (equivalent to be degenerated over F_p with 
	p>>0).
	For finite fields (``p`` is a given prime):
		- ``method = 'default'`` check the condition using evaluation over the (F_p^x)^n in the 
	system of equations
		- ``method = 'ideals'`` check the condition using ideals over the finite field.
	"""
	Q = f.newton_polytope()	#Global Newton Polyhedron of f
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
	in terms of symbolic variable ``s``:\n
		- ``local = True`` calculates the local Topological zeta function (at the origin).\n
		- ``weights`` -- a list of weights for the volume form.\n
		- ``info = True`` gives information of face ``tau``, cone of ``tau``, ``J_tau`` and
	``dim_tau!*Vol(tau)`` in the process.
		- ``check`` -- choose the method to check the non-degeneracy condition ('default' or \
	'ideals'). If ``check = 'no_check'``, degeneracy checking is omitted.
	"""
	s = var('s')
	P = newton_polyhedron(f)
	if check!='no_check':
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

def monodromy_zeta(f, weights = None, char = False, info = False):
	"""
	The Monodromy Zeta Function in the origin, in terms of variable ``s``.\n
	``weights`` is a list of weights if you want to considerate some weighting.\n
	``char = True`` prints the characteristic polynomial of the monodromy (only if ``f`` has 
	an isolated singularity in the origin).\n
	``info = True`` gives information of face ``tau``, cone of ``tau`` (all), ``J_tau`` and
	``dim_tau!*Vol(tau)`` in the process.
	"""
	n = f.parent().ngens()
	t = var('t')
	P = newton_polyhedron(f)
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



# -- INFORMATION ABOUT POLES IN TOPOLOGICAL ZETA
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
	Prints information about the candidate real poles for the topological zeta function of ``f``: 
	the orders and responsible faces of highest dimension.
	``local = True`` calculates the local (in the origin) Topological Zeta Function.\n
	``weights`` is a list of weights if you want to considerate some weighting.\n
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


