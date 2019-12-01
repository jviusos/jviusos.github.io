#  line_arrangements_vector_fields.sage   (Last update: 15-09-2015)
#  
#  Copyright 2015 Juan Viu-Sos <juan.viusos@univ-pau.fr>
#
#  For any bug or commentary, please contact the author.
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
# 

# TO DO LIST:
#   - generic_pols:
#       1) Extend to QQbar and AA
#       2) Only homogeneous components
#       3) pol_cond_to_matrix from any field and any pol

print "Suite of functions \"Filtration and dynamics of logarithmic vector fields of line arrangements in the plane\", coded by Juan Viu-Sos (<juan.viusos@univ-pau.fr>). To report any bug or commentary, please contact me.\n"

print "New functions for Line Arrangements declared: dual_point_to_eq, def_poly, eq2points, nlines, nslopes_and_parallels, singularities, nsing_multiplicities, info_arrangement, plot_arrangement.\n"

print "New functions for Log Vector Fields declared: matrix_eqs_vector_field, basis_vector_field_from_arrangement, line_is_invariant, curve_is_invariant, check_finiteness_base, check_finiteness_basis, check_finiteness_degree, singular_locus, linear_system, show_linear_system_by_points, plot_A_der.\n"

print "New (very) useful functions: basis_monomials, generic_pols, pol_cond_to_matrix, ordering_tuple, solve_to_points."




# -------------------------------------------------------------------------
# Arrangements
# -------------------------------------------------------------------------
def dual_point_to_eq(R,v,proj=False):
    """
    INPUT:
        - A polynomial ring ``R``.
        - A tuple `(a,b,c)` of elements of the base ring of ``R``.
        
    OUTPUT:
        - A linear form `ax+by+c` where `x,y` are the generators of ``R``.
        - If ``proj=True`` (by default, ``False``): A linear form `ax+by+cz` \
        where `x,y,z` are the generators of ``R``.
    """
    if proj: vgens=vector(R,list(R.gens()))
    else: vgens=vector(R,list(R.gens())+[1])
    eq=vgens.dot_product(vector(v))
    return eq


def def_poly(R,A,in_list=False,proj=False):
    """
    INPUT:
        - A polynomial ring ``R`` for defining polynomials.
        - A list `A=[(a_1,b_1,c_1),\ldots,(a_n,b_n,c_n)]` representing a \
        line arrangement `\mathcal{A}=\{\mathcal{l}_1,\ldots,\mathcal{l}_n\}`\
        whose elements are described in the dual projective plane \
        (`(a_i,b_i,c_i)\leftrightarrow L_i: a_ix+b_iy+c_i=0`).
        - If ``proj=True`` (by default, ``False``): `A` is consider as \
        a projective line arrangement.
        
    OUTPUT:
        - A product of affine/projective factors `\prod_{i=1}^nf_i` in ``R`` \
        representing the line arrangement \
        `\mathcal{A}=\{\mathcal{l}_1,\ldots,\mathcal{l}_n\}`.
        - If ``list=True`` (by default ``False``), returns a list with all \
        the equations.  
    """
    if proj: vgens=vector(R,list(R.gens()))
    else: vgens=vector(R,list(R.gens())+[1])
    eq_lines=map(lambda l: dual_point_to_eq(R,l,proj), A)
    if in_list==True: return eq_lines
    else: return prod(eq_lines)

def eq2points(R,P1,P2):
    """
    Enter 2 points to obtain the line which relates them in dual \
    coordinates, i.e. `(a:b:c)`.
    """
    return vector(R,(P2[1]-P1[1], P1[0]-P2[0],
                     (P2[0]-P1[0])*P1[1]-(P2[1]-P1[1])*P1[0]))

def nlines(A):
    """
    Returns the number of different lines contained in ``A``.
    """
    n=len(Set(A))
    if n is not len(A):
        print "Warning: There are repeated lines in the arrangement."
    return n

def nslopes_and_parallels(A,info=False):
    """
    Returns a tuple `(s,p)` where:
        - `s` is the number of different slopes contained in ``A``.\n
        - `p` is the maximum number of parallel lines contained in ``A``.\n
    If ``info=True``: detailed information is printed.
    """
    L=[]
    for l in A:
        if l[1]==0: L = L + [oo]
        else: L = L + [-l[0]/l[1]]
    slopes=Set(L)
    p=max(map(lambda s: L.count(s), slopes))
    if info: print 'Slopes and their lines:  ' + \
                    str(map(lambda s: (s,L.count(s)), slopes))
    return (len(slopes), p)

def singularities(R,A):
    """
    Returns a list with the singular points of the arrangement (with \
    coordinates in ``R``).
    """
    x=R.gens()
    Aeqs=def_poly(R,A,in_list=True)
    Sing=[]
    n=nlines(A)
    for i in range(n):
        for j in range(i,n):
            p=R.ideal([Aeqs[i],Aeqs[j]]).groebner_basis()
            if p[0].coefficient(x[0]) is not 0 and len(p)!=1:
                Sing.append((-p[0].constant_coefficient() / \
                            p[0].coefficient(x[0]),
                            -p[1].constant_coefficient()/ \
                            p[1].coefficient(x[1])))
            elif len(p)!=1:
                Sing.append((-p[0].constant_coefficient()/ \
                            p[0].coefficient(x[1]),
                            -p[1].constant_coefficient()/ \
                            p[1].coefficient(x[0])))
            elif len(p)!=1 and(p[0].degree()>1 or p[1].degree()>1):
                print "Warning: multiple solutions!"
    return list(Set(Sing))

def nsing_multiplicities(R,A):
    """
    Returns a list of tuples `[(1,n_1),\ldot,(k,n_k)]` where `k` is the \
    maximal multiplicity of the singular points of the arrangement and \
    `n_i` are the number of singular points of multiplicity `i`.
    """
    n = nlines(A); Aeqs=def_poly(R,A,in_list=True)
    sing = singularities(R,A)
    L=map(lambda i: vector((i,0)), range(2,n+1))
    for s in sing:
        k = len(filter(lambda p: p(s)==0, Aeqs))
        L[k-2] = L[k-2] + vector((0,1))
    while L[-1][1]==0: L.remove(L[-1])
    #L = filter(lambda v: v[1]!=0, L)
    return L


def info_arrangement(R,A,plot=True):
    """
    INPUT:
        - A ring ``R`` of defining polynomials.
        - A list `A=[(a_1,b_1,c_1),\ldots,(a_n,b_n,c_n)]` representing a \
        line arrangement `\mathcal{A}=\{\mathcal{l}_1,\ldots,\mathcal{l}_n\}`\
        whose elements are described in the dual projective plane \
        (`(a_i,b_i,c_i)\leftrightarrow L_i: a_ix+b_iy+c_i=0`).
        - A positive integer``d``.
        
    OUTPUT:
        - Prints a list of basic informations of ``A`` as the number of \
        lines, number of singularities, number of slopes, maximal number of \
        parallels.
        - Shows a plot of the arrangement containing of the singular points.
    """
    sing=singularities(R,A)
    print "Number of lines: " + str(nlines(A))
    print "Number of singularities: " + str(len(sing)) + '=' + \
                                         str(nsing_multiplicities(R,A))
    sandp=nslopes_and_parallels(A,info=True)
    print "Number of slopes: " + str(sandp[0])
    print "Maximal number of parallels: " + str(sandp[1])
    if plot:
        T=R.base_ring()
        x0=min(map(lambda s: T(s[0]), sing))-1/2
        x1=max(map(lambda s: T(s[0]), sing))+1/2
        y0=min(map(lambda s: T(s[1]), sing))-1/2
        y1=max(map(lambda s: T(s[1]), sing))+1/2
        show(plot_arrangement(R,A,(x0,x1),(y0,y1),color='red'))



# -------------------------------------------------------------------------
# Logarithmic vector fields
# -------------------------------------------------------------------------
def matrix_eqs_vector_field(R,A,d,info=False):
    """
    INPUT:
        - A ring ``R`` of defining polynomials.
        - A list `A=[(a_1,b_1,c_1),\ldots,(a_n,b_n,c_n)]` representing a \
        line arrangement `\mathcal{A}=\{\mathcal{l}_1,\ldots,\mathcal{l}_n\}`\
        whose elements are described in the dual projective plane \
        (`(a_i,b_i,c_i)\leftrightarrow L_i: a_ix+b_iy+c_i=0`).
        - A positive integer``d``.
        
    OUTPUT:
        - A matrix of conditions for with a polynomial vector fields in the \
        plane `\chi=P(x,y)\partial_x + Q(x,y)\partial_y`, where \
        `P,Q\in R` such that `\mathop{deg}P,\mathop{deg}Q\leq d`, and \
        `\mathcal{A}` is invariant by `\chi`. Each row of the matrix \
        represents a linear equation into with respect to coefficients \
        `\{a_{i,j}, b_{k,l}\}` of `P` and `Q` in generic form, i.e. \
        `P=\sum_{i,j=0}^{i+j\leq d}a_{i,j}\cdot x^i y^j` and \
        `Q=\sum_{k,l=0}^{k+l\leq d}b_{k,l}\cdot x^k y^l`. To obtain the \
        ordered set of basis coefficients, take ``factor=True`` \
        (by default ``False``).
    """
    [S1,P]=generic_pols(R,d,'a'); [S2,Q]=generic_pols(S1,d,'b',extension=True)
    eqs=[]
    for l in A:
        if l[1]==0:
            X=S2(-l[2]/l[0]); Y=S2(y)
        else:
            X=S2(-l[1]*y); Y=S2(l[0]*y-l[2]/l[1])
        P1 = l[0]*P(X,Y) + l[1]*Q(X,Y)
        eqs = eqs + P1.coefficients()
    S=S2.base_ring(); Sgens = (S2.base_ring()).gens()
    L=[]
    for e in eqs:
        L = L + map(lambda x: (S(e).numerator()).coefficient({x:1})/ \
                                S(e).denominator(), Sgens)
    M = matrix(S.base_ring(),len(eqs),(d+1)*(d+2),L)
    #matrix composed by eqs by rows and the 2x(d+1)*(d+2)/2 coefficients
    #of the polynomials
    if info==True: print S2.base_ring().gens()
    return M

def basis_vector_field_from_arrangement(R,A,d,factor=False):
    """
    INPUT:
        - A ring ``R`` of defining polynomials.
        - A list `A=[(a_1,b_1,c_1),\ldots,(a_n,b_n,c_n)]` representing a \
        line arrangement `\mathcal{A}=\{\mathcal{l}_1,\ldots,\mathcal{l}_n\}`\
        whose elements are described in the dual projective plane \
        (`(a_i,b_i,c_i)\leftrightarrow L_i: a_ix+b_iy+c_i=0`).
        - A positive integer``d``.
        
    OUTPUT:
        - A list contained a basis of polynomial vector fields in the plane \
        `\chi=P(x,y)\partial_x + Q(x,y)\partial_y`, where \
        `P,Q\in R` such that `\mathop{deg}P,\mathop{deg}Q\leq d`, and \
        `\mathcal{A}` is invariant by `\chi`.
        - If ``factor=True`` (by default ``False``): polynomials are given\
        in factorized form.

    EXAMPLES::

        sage: R.<x,y>=QQ[]; R
        Multivariate Polynomial Ring in x, y over Rational Field
        sage: A1=[(1,0,1),(1,0,-1),(0,1,1),(0,1,-1),(1,1,0),(1,-1,0)]
        sage: basis_vector_field_from_arrangement(R,A1,3)
        [[-x^3 + x, -y^3 + y], [-x^2*y + y, -x*y^2 + x]]
        sage: basis_vector_field_from_arrangement(R,A1,3,factor=True)
        [[(-1) * x * (x - 1) * (x + 1), (-1) * y * (y - 1) * (y + 1)],  \
        [(-1) * y * (x - 1) * (x + 1), (-1) * (y - 1) * (y + 1) * x]]
        sage: #Over a Number Field
        sage: var('t'); K.<sqrt5> = NumberField(t^2-5,embedding=2.23)
        sage: R.<x,y>=K[]; R
        Multivariate Polynomial Ring in x, y over Number Field in sqrt5 \
        with defining polynomial t^2 - 5
        sage: A2=[(1,0,sqrt5),(1,0,-1),(0,1,1),(0,1,-1),(1,1,0),(1,-1,0),\
        (sqrt5,0,1)]
        sage: basis_vector_field_from_arrangement(R,A2,4)
        [[0, -x^2*y^2 + y^4 + x^2 - y^2]]
        sage: #Over the Real (or Complex) Field
        sage: R.<x,y>=RR[]; R
        Multivariate Polynomial Ring in x, y over Real Field with 53 bits of \
        precision
        sage: A3=[(1,0,sqrt(5)),(1,0,-1),(0,1,1)]
        sage: basis_vector_field_from_arrangement(R,A3,5)
        [[-0.447213595499958*x^2 - 0.552786404500042*x + 1.00000000000000, 0]\
        ,  [0, -y^2 + 1.00000000000000], [0, x*y + x], [0, y^2 + y]]
    """
    M=matrix_eqs_vector_field(R,A,d)
    S=M.transpose().kernel(); Sb=S.basis()
    m=basis_monomials(R,d); k=(d+1)*(d+2)/2
    L=map(lambda v: vector(R,(v[0:k].dot_product(vector(R,m)),
                     v[k:2*k].dot_product(vector(R,m)))), Sb)
    #if len(L)==0: L = [[0,0]]
    if factor==True: L=map(lambda v: (v[0].factor(),v[1].factor()), L)
    return L

def line_is_invariant(R,L,V):
    """
    Returns if the line ``L`` written in dual coordinates `(a:b:c)` is \
    invariant by the vector field `V=[P(x,y),Q(x,y)]`. 
    """
    f=dual_point_to_eq(R,L)
    Vf= L[0]*V[0]+L[1]*V[1]
    return f.divides(Vf)

def curve_is_invariant(R,f,V):
    """
    Returns a boolean with the information of the invariance by ``V`` of the \
    algebraic curve defined by `\mathcal{C}=\{f=0\}`.
    """
    x=R.gens()
    Vf=V[0]*f.derivative(x[0])+V[1]*f.derivative(x[1])
    return f.divides(Vf)

def fixes_infinity_lines(R,V,(x0,y0)=(0,0)):
    """
    Test if a vector field `W` fixes a pencil of lines centered in `(x0,y0)` \
    (by default the origin (0,0)) or an infinity parallel lines.
    """
    Rgens= R.gens()
    if (Rgens[1]-y0)*V[0]-(Rgens[0]-x0)*V[1]==0: return True
    if V[1]==0: return True
    else:
        f=V[0]/V[1]
        if f.derivative(Rgens[0])==0 and f.derivative(Rgens[1])==0:
            return True
    return False
    

def check_infiniteness_base(W,S,(x0,y0)=(0,0)):
    """
    Checks if a subspace of vector fields `W` parametrized by coordinates \
    explicit in `S` fixes a pencil of lines centered in `(x0,y0)` \
    (by default the origin (0,0)) or an infinity parallel lines.
    """
    L=[]
    for i in S:
        for j in S:
            s='r'+str(j)+'=0'; exec s
        s='r'+str(i)+'=1'; exec s
        L = L + [(s,fixes_infinity_lines(R,V,(x0,y0)))]
    return L

def check_infiniteness_basis(R,B,points):
    """
    Checks if a subspace of vector fields given by a list of the basis ``B`` \
    fixes a pencil of lines centered in any point `(x0,y0)` contained in the \
    list ``points``, or an infinity parallel lines.\n
    Returns a list of booleans corresponding to each element of the basis.
    """
    L=[]
    for p in points:
        L = L + [ [p] + map(lambda V: fixes_infinity_lines(R,V,p),B) ]
    return L


def check_infiniteness_degree(R,A,d):
    """
    Checks if the basis of the subspace of vector fields of degree ``d`` \
    fixes a pencil of lines centered in the singularities of `A` or an \
    infinity parallel lines.\n
    Returns a list of booleans corresponding to each element of the basis.
    """
    B = basis_vector_field_from_arrangement(R,A,d)
    L = check_finiteness_basis(R,B,singularities(R,A))
    return L

def singular_locus(R,V):
    """
    [Should be improved] Returns a list with the singular REAL points of the \
    vector field (with coordinates in ``R``).
    """
    if R.base_ring()!=QQ:
        print
        "Warning: this function don't work correctly over number fields."
    x=R.gens()
    L=solve([SR(V[0]), SR(V[1])], SR(x[0]), SR(x[1]))
    J = solve_to_points(R,L)
    return J

def linear_system(R,V,(x0,y0),jordan=True):
    """
    Returns the linearized vector field of ``V`` defined in the ring of/ 
    polynomials ``R`` in the point ``(x0,y0)``.\\
    If ``jordan=True`` (by default), the linearized system is given \
    directly in the jordan form.
    """
    Rbase=R.base_ring(); Rgens=R.gens()
    L=V[0].gradient() + V[1].gradient(); L
    M=matrix(R,2,L)
    N=M.subs({Rgens[0]:x0, Rgens[1]:y0}).change_ring(Rbase)
    if jordan==True: N=N.jordan_form(subdivide=False)
    return N

def show_linear_system_by_points(R,V,S,jordan=True):
    """
    Prints (typed by LateX) a relation between the linearized vector field \
    of ``V`` defined in the ring of polynomials ``R`` in the points \
    contained in ``S``.\\
    If ``jordan=True`` (by default), the linearized system is given directly \
    in the Jordan form.
    """
    s='['
    for p in SV: s = s + latex(p) + '\\to' + \
                        latex(linear_system(R,V,p,jordan)) + ', '
    s = s + ']'
    show(s)

def homogenize_vf(S,V):
    """
    Returns the homogenization of the polynomial vector field ``V`` in the \
    ring of homogeneous coordinates ``S``.
    """
    Rgens=list(V[0].parent().gens())
    Sgens=list(S.gens())
    z=(Set(Sgens).difference(Set(Rgens)))[0]
    d=max([V[0].degree(),V[1].degree()])
    W=vector(S,(z^d*S(V[0])(Rgens[0]/z,Rgens[1]/z,1),
                z^d*S(V[1])(Rgens[0]/z,Rgens[1]/z,1),0))
    return W

def extended_vf(S,V,A):
    """
    For a polynomial vector field in the plane \
    `\chi=P(x,y)\partial_x + Q(x,y)\partial_y` (represented by ``V``) \
    fixing the arrangement ``\mathcal{A}``, returns the associated \
    derivation `\\bar{\chi}` fixing the cone of ``\mathcal{A}``, defined by 
    .. MATH::

        \\bar{\chi}=\chi^\\text{hom}-\\frac{1}{|\mathcal{A}|+1}K\\theta_E`
    
    where `\chi^\\text{hom}` is the homogenization of `\chi` in the ring \
    ``S``, `K` is the cofactor of the cone of `mathcal{A}` by \
    `\chi^\\text{hom}`, and `\\theta_E` is the Euler derivation in ``S``.
    """
    Sgens=S.gens(); x=Sgens[0]; y=Sgens[1]; z=Sgens[2]
    F=def_poly(R,A,proj=True)*z
    Fx=F.derivative(x); Fy=F.derivative(y); Fz=F.derivative(z)
    Vh=homogenize_vf(S,V)
    VhF=Vh[0]*Fx + Vh[1]*Fy + Vh[2]*Fz
    if F.divides(VhF)==False:
        print "The vector field in the input IS NOT fixing the cone of \
        arrangement.";
    return None
    K=S(VhF/F); nlines=nlines(A)+1
    W=Vh - K/nlines*vector(S,(x,y,z))
    return W



# -------------------------------------------------------------------------
# Plots
# -------------------------------------------------------------------------
def plot_arrangement(R,A,(x0,x1),(y0,y1),color='red'):
    """
    Plots a real arrangement ``A`` in the box ``(x0,x1),(y0,y1)``.
    """
    Rgens=R.gens()
    return implicit_plot(R(def_poly(R, A)),(Rgens[0],x0,x1),
            (Rgens[1],y0,y1),color=color,plot_points=1000)

def plot_A_der(R,A,V,(x0,x1),(y0,y1)):
    """
    Plots a real arrangement ``A`` and the vector field ``V`` in the box \
    ``(x0,x1),(y0,y1)``.
    """
    Rgens=R.gens()
    return plot_vector_field((V[0],V[1]), (Rgens[0],x0,x1),
        (Rgens[1],y0,y1),color='blue') + plot_arrangement(R,A,(x0,x1),(y0,y1))


#def subs_in_list(R,L,d):
#    J = []
#    for e in L:
#        svars = map(lambda v: str(v), e.variables())
#        S = PolynomialRing(R.base_ring(), svars)
#        et=S(e)
#        J = J+[R(et.subs(d))]
#    return J




# -------------------------------------------------------------------------
# Extra (very) Useful Functions
# -------------------------------------------------------------------------
def basis_monomials(R,deg):
    """
    Returns an ordered list of monomials up to degree ``deg``.
    """
    rbase=R.base_ring()
    rgens=R.gens(); ngens=len(rgens)    #variables of R
    exps=filter(lambda e: sum(e)<=deg,Tuples(range(deg+1),ngens).list())
    #Exponents of monomials
    exps=ordering_tuple(R,map(lambda e: tuple(e), exps))
    #Tuple is ordered with respect to the ordering of R
    mons=[]
    for e in exps:
        m = 1
        for i in range(ngens):
            m = m*rgens[i]^e[i]  #associated monomial
        mons.append(m)
    return mons

def generic_pols(R,deg,parm='a',extension=False):
    """
    This function allows us to work with generic polynomials in ``R`` \
    of degree ``deg`` parametrized by coefficients `a_{i_1,\ldots,i_n}`\
    (``parm``='a', by default).\n
    If ``extension=True``, we suppose that the coefficient field has already \
    parameters and we do a symbolic field extension with the new ones.
    Returns a list ``[S,p]`` where:
        - ``S`` is the new polynomial ring with same variables of ``R`` and \
        extended coefficients `a_{i_1,\ldots,i_n}`.
        - ``p`` is a generic polynomial \
        `p=\sum_{i_1,\ldots,i_n=0}^{i_1+\ldots+i_n\leq deg}a_{i_1,\ldots,i_n}\
        \cdot x_1^{i_1}\ldots x_n^{i_n}`.
    """
    rbase=R.base_ring()
    rgens=R.gens(); ngens=len(rgens)    #variables of R
    exps=filter(lambda e: sum(e)<=deg,Tuples(range(deg+1),ngens).list())
    #Exponents of monomials
    exps=ordering_tuple(R,map(lambda e: tuple(e), exps))
    #Tuple is ordered with respect to the ordering of R
    parms, mons=[],[]
    for e in exps:
        s = str(e[0]); m = rgens[0]^e[0]
        for i in range(1,ngens):
            s = s +'_' + str(e[i])    #parameter
            m = m*rgens[i]^e[i]  #associated monomial
        parms.append(parm + s); mons.append(m)
    if extension is True:
    #If the coefficient field has already parameters, we do a field extension
        parms_ext = map(lambda x: str(x), rbase.gens()) + parms
        S=PolynomialRing(
          FractionField(PolynomialRing(rbase.base_ring(), parms_ext)),rgens)
    else:
        S=PolynomialRing(FractionField(PolynomialRing(rbase, parms)),rgens)
    gens=map(lambda p: S(p), parms)
    p=S(vector(S,gens)*vector(S,mons))
    return [S,p]

def pol_cond_to_matrix(R,F,d):
    """
    INPUT:
        - A polynomial equation `F=F(P,Q)=0` in two generic polynomials \
        `P=\sum_{i_1,\ldots,i_n=0}^{i_1+\ldots+i_n\leq d}a_{i_1,\ldots,i_n}\
        \cdot x_1^{i_1}\ldots x_n^{i_n}` and \
        `Q=\sum_{i_1,\ldots,i_n=0}^{i_1+\ldots+i_n\leq d}b_{i_1,\ldots,i_n}\
        \cdot x_1^{i_1}\ldots x_n^{i_n}` of degree `d`.
    
    OUTPUT:
        - A matrix which entries in each row correspond to the coefficients \
        in `(a_{ij},b_{ij})` for each coefficient in \
        `x_1^{i_1}\ldots x_n^{i_n}` in `F`. 
    """
    Rgens = (R.base_ring()).gens()
    eqs=F.coefficients()
    L=[]
    for e in eqs:
        en=e.numerator()
        L = L + map(lambda v: en.coefficient({v:1}), Rgens)
    M = matrix(R.base_ring(),len(eqs),(d+1)*(d+2),L)
    return M

def solve_to_points(R,L):
    """
    Traduces an expression returned by the function``solve``\
    representing points into a list REAL of solutions in ``R``.
    """
    Sol=[]; iSol=[]; oSol=[]
    for e in L:
        if len(e)==2:
            s=[str(e[0].rhs()), str(e[1].rhs())]
            if 'I' in s[0] or 'I' in s[1]: iSol.append((e[0].rhs(),
                                                         e[0].rhs()))
            else:            
                try: Sol.append((R(s[0]), R(s[1])))
                except Exception as err: oSol.append((e[0].rhs(),
                                                        e[0].rhs()))
        else: oSol.append(str(e))
    if len(iSol)!=0: print "Warning! Complex solutions: " + str(iSol)
    if len(oSol)!=0: print "Warning! Other solutions: " + str(oSol)
    return Sol

def ordering_tuple(R,L):
    """
    INPUT:
        - A polynomial ring ``R`` of `n` variables.
        - A list of `n`-tuples of integers.
    
    OUTPUT:
        - An copy of ``L`` with the inherited ordering of ``R``.  
    """
    p=0
    for e in L: p = p + R({e:1})
    return p.exponents()
