#ifndef _QUADSOLVE_HPP_
#define _QUADSOLVE_HPP_


/*
 FILE uquadprog.hh
 
 NOTE: this is a modified of QuadProg++ package, originally developed by 
       Luca Di Gaspero, working with ublas data structures. 

 The quadprog_solve() function implements the algorithm of Goldfarb and Idnani 
 for the solution of a (convex) Quadratic Programming problem
by means of a dual method.
	 
The problem is in the form:

min 0.5 * x G x + g0 x
s.t.
    CE^T x + ce0 = 0
    CI^T x + ci0 >= 0
	 
 The matrix and vectors dimensions are as follows:
     G: n * n
		g0: n
				
		CE: n * p
	 ce0: p
				
	  CI: n * m
   ci0: m

     x: n
 
 The function will return the cost of the solution written in the x vector or
 std::numeric_limits::infinity() if the problem is infeasible. In the latter case
 the value of the x vector is not correct.
 
 References: D. Goldfarb, A. Idnani. A numerically stable dual method for solving
             strictly convex quadratic programs. Mathematical Programming 27 (1983) pp. 1-33.

 Notes:
  1. pay attention in setting up the vectors ce0 and ci0. 
	   If the constraints of your problem are specified in the form 
	   A^T x = b and C^T x >= d, then you should set ce0 = -b and ci0 = -d.
  2. The matrix G is modified within the function since it is used to compute
     the G = L^T L cholesky factorization for further computations inside the function. 
     If you need the original matrix G you should make a copy of it and pass the copy
     to the function.
    
 Author: Angelo Furfaro
  			 DEIS - University of Calabria, Italy
				 a.furfaro@dimes.unical.it
				 http://dev.dimes.unical.it/~angelo
 
 The author will be grateful if the researchers using this software will
 acknowledge the contribution of this modified function and of Di Gaspero's
 original version in their research papers.


LICENSE

Copyright (2016) Angelo Furfaro
Copyright (2008) Angelo Furfaro
Copyright (2006) Luca Di Gaspero


This file is a porting of QuadProg++ routine, originally developed
by Luca Di Gaspero, exploiting uBlas data structures for vectors and
matrices instead of native C++ array.

uquadprog is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

uquadprog is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with uquadprog; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/








#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

using namespace boost::numeric::ublas;


void cholesky_decomposition(matrix<double> & A);
void backward_elimination(matrix<double>& U,vector<double>& x, vector<double>& y);
void forward_elimination(matrix<double>& L, vector<double>& y, vector<double>& b);
void cholesky_solve(matrix<double>& L, vector<double> & x, vector<double> & b);

double distance(double a, double b);
void compute_d(vector<double> &d, matrix<double>& J, vector<double>& np);
void update_z(vector<double>& z, matrix<double>& J, vector<double>& d,  int iq);
void update_r(matrix<double>& R, vector<double> &r, vector<double> &d, int iq);
bool add_constraint(matrix<double>& R, matrix<double>& J, vector<double>& d, int& iq, double& R_norm);
void delete_constraint(matrix<double>& R, matrix<double>& J, vector<int>& A, vector<double>& u,  int p, int& iq, int l);


inline double solve_quadprog( matrix<double> & G,  vector<double> & g0,  
                      const matrix<double> & CE, const vector<double> & ce0,  
                      const matrix<double> & CI, const vector<double> & ci0, 
                      vector<double>& x)
{
  int i, j, k, l; /* indices */
  int ip, me, mi;
  int n=g0.size();  int p=ce0.size();  int m=ci0.size();  
  matrix<double> R(G.size1(),G.size2()), J(G.size1(),G.size2());
 
  vector<double> s(m+p), z(n), r(m + p), d(n),  np(n), u(m + p);
  vector<double> x_old(n), u_old(m + p);
  double f_value, psi, c1, c2, sum, ss, R_norm;
  const double inf = std::numeric_limits<double>::infinity();
  double t, t1, t2; /* t is the step lenght, which is the minimum of the partial step length t1 
    * and the full step length t2 */
  vector<int> A(m + p), A_old(m + p), iai(m + p);int q;
  int iq, iter = 0;
  bool iaexcl[m + p];
 	
  me = p; /* number of equality constraints */
  mi = m; /* number of inequality constraints */
  q = 0;  /* size of the active set A (containing the indices of the active constraints) */
  
  /*
   * Preprocessing phase
   */
	
  /* compute the trace of the original matrix G */
	c1 = 0.0;
	for (i = 0; i < n; i++)
	{
		c1 += G(i,i);
	}
	/* decompose the matrix G in the form L^T L */
	cholesky_decomposition(G);
 
/* initialize the matrix R */
	for (i = 0; i < n; i++)
	{
		d(i) = 0.0;
		for (j = 0; j < n; j++)
			R(i,j) = 0.0;
	}
	R_norm = 1.0; /* this variable will hold the norm of the matrix R */
  
	/* compute the inverse of the factorized matrix G^-1, this is the initial value for H */
	c2 = 0.0;
  for (i = 0; i < n; i++) 
	{
		d(i) = 1.0;
		forward_elimination(G, z, d);
		for (j = 0; j < n; j++)
			J(i,j) = z(j);
		c2 += z(i);
    d(i) = 0.0;
	}
#ifdef TRACE_SOLVER
 print_matrix("J", J, n);
#endif
  
	/* c1 * c2 is an estimate for cond(G) */
  
	/* 
   * Find the unconstrained minimizer of the quadratic form 0.5 * x G x + g0 x 
   * this is a feasible point in the dual space
	 * x = G^-1 * g0
   */
  cholesky_solve(G, x, g0);
  for (i = 0; i < n; i++)
		x(i) = -x(i);
	/* and compute the current solution value */ 
	f_value = 0.5 * inner_prod(g0, x);
#ifdef TRACE_SOLVER
  std::cerr << "Unconstrained solution: " << f_value << std::endl;
  print_vector("x", x, n);
#endif
  
	/* Add equality constraints to the working set A */
  iq = 0;
	for (i = 0; i < me; i++)
	{
		for (j = 0; j < n; j++)
			np(j) = CE(j,i);
    compute_d(d, J, np);
		update_z(z, J, d,  iq);
		update_r(R, r, d,  iq);
#ifdef TRACE_SOLVER
		print_matrix("R", R, iq);
		print_vector("z", z, n);
		print_vector("r", r, iq);
		print_vector("d", d, n);
#endif
    
    /* compute full step length t2: i.e., the minimum step in primal space s.t. the contraint 
      becomes feasible */
		t2 = 0.0;
		if (fabs(inner_prod(z, z)) > std::numeric_limits<double>::epsilon()) // i.e. z != 0
			t2 = (-inner_prod(np, x) - ce0(i)) / inner_prod(z, np);
    
    /* set x = x + t2 * z */
		for (k = 0; k < n; k++)
			x(k) += t2 * z(k);
			
    /* set u = u+ */
		u(iq) = t2;
		for (k = 0; k < iq; k++)
			u(k) -= t2 * r(k);
    
		/* compute the new solution value */
		f_value += 0.5 * (t2 * t2) * inner_prod(z, np);
		A(i) = -i - 1;
    
		if (!add_constraint(R, J, d, iq, R_norm))
		{
			// FIXME: it should raise an error
			// Equality constraints are linearly dependent
			return f_value;
		}
	}
  
	/* set iai = K \ A */
	for (i = 0; i < mi; i++)
		iai(i) = i;
  
l1:	iter++;
#ifdef TRACE_SOLVER
  print_vector("x", x, n);
#endif
  /* step 1: choose a violated constraint */
	for (i = me; i < iq; i++)
	{
	  ip = A(i);
		iai(ip) = -1;
	}
	
	/* compute s(x) = ci^T * x + ci0 for all elements of K \ A */
	ss = 0.0;
	psi = 0.0; /* this value will contain the sum of all infeasibilities */
	ip = 0; /* ip will be the index of the chosen violated constraint */
	for (i = 0; i < mi; i++)
	{
		iaexcl[i] = true;
		sum = 0.0;
		for (j = 0; j < n; j++)
			sum += CI(j,i) * x(j);
		sum += ci0(i);
		s(i) = sum;
		psi += std::min(0.0, sum);
	}
#ifdef TRACE_SOLVER
  print_vector("s", s, mi);
#endif

    
	if (fabs(psi) <= mi * std::numeric_limits<double>::epsilon() * c1 * c2* 100.0)
	{
    /* numerically there are not infeasibilities anymore */
    q = iq;
		return f_value;
  }
    
  /* save old values for u and A */
	for (i = 0; i < iq; i++)
	{
		u_old(i) = u(i);
		A_old(i) = A(i);
	}
   /* and for x */
	for (i = 0; i < n; i++)
		x_old(i) = x(i);
    
l2: /* Step 2: check for feasibility and determine a new S-pair */
	for (i = 0; i < mi; i++)
	{
		if (s(i) < ss && iai(i) != -1 && iaexcl[i])
		{
			ss = s(i);
			ip = i;
		}
	}
  if (ss >= 0.0)
  {
    q = iq;
    return f_value;
  }
    
  /* set np = n(ip) */
  for (i = 0; i < n; i++)
    np(i) = CI(i,ip);
  /* set u = (u 0)^T */
  u(iq) = 0.0;
  /* add ip to the active set A */
  A(iq) = ip;

#ifdef TRACE_SOLVER
	std::cerr << "Trying with constraint " << ip << std::endl;
	print_vector("np", np, n);
#endif
    
l2a:/* Step 2a: determine step direction */
  /* compute z = H np: the step direction in the primal space (through J, see the paper) */
  compute_d(d, J, np);
  update_z(z, J, d, iq);
  /* compute N* np (if q > 0): the negative of the step direction in the dual space */
  update_r(R, r, d, iq);
#ifdef TRACE_SOLVER
  std::cerr << "Step direction z" << std::endl;
		print_vector("z", z, n);
		print_vector("r", r, iq + 1);
    print_vector("u", u, iq + 1);
    print_vector("d", d, n);
    print_ivector("A", A, iq + 1);
#endif
    
  /* Step 2b: compute step length */
	l = 0;
  /* Compute t1: partial step length (maximum step in dual space without violating dual feasibility */
	t1 = inf; /* +inf */
  /* find the index l s.t. it reaches the minimum of u+(x) / r */
	for (k = me; k < iq; k++)
	{
		if (r(k) > 0.0)
		{
			if (u(k) / r(k) < t1)
			{
				t1 = u(k) / r(k);
				l = A(k);
       }
     }
  }
  /* Compute t2: full step length (minimum step in primal space such that the constraint ip becomes feasible */
	if (fabs(inner_prod(z, z))  > std::numeric_limits<double>::epsilon()) // i.e. z != 0
	  t2 = -s(ip) / inner_prod(z, np);
  else
    t2 = inf; /* +inf */
  
  /* the step is chosen as the minimum of t1 and t2 */
  t = std::min(t1, t2);
#ifdef TRACE_SOLVER
  std::cerr << "Step sizes: " << t << " (t1 = " << t1 << ", t2 = " << t2 << ") ";
#endif
  
  /* Step 2c: determine new S-pair and take step: */
  
  /* case (i): no step in primal or dual space */
  if (t >= inf)
  {
    /* QPP is infeasible */
    // FIXME: unbounded to raise
    q = iq;
    return inf;
  }
  /* case (ii): step in dual space */
  if (t2 >= inf)
  {
    /* set u = u +  t * [-r 1) and drop constraint l from the active set A */
    for (k = 0; k < iq; k++)
      u(k) -= t * r(k);
    u(iq) += t;
    iai(l) = l;
    delete_constraint(R, J, A, u, p, iq, l);
#ifdef TRACE_SOLVER
    std::cerr << " in dual space: " 
      << f_value << std::endl;
    print_vector("x", x, n);
    print_vector("z", z, n);
		print_ivector("A", A, iq + 1);
#endif
    goto l2a;
  }
  
  /* case (iii): step in primal and dual space */
  
  /* set x = x + t * z */
  for (k = 0; k < n; k++)
    x(k) += t * z(k);
  /* update the solution value */
  f_value += t * inner_prod(z, np) * (0.5 * t + u(iq));
  /* u = u + t * (-r 1) */
  for (k = 0; k < iq; k++)
		u(k) -= t * r(k);
  u(iq) += t;
#ifdef TRACE_SOLVER
  std::cerr << " in both spaces: " 
    << f_value << std::endl;
	print_vector("x", x, n);
	print_vector("u", u, iq + 1);
	print_vector("r", r, iq + 1);
	print_ivector("A", A, iq + 1);
#endif
  
  if (t == t2)
  {
#ifdef TRACE_SOLVER
    std::cerr << "Full step has taken " << t << std::endl;
    print_vector("x", x, n);
#endif
    /* full step has taken */
    /* add constraint ip to the active set*/
		if (!add_constraint(R, J, d, iq, R_norm))
		{
			iaexcl[ip] = false;
			delete_constraint(R, J, A, u, p, iq, ip);
#ifdef TRACE_SOLVER
      print_matrix("R", R, n);
      print_ivector("A", A, iq);
#endif
			for (i = 0; i < m; i++)
				iai(i) = i;
			for (i = 0; i < iq; i++)
			{
				A(i) = A_old(i);
				iai(A(i)) = -1;
				u(i) = u_old(i);
			}
			for (i = 0; i < n; i++)
				x(i) = x_old(i);
      goto l2; /* go to step 2 */
		}    
    else
      iai(ip) = -1;
#ifdef TRACE_SOLVER
    print_matrix("R", R, n);
    print_ivector("A", A, iq);
#endif
    goto l1;
  }
  
  /* a patial step has taken */
#ifdef TRACE_SOLVER
  std::cerr << "Partial step has taken " << t << std::endl;
  print_vector("x", x, n);
#endif
  /* drop constraint l */
	iai(l) = l;
	delete_constraint(R, J, A, u, p, iq, l);
#ifdef TRACE_SOLVER
  print_matrix("R", R, n);
  print_ivector("A", A, iq);
#endif
  
  /* update s[ip) = CI * x + ci0 */
	sum = 0.0;
	for (k = 0; k < n; k++)
		sum += CI(k,ip) * x(k);
	s(ip) = sum + ci0(ip);

#ifdef TRACE_SOLVER
  print_vector("s", s, mi);
#endif
  goto l2a;
}










inline void compute_d(vector<double> &d, matrix<double>& J, vector<double>& np)
{ int n=np.size();
  register int i, j;
  register double sum;
  
  /* compute d = H^T * np */
	for (i = 0; i < n; i++)
	{
		sum = 0.0;
		for (j = 0; j < n; j++)
			sum += J(j,i) * np(j);
		d(i) = sum;
	}
}

inline void update_z(vector<double>& z, matrix<double>& J, vector<double>& d,  int iq)
{
	register int i, j;
	int n=z.size();
	/* setting of z = H * d */
  for (i = 0; i < n; i++)
	{
		z(i) = 0.0;
		for (j = iq; j < n; j++)
			z(i) += J(i,j) * d(j);
	}
}

inline void update_r(matrix<double>& R, vector<double> &r, vector<double> &d, int iq) 
{
	register int i, j;
	register double sum;
  
  /* setting of r = R^-1 d */
	for (i = iq - 1; i >= 0; i--)
	{
		sum = 0.0;
		for (j = i + 1; j < iq; j++)
			sum += R(i,j) * r(j);
		r(i) = (d(i) - sum) / R(i,i);
	}
}




inline bool add_constraint(matrix<double>& R, matrix<double>& J, vector<double>& d, int& iq, double& R_norm)
{
 int n=J.size1();
#ifdef TRACE_SOLVER
  std::cerr << "Add constraint " << iq << '/';
#endif
	int i, j, k;
	double cc, ss, h, t1, t2, xny;
	
  /* we have to find the Givens rotation which will reduce the element
		d(j) to zero.
		if it is already zero we don't have to do anything, except of
		decreasing j */  
	for (j = n - 1; j >= iq + 1; j--)
	{
    /* The Givens rotation is done with the matrix (cc cs, cs -cc).
			 If cc is one, then element (j) of d is zero compared with element
			 (j - 1). Hence we don't have to do anything. 
			 If cc is zero, then we just have to switch column (j) and column (j - 1) 
			 of J. Since we only switch columns in J, we have to be careful how we
			 update d depending on the sign of gs.
			 Otherwise we have to apply the Givens rotation to these columns.
			 The i - 1 element of d has to be updated to h. */
		cc = d(j - 1);
		ss = d(j);
		h = distance(cc, ss);
		if (h == 0.0)
			continue;
		d(j) = 0.0;
		ss = ss / h;
		cc = cc / h;
		if (cc < 0.0)
		{
			cc = -cc;
			ss = -ss;
			d(j - 1) = -h;
		}
		else
			d(j - 1) = h;
		xny = ss / (1.0 + cc);
		for (k = 0; k < n; k++)
		{
			t1 = J(k,j - 1);
			t2 = J(k,j);
			J(k,j - 1) = t1 * cc + t2 * ss;
			J(k,j) = xny * (t1 + J(k,j - 1)) - t2;
		}
	}
  /* update the number of constraints added*/
	iq++;
  /* To update R we have to put the iq components of the d vector
    into column iq - 1 of R
    */
	for (i = 0; i < iq; i++)
		R(i,iq - 1) = d(i);
#ifdef TRACE_SOLVER
  std::cerr << iq << std::endl;
#endif
  
	if (fabs(d(iq - 1)) <= std::numeric_limits<double>::epsilon() * R_norm)
		// problem degenerate
		return false;
	R_norm = std::max<double>(R_norm, fabs(d(iq - 1)));
	return true;
}


inline void delete_constraint(matrix<double>& R, matrix<double>& J, vector<int>& A, vector<double>& u,  int p, int& iq, int l)
{

  int n=R.size1();
#ifdef TRACE_SOLVER
  std::cerr << "Delete constraint " << l << ' ' << iq;
#endif
	int i, j, k, qq;
	double cc, ss, h, xny, t1, t2;
  
	/* Find the index qq for active constraint l to be removed */
  for (i = p; i < iq; i++)
		if (A(i) == l)
		{
			qq = i;
			break;
		}
      
  /* remove the constraint from the active set and the duals */
  for (i = qq; i < iq - 1; i++)
  {
    A(i) = A(i + 1);
    u(i) = u(i + 1);
    for (j = 0; j < n; j++)
      R(j,i) = R(j,i + 1);
  }
      
  A(iq - 1) = A(iq);
  u(iq - 1) = u(iq);
  A(iq) = 0; 
  u(iq) = 0.0;
  for (j = 0; j < iq; j++)
    R(j,iq - 1) = 0.0;
  /* constraint has been fully removed */
  iq--;
#ifdef TRACE_SOLVER
  std::cerr << '/' << iq << std::endl;
#endif 
  
  if (iq == 0)
    return;
  
  for (j = qq; j < iq; j++)
  {
    cc = R(j,j);
    ss = R(j + 1,j);
    h = distance(cc, ss);
    if (h == 0.0)
      continue;
    cc = cc / h;
    ss = ss / h;
    R(j + 1,j) = 0.0;
    if (cc < 0.0)
    {
      R(j,j) = -h;
      cc = -cc;
      ss = -ss;
    }
    else
      R(j,j) = h;
    
    xny = ss / (1.0 + cc);
    for (k = j + 1; k < iq; k++)
    {
      t1 = R(j,k);
      t2 = R(j + 1,k);
      R(j,k) = t1 * cc + t2 * ss;
      R(j + 1,k) = xny * (t1 + R(j,k)) - t2;
    }
    for (k = 0; k < n; k++)
    {
      t1 = J(k,j);
      t2 = J(k,j + 1);
      J(k,j) = t1 * cc + t2 * ss;
      J(k,j + 1) = xny * (J(k,j) + t1) - t2;
    }
  }
}




inline double distance(double a, double b)
{
	register double a1, b1, t;
	a1 = fabs(a);
	b1 = fabs(b);
	if (a1 > b1) 
	{
		t = (b1 / a1);
		return a1 * sqrt(1.0 + t * t);
	}
	else
		if (b1 > a1)
		{
			t = (a1 / b1);
			return b1 * sqrt(1.0 + t * t);
		}
	return a1 * sqrt(2.0);
}


inline void cholesky_decomposition(matrix<double> & A) 
{
  register int i, j, k;
  register double sum;
  int n=A.size1();	
  for (i = 0; i < n; i++)
  {
    for (j = i; j < n; j++)
    {
      sum = A(i,j);
      for (k = i - 1; k >= 0; k--)
        sum -= A(i,k)*A(j,k);
      if (i == j) 
      {
        if (sum <= 0.0)
				{
          // raise error
			//		print_matrix("A", A, n);
          throw std::runtime_error("The matrix passed to the Cholesky A = L L^T decomposition is not positive definite");
				}
        A(i,i) = sqrt(sum);
      }
      else
        A(j,i) = sum / A(i,i);
    }
    for (k = i + 1; k < n; k++)
      A(i,k) = A(k,i);
  } 
}




inline void cholesky_solve(matrix<double>& L, vector<double> &x, vector<double> & b)
{
  
  vector<double> y(x.size());
	
	/* Solve L * y = b */
	forward_elimination(L, y, b);
  /* Solve L^T * x = y */
	backward_elimination(L, x, y);
}


inline void forward_elimination(matrix<double>& L, vector<double>& y, vector<double>& b)
{ 
	register int i, j;
	int n=y.size();
	y(0) = b(0) / L(0,0);
	for (i = 1; i < n; i++)
	{
		y(i) = b(i);
		for (j = 0; j < i; j++)
			y(i) -= L(i,j) * y(j);
		y(i) = y(i) / L(i,i);
	}
}





inline void backward_elimination(matrix<double>& U,vector<double>& x, vector<double>& y){
  int n=x.size();

  register int i, j;
	
  x(n - 1) = y(n - 1) / U(n - 1,n - 1);
  for (i = n - 2; i >= 0; i--){
    x(i) = y(i);
    for (j = i + 1; j < n; j++)
      x(i) -= U(i,j) * x(j);
    x(i) = x(i) / U(i,i);
  }
}


#endif
