% Figure 4-15 shows a three-node truss element of length L and a constant cross-section area A. It is made of a material of Young’s modulus E and density . The truss is subjected to a uniformly distributed force b. 
% a) Derive the stiffness matrix for the element.
% b) Write down the expression for the element mass matrix, and obtain  in terms of L, E, , and A. 
% c) Derive the external force vector. 
  
%% (a) Derive the stiffness matrix for the element
clear all;clc
syms L x 

% The natural coordinate  is defined as 
x_c = L/2;
l   = L;
xi  = 2 * (x-x_c)/l;

% Using the quadradic one-dimensional element with three nodes
N1 = -1/2 * xi *(1-xi);
N2 = (1+xi) * (1-xi);
N3 = 1/2 * xi *(1+xi);

N1 = simplify(N1) 
N2 = simplify(N2) 
N3 = simplify(N3) 

% The strain matrix  has the form which yeilds
N = [N1 N2 N3];

% Derivatives of N1, N2, N3 with respect to xi
dN1_dx = diff(N1);
dN1_dx = simplify(dN1_dx);
dN2_dx = diff(N2);
dN2_dx = simplify(dN2_dx);
dN3_dx = diff(N3);
dN3_dx = simplify(dN3_dx);

% Combine the results into a matrix
dN_dx = [dN1_dx, dN2_dx, dN3_dx];

% Display the result
disp(dN_dx);

% Now that the stain matrix has been obtained, using the stiffness matrix for the truss element can be obtained which yeilds
B   = dN_dx;
B_T = B.';
BB  = B_T * B;

% Size of BB
[n,m] = size(BB);

% Initialize a symbolic matrix for the integrated values
BB_integrated = sym(zeros(n, m));

% Integrate each element of BB
for i = 1:n
    for j = 1:m
        BB_integrated(i,j) = int(BB(i,j), x, 0, L); % From 0 to L
    end
end

BB_integrated = simplify(BB_integrated); 
%disp(BB_integrated);

%% b) Write down the expression for the element mass matrix, and obtain  in terms of L, E, , and A. 
syms A E
K_e=BB_integrated*A*E

% Integrating yeilds
m_e = [N1*N1 N1*N2 N1*N3;
       N2*N1 N2*N2 N2*N3;
       N3*N1 N3*N2 N3*N3;];

% Size of m_e
[n,m] = size(m_e);

% Initialize a symbolic matrix for the integrated values
m_e_integrated = sym(zeros(n, m));

% Integrate each element of m_e
for i = 1:n
    for j = 1:m
        m_e_integrated(i,j) = int(m_e(i,j), x, 0, L); % From 0 to L
    end
end

m_e_integrated = simplify(m_e_integrated); 
%disp(m_e_integrated);
syms A rho
m_e = m_e_integrated * A * rho
[row, col] = size(m_e);

for i = 1:row
    for j = 1:col
        varName = ['m' num2str(i) num2str(j)];
        eval([varName ' = m_e(i,j);']);
    end
end

disp(m11);

%% (c) Derive the external force vector
% The total nodal force vector for the three-node truss element can be obtained using the following equation


N_T = N.';

% Size of N_T
[n,m] = size(N_T);

% Initialize a symbolic matrix for the integrated values
N_T_integrated = sym(zeros(n, m));

% Integrate each element of N_T
for i = 1:n
    for j = 1:m
        N_T_integrated(i,j) = int(N_T(i,j), x, 0, L); % From 0 to L
    end
end

N_T_integrated = simplify(N_T_integrated);
%disp(N_T_integrated)
syms b A
f_e = b*A*N_T_integrated
