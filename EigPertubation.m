%B.W.
%6/13

%Nondegenerate eigenvalue pertubation method - following technique used in
%3rd edition of Griffiths Quantum Mechanics.
%%
clear
clc

%Generate random matrix
A0 = 10*(rand(6,6)-0.5);

%Find eigenvalues of unperturbed matrix A0
[V,eigs] = eig(A0);

%Perturb A0 with del_A
del_A = 0.01*(rand(6,6)-0.5);
A = A0 + del_A;

%Check to make sure no eigenvalues are shared between eigenvectors
if length(unique(diag(eigs))) ~= length(diag(eigs))
    fprintf('Degenerate Pertubation Theory Needed!')
end

%Find correction factors for the perturbed eigenvalues:
%First order
eigs_C1 = zeros(length(eigs(:,1)),1);
%Second order
eigs_C2 = zeros(length(eigs(:,1)),1);
for i = 1:length(eigs(:,1))
    eigs_C1(i) = conj(V(:,i).')*del_A*V(:,i);
    
    for j = 1:length(eigs(:,1))
        if j == i
            continue
        end
        eigs_C2(i) = eigs_C2(i) + (conj(V(:,j).')*del_A*V(:,i))*(conj(V(:,i).')*del_A*V(:,j))/(eigs(i,i)-eigs(j,j));
    end
end

eigs_perturbed = diag(eigs) + eigs_C1 + eigs_C2;

%Compare with actual eigenvalues
eigs_actual  = eig(A);

%Difference between estimated and actual eigenvalues
format long g
difference = eigs_perturbed - eigs_actual














