%% Main driver code
%Jonathan Engle
%Project 3
%Due 11/13/23
clear clc clear all 
format default
tic
%% Step 1
%Generate a random Symmetric positive definite matrix
global n 
%n=20
%n=200
n=3;
%trial used 3 
%n=3;

%Create symetric positive definite matrix
lil=tril(randi([1,1],[n,n]));
%Ensure that we are diagonally dominant
Spd=lil*transpose(lil)+(eye(n)*n);
%Checks imposed to verify that we have done this correctly
checkspd=rcond(Spd)>eps;
Anew=randi([1,1],[n,n]);
A2=Anew*transpose(Anew)+(eye(n)*n^4);

if checkspd==0;
    error('Singular Matrix, try again');
end
%Additional Check
if det(Spd)==0;
    error('Singular Matrix, try again')
end

%% Step 2
%%Step 2 Given a matrix A, we need to generate a vector b basaed
%on a chosen solution x. I.E for a given A choose a random x and
%generate b via A*x

%Generate random vector x
rng("default")
realx = rand(n,1);
%compute b
global b;
b=Spd*realx;

%% St;eepest Descent method
A=A2;
%[200 .01; 200 200];
eig(A);
%No pre conditioner
nP=eye(n);
%initializing our global tollerance
global tol
tol=10^(-6);
% Jacobi Pre conditioner i.e diags of A
JP=sqrt(diag(A)).*eye(n);

%initializing our max itteration number
%this ensures that it wont go off to infinity i.e finite run time
global intermax;
intermax=100000;

% Gauss Seidel pre conditioner i.e lower triangular
%matrix with a_ii/sqrt(a_ii) for i>=j
GSP= zeros(n,n);
for i = 1:n
    for j = 1:i
        GSP(i, j) = A(i, j) / sqrt(A(i, i));
    end
end
SGS=GSP*transpose(GSP);
%Preconditioner used
Precon=GSP;
%Running Steepest descent method
[xtild1,r,k,p,r_new,z] = psteep(A, Precon, b, intermax, tol);
xtild1=psteep(A,JP,b, intermax,tol)
k
% %Running conjugate gradient method
% %[conx, iter, conz_new,conr_new] = ConGrad(A, nP, b, intermax, tol);
% [xcon, iter, z_new,p,rcon] = PreConGrad2(A, Precon, b, intermax, tol);
% iter;
%Forward substitution for our given matrix
% y=forwardSubstitution(Precon,rcon);
%Backward substitution for our given matrix
% zed=backwardSubstitution(transpose(Precon),y);
%zed is used from now on
%% Error analysis
%Evation of x for Steepest Descent method
error_x= norm(xtild1-realx,2)/norm(realx,2);
%Evaluation of x for Conjugate Gradient method
error_conx= norm(xcon-realx,2)/norm(realx,2);

%Evaluation of r 
%solving rk Steepest Descent is done in the function
r_errorSteep=norm(z,2)/norm(b,2)
%
r_errorCon=norm(zed,2)/norm(b,2)
%% Plotting solutions
%Jacobi plot
%plot(abs(rcon))
iter;
%GS plot
% figure
% hold on
% plot(abs(r_new))
% plot(abs(rcon))
% legend('Gonjugate Gradient','Steepest Descent')
% hold off
%plot(abs(r));
% abs(zed);

%% Correctness test
%Setting up our big A
Atest=[5 7 6 5;
    7 10 8 7;
    6 8 10 9;
    5 7 9 10];
btest=[57;79;88;86];
%initiallizing our preconditioners
%No preconditioner
np_test=eye(4);
%Jacobi preconditioner
Jp_test=diag(Atest).*eye(4);
%Gauss preconditioner
GSPtest= zeros(4,4);
for i = 1:4
    for j = 1:i
        GSPtest(i, j) = Atest(i, j) / sqrt(Atest(i, i));
    end
end

GSPtest1=transpose(GSPtest)*GSPtest;
% Implement con grad
%No preconditioner
[x4,r4,ktestnp2,p4,r_new4,z4] = psteep(Atest, np_test, btest, intermax, tol);
%Jacobi Preconditioner
[x5,r5,ktestjp2,p5,r_new5,z5] = psteep(Atest, Jp_test, btest, intermax, tol);
%Gauss Seidel preconditioner
[x6,r6,ktestgs2,p6,r_new6,z6] = psteep(Atest, GSPtest1, btest, intermax, tol);
%Solving for z
y=forwardSubstitution(GSPtest,r_new6);
%Backward substitution for our given matrix
zedtest=backwardSubstitution(transpose(GSPtest),y);
r_errortestCon=norm(zedtest,2)/norm(btest,2);
%Output
ktestnp2;
ktestjp2;
ktestgs2;
eig(Atest);


% % Conjugate gradient method
% % Implement PSD
% [xcon4, iter4, z_new4,p4,rcon4] = PreConGrad2(Atest, np_test, btest, intermax, tol);
% iter4;
% %Jacobi
% [xcon5, iter5, z_new5,p5,rcon5] = PreConGrad2(Atest, Jp_test, btest, intermax, tol);
% iter5;
% %GS preconditioner
% [xcon6, iter6, z_new6,p6,rcon6] = PreConGrad2(Atest, GSPtest1, btest, intermax, tol);
% iter6;
% %Solving for z
% y1=forwardSubstitution(GSPtest,rcon6);
% %Backward substitution for our given matrix
% zedtest1=backwardSubstitution(transpose(GSPtest1),y1);
% r_errortestCon1=norm(zedtest1,2)/norm(btest,2);
% 

%No preconditioner
[x1,r1,ktestnp,p1,r_new1,z1] = psteep(Atest, np_test, btest, intermax, tol);
%Jacobi Preconditioner
[x2,r2,ktestjp,p2,r_new2,z2] = psteep(Atest, Jp_test, btest, intermax, tol);
%Gauss Seidel preconditioner
[x3,r3,ktestgs,p3,r_new3,z3] = psteep(Atest, GSPtest1, btest, intermax, tol);
%Solving for z
y=forwardSubstitution(GSPtest,r_new3);
%Backward substitution for our given matrix
zedtest=backwardSubstitution(transpose(GSPtest1),y);
r_errortestCon=norm(zedtest,2)/norm(btest,2);
%zed is used from now on
%Output
ktestnp;
ktestjp;
ktestgs;
% graphs
% figure
% hold on
% plot(abs(r1))
% plot(abs(r2))
% plot(abs(r3))
% hold off
%Report steps

%% Extra credit

%% Theorem 4.2
% If A is strictly diagonally dominant by rows the Jacobi and Gauss Seidel
%Are convergent

%Creating a diagonally dominant matrix by rows

%% Theorem 4.3 
% If A and 2D-A are positive definite matricies
%Then Jacobi method is convergent and \rho(B_j)=\|B_J\|_A=\|B_J\|_D
rng("default")
extrax = rand(n,1);
%Testing if 2D-A is SPD
extra=eye(n,n);


%% Theorem 4.5
%If A is spd 

toc