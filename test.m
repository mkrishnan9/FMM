N =10;
x=rand(N,1);
y=rand(N,1);
u=rand(N,1);
N_trials = [128,256,512,1024,2048,4096];
T_a = zeros(5,1); T_b = zeros(5,1);
T_c = zeros(5,1); T_d = zeros(5,1);
T_e = zeros(5,1);
for i=1:6
    N = N_trials(i);
    x=rand(N,1);
    y=rand(N,1);
    u=rand(N,1);
    tic
    soln_a = MatVecProd_a(x,y,u);
    T_a(i) = toc;
    tic
    soln_b = MatVecProd_b(x,y,u);
    T_b(i) = toc;
    tic
    soln_c = MatVecProd_c(x,y,u);
    T_c(i) = toc;
    tic
    soln_d = MatVecProd_d(x,y,u);
    T_d(i) = toc;
    tic
    soln_e = MatVecProd_e(x,y,u);
    T_e(i) = toc;
end
figure;
hold on;
loglog(N_trials,T_a,'-o')
loglog(N_trials,T_b,'-o')
loglog(N_trials,T_c,'-o')
loglog(N_trials,T_d,'-o')
loglog(N_trials,T_e,'-o')
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Size of Matrix')
ylabel('Time')
title('Time vs Size of Matrix')
legend('Colm Major', 'Row major', 'Colm Vector', 'Row Vector','Factorization','Location','northwest')

function v =  MatVecProd_a(x,y,u)
N = size(x,1);
A = zeros(N,N);
for j=1:N
    for i=1:N
        A(i,j)=(x(i)-y(j))^2;
    end
end
v = zeros(N,1);
for j=1:N
    for i=1:N
        v(i) = v(i)+ A(i,j)*u(j);
    end
end
end

function v =  MatVecProd_b(x,y,u)
N = size(x,1);
A = zeros(N,N);
for j=1:N
    for i=1:N
        A(i,j)=(x(i)-y(j))^2;
    end
end
v = zeros(N,1);
for i=1:N
    for j=1:N
        v(i) = v(i)+ A(i,j)*u(j);
    end
end
end

function v =  MatVecProd_c(x,y,u)
N = size(x,1);
A = zeros(N,N);
for j=1:N
    for i=1:N
        A(i,j)=(x(i)-y(j))^2;
    end
end
v = zeros(N,1);
for j = 1:N
    v = v + A(:,j)*u(j);
end
end

function v =  MatVecProd_d(x,y,u)
N = size(x,1);
A = zeros(N,N);
for j=1:N
    for i=1:N
        A(i,j)=(x(i)-y(j))^2;
    end
end
v = zeros(N,1);
for i = 1:N
    v(i) =  A(i,:)*u;
end
end

function v =  MatVecProd_e(x,y,u)
N = size(x,1);
v = zeros(N,1);
beta=0;
gamma = 0;
delta =0;
for i = 1:N
    beta = beta + u(i);
    gamma = gamma + u(i)*y(i)^2;
    delta = delta + u(i)*y(i);
end
for j = 1:N
    v(j) = beta*x(j)^2 + gamma - 2*x(j)*delta;
end
end