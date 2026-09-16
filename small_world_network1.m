%% convert quaternion to rotation matrix
tic
clear all
clc

n = 30; % number of nodes
k_neigh = 3; % order of ring lattice
p1 = 0.23; % rewiring probability

num_networks = 10; % number of small world network realizations
monte_carlo_sim_num = 1000; % number of initial conditions for each network
A = cell(num_networks,1); % initialize the adjacency matrix of all realizations as a single cell

T = 300; % simulation time
dt = 0.01; % time step size
t = 0 : dt : T; % simulation times

% adjacency matrices, uses the wattsstrogatz algorithm from matlab
for j = 1 : num_networks
    h = WattsStrogatz(n,k_neigh,p1);
    A1 = zeros(n,n);
    for i = 1 : size(h.Edges)
        edge = h.Edges(i,:);
        A1(edge.EndNodes(1),edge.EndNodes(2)) = 1;
        A1(edge.EndNodes(2),edge.EndNodes(1)) = 1;
    end
    A{j} = A1;
end

% total simulations to be done, necessary as we will parallelize across all
% simulations
total_simulations = num_networks * monte_carlo_sim_num;

number_of_failures = 0;
tol = 1e-2;
check_consensus = zeros(total_simulations,1); % stores consensus as success = 1, and non consensus as failure = 0
winding_num = zeros(total_simulations,1); % finds winding number for each simulation

[network_indices, mc_indices] = ndgrid(1:num_networks, 1:monte_carlo_sim_num); % grid of networks indices and mc indices

network_indices = network_indices(:);
mc_indices = mc_indices(:);

% this part displays how many simulations are completed
progress = 0;
q = parallel.pool.DataQueue;
afterEach(q,@updateProgress);

% parallel for loop
parfor sim_idx = 1:total_simulations
    network_idx = network_indices(sim_idx);
    mc_idx = mc_indices(sim_idx);
    A_current = A{network_idx};

    % Run one simulation
    [consensus, wind_num] = run_single_simulation(A_current, n, T);

    % Only index the output with the parfor loop variable
    check_consensus(sim_idx) = consensus;
    winding_num(sim_idx) = wind_num;
    %send(q,struct('increment',1,'total',total_simulations));
    send(q,1);
end

% reshape total simulations results by network realization num and monte
% carlo num
check_consensus = reshape(check_consensus, monte_carlo_sim_num, num_networks).';

function updateProgress(~)
    persistent count

    if isempty(count)
        count = 0;
    end
    count = count + 1;
    fprintf('\rCompleted %d simulations',count);
end
%% 
%
% Simulate SO3 kinematics using Runga Kutta Method and classify the results
% as consensus or non-consensus and find corresponding winding num

function [consensus, wind_num] = run_single_simulation(A, n, t)
R = zeros(3,3,n);
noedge = sum(sum(A,1))/2;
for i = 1 : n
        rand_vec = randomVectorWithNormBound(3,pi);
        R(:,:,i) = expm(skew(rand_vec));
end

    %% ODE 45
    ICs = zeros(9*n,1);
    for i=1:n
        ICs((i-1)*9+1:i*9) = reshape(R(:,:,i),[9,1]);
    end
    tolerance = 1e-9;
    options = odeset('RelTol',tolerance,'AbsTol',tolerance);
    [time, Rw] = ode45(@(time,Rw) f_Rw(Rw,A,n), [0,t], ICs, options);
    len=size(Rw,1);
    %% unwrap ODE soln
    praDiff=zeros(noedge,len);
    pairDiff = zeros(n,1);
    for k = 1:len
        ct = 1;
        for i=1:n
         Ri = reshape(Rw(k,1+(i-1)*9:9*i),[3,3]);
             for j=i+1:n
                  if A(i,j)~=0
                     Rj = reshape(Rw(k,1+(j-1)*9:9*j),[3,3]);
                     praDiff(ct,k) = phi(Ri*Rj');  
                      ct = ct+1;

                   if (k == len)
                       if (i == 1) && (j == n)
                           pairDiff(n) = phi(Ri*Rj');
                       elseif (i < n) && (j == i + 1)
                           pairDiff(i) = phi(Ri*Rj');
                       end
                   end
                  end

             end
        end
    end
    check_consensus_val = max(max(praDiff(:,2*floor(len/3):end)));
    if check_consensus_val < 1e-2
        consensus = 1;
        wind_num = 0;
    else
        consensus = 0;
        wind_num = round(sum(pairDiff)/360);
    end
    end
       
%%
% generates random prv vector for initial conditions
function v = randomVectorWithNormBound(n, R)
% Generates a single random vector of dimension n
% with a 2-norm less than or equal to R.

% 1. Generate an unconstrained sample with normally distributed elements
X = randn(n, 1);

% 2. Calculate the norm of the vector
x_norm = norm(X);

% 3. Generate a random magnitude (r) to ensure uniform volume distribution
% r needs to be scaled by rand()^(1/n)
r = R * rand()^(1/n);

% 4. Scale the vector to have the desired norm
v = X * (r / x_norm);
end

% Runga Kutta Method for SO3 simulation
function dRw=f_Rw(Rw,A,n)
    dRw=zeros(9*n,1);
    for i=1:n
        dR=zeros(3,3);
        Ri = reshape(Rw(1+(i-1)*9:9*i),[3,3]);
        det(Ri);
        wi = zeros(3,1);
        % taui = taui - ki*wi;
        for j=1:n
            if A(i,j)==0 || i==j
                continue;
            end
            Rj = reshape(Rw(1+(j-1)*9:9*j),[3,3]);
            Rij = Rj'*Ri;
            wi = wi + control(Rij);
        end
        % kin & dyn eqns
        dR = Ri*skew(wi);

        % formatting for ode45
        dRw(1+(i-1)*9:9*i) = reshape(dR,[9,1]);
    end
    end

    
  % pairwise control protocol
function vec_ = control(R)
kp = 10;
vec_ = -1/2*kp*unskew(R- R');
end

% find PRA of agents' attitudes
function val = phi(DCM)
if trace(DCM) > 2.992
    skew_mat = (DCM - DCM')/2;
    val = rad2deg(sqrt(skew_mat(2,3)^2 + skew_mat(3,1)^2 + skew_mat(1,2)^2));
else
    val = acosd((trace(DCM)-1)/2);
end
end

% cross operator
function [vecX] = skew(vec)
vec1 = vec(1); vec2 = vec(2); vec3 = vec(3);
vecX = [0 -vec3 vec2; vec3 0 -vec1; -vec2 vec1 0];
end

% vectorize operator which undoes skew operator
function [vec] = unskew(mat)
vec(1,1) = mat(3,2);
vec(2,1) = mat(1,3);
vec(3,1) = mat(2,1);
end