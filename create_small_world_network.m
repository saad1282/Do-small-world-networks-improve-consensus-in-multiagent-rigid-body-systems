%% Create small world network from a ring or lattice
function A_small_world = create_small_world_network(N,k,p)
% N = 100;
% k = 4;
% p = 0.2;
A = zeros(N,N); %adjacency matrix

%% create lattice first
for i = 1 : N
    for j = 1 : k
        ind1 = mod(i+j+N-1,N) + 1;
        ind2 = mod(i-j+N-1,N) + 1;
        A(i, ind1) = 1;
        A(i, ind2) = 1;
    end
end

%% create small world network
 A_small_world = A;

for j = 1 : k
    for i = 1 : N

        %select jth clockwise neighbor
        neigh_clock = mod(i+j+N-1,N) + 1;

        %random number which determines whether to rewire the edge
        rand_p = rand(1,1);
        if rand_p < p
            A_small_world(i,neigh_clock) = 0;
            A_small_world(neigh_clock,i) = 0;
        
            % new integer to connect node j to
            % new_ind = randi([1 N]);
            % while A_small_world(i, new_ind) == 1 || i == new_ind
            %     new_ind = randi([1 N]);
            % end
            possible_targets = find((1:N) ~= i & A_small_world(i,:) == 0 & (1:N)~= neigh_clock);
            new_ind = possible_targets(randi(length(possible_targets)));

            % add new node to small world adjacency matrix
            A_small_world(i,new_ind) = 1;
            A_small_world(new_ind,i) = 1;
        end
    end
end

% % Plot graphs
% figure;
% G1 = graph(A);
% subplot(1,2,1)
% plot(G1);
% title('Ring Lattice')
% 
% G2 = graph(A_small_world);
% subplot(1,2,2)
% plot(G2);
% title('Watts-Strogatz Small-World')
% A_small_world - A
% sum(sum(abs(A_small_world - A)))/2
% 
% save('A_100nodes_prob_p2.mat', "A_small_world") 

% %%
% k1 = 2;
% N1 = 8;
% for i = 1 : N1
%     for j = 1 : k1
%         ind1 = mod(i+j+N1-1,N1) + 1;
%         ind2 = mod(i-j+N1-1,N1) + 1;
%         A1(i, ind1) = 1;
%         A1(i, ind2) = 1;
%     end
% end
% L = diag(sum(A1,2)) - A1;
% k2 = 1;
% for i = 1 : N1
%     for j = 1 : k2
%         ind1 = mod(i+j+N1-1,N1) + 1;
%         ind2 = mod(i-j+N1-1,N1) + 1;
%         B(i, ind1) = 1;
%         B(i, ind2) = 1;
%     end
% end
% 
% C = A1 - B;
% L1 = diag(sum(B,2)) - B;
% L2 = diag(sum(C,2)) - C;
% 
% K2 = [0:0.01:1];
% ks = 0.309;
% 
% %eig(-cos(2*pi/N1) - K2(f)*((cos(6*pi/N1) - ks)*L1 - cos(4*pi/N1) - K2(f)*((cos(6*pi/N1) - ks)*L1 + cos(12*pi/N) - ks)*L2)
% eig(-cos(6*pi/N1)*L)
% eig(-cos(8*pi/N1)*L1 - cos(48*pi/N1)*L2)
% for f = 1 : length(K2)
%     %eig_val(f) = max(eig((-cos(2*pi/N1) - K2(f)*((cos(6*pi/N1) - ks)*L1 - cos(4*pi/N1) - K2(f)*((cos(6*pi/N1) - ks)*L1 + cos(12*pi/N) - ks)*L2));
%     eig_val(f) = max(eig((-cos(2*pi/N1) -K2(f)*(cos(6*pi/N1) -ks))*L1));
% end
% 
% plot(K2,eig_val)
% hold on;
% yline(0)