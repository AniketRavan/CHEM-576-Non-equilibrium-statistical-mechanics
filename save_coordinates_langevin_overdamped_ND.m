step = 0.25;
index = 1;
n_dim = 2;
initial_partition = 1;
count = 1;

for n_dim = 1
    clear entropy_rate
    for partitions = 10.^(initial_partition:step:initial_partition)
        X = [];
        mkdir(['coordinates_overdamped_0.8\X\',num2str(partitions)])
        mkdir(['coordinates_overdamped_0.8\V\',num2str(partitions)])
        [num2str(partitions),' ',num2str(n_dim)]
        partitions_idx = round((log(partitions)/log(10) - initial_partition)/step) + 1;
        x(1) = rand - 0.5;
        v = rand - 0.5;
        m = 1;
        gamma = 10;
        kbt = 0.8;
        dt = 0.001;
        tau = 0.8;
        interval = 0;
        side = (ceil(partitions^(1/n_dim)));
        n_states = side^n_dim;
        previous_state = floor((x(1) + 2 + 1e-15)*side/4) + 1;
        transition_matrix = zeros(n_states, n_states)*1e-15;
        eq_probability = zeros(n_states,1);
        x(2) = x(1);
        for i = 1:10000000
            x(1) = update_velocity_overdamped(x(1),gamma,kbt,dt,1);
            X = [X, x(1)];
            if (mod(i, 500000) == 0)
                save(['coordinates_overdamped_0.8\X\',num2str(partitions),'\x_',num2str(count),'.mat'],'X');
                X = [];
                count = count + 1;
            end
        end
    end
end