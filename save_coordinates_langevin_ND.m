step = 0.25;
index = 1;
n_dim = 2;
initial_partition = 1;
count = 1;
vid = VideoWriter('Langevin_projected.avi');
open(vid)
for n_dim = 1
    clear entropy_rate
    for partitions = 10.^(1)%10.^(initial_partition:step:4)
        X = [];
        V = [];
        mkdir(['coordinates_1.4\X\',num2str(partitions)])
        mkdir(['coordinates_1.4\V\',num2str(partitions)])
        [num2str(partitions),' ',num2str(n_dim)]
        partitions_idx = round((log(partitions)/log(10) - initial_partition)/step) + 1;
        x(1) = rand - 0.5;
        v = rand - 0.5;
        m = 1;
        gamma = 10;
        kbt = 0.8;
        dt = 0.02;
        tau = 0.8;
        interval = 0;
        side = (ceil(partitions^(1/n_dim)));
        n_states = side^n_dim;
        previous_state = floor((x(1) + 2 + 1e-15)*side/4) + 1;
        transition_matrix = zeros(n_states, n_states)*1e-15;
        eq_probability = zeros(n_states,1);
        x(2) = x(1);
        for i = 1:1000%10000000
            [x(1),v] = update_velocity(x(1),v,gamma,kbt,dt,1);
            X = [X, x(1)];
            V = [V, v];
            plot(X,V)
            %hold on 
            %plot(X(end),V(end),'o','MarkerFaceColor','r','LineStyle','-')
            %xlim([-2 2]); ylim([-3 3])
            %hold off
            set(gca,'FontSize',28)
            %xlabel('X'); ylabel('V')
            plot(1:length(X),X)
            hold on 
            plot(length(X),X(end),'o','MarkerFaceColor','r')
            hold off
            xlim([0 1000])
            ylim([-2 2])
            set(gca,'FontSize',28)
            drawnow
            f = getframe;
            writeVideo(vid,f)
            if (mod(i, 500000) == 0)
                save(['coordinates_1.4\X\',num2str(partitions),'\x_',num2str(count),'.mat'],'X');
                save(['coordinates_1.4\V\',num2str(partitions),'\v_',num2str(count),'.mat'],'V');
                X = [];
                V = [];
                count = count + 1;
            end
        end
    end
end
close(vid)