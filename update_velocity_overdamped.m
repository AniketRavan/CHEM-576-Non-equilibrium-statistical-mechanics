function [x_update] = update_velocity_overdamped(x,gamma,kbt,dt,nparticles)
    noise = normrnd(0,1,nparticles,1);
    x_dot = (- 8*x.*(x.^2 - 1) - sqrt(2*gamma*kbt/dt)*noise)/gamma;
    x_update = x + x_dot*dt;
end