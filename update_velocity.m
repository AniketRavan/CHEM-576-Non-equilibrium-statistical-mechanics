function [x_update,v_update] = update_velocity(x,v,gamma,kbt,dt,nparticles)
    noise = normrnd(0,1,nparticles,1);
    v_dot = -gamma*v - 8*x.*(x.^2 - 1) - sqrt(2*gamma*kbt/dt)*noise;
    x_dot = v;
    v_update = v + v_dot*dt;
    x_update = x + x_dot*dt;
end