%function for determing percent of time step present after tau time steps

function alph = alpha(tau, rise_time)
    if tau >= rise_time
        alph = 1.0;
    elseif tau <= rise_time - 1
        alph = (tau + .5) / rise_time;
    else
        alph = 1 - (rise_time - tau)*.5*(1 - tau / rise_time);
    end
    
    if alph < 0
        display([int2str(alph) ' error at tau = ' int2str(tau) ' and rt ' int2str(rise_time)])
    end
end
