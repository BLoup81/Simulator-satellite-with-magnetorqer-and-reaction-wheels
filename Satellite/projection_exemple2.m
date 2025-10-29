function proj = projection_exemple2(x,dx,step,x_max,x_min)

    proj = x + step*dx;

    if proj > x_max
        proj = (x_max-x)/step;
    elseif proj < x_min
        proj = (x_min-x)/step;
    end
end

