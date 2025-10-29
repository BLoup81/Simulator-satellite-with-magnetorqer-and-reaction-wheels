function u = command_exemple2(state)

    % This allocation is due to commandOption = 'cardan'
    % theta 1 = state(1)
    % theta 2 = state(2)
    % theta 3 = state(3)

    % omega 1 = state(4)
    % omega 2 = state(5)
    % omega 3 = state(6)

    % wheel 1 = state(7)
    % wheel 2 = state(8)
    % wheel 3 = state(9)

    % K_theta 1 = state(10)
    % K_theta 2 = state(11)
    % K_theta 3 = state(12)

    % K_omega 1 = state(13)
    % K_omega 2 = state(14)
    % K_omega 3 = state(15)

    u = zeros(6,1);

    u(1) = -state(10)*state(1) - state(13)*state(4);
    u(2) = -state(11)*state(2) - state(14)*state(5);
    u(3) = -state(12)*state(3) - state(15)*state(6);
    u(4) = 0;
    u(5) = 0;
    u(6) = 0;

end