function dx = derivativeAdditionState_exemple2(state)

    % In the derivative function for addition state, the commandOption do nothing
    % Quaternion 0 = state(1)
    % Quaternion 1 = state(2)
    % Quaternion 2 = state(3)
    % Quaternion 3 = state(4)

    % omega 1 = state(5)
    % omega 2 = state(6)
    % omega 3 = state(7)

    % wheel 1 = state(8)
    % wheel 2 = state(9)
    % wheel 3 = state(10)

    % K_theta 1 = state(11)
    % K_theta 2 = state(12)
    % K_theta 3 = state(13)

    % K_omega 1 = state(14)
    % K_omega 2 = state(15)
    % K_omega 3 = state(16)

    dx = zeros(6,1);

    % Definition of constant parameters
    G_theta = [32.92 32.92 32.92];
    G_omega = [798.4 798.4 798.4];

    c_theta = [1.0966 3.0462 1.0966];
    c_omega = [1.0966 3.0462 1.0966];

    K_p = [0.1 0.1 0.1];
    K_d = [2 2 2];

    % Conversion of the quaternion to Cardan's angle

    cardan = tools.quater2euler(state(1:4));

    % Derivative equation
    % dK_theta
    dx(1) = -G_theta(1)*(cardan(1)^2) - c_theta(1)*(state(11) - K_p(1));
    dx(2) = -G_theta(2)*(cardan(2)^2) - c_theta(2)*(state(12) - K_p(2));
    dx(3) = -G_theta(3)*(cardan(3)^2) - c_theta(3)*(state(13) - K_p(3));

    % dK_omega
    dx(4) = -G_omega(1)*(cardan(1)^2) - c_omega(1)*(state(14) - K_d(1));
    dx(5) = -G_omega(2)*(cardan(2)^2) - c_omega(2)*(state(15) - K_d(2));
    dx(6) = -G_omega(3)*(cardan(3)^2) - c_omega(3)*(state(16) - K_d(3));

end