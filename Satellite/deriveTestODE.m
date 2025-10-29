%{
function dx = deriveTestODE(time,x)

    % x = [ q0 q1 q2 q3  w1 w2 w3  ww  t  h  c1 c2 c3 c4 ]'
    % où c1..c4 = anciennement "command" et h = anciennement "th"

    % Décomposition des états (mêmes que toi)
    q = x(1:4);      % quaternion
    w = x(5:7);      % vitesse angulaire
    ww = x(8);       % vitesse roue
    t = x(9);        % état additionnel t
    h = x(10);       % état additionnel h = th
    command = x(11:14); % ancien global "command"

    % --- Modèle physique (inchangé) ---

    % Paramètres inertiels
    I = [10 0 0;
         0 5 0;
         0 0 15];

    Iw = 3.5e-5;

    % Champ magnétique (fixe pour test)
    B_brf = [25805;1315;13150]; % en BRF (puisque tools supprimé pour le test)

    % Quaternion dérivée
    dx = zeros(size(x));
    dx(1) = 0.5*(-q(2)*w(1) - q(3)*w(2) - q(4)*w(3));
    dx(2) = 0.5*(q(1)*w(1) - q(4)*w(2) + q(3)*w(3));
    dx(3) = 0.5*(q(4)*w(1) + q(1)*w(2) - q(2)*w(3));
    dx(4) = 0.5*(-q(3)*w(1) + q(2)*w(2) + q(1)*w(3));

    % Mise à jour de command (exactement comme toi)
    command(2:4) = -cross(B_brf,cross(B_brf,command(2:4)))/(norm(B_brf)^2);

    % Dynamique angulaire
    dx(5) = (command(1) + command(2) - w(2)*w(3)*(I(3,3)-I(2,2)))/I(1,1);
    dx(6) = (command(3) - w(1)*w(3)*(I(1,1)-I(3,3)))/I(2,2);
    dx(7) = (command(4) - w(1)*w(2)*(I(2,2)-I(1,1)))/I(3,3);

    % Roue
    dx(8) = -command(1)/Iw;

    % Etats additionnels (inchangés)
    dx(9) = 1;
    dx(10) = w(1) - t + 1;  % correspond à "w(1)-th"

    % Mise à jour interne des commandes (tu la stockais en global)
    new_command = [-w(1); w(2); w(3); ww];
    dx(11:14) = (new_command - command);  % dérivée = transition vers la prochaine valeur

end
%}

function dx = deriveTestODE(time,x)

%{
    earth = planet();
    earth = earth.setGravitationalParameter(398600.44);
    earth = earth.setRotation(2*pi/84600);
    earth = earth.setRadius(6378);

    orb = orbit();
    orb = orb.setRAAN(30);
    orb = orb.setInclination(87);
    orb = orb.setArgumentOfPeriapsis(20);
    orb = orb.setInitialTrueAnomaly(0);
    orb = orb.setMajorSemiAxis(6978);
    orb = orb.setEccentricity(0.001);
    orb = orb.setGHA(0);
    orb.attractorBody = earth;

    [lat,lon,alt] = orb.geocentric(time);

    B_local = igrf([2025 09 26 10 09 30],lat,lon,alt,'geocentric')';

    %B_local = [25805;1315;13150];

    B_eci = tools.conversionCoordinates(B_local,[lat;lon],orb,time,'ECI');

    B_brf = tools.frame2brf(B_eci,x(1:4));

    dx = zeros(10,1);
%}

    dx = zeros(25,1);

    q = x(1:4);
    w = x(5:7);
    ww = x(8:10);

    I = [10 0 0;
         0 5 0;
         0 0 15];
    
    Iw = [0.01 0 0;
          0 0.01 0;
          0 0 0.01];

    

    B_local = [1;2;3];
    B_brf = tools.frame2brf(B_local,q);

    dx(1) = 0.5*(-q(2)*w(1) - q(3)*w(2) - q(4)*w(3));
    dx(2) = 0.5*(q(1)*w(1) - q(4)*w(2) + q(3)*w(3));
    dx(3) = 0.5*(q(4)*w(1) + q(1)*w(2) - q(2)*w(3));
    dx(4) = 0.5*(-q(3)*w(1) + q(2)*w(2) + q(1)*w(3));    

    %trueCom = [command(1);0;0;0];

    %trueCom(2:4) = -cross(B_brf,cross(B_brf,[command(2);0;command(3)]))/(norm(B_brf)^2);

    command = [0.01*w(1);0.01*w(2);0.01*w(3);0.5];
    M_com = [0;command(4);0];
    M = cross(B_brf,M_com)/(norm(B_brf)^2);
    U_mgt = cross(M,B_brf);
    U_rw = [command(1);command(2);command(3)];
    %command = -cross(B_brf,cross(B_brf,command))/(norm(B_brf)^2);

    dx(5) = (U_rw(1) + U_mgt(1) - w(2)*w(3)*(I(3,3)-I(2,2)))/I(1,1);
    dx(6) = (U_rw(2) + U_mgt(2) - w(1)*w(3)*(I(1,1)-I(3,3)))/I(2,2);
    dx(7) = (U_rw(3) + U_mgt(3) - w(1)*w(2)*(I(2,2)-I(1,1)))/I(3,3);

    dx(8:10) = (-U_rw - cross(w,Iw*ww))/0.01;

    %dx(9) = 1;
    %dx(10) = w(1) - th;

    dx(11) = 1;
    dx(12) = w(1) + w(2) + w(3);

    dx(13:16) = command;
    dx(17:19) = M;
    dx(20:22) = U_rw;
    dx(23:25) = U_mgt;

end
%}