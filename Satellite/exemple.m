clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is an exemple to use the simulator, %
% show the diferents functions to interact or set %
% the diferents parameters to the satellite. It   %
% is also possible get these parameters.          %
%                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%
% Satellite dimensions:                           %
% - Length, width and height                      %
% - Mass                                          %
% - Center of mass                                %
%                                                 %
% State of the model:                             %
% - Current/initial quaternion                    %
% - Current/initial angular velocity              %
%                                                 %
% Actuators:                                      %
% - Reaction wheels:                              %
%    - {0, 1, 2, 3}                               %
%    - Axis of the wheel {1, 2, 3}                %
%    - Inertia of the wheel                       %
%    - Current/initial angular velocity           %
%    - Velocity bound                             %
%    - Torque bound                               %
% - Magnetorquers:                                %
%    - {0, 1, 2, 3}                               %
%    - Axis of the magnetorquer {1, 2, 3}         %
%    - Magnetic moment bound                      %
%                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Orbit %%%%%%%%%%%%%%%%%%%%%%%
% Attractor body:                                 %
% - Gravitational parameter                       %
% - Rotation rate                                 %
% - Radius                                        %
%                                                 %
% Orbit:                                          %
% - Right Ascension of Ascendent Node             %
% - Inclination                                   %
% - Argument of the periapsis                     %
% - Major semi axis                               %
% - Eccentricity                                  %
% - Initial true anomaly                          %
% - Initial Greenwich Hour Angle                  %
% - Initial date                                  %
% - Attractor body                                %
%                                                 %
% Orbit functions:                                %
% - Mean motion                                   %
% - Semi latus rectus ()                          %
% - Initial eccentric anomaly ()                  %
% - Mean anomaly (time)                           %
% - Eccentric anomaly (time)                      %
% - True anomaly (time)                           %
% - Radius (time, true anomaly)                   %
% - ECI coordinates (time, true anomaly, radius)  %
% - Calcul GHA (time)                             %
% - Geocentric coordinate (time)                  %
%                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Earth magnetic field %%%%%%%%%%%%%%%%
%                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Conversion %%%%%%%%%%%%%%%%%%%%%%
% Conversion:                                     %
% - Quaternion to Cardan's angles                 %
% - Cardan's angles to quaternion                 %
% - Cardan's angles to Cardan's angles            %
% - ECEF frame to ECI frame                       %
% - ECI frame to LVLH frame                       %
% - Fframe to BRF frame                           %
% - NED frame to ECEF frame                       %
%                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                 %
% Created by:                                     %
% Loup BESSARD                                    %
%                                                 %
% At LAAS-CNRS, Toulouse, France                  %
% The 26/09/2025                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% The satellite
sat = satellite();

fprintf("==============================================================\n");
fprintf("The satellite has the following caracteristics:\n");

% Dimensions of the satellite
sat = sat.setLength(0.1);       % length = 0.1m
sat = sat.setWidth(0.1);        % width = 0.1m
sat = sat.setHeight(0.3);       % height = 0.3m
sat = sat.setMass(4);           % mass = 4kg

fprintf("Length:\t %.4f m\n",sat.getLength());
fprintf("Width:\t %.4f m\n",sat.getWidth());
fprintf("Height:\t %.4f m\n",sat.getHeight());
fprintf("Mass:\t %.4f kg\n",sat.getMass());

% Center of mass
centerX = 0.05;                 % Center of mass following X axis [0, length]
centerY = 0.05;                 % Center of mass following Y axis [0, width]
centerZ = 0.15;                 % Center of mass following Z axis [0, height]
sat = sat.setCenterOfMass([centerX;centerY;centerZ]);

centerOfMass = sat.getCenterOfMass();

fprintf("\nThe center of mass is located in:\n");
fprintf("Axis X: %.4f m\n",centerOfMass(1));
fprintf("Axis Y: %.4f m\n",centerOfMass(2));
fprintf("Axis Z: %.4f m\n",centerOfMass(3));

% Inertia matrix
inertia = [10 0 0;0 2 0; 0 0 35];   % Inertia matrix
sat = sat.setInertia(inertia);

fprintf("\nThe inertia matrix of the satellite is:\n%.4f %.4f %.4f\n%.4f %.4f %.4f\n%.4f %.4f %.4f\n",sat.getInertia());

% It is possible to compute the inertia matrix and the center of mass with the dimensions of the
% satellite
sat = sat.computeCenterOfMass();
sat = sat.computeInertia();

fprintf("\nThe real inertia matrix of the satellite is:\n%.4f %.4f %.4f\n%.4f %.4f %.4f\n%.4f %.4f %.4f\n",sat.getInertia());


% Set the initial conditions
quaternion = [1;0;0;0];         % The initial quaternion
sat = sat.setQuaternion(quaternion);

angularVelocity = [0.56;0;0];      % The initial angular velocity
sat = sat.setAngularVelocity(angularVelocity);

fprintf("\nInital conditions are:\n");
fprintf("Quaternion: %.4f %.4f %.4f %.4f\n",sat.getQuaternion());
fprintf("Angular velocity: %.4f %.4f %.4f %.4f rad/s\n",sat.getAngularVelocity());

%% Reaction wheels
wheel = reactionWheel();

fprintf("\n\n==============================================================\n");
fprintf("The reaction wheel has the following caracteristics:\n");

% Inertia coefficient
inertiaWheel = 35;          % inertia of the wheel = 35 kg.m^2
wheel = wheel.setInertia(inertiaWheel);

fprintf("Inertia:\t %.f kg.m^2\n",wheel.getInertia());

% Rotation axis
axis = 1;                       % Axis of rotation of the wheel {1,2,3}
wheel = wheel.setAxis(axis);

fprintf("Axis number: %.f\n", wheel.getAxis());

%% Affection of reaction wheels to the satellite
sat.wheels = wheel;

% It is possible to change value of wheel's attributs from the satellite
sat.wheels = sat.wheels.setInertia(3.5e-5);

fprintf("New inertia: %.6f kg.m^2\n",sat.wheels.getInertia());


wheel2 = reactionWheel();
wheel2 = wheel2.setInertia(87);
wheel2 = wheel2.setAxis(2);

% Allocation of two reaction wheels to the satellite
sat.wheels = [wheel wheel2];

% Addition of one reaction wheel to the satellite
wheel3 = reactionWheel();
sat = sat.addWheel(wheel3);

sat.wheels(3) = sat.wheels(3).setAxis(3);
sat.wheels(3) = sat.wheels(3).setInertia(65);

fprintf("\nThere are %d reaction wheels:\n",length(sat.wheels));
fprintf("Axis:\t\t %d\t %d\t %d\n",sat.wheels(1).getAxis(),sat.wheels(2).getAxis(),sat.wheels(3).getAxis());
fprintf("Inertia:\t %.f\t %.f\t %.f\n",sat.wheels(1).getInertia(),sat.wheels(2).getInertia(),sat.wheels(3).getInertia());

% Reaction wheels are saturated
wheel = wheel.setVelocityBound(293);
wheel = wheel.setTorqueBound(0.5);

fprintf("\nBound saturation of the raction wheel:\n");
fprintf("Velocity bound:\t %.f rad/s\n",wheel.getVelocityBound());
fprintf("Torque bound:\t %.1f N.m\n",wheel.getTorqueBound());

%% Magnetorquer

mgt1 = magnetorquer();

fprintf("\n\n==============================================================\n");
fprintf("The magnetorquer has the following caracteristics:\n");

% Axis
axisMGT1 = 1;
mgt1 = mgt1.setAxis(axisMGT1);

fprintf("Axis number: %.f\n",mgt1.getAxis());

% Magnetic moment
saturationMGT1 = 0.55;              % Magnetic saturation 0.55 A.m^2 = 0.55 N.m/T
mgt1 = mgt1.setMagneticMomentBound(saturationMGT1);

fprintf("Magnetic saturation: %.2f A.m^2\n",mgt1.getMagneticMomentBound());

%% Affectation of magnetorquers to the satellite

mgt2 = magnetorquer();

% It is possible to set a vector of magnetorquers to the satellite
sat.magnetorquers = [mgt1 mgt2];

% Modify the attribut of the magnetorquer by the satellite
sat.magnetorquers(2) = sat.magnetorquers(2).setAxis(2);
sat.magnetorquers(2) = sat.magnetorquers(2).setMagneticMomentBound(0.005);

% It is also possible to add a magnetorquer
mgt3 = magnetorquer();
mgt3 = mgt3.setAxis(3);
mgt3 = mgt3.setMagneticMomentBound(0.32);

sat = sat.addMagnetorquer(mgt3);

fprintf("\nThere are %d magnetorquers:\n",length(sat.magnetorquers));
fprintf("Axis:\t\t %d\t\t %d\t\t %d\n",sat.magnetorquers(1).getAxis(),sat.magnetorquers(2).getAxis(),sat.magnetorquers(3).getAxis());
fprintf("Saturation:\t %.3f\t %.3f\t %.3f\n",sat.magnetorquers(1).getMagneticMomentBound(),sat.magnetorquers(2).getMagneticMomentBound(),sat.magnetorquers(3).getMagneticMomentBound());

%% Attractor body of the orbit, exemple of the Earth

earth = planet();

% Affectation of the parameters of the Earth
earth = earth.setGravitationalParameter(398600.44);         % Gravitational parameter of the Earth  [km^3/s^2]
earth = earth.setRotation(2*pi/84600);                      % Rotation rate of the Earth            [rad/s]
earth = earth.setRadius(6378);                              % Radius of the Earth                   [km]

%% Orbit

fprintf("\n\n==============================================================\n");
fprintf("The orbit has the following caracteristics\n");

orb = orbit();

% Affectation of the parameter's orbit
orb = orb.setRAAN(30);                              % Right Ascension of the Ascendent Node     [degree]
orb = orb.setInclination(90);                       % Inclination                               [degree]
orb = orb.setArgumentOfPeriapsis(20);               % Argument of periapsis (perigee on Earth)  [degree]
orb = orb.setMajorSemiAxis(6978);                   % Major semi axis                           [km]
orb = orb.setEccentricity(0.001);                   % Eccentricity
orb = orb.setGHA(0);                                % Greenwich Hour Angle                      [degree]
orb = orb.setInitialTrueAnomaly(45);                % Initial true anomaly                      [degree]
orb = orb.setDateInitial([2025 10 02 14 11 40]);    % Initial date: 2nd October 2025, 14h 11min 40s
orb.attractorBody = earth;                          % Affectation of the attractor body of the orbit, exemple of the Earth

fprintf("RAAN:\t\t\t\t\t %.f degree, %f rad\n",orb.getRAAN('deg'),orb.getRAAN('rad'));      % It is possible to get the attribut and specify the unit 'deg' or 'rad'
fprintf("i:\t\t\t\t\t\t %.f degree, %f rad\n",orb.getInclination(),orb.getInclination('rad'));  % If no unit is set, by default the attribut is in degree
fprintf("Argument of perigee:\t %.f degree, %f rad\n",orb.getArgumentOfPeriapsis(),orb.getArgumentOfPeriapsis('rad'));
fprintf("Major semi axis:\t\t %f km\n",orb.getMajorSemiAxis());
fprintf("Eccentricity:\t\t\t %f\n",orb.getEccentricity());
fprintf("GHA:\t\t\t\t\t %.f degree, %f rad\n",orb.getGHA('deg'),orb.getGHA('rad'));
fprintf("Initial true anomaly:\t %.f degree, %f rad\n",orb.getInitialTrueAnomaly(),orb.getInitialTrueAnomaly('rad'));
fprintf("Initial date;\t\t\t year: %d, month: %d, day: %d, %d h %d min %d s\n",orb.getDateInitial());

fprintf("\nParameters of the attractor body:\n");
fprintf("Gravitational parameter:\t %f km^3/s^2\n",orb.attractorBody.getGravitationalParameter());
fprintf("Radius:\t\t\t\t\t\t %.f km\n",orb.attractorBody.getRadius());
fprintf("Rotation rate:\t\t\t\t %f rad/s\n",orb.attractorBody.getRotation());

fprintf("\n\n==============================================================\n");
time = [0 10 100];
fprintf("Orbit functions for the time vector [%.f %.f %.f]:\n",time);
fprintf("Mean motion ():\t\t\t\t\t\t\t %f\n",orb.meanMotion());
fprintf("Semi latus rectus ():\t\t\t\t\t %f\n",orb.semi_latus_rectus());
fprintf("Mean anomaly (time):\t\t\t\t\t [%.3f %.3f %.3f]\n",orb.meanAnomaly(time));
fprintf("Eccentric anomaly (time):\t\t\t\t [%.3f %.3f %.3f]\n",orb.eccentricAnomaly(time));
trueAnomaly = orb.trueAnomaly(time);
fprintf("True anomaly (time):\t\t\t\t\t [%.3f %.3f %.3f]\n",trueAnomaly);
radius = orb.radius(time,orb.trueAnomaly(time));
fprintf("Raius (time, true anomaly):\t\t\t\t [%.3f %.3f %.3f]\n",radius);
fprintf("ECI coordinates (time, true anomaly, radius):\n");
fprintf("\t\t\t\t\t\t\t\t\t\t X: %.3f, Y: %.3f, Z: %.3f\n",orb.coordinateECI(time,trueAnomaly,radius));
fprintf("Calcul GHA (time):\t\t\t\t\t\t [%.3f %.3f %.3f]\n",orb.calculGHA(time));
fprintf("Geocentric coordinates (time):\n");
[latitude,longitude,height] = orb.geocentric(time);       % Compute the latitude and longitude in degree and the height in km
fprintf("\t\t\t\t\t\t\t\t\t\t Latitude = [%.3f %.3f %.3f]\n",latitude);
fprintf("\t\t\t\t\t\t\t\t\t\t Longitude = [%.3f %.3f %.3f]\n",longitude);
fprintf("\t\t\t\t\t\t\t\t\t\t Height = [%.3f %.3f %.3f]\n",height);

%% Earth magnetic field

fprintf("\n\n==============================================================\n");
fprintf("Earth magnetic field:\n");

addpath("earth_magnetic_field\m_IGRF-main\");

% To compute the earh magnetic field coordinates, we need the geocentric
% coordinates of the point we want know the magnetic field and the date.
time = 0;           % To compute geocentric coordinates
date = orb.getDateInitial();    % Get the date of the initial time of the simulation
[latitude,longitude,height] = orb.geocentric(time);     % Compute geocentric coordinates

magneticField = igrf(date,latitude,longitude,height,'geocentric');      % Magnetic field in the NED frame

fprintf("Magnetic field at the initial time:\t X: %.f, Y: %.f, Z: %.f\n",magneticField);
fprintf("WARNING: The magnetic field is in the North-East-Down frame (NED)\n");

%% Conversions

fprintf("\n\n==============================================================\n");
fprintf("Conversions:\n");

fprintf("Quaternion / Euler's angles");
% Quaternion -> Euler's angles
euler = [pi;0;pi/2];
quaternion = tools.euler2quater(euler);
fprintf("From Euler's angles [%.3f %.3f %.3f] it is possible get the associated quaternion [%.3f %.3f %.3f %.3f]\n",euler,quaternion);

% Euler's angles -> Quaternion
fprintf("From the quaternion [%.3f %.3f %.3f %.3f] it is possible go back to Euler's angles [%.3f %.3f %.3f]\n",quaternion,tools.quater2euler(quaternion));

% To use the magnetic field, we have to convert it in the body reference
% frame (BRF) of the satellite
fprintf("\nNED to BRF frame:\n");

magneticField = [1;1;1];
fprintf("Suppose that we have the magnetic field [%f %f %f] in the NED\n",magneticField);
latitude = 45;      % Exemple latitude
longitude = -90;    % Exemple longitude
time = 10000;          % Exemple time

% If the satellite's attitude is compared to the ECEF frame:
ecefCoordinates = tools.conversionCoordinates(magneticField,[latitude;longitude],orb,time,'ECEF');
% We can use: ECEF_coordinates = tools.ned2ecef(latitude,longitude,NED_coordinates)
fprintf("Magnetic field in the ECEF:\t [%.3f %.3f %.3f]\n",ecefCoordinates);

% If the satellite's attitude is compared to the ECI frame:
eciCoordinates = tools.conversionCoordinates(magneticField,[latitude;longitude],orb,time,'ECI');
% We can use: ECI_coordinates = tools.ecef2eci(ECEF_coordinates,earth.getRotation(),time);
fprintf("Magnetic field in the ECI:\t [%.3f %.3f %.3f]\n",eciCoordinates);

% If the satellite's attitude is compared to the LVLH frame:
lvlhCoordinates = tools.conversionCoordinates(magneticField,[latitude;longitude],orb,time,'LVLH');
% We can use: LVLH_coordinates = tools.eci2lvlh(ECI_coordinates,arguments_of_perigee,RAAN,inclination,true_anomaly);
fprintf("Magnetic field in the LVLH:\t [%.3f %.3f %.3f]\n",lvlhCoordinates);

% Conversion of the magnetic field from the ECI frame to the BRF
% We use the quaternion = [0 0.707 0.707 0] (before)
brfCoordinates = tools.frame2brf(eciCoordinates,quaternion);
fprintf("Magneic fielf in the BRF:\t [%.3f %.3f %.3f]\n",brfCoordinates);