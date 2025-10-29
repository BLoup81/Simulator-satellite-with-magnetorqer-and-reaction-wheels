clear all;
clc;

% This file is an exemple to use this simulator. We use the configuration
% of a satellite with 3 magnetorquers (MGT) and 3 reaction wheels (RW).


%% Earth

earth = planet();
earth = earth.setGravitationalParameter(398600.44);
earth = earth.setRadius(6378);
earth = earth.setRotation(2*pi/84600);

%% Orbit

% Parameters
orb = orbit();
orb = orb.setRAAN(30);
orb = orb.setInclination(87);           % Pass close to the poles
orb = orb.setArgumentOfPeriapsis(45);
orb = orb.setMajorSemiAxis(6978);       % 600 km altitude
orb = orb.setEccentricity(0.001);       % Approximatly circular

% Initial conditions
orb = orb.setGHA(0);
orb = orb.setInitialTrueAnomaly(0);
orb = orb.setDateInitial([2025 10 02 17 24 00]);

% Attractor body
orb.attractorBody = earth;

%% Reaction wheels

% Wheel 1
wheel1 = reactionWheel();
wheel1 = wheel1.setAxis(1);
wheel1 = wheel1.setInertia(3.5e-5);
wheel1 = wheel1.setVelocityBound(300);
wheel1 = wheel1.setTorqueBound(0.005);

% Whee 1: initial condition
wheel1 = wheel1.setVelocity(10);

% Wheel 2
wheel2 = reactionWheel();
wheel2 = wheel2.setAxis(2);
wheel2 = wheel2.setInertia(3.5e-5);
wheel2 = wheel2.setVelocityBound(300);
wheel2 = wheel2.setTorqueBound(0.005);

% Whee 2: initial condition
wheel2 = wheel2.setVelocity(-20);

% Wheel 3
wheel3 = reactionWheel();
wheel3 = wheel3.setAxis(3);
wheel3 = wheel3.setInertia(3.5e-5);
wheel3 = wheel3.setVelocityBound(300);
wheel3 = wheel3.setTorqueBound(0.005);

% Whee 3: initial condition
wheel3 = wheel3.setVelocity(0);

%% Magnetorquers

% Magnetorquer 1
mgt1 = magnetorquer();
mgt1 = mgt1.setAxis(1);
mgt1 = mgt1.setMagneticMomentBound(0.0005);

% Magnetorquer 2
mgt2 = magnetorquer();
mgt2 = mgt2.setAxis(2);
mgt2 = mgt2.setMagneticMomentBound(0.0005);

% Magnetorquer 3
mgt3 = magnetorquer();
mgt3 = mgt3.setAxis(3);
mgt3 = mgt3.setMagneticMomentBound(0.0005);

%% Satellite

sat = satellite();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dimension
sat = sat.setLength(0.9);       % 90 cm
sat = sat.setWidth(0.75);       % 75 cm
sat = sat.setHeight(1);         % 1 m

% Mass
sat = sat.setMass(10);          % 10 kg

% Center of mass
sat = sat.computeCenterOfMass();

% Inertia
sat = sat.computeInertia();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial conditions
% Quaternion:
% We want to start from Cardan's angles [0.1 -0.003 0.05]
cardan = [0.1;-0.003;0.05];

% Conversion from cardan to quaternion
initialQuaternion = tools.euler2quater(cardan);

% Affectation of the initial attitude
sat = sat.setQuaternion(initialQuaternion);

% Angular velocity
initialAngularVelocity = [0;0;0];

% Affectation of the initial angular velocity
sat = sat.setAngularVelocity(initialAngularVelocity);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Actuators
% Reaction wheels
sat.wheels = [wheel1 wheel2 wheel3];

% Magnetorquers
sat.magnetorquers = [mgt1 mgt2 mgt3];

% Time
step = 0.01;        % Step time
tf = 20;            % Final time

%% Simulation 1

% Time
time = 0:step:tf;

% Command
command = @(state) [-state(2)-state(5);     % Command RW 1
                    -state(3)-state(6);     % Command RW 2
                    -state(4)-state(7);     % Command RW 3
                    -state(2)-state(5);     % Command MGT 1
                    -state(3)-state(6);     % Command MGT 2
                    -state(4)-state(7)];    % Command MGT 3

% Command option: 'quaternion' or 'cardan'
commandOption = 'quaternion';

% Frame of reference to compare the attitude of the satellite
frame = 'ECI';

[T,Y,U,M] = sat.ODE(time,command,commandOption,orb,frame);

figure(1);
tools.plotAttitudeQuaternion(T,Y);   % Plot quaternion

figure(2);
tools.plotCardan(T,Y);      % Plot Cardan's angles

figure(3);
tools.plotSatelliteVelocity(T,Y);       % Plot angular velocity of the

figure(4);
tools.plotWheelVelocity(T,Y,sat);   % Plot wheel velocity
% Remark: The wheel on the yaw axis go to its bound (300 rad/s)

figure(5);
tools.plotWheelTorque(T,U,sat);     % Plot generated torque by reaction wheels
% Remark: The wheel on the yaw axis has a saturated torque

figure(6);
tools.plotMagnetorquerTorque(T,U,sat);      % Plot generated torque by magnetorquers

figure(7);
tools.plotMagneticMoment(T,M);       % Plot magnetic moment

input("Tap entry to see the seconde simulation","s");
close all;

%% Simulation 2

% We use the same configuration but we use an adaptative command only on reaction wheels.
% This simulation show that we can add state to the main model.

additionState = [zeros(6,1)];%;zeros(3,1)];      % Addition of 2*6 states equal to 0 at the initial time
% [K_theta;K_omega]

% Time
time = 0:step:tf;

% Command option: 'quaternion' or 'cardan'
commandOption = 'cardan';

% Frame of reference to compare the attitude of the satellite
frame = 'ECI';

[T,Y,U,M] = sat.ODE(time,@command_exemple2,commandOption,orb,frame,additionState,@derivativeAdditionState_exemple2);

figure(1);
tools.plotAttitudeQuaternion(T,Y);   % Plot quaternion

figure(2);
tools.plotCardan(T,Y);      % Plot Cardan's angles

figure(3);
tools.plotSatelliteVelocity(T,Y);       % Plot angular velocity of the satellite

figure(4);
tools.plotWheelVelocity(T,Y,sat);   % Plot wheel velocity
% Remark: The wheel on the yaw axis go to its bound (300 rad/s)

figure(5);
tools.plotWheelTorque(T,U,sat);     % Plot generated torque by reaction wheels
% Remark: The wheel on the yaw axis has a saturated torque

figure(6);
tools.plotMagnetorquerTorque(T,U,sat);      % Plot generated torque by magnetorquers

figure(7);
tools.plotMagneticMoment(T,M);       % Plot magnetic moment

figure(8);
plot(T,Y(11:end,:));            % Plot additional state
xlabel('Time (s)');
ylabel('Corrector gain');
legend('K_{\theta 1}','K_{\theta 2}','K_{\theta 3}','K_{\omega 1}','K_{\omega 2}','K_{\omega 3}');

input("Tap entry to see the third simulation","s");
close all;

%% Simulation 3

% The control law in the simulation 2 do not respect the control law on the
% paper and the attitude do not converge. So, we add a function called
% 'projection_exemple2' to respect the control law.
%
% This simulation show that we can add input parameters to the derivative
% function of the additional states

additionState = [zeros(6,1)];%;zeros(3,1)];      % Addition of 2*6 states equal to 0 at the initial time
% [K_theta;K_omega]

% Time
tf = 350;
time = 0:step:tf;

% Command option: 'quaternion' or 'cardan'
commandOption = 'cardan';

% Frame of reference to compare the attitude of the satellite
frame = 'ECI';

K_min = [0.002 1];
K_max = [0.5 2];

[T,Y,U,M] = sat.ODE(time,@command_exemple2,commandOption,orb,frame,additionState,@derivativeAdditionState_2_exemple2,@projection_exemple2,step,K_max,K_min);

figure(1);
tools.plotAttitudeQuaternion(T,Y);   % Plot quaternion

figure(2);
tools.plotCardan(T,Y);      % Plot Cardan's angles

figure(3);
tools.plotSatelliteVelocity(T,Y);       % Plot angular velocity of the satellite

figure(4);
tools.plotWheelVelocity(T,Y,sat);   % Plot wheel velocity
% Remark: The wheel on the yaw axis go to its bound (300 rad/s)

figure(5);
tools.plotWheelTorque(T,U,sat);     % Plot generated torque by reaction wheels
% Remark: The wheel on the yaw axis has a saturated torque

figure(6);
tools.plotMagnetorquerTorque(T,U,sat);      % Plot generated torque by magnetorquers

figure(7);
tools.plotMagneticMoment(T,M);       % Plot magnetic moment

figure(8);
plot(T,Y(11:end,:));            % Plot additional state
xlabel('Time (s)');
ylabel('Corrector gain');
legend('K_{\theta 1}','K_{\theta 2}','K_{\theta 3}','K_{\omega 1}','K_{\omega 2}','K_{\omega 3}');