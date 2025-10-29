% Test of the class satellite.m

% To run tests:
%       runtests('satelliteTest')

%% Test 0: Constructor

sat = satellite();

assert(isequal(sat.getInertia(),eye(3,3)));

%% Test 1.1: Getter/Setter length of the satellite

% Input data:
length = 3.65;

% Expected output data:
length_expect = 3.65;

% Test:
sat = satellite();

sat = sat.setLength(length);
length_test = sat.getLength();

assert(isequal(length_expect,length_test));

%% Test 1.2: Length positive exception

sat = satellite();

try
    sat.setLength(-1);
    error('Test failed: length must be positive');
catch ME
    assert(strcmp(ME.identifier,'satellite:BadInput'),'Length of satellite must be positive');
end

%% Test 2.1: Getter/Setter width of the satellite

% Input data:
width = 5.24;

% Expected output data:
width_expect = 5.24;

% Test:
sat = satellite();

sat = sat.setWidth(width);
width_test = sat.getWidth();

assert(isequal(width_expect,width_test));

%% Test 2.2: Width positive exception

sat = satellite();

try
    sat.setWidth(-1);
    error('Test failed: width must be positive');
catch ME
    assert(strcmp(ME.identifier,'satellite:BadInput'),'Width of satellite must be positive');
end

%% Test 3.1: Getter/Setter height of the satellite

% Input data:
height = 8.35;

% Expected output data:
height_expect = 8.35;

% Test:
sat = satellite();

sat = sat.setHeight(height);
height_test = sat.getHeight();

assert(isequal(height_expect,height_test));

%% Test 3.2: Height positive exception

sat = satellite();

try
    sat.setHeight(-1);
    error('Test failed: height must be positive');
catch ME
    assert(strcmp(ME.identifier,'satellite:BadInput'),'Height of satellite must be positive');
end

%% Test 4.1: Getter/Setter mass of the satellite

% Input data:
mass = 5.32;

% Expected output data:
mass_expect = 5.32;

% Test:
sat = satellite();

sat = sat.setMass(mass);
mass_test = sat.getMass();

assert(isequal(mass_expect,mass_test));

%% Test 4.2: Mass positive exception

sat = satellite();

try
    sat.setMass(-1);
    error('Test failed: mass must be positive');
catch ME
    assert(strcmp(ME.identifier,'satellite:BadInput'),'Mass of satellite must be positive');
end

%% Test 5.1: Getter/Setter center of mass

% Input data:
centerOfMass = [3.6; 5.6; 1.26];

% Expected output data:
centerOfMass_expect = [3.6; 5.6; 1.26];

% Test:
sat = satellite();

sat = sat.setCenterOfMass(centerOfMass);
centerOfMass_test = sat.getCenterOfMass();

assert(isequal(centerOfMass_expect, centerOfMass_test));

%% Test 5.2: Center of mass dimension exception

sat = satellite();

try
    sat.setCenterOfMass(1);
    error('Test failed: center of mass vector must be 3 by 1');
catch ME
    assert(strcmp(ME.identifier,'satellite:BadCenterOfMass'),'The center of mass vector must be 3 by 1');
end

%% Test 6.1: Getter/Setter inertia of satellite body

% Input data:
inertia = [5.6 3.15 0;
           6.36 25.6 3.17;
           98.4 6.45 79.1];

% Expected output data:
inertia_expect = [5.6 3.15 0;
                  6.36 25.6 3.17;
                  98.4 6.45 79.1];

% Test:
sat = satellite();

sat = sat.setInertia(inertia);
inertia_test = sat.getInertia();

assert(isequal(inertia_expect, inertia_test));

%% Test 6.2: Inertia definite positive exception

sat = satellite();

inertia2 = [-1 0 0;
           0 -2 0;
           0 0 3];

try
    sat.setInertia(inertia2);
    error('Test failed: inertia matrix could be negative definite');
catch ME
    assert(strcmp(ME.identifier,'satellite:BadInertia'),'The eigenvalues of inertia matrix must be positive');
end

%% Test 6.3: Inertia dimension exception

sat = satellite();

try
    sat.setInertia([1 2 3]);
    error('Test failed: inertia matrix could be different of 3 by 3');
catch ME
    assert(strcmp(ME.identifier,'satellite:BadInertia'),'The inertia matrix must be 3 by 3');
end


%% Test 7.1: Getter/Setter quaternion

% Input data:
quaternion = [0.7;-0.32;0.37;-0.52];

% Expected output data:
quaternion_expect = [0.7;-0.32;0.37;-0.52];

% Test:
sat = satellite();

sat = sat.setQuaternion(quaternion);
quaternion_test = sat.getQuaternion();

assert(isequal(quaternion_expect,quaternion_test));

%% Test 7.2: Queternion dimension exception

sat = satellite();

try
    sat.setQuaternion(1);
    error('Test failed: quaternion must be 4 by 1');
catch ME
    assert(strcmp(ME.identifier,'satellite:BadInput'),'Quaternion must be 4 by 1');
end

%% Test 7.3: Quaternion norm exception

sat = satellite();

try
    sat.setQuaternion([0.707;0.5;0;0]);
    error('Test failed: norm of quaternion must be equal to 1');
catch ME
    assert(strcmp(ME.identifier,'satellite:BadInput'),'Norm of quaternion must be equal to 1');
end

%% Test 8.1: Getter/Setter angular velocity

% Input data:
angularVelocity = [-0.036;1.28;54];

% Expected output data:
angularVelocity_expect = [-0.036;1.28;54];

% Test:
sat = satellite();

sat = sat.setAngularVelocity(angularVelocity);
angularVelocity_test = sat.getAngularVelocity();

assert(isequal(angularVelocity_expect,angularVelocity_test));

%% Test 8.2: Angular velocity dimension exception

sat = satellite();

try
    sat.setAngularVelocity(1);
    error('Test failed: angular velocity vector must be 3 by 1');
catch ME
    assert(strcmp(ME.identifier,'satellite:BadInput'),'Angular velocity vector must be 3 by 1');
end

%% Test 9: Getter/Setter and addWheel for several reaction Wheels

% Input data:
inertiaWheel_1 = 54;
axisWheel_1 = 1;

inertiaWheel_2 = 32;
axisWheel_2 = 2;

inertiaWheel_3 = 19;
axisWheel_3 = 3;

% Expected output data:
inertiaWheel_1_expect = 54;
axisWheel_1_expect = 1;

inertiaWheel_2_expect = 32;
axisWheel_2_expect = 2;

inertiaWheel_3_expect = 19;
axisWheel_3_expect = 3;

% Test:
sat = satellite();

wheel_1 = reactionWheel();
wheel_1 = wheel_1.setInertia(inertiaWheel_1);
wheel_1 = wheel_1.setAxis(axisWheel_1);

wheel_2 = reactionWheel();
wheel_2 = wheel_2.setInertia(inertiaWheel_2);
wheel_2 = wheel_2.setAxis(axisWheel_2);

wheel_3 = reactionWheel();
wheel_3 = wheel_3.setInertia(inertiaWheel_3);
wheel_3 = wheel_3.setAxis(axisWheel_3);

sat.wheels = [wheel_1 wheel_2];
sat = sat.addWheel(wheel_3);

wheels_test = sat.wheels;

assert(isequal(inertiaWheel_1_expect,wheels_test(1).getInertia()));
assert(isequal(axisWheel_1_expect,wheels_test(1).getAxis()));

assert(isequal(inertiaWheel_2_expect,wheels_test(2).getInertia()));
assert(isequal(axisWheel_2_expect,wheels_test(2).getAxis()));

assert(isequal(inertiaWheel_3_expect,wheels_test(3).getInertia()));
assert(isequal(axisWheel_3_expect,wheels_test(3).getAxis()));

%% Test 10: Getter/Setter and addMagnetorquer for several magnetorquers

% Input data:
axisMGT_1 = 1;
saturationMGT_1 = 328;

axisMGT_2 = 2;
saturationMGT_2 = 0.068;

axisMGT_3 = 3;
saturationMGT_3 = 54.6;

% Expected output data:
axisMGT_1_expect = 1;
saturationMGT_1_expect = 328;

axisMGT_2_expect = 2;
saturationMGT_2_expect = 0.068;

axisMGT_3_expect = 3;
saturationMGT_3_expect = 54.6;

% Test:
sat = satellite();

mgt1 = magnetorquer();
mgt1 = mgt1.setAxis(axisMGT_1);
mgt1 = mgt1.setMagneticMomentBound(saturationMGT_1);

mgt2 = magnetorquer();
mgt2 = mgt2.setAxis(axisMGT_2);
mgt2 = mgt2.setMagneticMomentBound(saturationMGT_2);

mgt3 = magnetorquer();
mgt3 = mgt3.setAxis(axisMGT_3);
mgt3 = mgt3.setMagneticMomentBound(saturationMGT_3);

sat.magnetorquers = [mgt1 mgt2];
sat = sat.addMagnetorquer(mgt3);

assert(isequal(axisMGT_1_expect,sat.magnetorquers(1).getAxis()));
assert(isequal(saturationMGT_1_expect,sat.magnetorquers(1).getMagneticMomentBound()));

assert(isequal(axisMGT_2_expect,sat.magnetorquers(2).getAxis()));
assert(isequal(saturationMGT_2_expect,sat.magnetorquers(2).getMagneticMomentBound()));

assert(isequal(axisMGT_3_expect,sat.magnetorquers(3).getAxis()));
assert(isequal(saturationMGT_3_expect,sat.magnetorquers(3).getMagneticMomentBound()));

%% Test 11.1: Compute the center of mass

% Input data:
length = 2;
width = 4;
height = 3.5;

% Expected output data:
centerOfMass_expect = [1;2;1.75];

% Test:
sat = satellite();

sat = sat.setLength(length);
sat = sat.setWidth(width);
sat = sat.setHeight(height);

sat = sat.computeCenterOfMass();

centerOfMass_test = sat.getCenterOfMass();

assert(isequal(centerOfMass_expect(1),centerOfMass_test(1)));
assert(isequal(centerOfMass_expect(2),centerOfMass_test(2)));
assert(isequal(centerOfMass_expect(3),centerOfMass_test(3)));

%% Test 11.2: Compute the center of mass, existing of dimensions

% Input data:
length = 2;
height = 3.5;

% Test:
sat = satellite();

sat = sat.setLength(length);
sat = sat.setHeight(height);

try
    sat.computeCenterOfMass();
    error('Test failed: Dimensions must be defined to compute the center of masse');
catch ME
    assert(strcmp(ME.identifier,'satellite:MissingParameter'),'Dimensions of the satellite must be defined to compute the center of mass');
end


%% Test 12.1: Compute the inertia matrix

% Input data:
length = 2;
width = 4;
height = 3;
mass = 4;

centerOfMass = [length/4;width/5;height/1.5];

% Expected output data:
inertia_expect = [1132/75 12/5 -1;
                  12/5 19/3 -12/5;
                  -1 -12/5 1007/75];

% Test:
sat = satellite();

sat = sat.setLength(length);
sat = sat.setWidth(width);
sat = sat.setHeight(height);
sat = sat.setMass(mass);
sat = sat.setCenterOfMass(centerOfMass);
sat = sat.computeInertia();

inertia_test = sat.getInertia();

for row=1:3
    for col=1:3
        assert(abs(inertia_expect(row,col) - inertia_test(row,col)) <= 1e-4);
    end
end

%% Test 12.2: Compute inertia, existing of mass

% Input data:
length = 2;
width = 4;
height = 3;

centerOfMass = [length/4;width/5;height/1.5];

% Test:
sat = satellite();

sat = sat.setLength(length);
sat = sat.setWidth(width);
sat = sat.setHeight(height);
sat = sat.setCenterOfMass(centerOfMass);

try
    sat.computeInertia();
    error('Test failed: Mass must be defined to compute the inertia matrix');
catch ME
    assert(strcmp(ME.identifier,'satellite:MissingParameter'),'Mass must be defined to compute the inertia matrix');
end

%% Test 12.3: Compute inertia, existing of dimensions

% Input data:
length = 2;
height = 3;
mass = 4;

centerOfMass = [length/4;3;height/1.5];

% Test:
sat = satellite();

sat = sat.setLength(length);
sat = sat.setHeight(height);
sat = sat.setMass(mass);
sat = sat.setCenterOfMass(centerOfMass);

try
    sat.computeInertia();
    error('Test failed: Dimensions must be defined to compute the inertia matrix');
catch ME
    assert(strcmp(ME.identifier,'satellite:MissingParameter'),'Dimensions must be defined to compute the inertia matrix');
end


%% Test 12.4: Compute inertia, existing of center of mass

% Input data:
length = 2;
width = 4;
height = 3;
mass = 4;

% Test:
sat = satellite();

sat = sat.setLength(length);
sat = sat.setWidth(width);
sat = sat.setHeight(height);
sat = sat.setMass(mass);

try
    sat.computeInertia();
    error('Test failed: Center of mass must be defined to compute the inertia matrix');
catch ME
    assert(strcmp(ME.identifier,'satellite:MissingParameter'),'Center of mass must be defined to compute the inertia matrix');
end

%% Test 13: Derivative of the quaternion

% Input data:
quaternion = [0.66;0.22;-0.58;sqrt(449)/50];
angularVelocity = [2;-1;3];

% Expected output data:
dq_expect = [-1.1457;0.0019;-0.2362;1.46];

% Test:
sat = satellite();

dq_test = sat.quaternionDerivative([quaternion;angularVelocity]);

for i=1:4
    assert(abs(dq_expect(i) - dq_test(i)) <= 1e-4);
end

%% Test 14: Derivative of the angular velocity of the satellite

% Input data:
command = [1;2.5;-3];
angularVelocity = [0.5;2;-1];
inertia = [2 0 0;
           0 3 0;
           0 0 5];

% Expected output data:
dw_expect = [2.5;1/3;-4/5];

% Test:
sat = satellite();

sat = sat.setInertia(inertia);

dw_test = sat.angularVelocityDerivative([zeros(4,1);angularVelocity],command);

assert(abs(dw_expect(1) - dw_test(1)) <= 1e-4);
assert(abs(dw_expect(2) - dw_test(2)) <= 1e-4);
assert(abs(dw_expect(3) - dw_test(3)) <= 1e-4);

%% Test 15.1: Rearrange 3 reaction wheels

% Input data:
wheel1 = reactionWheel();
wheel1 = wheel1.setAxis(1);

wheel2 = reactionWheel();
wheel2 = wheel2.setAxis(2);

wheel3 = reactionWheel();
wheel3 = wheel3.setAxis(3);

sat = satellite();

sat.wheels = [wheel3 wheel1 wheel2];

% Expected output data:
axis_expect = [1 2 3];

% Test:
sat = sat.rearrangeReactionWheel();

assert(isequal(sat.wheels(1).getAxis(),axis_expect(1)));
assert(isequal(sat.wheels(2).getAxis(),axis_expect(2)));
assert(isequal(sat.wheels(3).getAxis(),axis_expect(3)));

%% Test 15.2: Rearrange 2 reaction wheels

% Input data:
wheel2 = reactionWheel();
wheel2 = wheel2.setAxis(2);

wheel3 = reactionWheel();
wheel3 = wheel3.setAxis(3);

sat = satellite();

sat.wheels = [wheel3 wheel2];

% Expected output data:
axis_expect = [2 3];

% Test:
sat = sat.rearrangeReactionWheel();

assert(isequal(sat.wheels(1).getAxis(),axis_expect(1)));
assert(isequal(sat.wheels(2).getAxis(),axis_expect(2)));

%% Test 16.1: Derivative of the wheel velocity, 1 wheel

% Input data:
command = 1;

wheel = reactionWheel();
wheel = wheel.setAxis(2);
wheel = wheel.setInertia(0.5);

% Expected output data:
velocityDerivative_expect = -2;

% Test:
sat = satellite();

sat.wheels = wheel;
velocityDerivative_test = sat.wheelVelocityDerivative(ones(7,1),command);

assert(isequal(velocityDerivative_expect,velocityDerivative_test));

%% Test 16.2: Derivative of the wheel velocity, 2 wheels

% Input data:
command = [1;-1];

wheel1 = reactionWheel();
wheel1 = wheel1.setAxis(2);
wheel1 = wheel1.setInertia(0.5);

wheel2 = reactionWheel();
wheel2 = wheel2.setAxis(3);
wheel2 = wheel2.setInertia(0.2);

% Expected output data:
velocityDerivative_expect = [-2;5];

% Test:
sat = satellite();

sat.wheels = [wheel1 wheel2];
velocityDerivative_test = sat.wheelVelocityDerivative(ones(9,1),command);

assert(isequal(velocityDerivative_expect(1),velocityDerivative_test(1)));
assert(isequal(velocityDerivative_expect(2),velocityDerivative_test(2)));

%% Test 16.3: Derivative of the wheel velocity, 3 wheels

% Input data:
command = [1;-1;1];
angularVelocity = [3;1.5;-2];           % Angular velocity of the satellite

wheel1 = reactionWheel();
wheel1 = wheel1.setAxis(1);
wheel1 = wheel1.setInertia(0.5);
wheel1 = wheel1.setVelocity(20);

wheel2 = reactionWheel();
wheel2 = wheel2.setAxis(2);
wheel2 = wheel2.setInertia(0.2);
wheel2 = wheel2.setVelocity(-58);

wheel3 = reactionWheel();
wheel3 = wheel3.setAxis(3);
wheel3 = wheel3.setInertia(0.25);
wheel3 = wheel3.setVelocity(5);

state = [zeros(4,1);angularVelocity;wheel1.getVelocity();wheel2.getVelocity();wheel3.getVelocity()];

% Expected output data:
velocityDerivative_expect = [40.65;123.75;195.2];

% Test:
sat = satellite();

sat = sat.setAngularVelocity(angularVelocity);
sat.wheels = [wheel1 wheel2 wheel3];
velocityDerivative_test = sat.wheelVelocityDerivative(state,command);

assert(abs(velocityDerivative_expect(1) - velocityDerivative_test(1)) <= 1e-4);
assert(abs(velocityDerivative_expect(2) - velocityDerivative_test(2)) <= 1e-4);
assert(abs(velocityDerivative_expect(3) - velocityDerivative_test(3)) <= 1e-4);

%% Test 17.1: Rearrange 3 magnetorquers

% Input data:
mgt1 = magnetorquer();
mgt1 = mgt1.setAxis(1);

mgt2 = magnetorquer();
mgt2 = mgt2.setAxis(2);

mgt3 = magnetorquer();
mgt3 = mgt3.setAxis(3);

sat = satellite();

sat.magnetorquers = [mgt3 mgt1 mgt2];

% Expected output data:
axis_expect = [1 2 3];

% Test:
sat = sat.rearrangeMagnetorquers();

assert(isequal(sat.magnetorquers(1).getAxis(),axis_expect(1)));
assert(isequal(sat.magnetorquers(2).getAxis(),axis_expect(2)));
assert(isequal(sat.magnetorquers(3).getAxis(),axis_expect(3)));

%% Test 18.1: RK4, 0 RW + 3 MGT

% Input data:
quaternion = [0.66;0.22;-0.58;sqrt(449)/50];
angularVelocity = [2;-1;3];
step = 0.1;
command = @(state) [3;-1;0.5];
inertia = [10 0 0;0 5 0;0 0 15];
derivativeAdditionState = @(state) 0;
affectationCommand = eye(3,3);
magneticField = [1;2;3];

state = [quaternion;angularVelocity];

mgt = magnetorquer();

sat = satellite();
sat.magnetorquers = [mgt mgt mgt];
for i=1:3
    sat.magnetorquers(i) = sat.magnetorquers(i).setAxis(i);
end
sat = sat.setInertia(inertia);

% Expected output data:
quaternion_expect = [0.542759414523714;0.213597435654307;-0.580848354900013;0.567805901727814];
angularVelocity_expect = [2.213325257815706;-0.374083438307173;2.961834051966752];
state_expect = [quaternion_expect;angularVelocity_expect];

u_compute_expect = 10*[0.3;-0.1;0.05];
u_effective_expect = 10*[0.067070320786257;-0.029666081618935;0.157471173139955];
m_compute_expect = 10*[-0.014249519011756;-0.043856144830180;-0.002215175589825];
m_effective_expect = 10*[-0.014249519011756;-0.043856144830180;-0.002215175589825];

% Test:
[state_test,u_compute_test,u_effective_test,m_compute_test,m_effective_test] = sat.RK4(state,command,affectationCommand,"quaternion",step,magneticField,derivativeAdditionState);

for i=1:7
    assert(abs(state_expect(i) - state_test(i)) <= 1e-4);
end

for i=1:3
    assert(abs(u_compute_test(i) - u_compute_expect(i)) <= 1e-3);
    assert(abs(u_effective_test(i) - u_effective_expect(i)) <= 1e-3);
    assert(abs(m_compute_test(i) - m_compute_expect(i)) <= 1e-3);
    assert(abs(m_effective_test(i) - m_effective_expect(i)) <= 1e-3);
end

%% Test 18.2: RK4, 1 RW + 2 MGT

% Input data:
quaternion = [0.66;0.22;-0.58;sqrt(449)/50];
angularVelocity = [2;-1;3];
wheelVelocity = 0;
step = 0.1;
command = @(state) [3;-1;0.5];
inertia = [10 0 0;0 5 0;0 0 15];
derivativeAdditionState = @(state) 0;
affectationCommand = [0 1 0;1 0 0;0 0 1;1 0 0];
magneticField = [1;2;3];

state = [quaternion;angularVelocity;wheelVelocity];

wheel = reactionWheel();
wheel = wheel.setInertia(0.01);
wheel = wheel.setAxis(2);

mgt = magnetorquer();

sat = satellite();
sat = sat.setInertia(inertia);

sat.magnetorquers = [mgt mgt];
sat.magnetorquers(1) = sat.magnetorquers(1).setAxis(1);
sat.magnetorquers(2) = sat.magnetorquers(2).setAxis(3);

sat.wheels = wheel;

% Expected output data:
quaternion_expect = [0.543890536838764;0.212802294870440;-0.580028256771579;0.567860448156813];
angularVelocity_expect = [2.196390311063344;-0.316594851477514;2.954170588067389];
wheelVelocity_expect = -30;
state_expect = [quaternion_expect;angularVelocity_expect;wheelVelocity_expect];

u_compute_expect = 10*[0.3;-0.1;0.05];
u_effective_expect = 10*[0.3;-0.006496108858198;0.007011710572744];
m_compute_expect = 10*[-0.003507739472126;-0.000931673539393;-0.007015478944252];
m_effective_expect = 10*[-0.003507739472126;-0.000931673539393;-0.007015478944252];

% Test:
[state_test,u_compute_test,u_effective_test,m_compute_test,m_effective_test] = sat.RK4(state,command,affectationCommand,"quaternion",step,magneticField,derivativeAdditionState);

assert(abs(state_expect(1) - state_test(1)) <= 1e-3);
assert(abs(state_expect(2) - state_test(2)) <= 1e-3);
assert(abs(state_expect(3) - state_test(3)) <= 1e-3);
assert(abs(state_expect(4) - state_test(4)) <= 1e-3);
assert(abs(state_expect(5) - state_test(5)) <= 1e-3);
assert(abs(state_expect(6) - state_test(6)) <= 1e-2);
assert(abs(state_expect(7) - state_test(7)) <= 1e-2);
assert(abs(state_expect(8) - state_test(8)) <= 1e-2);

for i=1:3
    assert(abs(u_compute_test(i) - u_compute_expect(i)) <= 1e-3);
    assert(abs(u_effective_test(i) - u_effective_expect(i)) <= 1e-3);
    assert(abs(m_compute_test(i) - m_compute_expect(i)) <= 1e-3);
    assert(abs(m_effective_test(i) - m_effective_expect(i)) <= 1e-3);
end

%% Test 18.3: RK4, 2 RW + 2 MGT

% Input data:
quaternion = [0.66;0.22;-0.58;sqrt(449)/50];
angularVelocity = [2;-1;3];
wheelVelocity = [5;20];
step = 0.1;
command = @(state) [0.01*state(5);0.01*state(6);-1;0.5];
inertia = [10 0 0;0 5 0;0 0 15];
derivativeAdditionState = @(state) 0;
affectationCommand = [1 0 1 0;0 0 0 1;0 1 0 0;1 0 0 0;0 1 0 0];
magneticField = [1;2;3];

state = [quaternion;angularVelocity;wheelVelocity];

sat = satellite();
sat =sat.setInertia(inertia);

wheel = reactionWheel();
wheel = wheel.setInertia(0.01);

sat.wheels = [wheel wheel];
sat.wheels(1) = sat.wheels(1).setAxis(1);
sat.wheels(2) = sat.wheels(2).setAxis(3);

mgt = magnetorquer();

sat.magnetorquers = [mgt mgt];
sat.magnetorquers(1) = sat.magnetorquers(1).setAxis(1);
sat.magnetorquers(2) = sat.magnetorquers(2).setAxis(2);


% Expected output data:
quaternion_expect = [0.543181694627273;0.213517017074100;-0.580709486900449;0.567574936075765];
angularVelocity_expect = [2.203211868553483;-0.366207325417992;2.949152438871169];
wheelVelocity_expect = [4.788245352767515;20.068661095484917];

state_expect = [quaternion_expect;angularVelocity_expect;wheelVelocity_expect];

u_compute_expect = 10*[0.002117546472325;-0.000686610954849;-0.1;0.05];
u_effective_expect = 10*[0.002117546472325;-0.000686610954849;-0.012737324323312;0.023503001226750];
m_compute_expect = 10*[0.005366991485885;0.010733982971771;0.004627364356335];
m_effective_expect = 10*[0.005366991485885;0.010733982971771;0.004627364356335];

% Test:
[state_test,u_compute_test,u_effective_test,m_compute_test,m_effective_test] = sat.RK4(state,command,affectationCommand,"quaternion",step,magneticField,derivativeAdditionState);

assert(abs(state_expect(1) - state_test(1)) <= 1e-3);
assert(abs(state_expect(2) - state_test(2)) <= 1e-3);
assert(abs(state_expect(3) - state_test(3)) <= 1e-3);
assert(abs(state_expect(4) - state_test(4)) <= 1e-3);
assert(abs(state_expect(5) - state_test(5)) <= 1e-3);
assert(abs(state_expect(6) - state_test(6)) <= 1e-3);
assert(abs(state_expect(7) - state_test(7)) <= 3e-3);
assert(abs(state_expect(8) - state_test(8)) <= 1e-3);
assert(abs(state_expect(9) - state_test(9)) <= 1e-3);

for i=1:4
    assert(abs(u_compute_expect(i) - u_compute_test(i)) <= 1e-3);
    assert(abs(u_effective_expect(i) - u_effective_test(i)) <= 1e-3);
end

for i=1:3
    assert(abs(m_compute_expect(i) - m_compute_test(i)) <= 1e-3);
    assert(abs(m_effective_expect(i) - m_effective_test(i)) <= 1e-3);
end

%% Test 18.4: RK4, 3 RW + 1 MGT + 2 addition states

quaternion = [0.66;0.22;-0.58;sqrt(449)/50];
angularVelocity = [2;-1;3];
wheelVelocity = [5;20;-10];
additionState = [0;0];          % Initial conditions of the 2 new states
step = 0.1;
command = @(state) [0.01*state(5);0.01*state(6);0.01*state(7);0.5];
inertia = [10 0 0;0 5 0;0 0 15];
derivativeAdditionState = @(state) [1;state(5)+state(6)+state(7)];
affectationCommand = [1 0 0 0;0 1 0 1;0 0 1 0;1 0 0 0;0 1 0 0;0 0 1 0];
magneticField = [1;2;3];

state = [quaternion;angularVelocity;wheelVelocity;additionState];

sat = satellite();
sat =sat.setInertia(inertia);

wheel = reactionWheel();
wheel = wheel.setInertia(0.01);

sat.wheels = [wheel wheel wheel];
sat.wheels(1) = sat.wheels(1).setAxis(1);
sat.wheels(2) = sat.wheels(2).setAxis(2);
sat.wheels(3) = sat.wheels(3).setAxis(3);

mgt = magnetorquer();
mgt = mgt.setAxis(2);

sat.magnetorquers = mgt;

% Expected output data:
quaternion_expect = [0.543190631156058;0.213457873565334;-0.580639394231282;0.567659901655724];
angularVelocity_expect = [2.205122531272918;-0.361269223600233;2.951848444697170];
wheelVelocity_expect = [9.276796097843512;15.301250010959478;-14.545379186069193];
addState_expect = [0.1;0.440641865788355];

state_expect = [quaternion_expect;angularVelocity_expect;wheelVelocity_expect;addState_expect];

u_compute_expect = 10*[0.002118468694874;-0.000684689130766;0.002972639093775;0.05];
u_effective_expect = 10*[0.002118468694874;-0.000684689130766;0.002972639093775;0.04606627701666];
m_compute_expect = 10*[0.005366214924939;0;0.011642278003170];
m_effective_expect = 10*[0.005366214924939;0;0.011642278003170];

% Test:
[state_test,u_compute_test,u_effective_test,m_compute_test,m_effective_test] = sat.RK4(state,command,affectationCommand,"quaternion",step,magneticField,derivativeAdditionState);

for i=1:12
    assert(abs(state_expect(i) - state_test(i)) <= 2e-3);
end

for i=1:4
    assert(abs(u_compute_expect(i) - u_compute_test(i)) <= 1e-3);
    assert(abs(u_effective_expect(i) - u_effective_test(i)) <= 1e-3);
end

for i=1:3
    assert(abs(m_compute_expect(i) - m_compute_test(i)) <= 1e-3);
    assert(abs(m_effective_expect(i) - m_effective_test(i)) <= 1e-3);
end

%% Test 19.1: Affectation command

% Input data:
wheel1 = reactionWheel();
wheel1 = wheel1.setAxis(1);

wheel2 = reactionWheel();
wheel2 = wheel2.setAxis(2);

wheel3 = reactionWheel();
wheel3 = wheel3.setAxis(3);

sat = satellite();

sat.wheels = [wheel2 wheel1 wheel3];

mgt1 = magnetorquer();
mgt1 = mgt1.setAxis(1);

mgt2 = magnetorquer();
mgt2 = mgt2.setAxis(2);

mgt3 = magnetorquer();
mgt3 = mgt3.setAxis(3);

sat.magnetorquers = [mgt3 mgt2 mgt1];

% Expected output data:
matrix_expect = [0 1 0 0 0 1;
                 1 0 0 0 1 0;
                 0 0 1 1 0 0;
                 1 0 0 0 0 0;
                 0 1 0 0 0 0;
                 0 0 1 0 0 0];

% Test:
matrix_test = sat.affectationCommand();

for row=1:6
    for col=1:6
        assert(isequal(matrix_expect(row,col),matrix_test(row,col)));
    end
end

%% Test 19.2: Affectation command, 1 wheel (roll) and 3 magnetorquers

% Input data:
wheel = reactionWheel();
wheel = wheel.setAxis(3);

mgt = magnetorquer();

sat = satellite();

sat.wheels = wheel;
sat.magnetorquers = [mgt mgt mgt];
sat.magnetorquers(1) = sat.magnetorquers(1).setAxis(1);
sat.magnetorquers(2) = sat.magnetorquers(2).setAxis(2);
sat.magnetorquers(3) = sat.magnetorquers(3).setAxis(3);

% Expected output data:
matrix_expect = [0 1 0 0;
                 0 0 1 0;
                 1 0 0 1;
                 1 0 0 0];

matrix_test = sat.affectationCommand();

size_matrix = size(matrix_test);

for row=1:size_matrix(1)
    for col=1:size_matrix(2)
        assert(isequal(matrix_expect(row,col),matrix_test(row,col)));
    end
end

%% Test 21: Wheel saturation

% Input data:
wheel = reactionWheel();

sat = satellite();

sat.wheels = [wheel wheel wheel];

velocityBound = [150;150;150];
torqueBound = [0.0005;0.0005;0.0005];
velocity = [-145;0;64];
command = [0.0004;0.05;0.00023];        % 1) sature vitesse, 2) sature vitesse puis couple, 3) ne sature pas
inertia = [3.5e-5;3.5e-5;3.5e-5];
step = 1;

% Expected output data:
command_expect = [1.75e-4;0.0005;0.00023];

% Test:
for i=1:3
    sat.wheels(i) = sat.wheels(i).setVelocityBound(velocityBound(i));
    sat.wheels(i) = sat.wheels(i).setTorqueBound(torqueBound(i));
    sat.wheels(i) = sat.wheels(i).setInertia(inertia(i));
end

command_test = sat.wheelSaturation(velocity,command,step);

for i=1:3
    assert(abs(command_expect(i) - command_test(i)) <= 1e-8);
end

%% Test 22: Magnetorquer saturation

% Input data:
mgt = magnetorquer();

sat = satellite();

sat.magnetorquers = [mgt mgt mgt];
sat.magnetorquers(1) = sat.magnetorquers(1).setAxis(1);
sat.magnetorquers(2) = sat.magnetorquers(2).setAxis(2);
sat.magnetorquers(3) = sat.magnetorquers(3).setAxis(3);

magneticSaturation = [0.005;0.005;0.005];
command = [-1;1;0.0025];

% Expected output data:
command_expect = [-0.005;0.005;0.0025];

% Test:
for i=1:3
    sat.magnetorquers(i) = sat.magnetorquers(i).setMagneticMomentBound(magneticSaturation(i));
end

command_test = sat.magnetorquerSaturation(command);

for i=1:3
    assert(isequal(command_expect(i),command_test(i)));
end

%% Test 23.1: Allocation magnetic moment, 3 magnetorquers

% Input data:
magneticField = [1;-2;3];
requiredTorque = [-4;5;6];

mgt = magnetorquer();

sat = satellite();
sat.magnetorquers = [mgt mgt mgt];
sat.magnetorquers(1) = sat.magnetorquers(1).setAxis(1);
sat.magnetorquers(2) = sat.magnetorquers(2).setAxis(2);
sat.magnetorquers(3) = sat.magnetorquers(3).setAxis(3);

% Expected output data:
magneticMoment_expect = [-27/14;-18/14;-3/14];

% Test:
magneticField_test = sat.allocationMagneticMoment(magneticField,requiredTorque);

assert(norm(magneticMoment_expect - magneticField_test) <= 1e-8);

%% Test 23.2: Allocation magnetic moment, 2 magnetorquers

% Input data:
magneticField = [1;-2;3];
requiredTorque = [-4;0;5];

sat = satellite();

% Expected output data:
magneticMoment_expect = [-10/14;-17/14;-8/14];

% Test:
magneticField_test = sat.allocationMagneticMoment(magneticField,requiredTorque);

assert(norm(magneticMoment_expect - magneticField_test) <= 1e-8);

%% Test 23.3: Allocation magnetic moment, 1 magnetorquer

% Input data:
magneticField = [1;-2;3];
requiredTorque = [0;-4;0];

sat = satellite();

% Expected output data:
magneticMoment_expect = [12/14;0;-4/14];

% Test:
magneticField_test = sat.allocationMagneticMoment(magneticField,requiredTorque);

assert(norm(magneticMoment_expect - magneticField_test) <= 1e-8);

%% Test 24.1: Application magnetic moment, 3 magnetorquers

% Input data:
magneticMoment = [1;2;3];
magneticField = [4;5;6];

mgt = magnetorquer();

sat = satellite();

sat.magnetorquers = [mgt mgt mgt];

nb_magnetorquer = numel(sat.magnetorquers);

for i=1:nb_magnetorquer
    sat.magnetorquers(i) = sat.magnetorquers(i).setAxis(i);
end

% Expected output data:
torque_expect = [-3;6;-3];

% Test:
torque_test = sat.applicationMagneticMoment(magneticField,magneticMoment);

for i=1:nb_magnetorquer
    assert(isequal(torque_expect(i),torque_test(i)));
end

%% Test 24.2: Application magnetic moment, 2 magnetorquers

% Input data:
magneticMoment = [1;2;3];
magneticField = [4;5;6];

mgt = magnetorquer();

sat = satellite();

sat.magnetorquers = [mgt mgt];

nb_magnetorquer = numel(sat.magnetorquers);

sat.magnetorquers(1) = sat.magnetorquers(1).setAxis(1);
sat.magnetorquers(2) = sat.magnetorquers(2).setAxis(3);

% Expected output data:
torque_expect = [-3;-3];

% Test:
torque_test = sat.applicationMagneticMoment(magneticField,magneticMoment);

for i=1:nb_magnetorquer
    assert(isequal(torque_expect(i),torque_test(i)));
end

%% Test 24.3: Application magnetic moment, 1 magnetorquer

% Input data:
magneticMoment = [1;2;3];
magneticField = [4;5;6];

mgt = magnetorquer();

sat = satellite();

sat.magnetorquers = mgt;

nb_magnetorquer = numel(sat.magnetorquers);

sat.magnetorquers(1) = sat.magnetorquers(1).setAxis(2);

% Expected output data:
torque_expect = 6;

% Test:
torque_test = sat.applicationMagneticMoment(magneticField,magneticMoment);

assert(isequal(torque_expect,torque_test));

%% Test 25.1: ODException, initial quaternion

% Input data:
sat = satellite();
time = 0;
command = @(state) 0;

try
    sat.ODException(time,command)
    error('Test failed: Initial quaternion must be set');
catch ME
    assert(strcmp(ME.identifier,'ODE:MissingInput'),'Initial quaternion did not set');
end

%% Test 25.2: ODException, initial angular velocity

% Input data:
sat = satellite();
sat = sat.setQuaternion([1;0;0;0]);
time = 0;
command = @(state) 0;

try
    sat.ODException(time,command);
    error('Test failed: Initial angular velocity must be set');
catch ME
    assert(strcmp(ME.identifier,'ODE:MissingInput'),'Initial angular velocity of the satellite did not set');
end

%% Test 25.3: ODException, more than 3 wheels

% Input data:
wheel = reactionWheel();
sat = satellite();
sat = sat.setQuaternion([1;0;0;0]);
sat = sat.setAngularVelocity([0;0;0]);
sat.wheels = [wheel wheel wheel wheel];
time = 0;
command = @(state) 0;

try
    sat.ODException(time,command);
    error('Test failed: The satellite can have 3 reaction wheels maximum');
catch ME
    assert(strcmp(ME.identifier,'ODE:ReactionWheel'),'The number of reaction wheels must not exceed 3');
end

%% Test 25.4: ODException, wheels with the same axis

% Input data:
wheel = reactionWheel();
wheel = wheel.setAxis(2);
sat = satellite();
sat = sat.setQuaternion([1;0;0;0]);
sat = sat.setAngularVelocity([0;0;0]);
sat.wheels = [wheel wheel];
time = 0;
command = @(state) 0;

try
    sat.ODException(time,command);
    error('Test failed: Reaction wheels must have differents axis number');
catch ME
    assert(strcmp(ME.identifier,'ODE:ReactionWheel'),'Reaction wheels must have different axis number');
end

sat.wheels = [wheel wheel wheel];

try
    sat.ODException(time,command);
    error('Test failed: Reaction wheels must have differents axis number');
catch ME
    assert(strcmp(ME.identifier,'ODE:ReactionWheel'),'Reaction wheels must have different axis number');
end

%% Test 25.5: ODException, more than 3 magnetorquers

% Input data:
mgt = magnetorquer();
sat = satellite();
sat = sat.setQuaternion([1;0;0;0]);
sat = sat.setAngularVelocity([0;0;0]);
sat.magnetorquers = [mgt mgt mgt mgt];
time = 0;
command = @(state) 0;

try
    sat.ODException(time,command);
    error('Test failed: The satellite can have 3 magnetorquers maximum');
catch ME
    assert(strcmp(ME.identifier,'ODE:Magnetorquer'),'The number of magnetorquer must not exceed 3');
end

%% Test 25.6: ODExecption, magnetorquers with the same axis

% Input data:
mgt = magnetorquer();
mgt = mgt.setAxis(3);
sat = satellite();
sat = sat.setQuaternion([1;0;0;0]);
sat = sat.setAngularVelocity([0;0;0]);
sat.magnetorquers = [mgt mgt];
time = 0;
command = @(state) 0;

try
    sat.ODException(time,command);
    error('Test failed: Magnetorquers must have differents axis number');
catch ME
    assert(strcmp(ME.identifier,'ODE:Magnetorquer'),'Magnetorquers must have different axis number');
end

sat.magnetorquers = [mgt mgt mgt];
try
    sat.ODException(time,command);
    error('Test failed: Magnetorquers must have differents axis number');
catch ME
    assert(strcmp(ME.identifier,'ODE:Magnetorquer'),'Magnetorquers must have different axis number');
end

%% Test 25.7: ODException, existing command

% Input data:
sat = satellite();
sat = sat.setQuaternion([1;0;0;0]);
sat = sat.setAngularVelocity([0;0;0]);
time = 0;

try
    sat.ODException(time,1);
    error('Test failed: ODE must take a command function');
catch ME
    assert(strcmp(ME.identifier,'ODE:ExistingCommand'),'ODE function must take a command function');
end

%% Test 25.8: ODException, command option is a char or a string

% Input data:
sat = satellite();
sat = sat.setQuaternion([1;0;0;0]);
sat = sat.setAngularVelocity([0;0;0]);
time = 0;
command = @(state) 0;

try
    sat.ODException(time,command,1);
    error('Test failed: Command option must be a char or a string');
catch ME
    assert(strcmp(ME.identifier,'ODE:BadInput'),'The command option must be a char or a string');
end

%% Test 25.9: ODException, command option is "quaternion" or "cardan"

% Input data:
sat = satellite();
sat = sat.setQuaternion([1;0;0;0]);
sat = sat.setAngularVelocity([0;0;0]);
time = 0;
command = @(state) 0;

try
    sat.ODException(time,command,"not_quaternion");
    error('Test failed: Command option must be "quaternion" or "cardan"');
catch ME
    assert(strcmp(ME.identifier,'ODE:BadInput'),'Command option must be "quaternion" or "cardan"');
end

%% Test 25.10: ODException, existing orbit

% Input data:
mgt = magnetorquer();
sat = satellite();
sat = sat.setQuaternion([1;0;0;0]);
sat = sat.setAngularVelocity([0;0;0]);
sat.magnetorquers = mgt;
time = 0;
command = @(state) 0;

try
    sat.ODException(time,command,"quaternion");
    error('Test failed: There are some magnetorquers, an orbit instance must be give to ODE function');
catch ME
    assert(strcmp(ME.identifier,'ODE:MissingOrbit'),'There are some magnetorquers, an orbit instance must be give to ODE function');
end

try
    sat.ODException(time,command,"quaternion",1);
    disp('Apres');
    error('Test failed: There are some magnetorquers, an orbit instance must be give to ODE function');
catch ME
    assert(strcmp(ME.identifier,'ODE:MissingOrbit'),'There are some magnetorquers, an orbit instance must be give to ODE function');
end

%% Test 25.11: ODException, existing planet istance

% Input data:
sat = satellite();
sat = sat.setQuaternion([1;0;0;0]);
sat = sat.setAngularVelocity([0;0;0]);
time = 0;
command = @(state) 0;
orb = orbit();

try
    sat.ODException(time,command,"quaternion",orb);
    error('Test failded: Orbit must have a planet instance');
catch ME
    assert(strcmp(ME.identifier,'ODE:MissingPlanet'),'Orbit must have a planet instance');
end

%% Test 25.12: ODException, frame is a char or a string

% Input data:
sat = satellite();
sat = sat.setQuaternion([1;0;0;0]);
sat = sat.setAngularVelocity([0;0;0]);
time = 0;
command = @(state) 0;
earth = planet();
orb = orbit();
orb.attractorBody = earth;

try
    sat.ODException(time,command,"quaternion",orb,1);
    error('Test failed: Frame must be a char or string');
catch ME
    assert(strcmp(ME.identifier,'ODE:BadInput'),'Frame must be a char or string');
end

%% Test 25.13: ODException, derivative addition state

% Input data:
sat = satellite();
sat = sat.setQuaternion([1;0;0;0]);
sat = sat.setAngularVelocity([0;0;0]);
time = 0;
command = @(state) 0;
earth = planet();
orb = orbit();
orb.attractorBody = earth;

try
    sat.ODException(time,command,"quaternion",orb,"ECI",1,1);
    disp('apres');
    error('Test failed: There some addition state, ODE must take a function to derivate them');
catch ME
    assert(strcmp(ME.identifier,'ODE:ExistingDerivativeAdditionState'),'There some addition state, ODE must take a function to derivate them');
end

%% Test 26: Dynamic model

% Input data:
quaternion = [0.66;0.22;-0.58;sqrt(449)/50];
angularVelocity = [2;-1;3];
wheelVelocity = 5;
wheelInertia = 0.01;
inertia = [1 0 0;0 2 0;0 0 3];
command = [1;2;3;0.2];
addState = 0;

state = [quaternion;angularVelocity;wheelVelocity;addState];

deriveAddState = @(state) state(5) + state(8);

wheel = reactionWheel();
wheel = wheel.setAxis(3);
wheel = wheel.setInertia(wheelInertia);

sat = satellite();
sat = sat.setInertia(inertia);
sat.wheels = wheel;

% Expected output data:
state_expect = [-1.1457;0.0019;-0.2362;1.46;4;7;1.6667;-20;7];

% Test:
state_test = sat.dynamicModel(state,command,deriveAddState);

for i=1:9
    assert(abs(state_expect(i) - state_test(i)) <= 1e-4);
end

%% Test 27.1: Applied command, option quaternion

% Input data:
quaternion = [0;1;0;0];
angularVelocity = [0;0;0];
wheelVelocity = 20;
command = @(state) [state(2);1;2;3];        % 1 RW + 3 MGT
wheelInertia = 3.5e-5;
step = 1;
commandOption = "quaternion";
magneticField = [2;1;0];

state = [quaternion;angularVelocity;wheelVelocity];

wheel = reactionWheel();
wheel = wheel.setInertia(wheelInertia);
wheel = wheel.setTorqueBound(0.5);
wheel = wheel.setVelocityBound(100);

mgt = magnetorquer();
mgt = mgt.setMagneticMomentBound(1);

sat = satellite();
sat.wheels = wheel;
sat.magnetorquers = [mgt mgt mgt];
sat.magnetorquers(1) = sat.magnetorquers(1).setAxis(1);
sat.magnetorquers(2) = sat.magnetorquers(2).setAxis(2);
sat.magnetorquers(3) = sat.magnetorquers(3).setAxis(3);

% Expected output data:
U_compute_expect = [1;1;2;3];
U_effective_expect = [4.2e-3;1;2;13/5];
M_compute_expect = [-3/5;-6/5;1];
M_effective_expect = [-3/5;-1;1];

% Test:
[U_compute_test,U_effective_test,M_compute_test,M_effective_test] = sat.appliedCommand(state,command,commandOption,step,magneticField);

for i=1:3
    assert(abs(M_compute_test(i) - M_compute_expect(i)) <= 1e-4);
    assert(abs(M_effective_test(i) - M_effective_expect(i)) <= 1e-4);
end

for i=1:4
    assert(abs(U_compute_expect(i) - U_compute_test(i)) <= 1e-4);
    assert(abs(U_effective_expect(i) - U_effective_test(i)) <= 1e-4);
end

%% Test 27.2: Applied command, option cardan

% Input data:
quaternion = [0;1;0;0];
angularVelocity = [0;0;0];
wheelVelocity = 20;
command = @(state) state(1);        % 1 RW
wheelInertia = 3.5e-5;
step = 1;
commandOption = "cardan";
magneticField = [2;1;0];

state = [quaternion;angularVelocity;wheelVelocity];

wheel = reactionWheel();
wheel = wheel.setInertia(wheelInertia);
wheel = wheel.setTorqueBound(0.5);

sat = satellite();
sat.wheels = wheel;

% Expected output data:
U_compute_expect = pi;
U_effective_expect = 0.5;
M_effective_expect = [0;0;0];
M_compute_expect = [0;0;0];

% Test:
[U_compute_test,U_effective_test,M_compute_test,M_effective_test] = sat.appliedCommand(state,command,commandOption,step,magneticField);

for i=1:3
    assert(abs(M_compute_test(i) - M_compute_expect(i)) <= 1e-4);
    assert(abs(M_effective_test(i) - M_effective_expect(i)) <= 1e-4);
end

for i=1:1
    assert(abs(U_compute_expect(i) - U_compute_test(i)) <= 1e-4);
    assert(abs(U_effective_expect(i) - U_effective_test(i)) <= 1e-4);
end

%% Test 28.1: ODE, exception addition state without derivative for them

% Input data:
time = 0:1:10;
command = @(state) 0;
commandOption = 'quaternion';
orb = orbit();
frame = 'ECI';
additionState = 0;

sat = satellite();

% Test
try
    sat.ODE(time,command,commandOption,orb,frame,additionState);
    error('Test failed: There are 1 or more addition state, ODE must take a function to derivate them');
catch ME
    assert(strcmp(ME.identifier,'ODE:MissingInput'),'There are 1 or more addition state, ODE must take a function to derivate them');
end

%% Test 28.2: ODE with 1 reaction wheel and 3 magnetorquers without saturations

% Input data:
time = 0:0.0001:1;
quaternion = [0;1;0;0];
angularVelocity = [0.01;0;0];
wheelVelocity = 20;
inertia = [10 0 0;0 5 0;0 0 15];
wheelInertia = 3.5e-5;

sat = satellite();

wheel = reactionWheel();
sat.wheels = wheel;

mgt = magnetorquer();
sat.magnetorquers = [mgt mgt mgt];

sat.wheels(1) = sat.wheels(1).setAxis(3);
sat.wheels(1) = sat.wheels(1).setInertia(wheelInertia);
sat.wheels(1) = sat.wheels(1).setVelocity(wheelVelocity);

for i = 1:3    
    sat.magnetorquers(i) = sat.magnetorquers(i).setAxis(i);
end

sat = sat.setQuaternion(quaternion);
sat = sat.setAngularVelocity(angularVelocity);
sat = sat.setInertia(inertia);

com = @(state) [-state(5);
                state(6);
                state(7);
                state(8)];

earth = planet();
earth = earth.setRadius(6378);
earth = earth.setGravitationalParameter(398600.44);
earth = earth.setRotation(2*pi/84600);

orb = orbit();
orb = orb.setInclination(87);
orb = orb.setArgumentOfPeriapsis(20);
orb = orb.setRAAN(30);
orb = orb.setMajorSemiAxis(6978);
orb = orb.setEccentricity(0.001);
orb = orb.setInitialTrueAnomaly(0);
orb = orb.setGHA(0);
orb = orb.setDateInitial([2025 09 26 10 09 30]);
orb.attractorBody = earth;

Th = @(state) state(9) - 1;
additionState = [0;0];
derivative = @(state,Th) [1;state(5)-Th(state)];

% Expected out data:
state_expect = (1e2)*[0.000213006890426;-0.000837147873711;0.009962613578231;-0.000075802581904;-0.003372602082310;0.000123850541733;0.213787488026439;5.681522954725107;0.01;0.005191853303415];

% Test:

[T,Y,U,M] = sat.ODE(time,com,'quaternion',orb,'ECI',additionState,derivative,Th);

% Quaternion
assert(abs(state_expect(1) - Y(1,end)) <= 1e-3);
assert(abs(state_expect(2) - Y(2,end)) <= 1e-3);
assert(abs(state_expect(3) - Y(3,end)) <= 1e-3);
assert(abs(state_expect(4) - Y(4,end)) <= 1e-3);

% Angular velocity
assert(abs(state_expect(5) - Y(5,end)) <= 1e-3);
assert(abs(state_expect(6) - Y(6,end)) <= 1e-3);
assert(abs(state_expect(7) - Y(7,end)) <= 2e-3);

% Wheel + 2 addition state
assert(abs(state_expect(8) - Y(8,end)) <= 1);
assert(abs(state_expect(9) - Y(9,end)) <= 1e-3);
assert(abs(state_expect(10) - Y(10,end)) <= 1e-3);