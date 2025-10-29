% Test of the class reactionWheel.m

% To run tests:
%       runtests('reactionWheelTest')

%% Test 0: Constructor

wheel = reactionWheel();

assert(isequal(wheel.getInertia(),1));
assert(isequal(wheel.getVelocityBound(),1e30));
assert(isequal(wheel.getTorqueBound(),1e30));
assert(isequal(wheel.getVelocity(),0));
assert(isequal(wheel.getAxis(),1));

%% Test 1.1: Getter/Setter inertia

% Input data:
inertia = 5.36;

% Expected output data:
inertia_expect = 5.36;

% Test:
wheel = reactionWheel();

wheel = wheel.setInertia(inertia);
inertia_test = wheel.getInertia();

assert(inertia_expect,inertia_test);

%% Test 1.2: Inertia positive exception

wheel = reactionWheel();

try
    wheel.setInertia(-1);
    error('Test failed: The reaction wheel inertia must be positive');
catch ME
    assert(strcmp(ME.identifier,'reactionWheel:BadInput'),'Reaction wheel inertia must be positive');
end

%% Test 2.1: Getter/Setter axis

% Input data:
axis = 3;

% Expected output data:
axis_expect = 3;

% Test:
wheel = reactionWheel();

wheel = wheel.setAxis(axis);
axis_test = wheel.getAxis();

assert(axis_test,axis_expect);

%% Test 2.2: Axis equal to 1, 2 or 3 exception

wheel = reactionWheel();

try
    wheel.setAxis(1.5);
    error('Test failed: Wheel axis must be equal to 1, 2 or 3');
catch ME
    assert(strcmp(ME.identifier,'reactionWheel:BadInput'),'Wheel axis must be equal to 1, 2 or 3');
end

%% Test 3.1: Getter/Setter velocity bound

% Input data:
velocityBound = 293;

% Expected output data:
velocityBound_expect = 293;

% Test:
wheel = reactionWheel();

wheel = wheel.setVelocityBound(velocityBound);
velocityBound_test = wheel.getVelocityBound();

assert(isequal(velocityBound_expect,velocityBound_test));

%% Test 3.2: Velocity bound exception

wheel = reactionWheel();

try
    wheel.setVelocityBound(-1);
    error('Test failed: Velocity bound must be positive');
catch ME
    assert(strcmp(ME.identifier,'reactionWheel:BadInput'),'Velocity bound must be positive');
end

%% Test 4.1: Getter/Setter torque bound

% Input data:
torqueBound = 0.0005;

% Expected output data:
torqueBound_expect = 0.0005;

% Test:
wheel = reactionWheel();

wheel = wheel.setTorqueBound(torqueBound);
torqueBound_test = wheel.getTorqueBound();

assert(isequal(torqueBound_expect,torqueBound_test));

%% Test 4.2: Torque bound exception

wheel = reactionWheel();

try
    wheel.setTorqueBound(-1);
    error('Test failed: Torque bound must be positive');
catch ME
    assert(strcmp(ME.identifier,'reactionWheel:BadInput'),'Torque bound must be positive');
end

%% Test 5: Getter/Setter velocity

% Input data:
velocity = 6.35;

% Expected output data:
velocity_expect = 6.35;

% Test:
wheel = reactionWheel();

wheel = wheel.setVelocity(velocity);
velocity_test = wheel.getVelocity();

assert(isequal(velocity_expect,velocity_test));

%% Test 6.1: Velocity saturation, with saturation

% Input data:
velocity = -140;
inertia = 3.5e-5;
velocityBound = 150;
step = 2;
command = 1;

wheel = reactionWheel();

% Expected output data:
command_expect = 1.75e-4;

% Test:
wheel = wheel.setInertia(inertia);
wheel = wheel.setVelocityBound(velocityBound);

command_test = wheel.velocitySaturation(velocity,command,step);

assert(abs(command_expect - command_test) <= 1e-5);

%% Test 6.2: Velocity saturation, without saturation

% Input data:
velocity = -140;
inertia = 3.5e-5;
velocityBound = 150;
step = 2;
command = -1e-3;

wheel = reactionWheel();

% Expected output data:
command_expect = -1e-3;

% Test:
wheel = wheel.setInertia(inertia);
wheel = wheel.setVelocityBound(velocityBound);

command_test = wheel.velocitySaturation(velocity,command,step);

assert(abs(command_expect - command_test) <= 1e-5);

%% Test 7.1: Torque saturation, with saturation

% Input data:
torqueBound = 0.005;
command = 1;

wheel = reactionWheel();

% Expected output data:
command_expect = 0.005;

% Test:
wheel = wheel.setTorqueBound(torqueBound);

command_test = wheel.torqueSaturation(command);

assert(isequal(command_expect,command_test));

%% Test 7.2: Torque saturation, without saturation

% Input data:
torqueBound = 0.005;
command = -0.002;

wheel = reactionWheel();

% Expected output data:
command_expect = -0.002;

% Test:
wheel = wheel.setTorqueBound(torqueBound);

command_test = wheel.torqueSaturation(command);

assert(isequal(command_expect,command_test));