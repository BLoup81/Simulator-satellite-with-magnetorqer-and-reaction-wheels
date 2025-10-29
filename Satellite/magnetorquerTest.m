% Test of the class magnetorquer

% To run test:
%       runtests('magnetorquerTest')

%% Test 1: Constructor

mgt = magnetorquer();

assert(isequal(mgt.getMagneticMomentBound(),1e30));
assert(isequal(mgt.getAxis(),1));

%% Test 2.1: Getter/Setter axis

% Input data:
axis = 2;

% Expected output data:
axis_expect = 2;

% Test:
mgt = magnetorquer();

mgt = mgt.setAxis(axis);
axis_test = mgt.getAxis();

assert(isequal(axis_expect,axis_test));

%% Test 2.2: Axis equal to 1, 2 or 3 exception

mgt = magnetorquer();

try
    mgt.setAxis(1.5);
    error('Test failed: Magnetorquer must be equal to 1, 2 or 3');
catch ME
    assert(strcmp(ME.identifier,'magnetorquer:BadInput'),'Magnetorquer axis must be equal to 1, 2 or 3');
end

%% Test 3.1: Getter/Setter magnetic moment saturation

% Input data:
saturation = 0.0056;

% Expected output data:
saturation_expect = 0.0056;

% Test:
mgt = magnetorquer();

mgt = mgt.setMagneticMomentBound(saturation);
saturation_test = mgt.getMagneticMomentBound();

assert(isequal(saturation_expect,saturation_test));

%% Test 3.2: Magnetic saturation positive exception

mgt = magnetorquer();

try
    mgt.setMagneticMomentBound(-1);
    error('Test failed: Magnetic saturation must be positive');
catch ME
    assert(strcmp(ME.identifier,'magnetorquer:BadInput'),'Magnetic saturation must be positive');
end

%% Test 4.1: Magnetic moment saturation, with saturation

% Input data:
saturation = 0.002;
command = -1;

mgt = magnetorquer();

% Expected output data:
command_expect = -0.002;

% Test:
mgt = mgt.setMagneticMomentBound(saturation);

command_test = mgt.magneticMomentSaturation(command);

assert(isequal(command_expect,command_test));