% Test of the class planet

% To run tests:
%       runtests('planetTest')

%% Test 0: Constructor

earth = planet();

assert(isequal(earth.getRadius(),6378));
assert(isequal(earth.getRotation(),2*pi/86400));
assert(isequal(earth.getGravitationalParameter(),398600.44));

%% Test 1.1: Getter/Setter radius

% Input data:
radius = 6378;

% Expected output data:
radius_expect = 6378;

% Test:
earth = planet();

earth = earth.setRadius(radius);

assert(isequal(earth.getRadius(),radius_expect));

%% Test 1.2: Set radius, positivity exception

earth = planet();

try
    earth.setRadius(-1);
    error('Test failed: Radius of the planet must be positive');
catch ME
    assert(strcmp(ME.identifier,'planet:BadInput'),'Radius of the planet must be positive');
end


%% Test 2.1: Getter/Setter gravitational parameter

% Input data:
gravitationalParameter = 398600.44;

% Expected output data:
gravitationalParameter_expect = 398600.44;

% Test:
earth = planet();

earth = earth.setGravitationalParameter(gravitationalParameter);
assert(isequal(earth.getGravitationalParameter(),gravitationalParameter_expect));

%% Test 2.2: Set gravitational parameter, positivity exception

earth = planet();

try
    earth.setGravitationalParameter(-1);
    error('Test failed: Gravitational parameter must be positive');
catch ME
    assert(strcmp(ME.identifier,'planet:BadInput'),'Gravitational parameter must be positive');
end

%% Test 3: Getter/Setter rotation speed

% Input data:
rotation = 7.292e-5;

% Expected output data:
rotation_expect = 7.292e-5;

% Test:
earth = planet();

earth = earth.setRotation(rotation);
assert(isequal(earth.getRotation(),rotation_expect));