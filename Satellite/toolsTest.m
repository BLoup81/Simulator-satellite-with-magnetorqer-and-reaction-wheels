% Test of the file tools.m

% To run tests:
%       runtests('toolsTest')

%% Test 1: Euler to quaternion

% Input data:
euler = [0 pi -pi/2;
          0 0 -pi;
          0 -pi/2 pi];

% Expected data output:
quat_expect = [1 0 sqrt(2)/2;
               0 sqrt(2)/2 sqrt(2)/2;
               0 -sqrt(2)/2 0;
               0 0 0];

err = 1e-8;

quat = tools.euler2quater(euler);

assert(norm(quat(:,1) - quat_expect(:,1)) <= err);
assert(norm(quat(:,2) - quat_expect(:,2)) <= err);
assert(norm(quat(:,3) - quat_expect(:,3)) <= err);


%% Test 2: Quaternion to euler

% Input data:
quat = [1 0 sqrt(2)/2;
        0 sqrt(2)/2 sqrt(2)/2;
        0 -sqrt(2)/2 0;
        0 0 0];

% Expected output data:
euler_expect = [0 pi pi/2;
         0 0 0;
         0 -pi/2 0];

err = 1e-8;

euler = tools.quater2euler(quat);

assert(norm(euler(:,1) - euler_expect(:,1)) <= err);
assert(norm(euler(:,2) - euler_expect(:,2)) <= err);
assert(norm(euler(:,3) - euler_expect(:,3)) <= err);

%% Test 3.1: ECEF to ECI

% Input data:
ecef = [1 1 1;
        1 1 1;
        1 1 1];
earthRotation = 2*pi/86164;
time = [0;10000;86164];

% Expected output data:
eci_expect = [1 0.079416937023335 1;
              1 1.411981922729123 1;
              1 1 1];

% Test:
eci = tools.ecef2eci(ecef,earthRotation,time);

err = 1e-8;
assert(norm(eci_expect(:,1) - eci(:,1)) <= err);
assert(norm(eci_expect(:,2) - eci(:,2)) <= err);
assert(norm(eci_expect(:,3) - eci(:,3)) <= err);


%% Test 3.2: ECEF to ECI, number element exception

% Input data:
ecef = [1 1 1;
        1 1 1;
        1 1 1];
earthRotation = 2*pi/86164;
time = [0;3056];

try
    tools.ecef2eci(ecef,earthRotation,time);
    error('Test failed: Number of time and ECEF coordinates must be the same');
catch ME
    assert(strcmp(ME.identifier,'ecef2eci:DimensionError'),'Number of time and ECEF coordinates must be the same');
end

%% Test 4.1: ECI to LVLH

% Input data:
inclination = pi/6;
argumentOfPerigee = pi/9;
raan = pi/4;
trueAnomaly = [0;pi/2;pi/4;17*pi/20];
eci = [1 2 1 -3;
       1 3 1 0.5;
       1 4 1 2];

% Expected output data:
lvlh_expect = [-0.01384321 -4.21579912 -1.07040363 -2.90443722;
               -0.86602540 -3.11054822 -0.86602540 -0.49461394;
               -1.49993612 -1.24560329 -1.05082637 -2.13766262];

% Test:
lvlh = tools.eci2lvlh(eci,argumentOfPerigee,raan,inclination,trueAnomaly);

err = 1e-8;
assert(norm(lvlh_expect(:,1) - lvlh(:,1)) <= err);
assert(norm(lvlh_expect(:,2) - lvlh(:,2)) <= err);
assert(norm(lvlh_expect(:,3) - lvlh(:,3)) <= err);
assert(norm(lvlh_expect(:,4) - lvlh(:,4)) <= err);

%% Test 4.2: ECI to LVLH, number element exception

% Input data:
inclination = pi/6;
argumentOfPerigee = pi/9;
raan = pi/4;
trueAnomaly = [0;0.5;1];
eci = [1 2 3 4;
       1 2 3 4;
       1 2 3 4];

try
    tools.eci2lvlh(eci,argumentOfPerigee,raan,inclination,trueAnomaly);
    error('Test failed: True anomaly and ECI coordinates must have the same number of elements');
catch ME
    assert(strcmp(ME.identifier,'eci2lvlh:DimensionError'),'True anomaly and ECI coordinates must have the same number of elements');
end

%% Test 4.3: ECI to LVLH, inclination exception

% Input data:
inclination = [pi -pi];
argumentOfPerigee = pi/9;
raan = pi/4;
trueAnomaly = 1;
eci = [1;1;1];

try
    tools.eci2lvlh(eci,argumentOfPerigee,raan,inclination(1),trueAnomaly);
    error('Test failed: Inclination must be in radian between [-pi/2 pi/2]');
catch ME
    assert(strcmp(ME.identifier,'eci2lvlh:BadInput'),'Inclination must be in radian between [-pi/2 pi/2]');
end

try
    tools.eci2lvlh(eci,argumentOfPerigee,raan,inclination(2),trueAnomaly);
    error('Test failed: Inclination must be in radian between [-pi/2 pi/2]');
catch ME
    assert(strcmp(ME.identifier,'eci2lvlh:BadInput'),'Inclination must be in radian between [-pi/2 pi/2]');
end

%% Test 4.4: ECI to LVLH, argument of perigee exception

% Input data:
inclination = 0;
argumentOfPerigee = [-1 7];
raan = pi/4;
trueAnomaly = 1;
eci = [1;1;1];

try
    tools.eci2lvlh(eci,argumentOfPerigee(1),raan,inclination,trueAnomaly);
    error('Test failed: Argument of perigee must be in radian between [0 2*pi]');
catch ME
    assert(strcmp(ME.identifier,'eci2lvlh:BadInput'),'Argument of perigee must be in radian between [0 2*pi]');
end

try
    tools.eci2lvlh(eci,argumentOfPerigee(2),raan,inclination,trueAnomaly);
    error('Test failed: Argument of perigee must be in radian between [0 2*pi]');
catch ME
    assert(strcmp(ME.identifier,'eci2lvlh:BadInput'),'Argument of perigee must be in radian between [0 2*pi]');
end

%% Test 4.5: ECI to LVLH, RAAN exception

% Input data:
inclination = 0;
argumentOfPerigee = pi/9;
raan = [-1 7];
trueAnomaly = 1;
eci = [1;1;1];

try
    tools.eci2lvlh(eci,argumentOfPerigee,raan(1),inclination,trueAnomaly);
    error('Test failed: RAAN must be in radian between [0 2*pi]');
catch ME
    assert(strcmp(ME.identifier,'eci2lvlh:BadInput'),'RAAN must be in radian between [0 2*pi]');
end

try
    tools.eci2lvlh(eci,argumentOfPerigee,raan(2),inclination,trueAnomaly);
    error('Test failed: RAAN must be in radian between [0 2*pi]');
catch ME
    assert(strcmp(ME.identifier,'eci2lvlh:BadInput'),'RAAN must be in radian between [0 2*pi]');
end


%% Test 5.1: Frame to BRF, multiple rotation

% rotation x: pi/3, y: pi/6, z: pi/4

% Input data:
quat = [0.822363171905999;0.360423405650356;0.391903837329120;0.200562121146575];

inputCoordinates = [1;-0.5;3];

% Expected output data:
brf_expect = [-1.193813782152102;1.872763023034038;2.305985106868668];

% Test:
brf_test = tools.frame2brf(inputCoordinates,quat);

err = 1e-8;
assert(norm(brf_test - brf_expect) <= err);

%% Test 5.2: Frame to BRF, several coordinates

% Input data:
quat = [0.965925826 sqrt(2)/2 sqrt(2)/2 sqrt(2)/2;
        0.258819045 sqrt(2)/2 0 0;
        0 0 sqrt(2)/2 0;
        0 0 0 sqrt(2)/2];

inputCoordinates = [1 1 1 1;
                    1 1 1 1;
                    1 1 1 1];

% Expected output data:
brf_expect = [1 1 -1 1;
              1.36602540 1 1 -1;
              0.36602540 -1 1 1];

% Test:
brf_test = tools.frame2brf(inputCoordinates,quat);

err = 1e-8;
for index = 1:length(inputCoordinates(1,:))
    assert(norm(brf_expect(:,index) - brf_test(:,index)) <= err);
end

%% Test 5.3: Frame to BRF, dimension exception

% Input data:
quat = [1 0 0;
        0 1 0;
        0 0 1;
        0 0 0];

inputCoordinates = [1 1 1 1;
                    1 1 1 1;
                    1 1 1 1];

try
    tools.frame2brf(inputCoordinates,quat);
    error('Test failed: Input coordinates and quaternions must have the same number of elements');
catch ME
    assert(strcmp(ME.identifier,'frame2brf:BadInput'),'Input coordinates and quaternions must have the same number of elements');
end

%% Test 6.1: NED to ECEF

% Input data:
lat = [0;0;pi/2];
lon = [0;pi/2;pi/2];
ned = [1 1 1;
       1 1 1;
       1 1 1];

% Expected output data:
ecef_expect = [-1 -1 -1;
               1 -1 -1;
               1 1 -1];

% Test:
ecef = tools.ned2ecef(lat,lon,ned);
err = 1e-8;

assert(norm(ecef_expect(:,1) - ecef(:,1)) <= err);
assert(norm(ecef_expect(:,2) - ecef(:,2)) <= err);
assert(norm(ecef_expect(:,3) - ecef(:,3)) <= err);

%% Test 6.2: NED to ECEF, dimension exception

% Input data:
lat = 0;
lon = 0;
ned = [1 2;
       1 2;
       1 2];

try
    tools.ned2ecef(lat,lon,ned);
    error('Test failed: Latitude, longitude and NED coordinates must have the same size');
catch ME
    assert(strcmp(ME.identifier,'ned2ecef:DimensionError'),'Latitude, longitude and NED coordinates must have the same size');
end

%% Test 7.1: Conversion coordinates, NED coordinates

% Input data:
earth = planet();
earth = earth.setRotation(2*pi/86400);
earth = earth.setGravitationalParameter(398600.44);

orb = orbit();
orb = orb.setInclination(87);
orb = orb.setRAAN(30);
orb = orb.setArgumentOfPeriapsis(20);
orb = orb.setEccentricity(0.001);
orb = orb.setMajorSemiAxis(6978);
orb = orb.setInitialTrueAnomaly(0);
orb = orb.setGHA(0);
orb.attractorBody = earth;

time = 0:1000:1000;
nb = numel(time);

inputCoordinates = [ones(1,nb);2*ones(1,nb);3*ones(1,nb)];

% Expected output data:
outputCoordinates_expect = [1 1;
                            2 2;
                            3 3];

% Test:

[lat,lon,alt] = orb.geocentric(time);

geocentric = [lat';lon'];

outputCoordinates_test = tools.conversionCoordinates(inputCoordinates,geocentric,orb,time,'NED');

assert(norm(outputCoordinates_expect - outputCoordinates_test) <= 1e-8);

%% Test 7.2: Conversion coordinates, ECEF coordinates

% Input data:
earth = planet();
earth = earth.setRotation(2*pi/86400);
earth = earth.setGravitationalParameter(398600.44);

orb = orbit();
orb = orb.setInclination(87);
orb = orb.setRAAN(30);
orb = orb.setArgumentOfPeriapsis(20);
orb = orb.setEccentricity(0.001);
orb = orb.setMajorSemiAxis(6978);
orb = orb.setInitialTrueAnomaly(0);
orb = orb.setGHA(0);
orb.attractorBody = earth;

time = 0:1000:1000;
nb = numel(time);

inputCoordinates = [ones(1,nb);2*ones(1,nb);3*ones(1,nb)];


% Expected output data:
outputCoordinates_expect = [-1.449895 1.65836;
                            2.5982 -3.125263;
                            -2.268736 1.217610];

% Test:
[lat,lon,alt] = orb.geocentric(time);

geocentric = [lat;lon];

outputCoordinates_test = tools.conversionCoordinates(inputCoordinates,geocentric,orb,time,'ECEF');

assert(norm(outputCoordinates_expect - outputCoordinates_test) <= 1e-5);

%% Test 7.3: Conversion coordinates, ECI coordinates

% Input data:
earth = planet();
earth = earth.setRotation(2*pi/86400);
earth = earth.setGravitationalParameter(398600.44);

orb = orbit();
orb = orb.setInclination(87);
orb = orb.setRAAN(30);
orb = orb.setArgumentOfPeriapsis(20);
orb = orb.setEccentricity(0.001);
orb = orb.setMajorSemiAxis(6978);
orb = orb.setInitialTrueAnomaly(0);
orb = orb.setGHA(0);
orb.attractorBody = earth;

time = 0:1000:1000;
nb = numel(time);

inputCoordinates = [ones(1,nb);2*ones(1,nb);3*ones(1,nb)];

% Expected output data:
outputCoordinates_expect = [-1.449895 1.881052;
                            2.5982 -2.996509
                            -2.268736 1.217610];

% Test:
[lat,lon,alt] = orb.geocentric(time);

geocentric = [lat;lon];

outputCoordinates_test = tools.conversionCoordinates(inputCoordinates,geocentric,orb,time,'ECI');

assert(norm(outputCoordinates_expect - outputCoordinates_test) <= 1e-5);

%% Test 7.4: Conversion coordinates, LVLH coordinates

% Input data:
earth = planet();
earth = earth.setRotation(2*pi/86400);
earth = earth.setGravitationalParameter(398600.44);

orb = orbit();
orb = orb.setInclination(87);
orb = orb.setRAAN(30);
orb = orb.setArgumentOfPeriapsis(20);
orb = orb.setEccentricity(0.001);
orb = orb.setMajorSemiAxis(6978);
orb = orb.setInitialTrueAnomaly(0);
orb = orb.setGHA(0);
orb.attractorBody = earth;

time = 0:1000:1000;
nb = numel(time);

inputCoordinates = [ones(1,nb);2*ones(1,nb);3*ones(1,nb)];

% Expected output data:
outputCoordinates_expect = [-1.997542 0.011082;
                            3.089714 -3.594459;
                            0.680803 -1.039107];

% Test:
[lat,lon,alt] = orb.geocentric(time);

geocentric = [lat;lon];

outputCoordinates_test = tools.conversionCoordinates(inputCoordinates,geocentric,orb,time,'LVLH');

assert(norm(outputCoordinates_expect - outputCoordinates_test) <= 1e-5);

%% Test 7.5: Conversion coordinates, exception frame unknow

% Input data:
time = 0;
geocentric = [1;2];
orb = orbit();
inputCoordinates = [4;5;6];

try
    tools.conversionCoordinates(inputCoordinates,geocentric,orb,time,'frame test');
    error('Test failed: The given frame test is unknow');
catch ME
    assert(strcmp(ME.identifier,'conversionCoordinates:UnknowFrame'), 'The given frame test is unknow');
end

%% Test 8.1: Set legend cardan, 1 actuator

% Input data:
wheel = reactionWheel();

sat = satellite();

sat.wheels = wheel;

% Test:
% Actuator on yaw angle
sat.wheels = sat.wheels.setAxis(1);
f = figure('Visible','off');
plot(1,1);
tools.setLegendCardan(sat.wheels);
lgd = legend;

assert(isequal(lgd.String,{'Yaw'}));
close(f);

% Actuator on pitch angle
sat.wheels = sat.wheels.setAxis(2);
f = figure('Visible','off');
plot(1,1);
tools.setLegendCardan(sat.wheels);
lgd = legend;

assert(isequal(lgd.String,{'Pitch'}));
close(f);

% Actuator on roll angle
sat.wheels = sat.wheels.setAxis(3);
f = figure('Visible','off');
plot(1,1);
tools.setLegendCardan(sat.wheels);
lgd = legend;

assert(isequal(lgd.String,{'Roll'}));
close(f);

%% Test 8.2:Set legend Cardan, 2 actuators

% Input data:
wheel = reactionWheel();

sat = satellite();

sat.wheels = [wheel wheel];

% Test:
% Actuators on yaw and pitch angles
sat.wheels(1) = sat.wheels(1).setAxis(1);
sat.wheels(2) = sat.wheels(2).setAxis(2);
f = figure('Visible','off');
plot(0,[0;1]);
tools.setLegendCardan(sat.wheels);
lgd = legend;

label = cellstr(lgd.String);
assert(isequal(label,{'Yaw','Pitch'}));
close(f);

% Actuators on yaw and roll angles
sat.wheels(1) = sat.wheels(1).setAxis(1);
sat.wheels(2) = sat.wheels(2).setAxis(3);
f = figure('Visible','off');
plot(0,[0;1]);
tools.setLegendCardan(sat.wheels);
lgd = legend;

assert(isequal(cellstr(lgd.String),{'Yaw','Roll'}));
close(f);

% Actuators on pitch and roll angles
sat.wheels(1) = sat.wheels(1).setAxis(2);
sat.wheels(2) = sat.wheels(2).setAxis(3);
f = figure('Visible','off');
plot(0,[0;1]);
tools.setLegendCardan(sat.wheels);
lgd = legend;

assert(isequal(cellstr(lgd.String),{'Pitch','Roll'}));
close(f);

%% Test 8.3: Set legend Cardan, 3 actuators

% Input data:
wheel = reactionWheel();

sat = satellite();

sat.wheels = [wheel wheel wheel];

sat.wheels(1) = sat.wheels(1).setAxis(1);
sat.wheels(2) = sat.wheels(2).setAxis(2);
sat.wheels(3) = sat.wheels(3).setAxis(3);

% Test:
f = figure('Visible','off');
plot(0,[0;1;2]);
tools.setLegendCardan(sat.wheels);
lgd = legend;

assert(isequal(cellstr(lgd.String),{'Yaw','Pitch','Roll'}));
close(f);