% Tests of the class orbit

% To run tests:
%       runtests('orbitTest')

%% Test 0: Constructor

orb = orbit();

assert(isequal(orb.getInclination(),0));
assert(isequal(orb.getRAAN(),0));
assert(isequal(orb.getArgumentOfPeriapsis(),0));
assert(isequal(orb.getMajorSemiAxis(),0));
assert(isequal(orb.getEccentricity(),0));
assert(isequal(orb.getDateInitial(),[2000 01 01 00 00 00]));
assert(isequal(orb.getGHA(),0));
assert(isequal(orb.getInitialTrueAnomaly(),0));

%% Test 1.1: Getter/Setter inclination

% Input data:
inclination = 45;

% Expected output data:
inclination_expect = 45;
inclination_rad_expect = pi*inclination_expect/180;

% Test:
orb = orbit();

orb = orb.setInclination(inclination);
assert(isequal(orb.getInclination(),inclination_expect));
assert(isequal(orb.getInclination('deg'),inclination_expect));
assert(isequal(orb.getInclination('rad'),inclination_rad_expect));

%% Test 1.2: Getter/Setter inclination exception

orb = orbit();

try
    orb.setInclination(-91);
    error('Test failed: Inclination must be in degree between [-90,90]');
catch ME
    assert(strcmp(ME.identifier,'orbit:BadInput'),'Inclination must be in degree between [-90,90]');
end

try
    orb.setInclination(91);
    error('Test failed: Inclination must be in degree between [-90,90]');
catch ME
    assert(strcmp(ME.identifier,'orbit:BadInput'),'Inclination must be in degree between [-90,90]');
end

%% Test 2.1: Getter/Setter RAAN

% Input data:
raan = 90;

% Expected output data:
raan_expect = 90;
raan_rad_expect = pi*raan_expect/180;

% Test:
orb = orbit();

orb = orb.setRAAN(raan);
assert(isequal(orb.getRAAN(),raan_expect));
assert(isequal(orb.getRAAN('deg'),raan_expect));
assert(isequal(orb.getRAAN('rad'),raan_rad_expect));

%% Test 2.2: Getter/Setter RAAN exception

orb = orbit();

try
    orb.setRAAN(-1);
    error('Test failed: RAAN must be in degree between [0,360[');
catch ME
    assert(strcmp(ME.identifier,'orbit:BadInput'),'RAAN must be in degree between [0,360[');
end

try
    orb.setRAAN(361);
    error('Test failed: RAAN must be in degree between [0,360[');
catch ME
    assert(strcmp(ME.identifier,'orbit:BadInput'),'RAAN must be in degree between [0,360[');
end

%% Test 3.1: Getter/Setter periapsis argument

% Input data:
argumentOfPerigee = 90;

% Expected output data:
argumentOfPerigee_expect = 90;
argumentOfPerigee_rad_expect = pi*argumentOfPerigee_expect/180;

% Test:
orb = orbit();

orb = orb.setArgumentOfPeriapsis(argumentOfPerigee);
assert(isequal(orb.getArgumentOfPeriapsis(),argumentOfPerigee_expect));
assert(isequal(orb.getArgumentOfPeriapsis('deg'),argumentOfPerigee_expect));
assert(isequal(orb.getArgumentOfPeriapsis('rad'),argumentOfPerigee_rad_expect));

%% Test 3.2: Getter/Setter, Argument of periapsis exception

orb = orbit();

try
    orb.setArgumentOfPeriapsis(-1);
    error('Test failed: Argument of periapsis must be in degree between [0,360[');
catch ME
    assert(strcmp(ME.identifier,'orbit:BadInput'),'Argument of periapsis must be in degree between [0,360[');
end

try
    orb.setArgumentOfPeriapsis(361);
    error('Test failed: Argument of periapsis must be in degree between [0,360[');
catch ME
    assert(strcmp(ME.identifier,'orbit:BadInput'),'Argument of periapsis must be in degree between [0,360[');
end

%% Test 4: Getter/Setter major semi-axis

% Input data:
majorSemiAxis = 9856;

% Expected output data:
majorSemiAxis_expect = 9856;

% Test:
orb = orbit();

orb = orb.setMajorSemiAxis(majorSemiAxis);
assert(isequal(orb.getMajorSemiAxis(),majorSemiAxis_expect));

%% Test 5: Getter/Setter eccentricity

% Input data:
eccentricity = 0.05;

% Expected output data:
eccentricity_expect = 0.05;

% Test:
orb = orbit();

orb = orb.setEccentricity(eccentricity);

eccentricity_test = orb.getEccentricity();

assert(isequal(eccentricity_test,eccentricity_expect));

%% Test 6.1: Getter/Setter GHA

% Input data:
gha = 64;

% Expected output data:
gha_expect = 64;
gha_rad_expect = pi*gha_expect/180;

% Test:
orb = orbit();

orb = orb.setGHA(gha);
assert(isequal(orb.getGHA(),gha_expect));
assert(isequal(orb.getGHA('deg'),gha_expect));
assert(isequal(orb.getGHA('rad'),gha_rad_expect));

%% Test 6.2: Getter/Setter GHA exception

orb = orbit();

try
    orb.setInitialTrueAnomaly(-1);
    error('Test failed: GHA must be in degree between [0,360[');
catch ME
    assert(strcmp(ME.identifier,'orbit:BadInput'),'GHA must be in degree between [0,360[');
end

try
    orb.setGHA(361);
    error('Test failed: GHA must be in degree between [0,360[');
catch ME
    assert(strcmp(ME.identifier,'orbit:BadInput'),'GHA must be in degree between [0,360[');
end

%% Test 7.1: Getter/Setter initial true anomaly

% Input data:
initialTrueAnomaly = 90;

% Expected output data:
initialTrueAnomaly_expect = 90;
initialTrueAnomaly_rad_expect = pi*initialTrueAnomaly_expect/180;

% Test:
orb = orbit();

orb = orb.setInitialTrueAnomaly(initialTrueAnomaly);
assert(isequal(orb.getInitialTrueAnomaly(),initialTrueAnomaly_expect));
assert(isequal(orb.getInitialTrueAnomaly('deg'),initialTrueAnomaly_expect));
assert(isequal(orb.getInitialTrueAnomaly('rad'),initialTrueAnomaly_rad_expect));

%% Test 7.2: Getter/Setter, Initial true anomaly exception

orb = orbit();

try
    orb.setRAAN(-1);
    error('Test failed: Initial true anomaly must be in degree between [0,360[');
catch ME
    assert(strcmp(ME.identifier,'orbit:BadInput'),'Initial true anomaly must be in degree between [0,360[');
end

try
    orb.setRAAN(361);
    error('Test failed: Initial true anomaly must be in degree between [0,360[');
catch ME
    assert(strcmp(ME.identifier,'orbit:BadInput'),'Initial true anomaly must be in degree between [0,360[');
end

%% Test 8: Getter/Setter initial date

% Input data:
initialDate = [2025 05 02 14 21 32];
orb = orbit();

% Expect output data:
initialDate_expect = [2025 05 02 14 21 32];

% Test:
orb = orb.setDateInitial(initialDate);
err = 1e-8;

assert(norm(orb.getDateInitial() - initialDate_expect) <= err);

%% Test 9: Getter/Setter attractor body

% Input data:
gravitationnalParameter = 398600.44;
radius = 6378;
rotationRate = 0.001;

earth = planet();

orb = orbit();
orb.attractorBody = earth;

orb.attractorBody = orb.attractorBody.setRadius(radius);
orb.attractorBody = orb.attractorBody.setGravitationalParameter(gravitationnalParameter);
orb.attractorBody = orb.attractorBody.setRotation(rotationRate);

% Expected output data:
gravitationnalParameter_expect = 398600.44;
radius_expect = 6378;
rotationRate_expect = 0.001;

% Test:
earth_test = orb.attractorBody;
gravitationnalParameter_test = earth_test.getGravitationalParameter();
radius_test = earth_test.getRadius();
rotationRate_test = earth_test.getRotation();

assert(isequal(rotationRate_expect,rotationRate_test));
assert(isequal(radius_expect,radius_test));
assert(isequal(gravitationnalParameter_expect,gravitationnalParameter_test));

%% Test 10: Mean motion

% Input data:
gravitationalParameter = 398600.44;
majorSemiAxis = 6978;

earth = planet();
earth = earth.setGravitationalParameter(gravitationalParameter);

orb = orbit();
orb = orb.setMajorSemiAxis(majorSemiAxis);
orb.attractorBody = earth;

% Expected output data:
meanMotion_expect = 1.083109685e-3;

% Test:
meanMotion_test = orb.meanMotion();
err = 1e-10;

assert(abs(meanMotion_test - meanMotion_expect) <= err);

%% Test 11: Semi-latus rectus

% Input data:
majorSemiAxis = 6978;
eccentricity = 0.08918;

orb = orbit();
orb = orb.setMajorSemiAxis(majorSemiAxis);
orb = orb.setEccentricity(eccentricity);

% Expected output data:
semiLatusRectus_expect = 6922.503461;

% Test:
semiLatusRectus_test = orb.semi_latus_rectus();
err = 1e-6;

assert(abs(semiLatusRectus_test - semiLatusRectus_expect) <= err);

%% Test 12: Initial eccentric anomaly

% Input data:
eccentricity = 0.0818;
initialTrueAnomaly = 30;

orb = orbit();
orb = orb.setInitialTrueAnomaly(initialTrueAnomaly);
orb = orb.setEccentricity(eccentricity);

% Expected output data:
initialEccentricAnomaly_expect = 0.48403905;

% Test:
initialEccentricAnomaly_test = orb.initialEccentricAnomaly();
err = 1e-7;

assert(abs(initialEccentricAnomaly_test - initialEccentricAnomaly_expect) <= err);

%% Test 13: Mean anomaly

% Input data:
gravitationalParameter = 398600.44;

majorSemiAxis = 6978;
eccentricity = 0.0818;
initialTrueAnomaly = 30;
time = [0;0.5;10;1000];

earth = planet();
earth = earth.setGravitationalParameter(gravitationalParameter);

orb = orbit();
orb = orb.setInitialTrueAnomaly(initialTrueAnomaly);
orb = orb.setMajorSemiAxis(majorSemiAxis);
orb = orb.setEccentricity(eccentricity);
orb.attractorBody = earth;

% Expected output data:
meanAnomaly_expect = [0.44597276;0.44651432;0.45680386;1.529082449];

% Test:
meanAnomaly_test = orb.meanAnomaly(time);
err = 1e-8;

assert(norm(meanAnomaly_test - meanAnomaly_expect) <= err);

%% Test 14: Eccentric anomaly

% Input data:
gravitationalParameter = 398600.44;
majorSemiAxis = 6978;
eccentricity = 0.0818;
initialTrueAnomaly = 30;
time = [0;0.5;10;1000];

earth = planet();
earth = earth.setGravitationalParameter(gravitationalParameter);

orb = orbit();
orb = orb.setEccentricity(eccentricity);
orb = orb.setMajorSemiAxis(majorSemiAxis);
orb = orb.setInitialTrueAnomaly(initialTrueAnomaly);

orb.attractorBody = earth;

% Expected output data:
eccentricAnomaly_expect = [0.484039045;0.484622869;0.49571274;1.61081695];

% Test:
eccentricAnomaly_test = orb.eccentricAnomaly(time);

err = 1e-4;
assert(norm(eccentricAnomaly_test - eccentricAnomaly_expect) <= err);

%% Test 15: True anomaly

% Input data:
gravitationalParameter = 398600.44;
majorSemiAxis = 6978;
eccentricity = 0.0818;
initialTrueAnomaly = 30;
time = [0;0.5;10;1000];

earth = planet();
earth = earth.setGravitationalParameter(gravitationalParameter);

orb = orbit();
orb = orb.setEccentricity(eccentricity);
orb = orb.setMajorSemiAxis(majorSemiAxis);
orb = orb.setInitialTrueAnomaly(initialTrueAnomaly);
orb.attractorBody = earth;

% Expected output data:
trueAnomaly_expect = [0.52359877;0.52422605;0.53613845;1.69250919];

% Test:
trueAnomaly_test = orb.trueAnomaly(time);

err = 1e-4;
assert(norm(trueAnomaly_test - trueAnomaly_expect) <= err);

%% Test 16: Radius

% Input data:
majorSemiAxis = 6978;
eccentricity = 0.0818;
initialTrueAnomaly = 30;
time = [0;0.5;10;1000];
trueAnomaly = [0.523598775514615;0.524226046962115;0.536138453016618;1.692509190590019];

orb = orbit();
orb = orb.setEccentricity(eccentricity);
orb = orb.setMajorSemiAxis(majorSemiAxis);
orb = orb.setInitialTrueAnomaly(initialTrueAnomaly);

% Expected output data:
radius_expect = [6472.771696;6472.926862;6475.906894;7000.837688];

% Test:
radius_test = orb.radius(time,trueAnomaly);

err = 1e-5;
assert(norm(radius_test - radius_expect) <= err);

%% Test 17: ECI coordinates

% Input data:
gravitationalParameter = 398600.44;
majorSemiAxis = 6978;
eccentricity = 0.0818;
initialTrueAnomaly = 30;
inclination = 60;
raan = 45;
argumentOfPerigee = 28.64788976;
time = [0;0.5;10;1000];
trueAnomaly = [0.523598775514615;0.524226046962115;0.536138453016618;1.692509190590019];
radius = [6472.771696;6472.926862;6475.906894;7000.837688];

earth = planet();
earth = earth.setGravitationalParameter(gravitationalParameter);

orb = orbit();
orb = orb.setEccentricity(eccentricity);
orb = orb.setMajorSemiAxis(majorSemiAxis);
orb = orb.setInitialTrueAnomaly(initialTrueAnomaly);
orb = orb.setArgumentOfPeriapsis(argumentOfPerigee);
orb = orb.setInclination(inclination);
orb = orb.setRAAN(raan);

orb.attractorBody = earth;

% Expected output data:
x_eci_expect = [427.042492 423.853864 363.242544 -4895.23948];
y_eci_expect = [4335.68570 4334.08382 4303.34675 -871.197209];
z_eci_expect = [4787.09073 4789.03409 4825.62242 4928.42513];

% Test:
eci_test = orb.coordinateECI(time,trueAnomaly,radius);

err = 1e-4;
assert(norm(x_eci_expect - eci_test(1,:)) <= err);
assert(norm(y_eci_expect - eci_test(2,:)) <= err);
assert(norm(z_eci_expect - eci_test(3,:)) <= err);

%% Test 18: Calcul GHA

% Input data:
rotationRate = 2*pi/86164;
initialGHA = 50;
time = [-1e5;0;10;1e5];

earth = planet();
earth = earth.setRotation(rotationRate);

orb = orbit();
orb = orb.setGHA(initialGHA);
orb.attractorBody = earth;

% Expected output data:
gha_expect = [6.14691172;5*pi/18;0.87339384;1.881602833];

% Test:
gha_test = orb.calculGHA(time);

err = 1e-8;
assert(norm(gha_expect - gha_test) <= err);

%% Test 19: Geocentric coordinates

% Input data:
gravitationalParameter = 398600.44;
majorSemiAxis = 6978;
eccentricity = 0.0818;
initialTrueAnomaly = 30;
inclination = 60;
raan = 45;
argumentOfPerigee = 28.64788976;
initialGHA = 50;
rotationRate = 2*pi/86164;
time = [0;0.5;10;1000];

orb = orbit();
orb = orb.setGHA(initialGHA);
orb = orb.setInclination(inclination);
orb = orb.setMajorSemiAxis(majorSemiAxis);
orb = orb.setEccentricity(eccentricity);
orb = orb.setRAAN(raan);
orb = orb.setInitialTrueAnomaly(initialTrueAnomaly);
orb = orb.setArgumentOfPeriapsis(argumentOfPerigee);

earth = planet();
earth = earth.setRotation(rotationRate);
earth = earth.setGravitationalParameter(gravitationalParameter);

orb.attractorBody = earth;

% Expected output data:
lat_expect = [47.695029750951214 47.719082959391052 48.173246265811919 44.746840375210070];
lon_expect = [34.3746897828045 34.4122851507153 35.1332399841788 -224.0867978736240];
alt_expect = [6472.771236103315 6472.926399526058 6475.906426237519 7000.838744468322];

% Test:
[lat_test,lon_test,alt_test] = orb.geocentric(time);

err = 1e-8;
assert(norm(lat_expect - lat_test) <= err);
assert(norm(lon_expect - lon_test) <= err);
assert(norm(alt_expect - alt_test) <= err);