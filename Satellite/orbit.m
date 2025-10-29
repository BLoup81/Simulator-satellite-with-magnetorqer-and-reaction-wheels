classdef orbit
    properties (Access = private)
        inclination {mustBeNumeric};                % Inclination                           [-90,90]    [deg]
        raan {mustBeNumeric};                       % Right Ascension of Ascendent Node     [0,360[     [deg]
        argumentOfPerigee {mustBeNumeric};          % Argument of periapsis                 [0,360[     [deg]
        majorSemiAxis {mustBeNumeric};               % Major semi-axis                       (> 6378)    [km]
        eccentricity {mustBeNumeric};               % Eccentricity                          (> 1-radius/majorSemiAxis)
        initialGHA {mustBeNumeric};                 % Initial Greenwich Hour Angle          [0,360[     [deg]
        initialTrueAnomaly {mustBeNumeric};         % Initial true anomaly                  [0,360[     [deg]
        initialDate;                                % Initial absolute date                 initialDate = [yyyy mm dd hh mm ss]
        attractorBody_ planet = planet.empty;       % Massive attractor body like the Earth
    end

    properties (Dependent)
        attractorBody;
    end

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%% Constructor %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = orbit()
            obj.inclination = 0;
            obj.raan = 0;
            obj.argumentOfPerigee = 0;
            obj.majorSemiAxis = 0;
            obj.eccentricity = 0;
            obj.initialGHA = 0;
            obj.initialTrueAnomaly = 0;
            obj.initialDate = [2000 01 01 00 00 00];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%% Setter %%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = setInclination(obj,inclination)
            if (inclination < -90) || (inclination > 90)
                throw(MException('orbit:BadInput','Inclination must be in degree between [-90,90]'));
            end
            obj.inclination = pi*inclination/180;
        end

        function obj = setRAAN(obj,raan)
            if (raan < 0) || (raan >= 360)
                throw(MException('orbit:BadInput','RAAN must be in degree between [0,360['));
            end
            obj.raan = pi*raan/180;
        end

        function obj = setArgumentOfPeriapsis(obj,periapsisArgument)
            if (periapsisArgument < 0) || (periapsisArgument >= 360)
                throw(MException('orbit:BadInput','Argument of the periapsis must be in degree between [0,360['));
            end
            obj.argumentOfPerigee = pi*periapsisArgument/180;
        end

        function obj = setMajorSemiAxis(obj,majorSemiAxis)
            obj.majorSemiAxis = majorSemiAxis;
        end

        function obj = setEccentricity(obj,eccentricity)
            obj.eccentricity = eccentricity;
        end

        function obj = setGHA(obj,gha)
            if (gha < 0) || (gha >= 360)
                throw(MException('orbit:BadInput','GHA must be in degree between [0,360['));
            end
            obj.initialGHA = pi*gha/180;
        end

        function obj = setInitialTrueAnomaly(obj,anomaly)
            if (anomaly < 0) || (anomaly >= 360)
                throw(MException('orbit:BadInput','Initial true anomaly must be in degree between [0,360['));
            end
            obj.initialTrueAnomaly = pi*anomaly/180;
        end

        function obj = setDateInitial(obj,date)
            obj.initialDate = [date(1) date(2) date(3) date(4) date(5) date(6)];
        end

        function obj = set.attractorBody(obj,planet)
            if isa(planet, "planet")
                obj.attractorBody_ = planet;
            else
                error("Attractor body must be a planet instance");
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%% Getter %%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function inclination = getInclination(obj,unit)
            if ~exist("unit")
                inclination = 180*obj.inclination/pi;
            else if exist("unit")
                
                if (unit == 'rad')
                    inclination = obj.inclination;
                else if (unit == 'deg')
                    inclination = 180*obj.inclination/pi;
                else
                    throw(MException('orbit:BadInput','Unit angle in getInclination is not valid'));
                end
                end
            end
            end
        end

        function raan = getRAAN(obj,unit)
            if ~exist("unit")
                raan = 180*obj.raan/pi;
            else if exist("unit")
                
                if (unit == 'rad')
                    raan = obj.raan;
                else if (unit == 'deg')
                        raan = 180*obj.raan/pi;
                    else
                        throw(MException('orbit:BadInput','Unit angle in getRAAN is not valid'));
                    end
                end
                end
            end
        end

        function periapsisArgument = getArgumentOfPeriapsis(obj,unit)
            if ~exist("unit")
                periapsisArgument = 180*obj.argumentOfPerigee/pi;
            else if exist("unit")
                
                if (unit == 'rad')
                    periapsisArgument = obj.argumentOfPerigee;
                else if (unit == 'deg')
                    periapsisArgument = 180*obj.argumentOfPerigee/pi;
                    else
                    throw('orbit:BadInput','Unit angle in getArgumentOfPeriapsis is not valid');
                    end
                end
                end
            end
        end

        function majorSemiAxis = getMajorSemiAxis(obj)
            majorSemiAxis = obj.majorSemiAxis;
        end

        function eccentricity = getEccentricity(obj)
            eccentricity = obj.eccentricity;
        end

        function gha = getGHA(obj,unit)
            if ~exist("unit")
                gha = 180*obj.initialGHA/pi;
            else if exist("unit")
                
                if (unit == 'rad')
                    gha = obj.initialGHA;
                else if (unit == 'deg')
                    gha = 180*obj.initialGHA/pi;
                    else
                    throw(MException('orbit:BadInput','Unit angle in getGHA is not valid'));
                    end
                end
                end
            end
        end

        function anomaly = getInitialTrueAnomaly(obj,unit)
            if ~exist("unit")
                anomaly = 180*obj.initialTrueAnomaly/pi;
            else if exist("unit")
                
                if (unit == 'rad')
                    anomaly = obj.initialTrueAnomaly;
                else if (unit == 'deg')
                    anomaly = 180*obj.initialTrueAnomaly/pi;
                    else
                    throw(MException('orbit:BadInput','Unit angle in getInitialTrueAnomaly is not valid'));
                    end
                end
                end
            end
        end

        function date = getDateInitial(obj)
            date = obj.initialDate;
        end

        function planet = get.attractorBody(obj)
            planet = obj.attractorBody_;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function n = meanMotion(obj)
            n = sqrt(obj.attractorBody.getGravitationalParameter()/(obj.majorSemiAxis^3));
        end

        function p = semi_latus_rectus(obj)
            p = obj.majorSemiAxis*(1 - obj.eccentricity^2);
        end

        function ecc_an0 = initialEccentricAnomaly(obj)
            cosE0 = (obj.eccentricity + cos(obj.initialTrueAnomaly))/(1 + obj.eccentricity*cos(obj.initialTrueAnomaly));
            sinE0 = (sqrt(1 - (obj.eccentricity^2))*sin(obj.initialTrueAnomaly))/(1 + obj.eccentricity*cos(obj.initialTrueAnomaly));
            ecc_an0 = atan2(sinE0,cosE0);

            if ecc_an0 < 0
                ecc_an0 = ecc_an0 + 2*pi;
            end
        end

        function meanAnomaly = meanAnomaly(obj,time)
            
            nb = length(time);

            meanAnomaly = zeros(nb,1);

            meanMotion = obj.meanMotion();

            initialEccentricAnomaly = obj.initialEccentricAnomaly();

            % Calcul of the passed time at the periapsis
            perigeeTime = (-initialEccentricAnomaly + obj.eccentricity*sin(initialEccentricAnomaly))/meanMotion;
            for index = 1:nb
                meanAnomaly(index) = meanMotion*(time(index) - perigeeTime);
            end
        end

        function eccentricAnomaly = eccentricAnomaly(obj,time)

            meanAnomaly = obj.meanAnomaly(time);
            err = 1e-5;
            eccentricAnomaly = meanAnomaly;
            
            % Use of Newton optimisation method
            for index = 1:length(time)
                f = eccentricAnomaly(index) - obj.eccentricity*sin(eccentricAnomaly(index)) - meanAnomaly(index);       
                while abs(f) > err
                    eccentricAnomaly(index) = eccentricAnomaly(index) - f/(1 + obj.eccentricity*cos(eccentricAnomaly(index)));
                    f = eccentricAnomaly(index) - obj.eccentricity*sin(eccentricAnomaly(index)) - meanAnomaly(index);
                end
            end
        end

        function trueAnomaly = trueAnomaly(obj,time)
            
            nb = length(time);

            trueAnomaly = zeros(nb,1);

            eccentricAnomaly = obj.eccentricAnomaly(time);

            for index = 1:nb
                sin_trueAnomaly = (sqrt(1 - obj.eccentricity^2)*sin(eccentricAnomaly(index)))/(1 - obj.eccentricity*cos(eccentricAnomaly(index)));
                cos_trueAnomaly = (cos(eccentricAnomaly(index)) - obj.eccentricity)/(1 - obj.eccentricity*cos(eccentricAnomaly(index)));
                trueAnomaly(index,1) = atan2(sin_trueAnomaly,cos_trueAnomaly);
            end
        end

        function r = radius(obj,time,trueAnomaly)

            nb = length(time);

            r = zeros(nb,1);

            for index = 1:nb
                r(index) = obj.majorSemiAxis*(1 - obj.eccentricity^2)/(1 + obj.eccentricity*cos(trueAnomaly(index)));
            end
        end

        function eci = coordinateECI(obj,time,trueAnomaly,radius)
            
            nb = length(time);

            eci = zeros(3,nb);

            for index = 1:nb
                theta = trueAnomaly(index) + obj.argumentOfPerigee;

                xp = radius(index)*cos(theta);
                yp = radius(index)*sin(theta);
    
                eci(1,index) = xp*cos(obj.raan) - yp*cos(obj.inclination)*sin(obj.raan);
                eci(2,index) = xp*sin(obj.raan) + yp*cos(obj.inclination)*cos(obj.raan);
                eci(3,index) = yp*sin(obj.inclination);
            end
        end

        function gha = calculGHA(obj,time)

            nb = length(time);

            gha = zeros(nb,1);

            for index = 1:nb
                gha(index) = obj.attractorBody.getRotation()*time(index) + obj.initialGHA;

                if gha(index) < -pi
                    gha(index) = gha(index) + ceil(gha(index)/(-2*pi))*2*pi;
                elseif gha(index) > pi
                        gha(index) = gha(index) - fix(gha(index)/(2*pi))*2*pi;
                end
            end
        end

        function [lat,lon,alt] = geocentric(obj,time)
            nb = length(time);

            lat = zeros(1,nb);
            lon = zeros(1,nb);
            alt = zeros(1,nb);

            trueAnomaly = obj.trueAnomaly(time);
            radius = obj.radius(time,trueAnomaly);
            eci = obj.coordinateECI(time,trueAnomaly,radius);
            gha = obj.calculGHA(time);
            
            
            for index = 1:nb
                theta = trueAnomaly(index) + obj.argumentOfPerigee;
                lat(index) = asin(sin(obj.inclination)*sin(theta))*180/pi;
                lon(index) = (atan2(eci(2,index)/radius(index),eci(1,index)/radius(index)) - gha(index))*180/pi;
                alt(index) = radius(index);
            end
        end
    end
end