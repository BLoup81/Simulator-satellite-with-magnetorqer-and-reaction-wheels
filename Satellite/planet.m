classdef planet
    properties (Access = private)
        radius {mustBeNumeric};                      % Radius of the planet      [km]
        gravitationalParameter {mustBeNumeric};      % Gravitational parameter   [km^3/s^2]
        angularVelocity {mustBeNumeric};             % Speed of rotation         [rad/s]
    end

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%% Constructor %%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = planet()
            obj.radius = 6378;
            obj.gravitationalParameter = 398600.44;
            obj.angularVelocity = 2*pi/86400;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% Setter %%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = setRadius(obj,radius)
            if radius <= 0
                throw(MException('planet:BadInput','Radius of the planet must be positive'));
            end
            obj.radius = radius;
        end

        function obj = setGravitationalParameter(obj,gravitationalParameter)
            if gravitationalParameter <= 0
                throw(MException('planet:BadInput','Gravitational parameter must be positive'));
            end
            obj.gravitationalParameter = gravitationalParameter;
        end

        function obj = setRotation(obj,rotationSpeed)
            obj.angularVelocity = rotationSpeed;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%% Getter %%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function radius = getRadius(obj)
            radius = obj.radius;
        end

        function gravitationalParameter = getGravitationalParameter(obj)
            gravitationalParameter = obj.gravitationalParameter;
        end

        function rotationSpeed = getRotation(obj)
            rotationSpeed = obj.angularVelocity;
        end
    end
end