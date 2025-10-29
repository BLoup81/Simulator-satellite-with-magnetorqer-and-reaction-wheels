classdef magnetorquer
    properties (Access = private)
        axis {mustBeNumeric};       % Axis of the magnetorquer
        magneticMomentBound {mustBeNumeric};   % Saturation of magnetic moment [A.m^2]
    end

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% Constructor %%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = magnetorquer()
            obj.magneticMomentBound = 1e30;
            obj.axis = 1;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%% Setter %%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = setAxis(obj,axis)
            if ((axis ~= 1) && (axis ~= 2) && (axis ~= 3))
                throw(MException('magnetorquer:BadInput','Magnetorquer axis must be equal to 1, 2 or 3'));
            end
            obj.axis = axis;
        end

        function obj = setMagneticMomentBound(obj,saturation)
            if saturation <= 0
                throw(MException('magnetorquer:BadInput','Magnetic saturation must be positive'));
            end
            obj.magneticMomentBound = saturation;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%% Getter %%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function axis = getAxis(obj)
            axis = obj.axis;
        end

        function magneticMomentBound = getMagneticMomentBound(obj)
            magneticMomentBound = obj.magneticMomentBound;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% Saturation %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function torque = magneticMomentSaturation(obj,command)
            if abs(command) > obj.magneticMomentBound
                torque = sign(command)*obj.magneticMomentBound;
            else
                torque = command;
            end
        end
    end
end