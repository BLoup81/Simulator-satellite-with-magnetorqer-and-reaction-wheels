classdef reactionWheel
    properties (Access = private)
        inertia {mustBeNumeric};    % Inertia of the wheel                          [kg.m^2]
        axis {mustBeNumeric};       % Rotation axis of the wheel                    [1, 2, 3]
        velocityBound {mustBeNumeric}   % Maximum velocity of the reaction wheel    [rad/s]
        torqueBound {mustBeNumeric} % Maximum torque created by reaction wheel      [N.m]
        velocity {mustBeNumeric}    % Angular velocity of the reaction wheel        [rad/s]
    end

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%% Constructor %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = reactionWheel()
            obj.inertia = 1;
            obj.velocityBound = 1e30;
            obj.torqueBound = 1e30;
            obj.velocity = 0;
            obj.axis = 1;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%% Setter %%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = setInertia(obj,inertia)
            if inertia <= 0
                throw(MException('reactionWheel:BadInput','Reaction wheel inertia must be positive'));
            end
            obj.inertia = inertia;
        end

        function obj = setAxis(obj,axis)
            if ((axis ~= 1) && (axis ~= 2) && (axis ~= 3))
                throw(MException('reactionWheel:BadInput','Reaction wheel axis must be equal to 1, 2 or 3'));
            end
            obj.axis = axis;
        end

        function obj = setVelocityBound(obj,velocityBound)
            if velocityBound <= 0
                throw(MException('reactionWheel:BadInput','Velocity bound of the reaction wheel must be positive'));
            end
            obj.velocityBound = velocityBound;
        end

        function obj = setTorqueBound(obj,torqueBound)
            if torqueBound <= 0
                throw(MException('reactionWheel:BadInput','Torque bound of the reaction wheel must be positive'));
            end
            obj.torqueBound = torqueBound;
        end

        function obj = setVelocity(obj,velocity)
            obj.velocity = velocity;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%% Getter %%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function inertia = getInertia(obj)
            inertia = obj.inertia;
        end

        function axis = getAxis(obj)
            axis = obj.axis;
        end

        function velocityBound = getVelocityBound(obj)
            velocityBound = obj.velocityBound;
        end

        function torqueBound = getTorqueBound(obj)
            torqueBound = obj.torqueBound;
        end

        function velocity = getVelocity(obj)
            velocity = obj.velocity;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%% Saturation %%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function torque = velocitySaturation(obj,wheelVelocity,command,step)

            if abs(wheelVelocity - step*command/obj.inertia) > obj.velocityBound
                torque = -obj.inertia*(sign(-command)*obj.velocityBound - wheelVelocity)/step;
            else
                torque = command;
            end
        end

        function torque = torqueSaturation(obj,command)
            if abs(command) > obj.getTorqueBound
                torque = sign(command)*obj.torqueBound;
            else
                torque = command;
            end
        end
    end
end