classdef satellite
    properties (Access = private)
        length {mustBeNumeric};             % The length of the satellite               [m]
        width {mustBeNumeric};              % The width of the satellite                [m]
        height {mustBeNumeric};             % The height of the satellite               [m]
        mass {mustBeNumeric};               % The mass of the satellite                 [m]
        centerOfMass (3,1) {mustBeVector};  % The vector 3 of center of mass            [m]
        inertia (3,3) {mustBeMatrix};       % Inertia matrix 3 by 3 of the satellite    [kg.m^2]
        invInertia (3,3) {mustBeMatrix}     % Inverse of the inertia matrix
        quaternion (4,1) {mustBeVector};     % The vector 4 of the quaternion
        angularVelocity (3,1) {mustBeVector};   % The vector 3 of the angular velocity of the satellite [rad/s]
        wheels_ reactionWheel = reactionWheel.empty;                % Reaction wheels of the satellite
        magnetorquers_ magnetorquer = magnetorquer.empty;           % Magnetorquers of the satellite
    end

    properties (Dependent)
        wheels;
        magnetorquers;
    end

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% Constructor %%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = satellite()
            obj.inertia = eye(3,3);
            obj.invInertia = eye(3,3);
            obj.length = -1;
            obj.width = -1;
            obj.height = -1;
            obj.mass = -1;
            obj.centerOfMass = [-1;-1;-1];
            obj.quaternion = [-1;-1;-1;-1];
            obj.angularVelocity = [1e30;1e30;1e30];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%% Setter %%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = setLength(obj,length)
            if length <= 0      % Check if the length is positive
                throw(MException('satellite:BadInput','Length of satellite must be positive'));
            end
            obj.length = length;
        end

        function obj = setWidth(obj,width)
            if width <= 0       % Check if the width is positive
                throw(MException('satellite:BadInput','Width of satellite must be positive'));
            end
            obj.width = width;
        end

        function obj = setHeight(obj,height)
            if height <= 0      % Check if the height is positive
                throw(MException('satellite:BadInput','Height of satellite must be positive'));
            end
            obj.height = height;
        end

        function obj = setMass(obj,mass)
            if mass <= 0        % Check if the mass is positive
                throw(MException('satellite:BadInput','Mass of satellite must be positive'));
            end
            obj.mass = mass;
        end

        function obj = setCenterOfMass(obj,centerOfMass)
            dim = size(centerOfMass);
            if ~isequal(dim,[3 1])      % Check if dimensions of the center of mass vector are 3 by 1
                throw(MException('satellite:BadCenterOfMass','The center of mass vector must be 3 by 1'));
            end
            obj.centerOfMass = centerOfMass;
        end

        function obj = setInertia(obj,inertia)
            dim = size(inertia);
            if ~isequal(dim,[3 3])      % Check if dimensions of the inertia matrix are 3 by 3
                throw(MException('satellite:BadInertia', 'The inertia matrix must be 3 by 3'));
            end
            vp = eig(inertia);
            if (vp(1) <= 0) || (vp(2) <= 0) || (vp(3) <= 0)     % Check if the eigenvalues of the inertia matrix are positive
                throw(MException('satellite:BadInertia','The eigenvalues of inertia matrix must be positive'));
            end
            obj.inertia = inertia;
            obj.invInertia = inv(inertia);
        end

        function obj = setQuaternion(obj,quaternion)
            dim = size(quaternion);
            if ~isequal(dim,[4 1])      % Check if dimensions of the quaternion vector are 4 by 1
                throw(MException('satellite:BadInput','Quaternion must be 4 by 1'));
            end

            if abs(norm(quaternion) - 1) >= 1e-2        % Check if the norm of the quaternion is equal to 1
                throw(MException('satellite:BadInput','Norm of the quaternion must be equal to 1'));
            end
            obj.quaternion = quaternion;
        end

        function obj = setAngularVelocity(obj,angularVelocity)
            dim = size(angularVelocity);
            if ~isequal(dim,[3 1])      % Check if dimension of the angular velocity vector are 3 by 1
                throw(MException('satellite:BadInput','Angular velocity vector must be 3 by 1'));
            end
            obj.angularVelocity = angularVelocity;
        end

        function obj = set.wheels(obj,wheel)
            if isa(wheel, "reactionWheel")
                obj.wheels_ = wheel;
            else
                error("wheel must be a reactionWheel instance");
            end
        end

        function obj = addWheel(obj,wheel)
            if isa(wheel, "reactionWheel")
                obj.wheels_(end+1) = wheel;
            else
                error("addWheel: argument must be a reactionWheel instance");
            end
        end

        function obj = set.magnetorquers(obj,mgt)
            if isa(mgt, "magnetorquer")
                obj.magnetorquers_ = mgt;
            else
                error("mgt must be a magnetorquer instance");
            end
        end

        function obj = addMagnetorquer(obj,mgt)
            if isa(mgt, "magnetorquer")
                obj.magnetorquers_(end+1) = mgt;
            else
                error("addMagnetorquer: argument must be a magnetorquer instance");
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% Getter %%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function length = getLength(obj)
            length = obj.length;
        end

        function width = getWidth(obj)
            width = obj.width;
        end

        function height = getHeight(obj)
            height = obj.height;
        end

        function mass = getMass(obj)
            mass = obj.mass;
        end

        function centerOfMass = getCenterOfMass(obj)
            centerOfMass = obj.centerOfMass;
        end

        function inertia = getInertia(obj)
            inertia = obj.inertia;
        end

        function quaternion = getQuaternion(obj)
            quaternion = obj.quaternion;
        end

        function anglarVelocity = getAngularVelocity(obj)
            anglarVelocity = obj.angularVelocity;
        end

        function wheel = get.wheels(obj)
            wheel = obj.wheels_;
        end

        function mgt = get.magnetorquers(obj)
            mgt = obj.magnetorquers_;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%% Inertia matrix %%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = computeCenterOfMass(obj)
            if isequal(obj.length,-1) || isequal(obj.width,-1) || isequal(obj.height,-1)
                throw(MException('satellite:MissingParameter','Dimensions of the satellite must be defined to compute the center of mass'));
            end
            
            obj.centerOfMass = [obj.length/2; obj.width/2; obj.height/2];
        end

        function obj = computeInertia(obj)
            if isequal(obj.length,-1) || isequal(obj.width,-1) || isequal(obj.height,-1)
                throw(MException('satellite:MissingParameter','Dimensions of the satellite must be defined to compute the inertia matrix'));
            end

            if isequal(obj.centerOfMass(1),-1) || isequal(obj.centerOfMass(2),-1) || isequal(obj.centerOfMass(3),-1)
                throw(MException('satellite:MissingParameter','Center of mass must be defined'));
            end

            if isequal(obj.mass,-1)
                throw(MException('satellite:MissingParameter','Mass must be defined to compute the inertia matrix'));
            end

            a = obj.centerOfMass(1)/obj.length;
            b = obj.centerOfMass(2)/obj.width;
            c = obj.centerOfMass(3)/obj.height;

            % Calcul of diagonal elements
            P_a = 1 - 3*a + 3*a^2;
            P_b = 1 - 3*b + 3*b^2;
            P_c = 1 - 3*c + 3*c^2;

            obj.inertia(1,1) = obj.mass*(P_c*obj.height^2 + P_b*obj.width^2)/3;
            obj.inertia(2,2) = obj.mass*(P_a*obj.length^2 + P_c*obj.height^2)/3;
            obj.inertia(3,3) = obj.mass*(P_a*obj.length^2 + P_b*obj.width^2)/3;

            % Calcul of no diagonal elements
            obj.inertia(1,2) = obj.mass*obj.length*obj.width*(1-2*a)*(1-2*b)/4;
            obj.inertia(2,1) = obj.inertia(1,2);

            obj.inertia(1,3) = obj.mass*obj.length*obj.height*(1-2*a)*(1-2*c)/4;
            obj.inertia(3,1) = obj.inertia(1,3);

            obj.inertia(2,3) = obj.mass*obj.width*obj.height*(1-2*b)*(1-2*c)/4;
            obj.inertia(3,2) = obj.inertia(2,3);

            % Compute the inverse of the inertia
            obj.invInertia = inv(obj.inertia);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Derivative of states of the model %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function dq = quaternionDerivative(obj,state)

            q = state(1:4);         % Variable of the satellite quaternions
            w = state(5:7);         % Variable of the satellite angular velocity

            dq = zeros(4,1);        % Variable of the derivative of the quaternion
        
            dq(1) = -0.5*(q(2)*w(1) + q(3)*w(2) + q(4)*w(3));
            dq(2) = 0.5*(q(1)*w(1) + q(3)*w(3) - q(4)*w(2));
            dq(3) = 0.5*(q(1)*w(2) + q(4)*w(1) - q(2)*w(3));
            dq(4) = 0.5*(q(1)*w(3) + q(2)*w(2) - q(3)*w(1));
        end

        function dw = angularVelocityDerivative(obj,state,command)
            dw = obj.invInertia*(command - cross(state(5:7),obj.inertia*state(5:7)));
        end

        function obj = rearrangeReactionWheel(obj)
            % This rearrange the vector "wheels" to have the wheel on the
            % axis 1 on the first place of "wheel", the wheel on the axis 2
            % in second place and the wheel on the axis 3 in the third
            % place.
            vectWheels = obj.wheels;

            nb_wheels = numel(vectWheels);

            if isequal(nb_wheels,2)
                if vectWheels(1).getAxis() > vectWheels(2).getAxis()
                    obj.wheels(1) = vectWheels(2);
                    obj.wheels(2) = vectWheels(1);
                end
            elseif isequal(nb_wheels,3)
                for i=1:nb_wheels
                    if isequal(vectWheels(i).getAxis(),1)
                        obj.wheels(1) = vectWheels(i);
                    elseif isequal(vectWheels(i).getAxis(),2)
                        obj.wheels(2) = vectWheels(i);
                    else
                        obj.wheels(3) = vectWheels(i);
                    end
                end
            end
        end

        function dw = wheelVelocityDerivative(obj,state,command)

            nb_wheels = numel(obj.wheels);
            dw = zeros(nb_wheels,1);

            if isequal(nb_wheels,2) || isequal(nb_wheels,1)
                for i=1:nb_wheels
                    dw(i) = -command(i)/obj.wheels(i).getInertia();
                end
            elseif isequal(nb_wheels,3)
                dw(1) = (-command(1) - state(6)*obj.wheels(3).getInertia()*state(10) + state(7)*obj.wheels(2).getInertia()*state(9))/obj.wheels(1).getInertia();
                dw(2) = (-command(2) - state(7)*obj.wheels(1).getInertia()*state(8) + state(5)*obj.wheels(3).getInertia()*state(10))/obj.wheels(2).getInertia();
                dw(3) = (-command(3) - state(5)*obj.wheels(2).getInertia()*state(9) + state(6)*obj.wheels(1).getInertia()*state(8))/obj.wheels(3).getInertia(); 
            end
        end

        function obj = rearrangeMagnetorquers(obj)
            % This rearrange the vector "magnetorquers" to have the magnetorquer on the
            % axis 1 on the first place of "magnetorquers", the magnetorquer on the axis 2
            % in second place and the magnetorquer on the axis 3 in the third
            % place.
            vectMagnetorquers = obj.magnetorquers;

            nb_magnetorquers = numel(vectMagnetorquers);

            if isequal(nb_magnetorquers,2)
                if vectMagnetorquers(1).getAxis() > vectMagnetorquers(2).getAxis()
                    obj.magnetorquers(1) = vectMagnetorquers(2);
                    obj.magnetorquers(2) = vectMagnetorquers(1);
                end
            elseif isequal(nb_magnetorquers,3)
                for i=1:nb_magnetorquers
                    if isequal(vectMagnetorquers(i).getAxis(),1)
                        obj.magnetorquers(1) = vectMagnetorquers(i);
                    elseif isequal(vectMagnetorquers(i).getAxis(),2)
                        obj.magnetorquers(2) = vectMagnetorquers(i);
                    else
                        obj.magnetorquers(3) = vectMagnetorquers(i);
                    end
                end
            end
        end

        function derivative = dynamicModel(obj,state,command,derivativeAdditionState,varargin)
            nb_state = numel(state);
            nb_wheel = numel(obj.wheels);

            derivative = zeros(nb_state,1);

            derivative(1:4) = obj.quaternionDerivative(state);
            derivative(5:7) = obj.angularVelocityDerivative(state,command(1:3));
            if nb_wheel > 0         % There are 1 or more reaction wheels
                derivative(8:7+nb_wheel) = obj.wheelVelocityDerivative(state,command(4:end));
            end
            if nb_state - nb_wheel > 7
                derivative(8+nb_wheel:end) = derivativeAdditionState(state,varargin{:});
            end
        end



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%% RK4 %%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [state_next,U_compute,U_effective,M_compute,M_effective] = RK4(obj,state,command,affectationCommand,commandOption,step,magneticField,derivativeAddState,varargin)

            nb_wheel = numel(obj.wheels);       % {1,2,3} reaction wheels
            n = 7 + nb_wheel;                   % 4 (quaternion) + 3 (angular velocity) + nb_wheel
            nb_addState = numel(state) - n;
        
            % Creation of the state at the next time
            state_next = zeros(nb_addState+n,1);
            
            % Creation of the method's variables
            k = zeros(nb_addState+n,4);
        
            % Compute of the k1 for each state
                [U_compute1,U_effective1,M_compute1,M_effective1] = obj.appliedCommand(state,command,commandOption,step,magneticField);
                k(:,1) = obj.dynamicModel(state,affectationCommand*U_effective1,derivativeAddState,varargin{:});

            % Compute of the k2 for each state
                state2 = state + 0.5*step*k(:,1);
                [U_compute2,U_effective2,M_compute2,M_effective2] = obj.appliedCommand(state2,command,commandOption,step,magneticField);
                k(:,2) = obj.dynamicModel(state2,affectationCommand*U_effective2,derivativeAddState,varargin{:});

            % Compute of the k3 for each state
                state2 = state + 0.5*step*k(:,2);
                [U_compute3,U_effective3,M_compute3,M_effective3] = obj.appliedCommand(state2,command,commandOption,step,magneticField);
                k(:,3) = obj.dynamicModel(state2,affectationCommand*U_effective3,derivativeAddState,varargin{:});
        
            % Compute of the k4 for each state
                state2 = state + step*k(:,3);
                [U_compute4,U_effective4,M_compute4,M_effective4] = obj.appliedCommand(state2,command,commandOption,step,magneticField);
                k(:,4) = obj.dynamicModel(state2,affectationCommand*U_effective4,derivativeAddState,varargin{:});

            % Compute of the next state
                state_next(:,1) = state + step*(k(:,1) + 2*k(:,2) + 2*k(:,3) + k(:,4))/6;
                state_next(1:4,1) = state_next(1:4,1)/norm(state_next(1:4,1));
                U_compute = (U_compute1 + 2*U_compute2 + 2*U_compute3 + U_compute4)/6;
                U_effective = (U_effective1 + 2*U_effective2 + 2*U_effective3 + U_effective4)/6;
                M_compute = (M_compute1 + 2*M_compute2 + 2*M_compute3 + M_compute4)/6;
                M_effective = (M_effective1 + 2*M_effective2 + 2*M_effective3 + M_effective4)/6;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%% Saturation %%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function saturatedCommand = wheelSaturation(obj,wheelVelocity,command,step)
            nb_wheel = numel(obj.wheels);

            saturatedCommand = zeros(nb_wheel,1);

            for i=1:nb_wheel
                saturatedCommand(i) = obj.wheels(i).velocitySaturation(wheelVelocity(i),command(i),step);
                saturatedCommand(i) = obj.wheels(i).torqueSaturation(saturatedCommand(i));
            end
        end

        function saturatedCommand = magnetorquerSaturation(obj,command)
            nb_magnetorquer = numel(obj.magnetorquers);

            saturatedCommand = command;

            for i=1:nb_magnetorquer
                saturatedCommand(obj.magnetorquers(i).getAxis()) = obj.magnetorquers(i).magneticMomentSaturation(command(obj.magnetorquers(i).getAxis()));
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%% Magnetic moment %%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function magneticMoment = allocationMagneticMoment(obj,magneticField,requiredTorque)

            nrmMagField = norm(magneticField);

            magneticMoment = cross(magneticField,requiredTorque);

            magneticMoment = magneticMoment/(nrmMagField^2);
        end

        function torque = applicationMagneticMoment(obj,magneticField,magneticMoment)
            magneticTorque = cross(magneticMoment,magneticField);

            nb_magnetorquer = numel(obj.magnetorquers);
            torque = zeros(nb_magnetorquer,1);

            if nb_magnetorquer == 1
                torque = magneticTorque(obj.magnetorquers.getAxis());
            elseif nb_magnetorquer == 2
                for i=1:2
                    torque(i) = magneticTorque(obj.magnetorquers(i).getAxis());
                end
            else
                torque = magneticTorque;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% Applied command %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [U_compute,U_effective,M_compute,M_saturated] = appliedCommand(obj,state,command,commandOption,step,magneticField)

            nb_wheel = numel(obj.wheels);
            nb_magnetorquer = numel(obj.magnetorquers);

            if nb_wheel + nb_magnetorquer > 0
                U_effective = zeros(nb_magnetorquer+nb_wheel,1);
            else
                U_effective = 0;
            end

            M_compute = zeros(3,1);
            M_saturated = zeros(3,1);

            % If commandOption is "cardan", conversion of the quaternion to Cardan's angles
            if commandOption == "cardan"
                state4command = [tools.quater2euler(state(1:4));state(5:end)];      % Conversion of quaternion to Cardan's angles
            else
                state4command = state;
            end

            % Control compute by the function command
            U_compute = command(state4command);

            % Saturation of the angular velocity of the wheel if 1 is there or more
            if nb_wheel > 0
                U_effective(1:nb_wheel) = obj.wheelSaturation(state(7+1:7+nb_wheel),U_compute(1:nb_wheel),step);
            end

            U_magnetorquer = zeros(3,1);
            for i=1:nb_magnetorquer
                U_magnetorquer(obj.magnetorquers(i).getAxis()) = U_compute(nb_wheel+i);
            end

            % Allocation, saturation and application of the magnetic moment
            if nb_magnetorquer > 0
                magneticField = tools.frame2brf(magneticField,state(1:4));      % Conversion of the magnetic fielf from the referenced frame to the BRF
                M_compute = obj.allocationMagneticMoment(magneticField,U_magnetorquer);  % Allocation
                M_saturated = obj.magnetorquerSaturation(M_compute);    % Saturation
                U_effective(nb_wheel+1:end) = obj.applicationMagneticMoment(magneticField,M_saturated);    % Application
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%% ODE %%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function mat = affectationCommand(obj)
            % This function create a matrix to sum the differente command.
            % Exemple: if there are 1 reaction wheel on each axis and 1
            % magnetorquer on the second axis, the matrix is:
            % angular velocity 1 [1 0 0 0]
            % angular velocity 2 [0 1 0 1]
            % angular velocity 3 [0 0 1 0]
            %   wheel velocity 1 [1 0 0 0] * [RW_1 RW_2 RW_3 MGT_2]^T
            %   wheel velocity 2 [0 1 0 0]
            %   wheel velocity 3 [0 0 1 0]

            nb_wheels = numel(obj.wheels);

            nb_magnetorquers = numel(obj.magnetorquers);

            mat = zeros(3 + nb_wheels,nb_magnetorquers + nb_wheels);

            for i=1:nb_wheels
                mat(obj.wheels(i).getAxis(),i) = 1;
                mat(3 + i,i) = 1;
            end

            for i=1:nb_magnetorquers
                mat(obj.magnetorquers(i).getAxis(),nb_wheels + i) = 1;
            end
        end

        function ODException(obj,time,command,commandOption,orbit,frame,additionState,derivativeAdditionState,varargin)

            % Exception satellite
                    % Initial quaternion
            if isequal(obj.quaternion(1),-1) && isequal(obj.quaternion(2),-1)
                throw(MException('ODE:MissingInput','Initial quaternion did not set'));
            end

                    % Initial angular velocity
            if isequal(obj.angularVelocity(1),1e30) && isequal(obj.angularVelocity(2),1e30) && isequal(obj.angularVelocity(3),1e30)
                throw(MException('ODE:MissingInput','Initial angular velocity of the satellite did not set'));
            end

            % Exception reaction wheel
            nb_wheel = numel(obj.wheels);
                    % Number of reaction wheel
            if nb_wheel > 3
                throw(MException('ODE:ReactionWheel','The number of reaction wheel must not exceed 3'));
            end

                    % Differents axis
            if nb_wheel == 2
                if obj.wheels(1).getAxis() == obj.wheels(2).getAxis()
                    throw(MException('ODE:ReactionWheel','Reaction wheels must have different axis number'));
                end
            elseif nb_wheel == 3
                if (obj.wheels(1).getAxis() == obj.wheels(2).getAxis()) ||(obj.wheels(1).getAxis() == obj.wheels(3).getAxis()) || (obj.wheels(2).getAxis() == obj.wheels(3).getAxis())
                    throw(MException('ODE:ReactionWheel','Reaction wheels must have different axis number'));
                end
            end

            % Exception magnetorquer
            nb_magnetorquer = numel(obj.magnetorquers);
                    % Number of magnetorquer
            if nb_magnetorquer > 3
                throw(MException('ODE:Magnetorquer','The number of magnetorquer must not exceed 3'));
            end

                    % Differents axis
            if nb_magnetorquer == 2
                if obj.magnetorquers(1).getAxis() == obj.magnetorquers(2).getAxis()
                    throw(MException('ODE:Magnetorquer','Magnetorquers must have different axis number'));
                end
            elseif nb_magnetorquer == 3
                if (obj.magnetorquers(1).getAxis() == obj.magnetorquers(2).getAxis()) ||(obj.magnetorquers(1).getAxis() == obj.magnetorquers(3).getAxis()) || (obj.magnetorquers(2).getAxis() == obj.magnetorquers(3).getAxis())
                    throw(MException('ODE:Magnetorquer','Magnetorquers must have different axis number'));
                end
            end

            % Exception command
            if isa(command,'function_handle')
                % nothing to do
            elseif ischar(command) || isstring(command)
                % Check if the name match witch a file function
                if ~(exist(command,'file') || exist(command,'builtin'))
                    throw(MException('ODE:ExistingCommand','ODE function must be an existing .m file or builtin'))
                end
            else
                throw(MException('ODE:ExistingCommand','ODE function must be a handle or valid function name'))
            end

            % Exception commandOption
            if ~(isstring(commandOption) || ischar(commandOption))
                throw(MException('ODE:BadInput','The commandOption must be a char or a string'));
            else
                commandOption = string(commandOption);
            end

            if (commandOption ~= "quaternion") && (commandOption ~= "cardan")
                throw(MException('ODE:BadInput','The commandOption must be "quaternion" or "cardan"'));
            end

            % Exception orbit
            if (nb_magnetorquer > 0) && (~exist("orbit","var") || ~isa(orbit,"orbit"))
                throw(MException('ODE:MissingOrbit','There are some magnetorquers, an orbit instance must be give to ODE function'));
            end

            if exist("orbit","var") && isempty(orbit.attractorBody)
                throw(MException('ODE:MissingPlanet','Orbit must have a planet instance'));
            end

            % Exception frame
            if exist("frame","var") && ~(isstring(frame) || ischar(frame))
                throw(MException('ODE:BadInput','Frame must be a char or string'));
            end

            % Exception addition state
            if isa(derivativeAdditionState,'function_handle')
                % nothing to do
            elseif ischar(derivativeAdditionState) || isstring(derivativeAdditionState)
                % Check if the name match witch a file function
                if ~(exist("derivativeAdditionState",'file') || exist("derivativeAdditionState",'builtin'))
                    throw(MException('ODE:ExistingDerivativeAdditionState',char(derivativeAdditionState) + ' is not a valid name'));
                end
            else
                throw(MException('ODE:ExistingDerivativeAdditionState','derivativeAdditionState must be a function'));
            end
        end

        function [T,Y,U,M] = ODE(obj,time,command,commandOption,orbit,frame,additionState,derivativeAdditionState,varargin)

            % Exception
            if nargin < 7
                additionState = [];
                derivativeAdditionState = @(state) 0;
            elseif nargin < 8
                throw(MException('ODE:MissingInput','There are 1 or more addition state, ODE must take a function to derivate them'));
            end

            obj.ODException(time,command,commandOption,orbit,frame,additionState,derivativeAdditionState,varargin);

            commandOption = string(commandOption);

            % Get the number of step time
            nb = numel(time);
            
            % Magnetic moment
            % Creation of magnetic moment and magnetic field vector
            magneticMoment_compute = zeros(3,nb);
            magneticMoment_saturated = zeros(3,nb);
            magneticField = zeros(3,nb);

            % Magnetic field
            if exist("orbit","var")       % Compute the orbit trajectory and the earth magnetic field
                [lat,lon,alt] = orbit.geocentric(time);
                addpath("earth_magnetic_field\m_IGRF-main\");
                resultIGRF = igrf(orbit.getDateInitial(),lat,lon,alt,'geocentric');
                for i=1:nb
                    magneticField(1,i) = resultIGRF(i,1);
                    magneticField(2,i) = resultIGRF(i,2);
                    magneticField(3,i) = resultIGRF(i,3);
                end
                if ~exist("frame","var")
                    frame = 'ECI';
                end
                magneticField = tools.conversionCoordinates(magneticField,[lat;lon],orbit,time,frame);
            end

            % Creation of state vectors
            nb_wheel = numel(obj.wheels);       % Number of reaction wheel

            nb_magnetorquer = numel(obj.magnetorquers);      % Number of magnetorquer

            if exist("additionState","var")
                nb_additionState = numel(additionState);
            else
                nb_additionState = 0;
            end

            nb_state = 7 + nb_wheel + nb_additionState;

            T = time;
            Y = zeros(nb_state,nb);         % State vector: 4 quaternion + 3 angular velocity satellite + number of wheel velocity + number of command state
            U_compute = zeros(nb_wheel + nb_magnetorquer,nb);   % Command vector: number of wheel + number of magnetorquer, before saturation
            U_effective = zeros(nb_wheel + nb_magnetorquer,nb);   % Command vector: number of wheel + number of magnetorquer, after saturation

            % Initial conditions
            Y(1:7,1) = [obj.quaternion;obj.angularVelocity];
            for i = 1:nb_wheel
                Y(7+i,1) = obj.wheels(i).getVelocity();
            end

            for i=1:nb_additionState
                Y(7+nb_wheel+i,1) = additionState(i,1);
            end

            % Rearrange "wheels" and "magnetorquers" vectors
            obj = obj.rearrangeReactionWheel();
            obj = obj.rearrangeMagnetorquers();

            affectationCommand = obj.affectationCommand();

            for i=2:nb
                % Step time
                step = time(i) - time(i-1);

                % Evolution of state, computation of commands and
                % saturations
                [Y(:,i),U_compute(:,i),U_effective(:,i),magneticMoment_compute(:,i),magneticMoment_saturated(:,i)] = obj.RK4(Y(:,i-1),command,affectationCommand,commandOption,step,magneticField(:,i),derivativeAdditionState,varargin{:});
            end
            % Get the final value of each state
            obj.quaternion = Y(1:4,end);
            obj.angularVelocity = Y(5:7,end);
            if nb_wheel > 0
                for i = 1:nb_wheel
                    obj.wheels(i) = obj.wheels(i).setVelocity(Y(7+i,end));
                end
            end

            U = [U_compute;U_effective];
            M = [magneticMoment_compute;magneticMoment_saturated];
        end
    end
end