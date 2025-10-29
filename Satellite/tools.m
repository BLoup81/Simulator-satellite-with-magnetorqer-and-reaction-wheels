classdef tools
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%% Conversion %%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function q = euler2quater(euler)
            
            dim = size(euler);
            q = zeros(4,dim(2));
            
            for i=1:dim(2)
                roll = euler(1,i);
                pitch = euler(2,i);
                yaw = euler(3,i);
            
                cr = cos(0.5*roll);
                sr = sin(0.5*roll);
                cp = cos(0.5*pitch);
                sp = sin(0.5*pitch);
                cy = cos(0.5*yaw);
                sy = sin(0.5*yaw);
            
                q(1,i) = cr*cp*cy + sr*sp*sy;
                q(2,i) = sr*cp*cy - cr*sp*sy;
                q(3,i) = cr*sp*cy + sr*cp*sy;
                q(4,i) = cr*cp*sy - sr*sp*cy;
            end
        end
        
        function euler = quater2euler(q)
        
            dim = size(q);
        
            euler = zeros(3,dim(2));
        
            for i=1:dim(2)
                w = q(1,i);
                x = q(2,i);
                y = q(3,i);
                z = q(4,i);
            
                euler(1,i) = atan2(2*(w*x + y*z),1 - 2*((x^2) + (y^2)));
                euler(2,i) = -(pi/2) + 2*atan2(sqrt(abs(1+2*(w*y - x*z))) , sqrt(abs(1-2*(w*y - x*z))));
                euler(3,i) = atan2(2*(w*z + x*y),1 - 2*((y^2) + (z^2)));
            end
        end
        
        
        function eci = ecef2eci(ecef,earthRotation,time)
        
            nb = length(time);

            size_ecef = size(ecef);

            if ~isequal(size_ecef(2),nb)
                throw(MException('ecef2eci:DimensionError','Number of time and ECEF coordinates must be the same'));
            end
        
            eci = zeros(3,nb);
        
            for index = 1:nb
                angle = earthRotation*time(index);
                c = cos(angle);
                s = sin(angle);
                rot = [c -s 0;
                       s c 0;
                       0 0 1];
                
                eci(:,index) = rot*ecef(:,index);
            end
        end
        
        function lvlh = eci2lvlh(eci,argumentOfPerigee,raan,inclination,trueAnomaly)
            % ArgumentOfPerigee, raan, inclination and trueAnomaly must be
            % in radian
        
            nb = length(trueAnomaly);

            size_eci = size(eci);

            if ~isequal(size_eci(2),nb)
                throw(MException('eci2lvlh:DimensionError','True anomaly and ECI coordinates must have the same number of elements'));
            end

            if (inclination < -pi/2) || (inclination > pi/2)
                throw(MException('eci2lvlh:BadInput','Inclination must be in radian between [-pi/2 pi/2]'));
            end

            if (argumentOfPerigee < 0) || (argumentOfPerigee > 2*pi)
                throw(MException('eci2lvlh:BadInput','Argument of perigee must be in radian between [0 2*pi]'));
            end

            if (raan < 0) || (raan > 2*pi)
                throw(MException('eci2lvlh:BadInput','RAAN must be in radian between [0 2*pi]'));
            end
        
            lvlh = zeros(3,nb);
        
            c_i = cos(inclination);
            s_i = sin(inclination);
            c_raan = cos(raan);
            s_raan = sin(raan);
        
            r1 = [0 1 0;
                  0 0 -1;
                  -1 0 0];
        
            r3 = [1 0 0;
                  0 c_i s_i;
                  0 -s_i c_i];
        
            r4 = [c_raan s_raan 0;
                  -s_raan c_raan 0;
                  0 0 1];
        
            for index = 1:nb
                c_t = cos(trueAnomaly(index) + argumentOfPerigee);
                s_t = sin(trueAnomaly(index) + argumentOfPerigee);
                r2 = [c_t s_t 0;
                      -s_t c_t 0;
                      0 0 1];
        
                lvlh(:,index) = r1*r2*r3*r4*eci(:,index);
            end
        end
        
        function brf = frame2brf(inputCoordinates,att)
            size_coor = size(inputCoordinates);
            nb = size_coor(2);

            size_att = size(att);

            if ~isequal(size_att(2),nb)
                throw(MException('frame2brf:BadInput','Input coordinates and quaternions must have the same number of elements'));
            end
        
            brf = zeros(3,nb);
        
            for index = 1:nb
                rot = [1-2*((att(3,index)^2)+(att(4,index)^2)) 2*(att(2,index)*att(3,index)+att(1,index)*att(4,index)) 2*(att(2,index)*att(4,index)-att(1,index)*att(3,index));
                       2*(att(3,index)*att(2,index)-att(1,index)*att(4,index)) 1-2*((att(2,index)^2)+(att(4,index)^2)) 2*(att(3,index)*att(4,index)+att(1,index)*att(2,index));
                       2*(att(4,index)*att(2,index)+att(1,index)*att(3,index)) 2*(att(4,index)*att(3,index)-att(1,index)*att(2,index)) 1-2*((att(2,index)^2)+(att(3,index)^2))];
            
                brf(:,index) = rot*inputCoordinates(:,index);
            end
        end
        
        function ecef = ned2ecef(lat,lon,ned)
            % Input:
            %   lat = latitude  [rad]
            %   lon = longitude [rad]
            %   ned = north, east, down (NED) coordinates
            %
            % Output:
            %   ecef = Earth-Centered Earth-Fixed frame (ECEF) coordinates
        
            nb_lat = length(lat);
            nb_lon = length(lon);

            size_ned = size(ned);
            nb_ned = size_ned(2);

            if (~isequal(nb_lat,nb_ned)) || (~isequal(nb_lat,nb_lon))
                throw(MException('ned2ecef:DimensionError','Latitude, longitude and NED coordinates must have the same size'));
            end
        
            ecef = zeros(3,nb_lat);
        
            for index = 1:nb_lat
                slat = sin(lat(index));
                clat = cos(lat(index));
                slon = sin(lon(index));
                clon = cos(lon(index));
            
                rot = [-clon*slat -slon -clat*clon;
                       -slat*slon clon -slon*clat;
                       clat 0 -slat];
            
                ecef(:,index) = rot*ned(:,index);
            end
        end

        function outputCoordinates = conversionCoordinates(inputCoordinate,geocentric,orbit,time,frame)

            frame = string(frame);

            if frame == "NED"
                outputCoordinates = inputCoordinate;
            elseif frame == "ECEF"
                outputCoordinates = tools.ned2ecef(geocentric(1,:),geocentric(2,:),inputCoordinate);
            elseif frame == "ECI"
                ecef = tools.ned2ecef(geocentric(1,:),geocentric(2,:),inputCoordinate);
                outputCoordinates = tools.ecef2eci(ecef,orbit.attractorBody.getRotation(),time);
            elseif frame == "LVLH"
                ecef = tools.ned2ecef(geocentric(1,:),geocentric(2,:),inputCoordinate);
                eci = tools.ecef2eci(ecef,orbit.attractorBody.getRotation(),time);
                outputCoordinates = tools.eci2lvlh(eci,orbit.getArgumentOfPeriapsis('rad'),orbit.getRAAN('rad'),orbit.getInclination('rad'),orbit.trueAnomaly(time));
            else
                msg = "The frame " + frame + " is unknow";
                throw(MException('conversionCoordinates:UnknowFrame',msg));
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function plotAttitudeQuaternion(time,state)
            subplot(4,1,1);
            plot(time,state(1,:));
            ylim([-1 1]);
            ylabel('q_0');

            subplot(4,1,2);
            plot(time,state(2,:));
            ylim([-1 1]);
            ylabel('q_1');
            
            subplot(4,1,3);
            plot(time,state(3,:));
            ylim([-1 1]);
            ylabel('q_2');

            subplot(4,1,4);
            plot(time,state(4,:));
            ylim([-1 1]);
            ylabel('q_3');

            xlabel('Time (s)');
            
        end

        function plotSatelliteVelocity(time,state)
            plot(time,state(5:7,:));
            ylabel('Satellite angular velocity (rad/s)');
            xlabel('Time (s)');
            legend('Yaw','Pitch','Roll');
        end

        function setLegendCardan(actuators)
            nb_actuator = numel(actuators);

            if isequal(nb_actuator,3)
                legend('Yaw','Pitch','Roll');
            elseif isequal(nb_actuator,2)
                if isequal(actuators(1).getAxis(),1)
                    if isequal(actuators(2).getAxis(),2)
                        legend('Yaw','Pitch');
                    else
                        legend('Yaw','Roll');
                    end
                else
                    legend('Pitch','Roll');
                end
            else
                if isequal(actuators(1).getAxis(),1)
                    legend('Yaw');
                elseif isequal(actuators(1).getAxis(),2)
                    legend('Pitch');
                else
                    legend('Roll');
                end
            end
        end

        function plotWheelVelocity(time,state,satellite)

            nb_wheel = numel(satellite.wheels);
            plot(time,state(8:7+nb_wheel,:));
            ylabel('Wheel angular velocity (rad/s)');
            xlabel('Time (s)');

            tools.setLegendCardan(satellite.wheels);
        end

        function plotWheelTorque(time,torque,satellite)
            nb_wheel = numel(satellite.wheels);
            nb_magnetorquer = numel(satellite.magnetorquers);

            if nb_wheel < 1
                disp('Plot of wheel torques is impossible, there are no wheels');
            else
                subplot(2,1,1);
                plot(time,torque(1:nb_wheel,:));
                ylabel('Compute torque RW (N.m)');
                tools.setLegendCardan(satellite.wheels);
    
                subplot(2,1,2);
                plot(time,torque(nb_wheel + nb_magnetorquer + 1:2*nb_wheel + nb_magnetorquer,:));
                ylabel('Effective torque RW (N.m)');
                xlabel('Time (s)');
                tools.setLegendCardan(satellite.wheels);
            end
        end

        function plotMagnetorquerTorque(time,torque,satellite)
            nb_wheel = numel(satellite.wheels);
            nb_magnetorquer = numel(satellite.magnetorquers);

            if nb_magnetorquer < 1
                disp('Plot of magnetic torque is impossible, there are no magnetorquers');
            else
                subplot(2,1,1);
                plot(time,torque(nb_wheel + 1:nb_wheel + nb_magnetorquer,:));
                ylabel('Compute torque MGT (N.m)');
                tools.setLegendCardan(satellite.magnetorquers);
    
                subplot(2,1,2);
                plot(time,torque(2*nb_wheel + nb_magnetorquer + 1:2*nb_wheel + 2*nb_magnetorquer,:));
                ylabel('Effective torque MGT (N.m)');
                xlabel('Time (s)');
                tools.setLegendCardan(satellite.magnetorquers);
            end
        end

        function plotCardan(time,state)
            cardan = tools.quater2euler(state(1:4,:));
            plot(time,cardan);
            xlabel('Time (s)');
            ylabel('Angles of Cardan');
            legend('Yaw','Pitch','Roll');
        end

        function plotMagneticMoment(time,magneticMoment)
            subplot(2,1,1);
            plot(time,magneticMoment(1:3,:));
            ylabel('Magnetic moment compute (N.m.T^{-1})');
            legend('Mx','My','Mz');

            subplot(2,1,2);
            plot(time,magneticMoment(4:6,:));
            ylabel('Magnetic moment saturated (N.m.T^{-1}');
            xlabel('Time (s)');
            legend('Mx','My','Mz');
        end
    end
end