function [return_val] = rk4_PowerSupply(target_val, target_coeff, constant_val, time_resolution)
% [V_s] = rk4_PowerSupply(target_val, target_coeff, constant_val, time_resolution)
% Function takes in a voltage than is to be computed <target_val>, the
% target_coefficient of this value <target_coeff>, the rest of the equation
% is entered as <otherval>. <time_resolution> is the timestep. This
% function requires all inputs or it will throw an error. Can solve both one equation
% or a pair of equations
% Version 2 (05April2017, Andrew Mannion & Conor Dooley)

% Input checking
if(nargin ~= 4) % first check that correct number of inputs
    error('Incorrect number of input arguments. rk4_Power_supply() requires 4 inputs and you gave %d.', nargin);
else
    if( isequal(size(target_val),[1,1]) && isequal(size(target_coeff),[1,1]) && isequal(size(constant_val),[1,1]) )
        %rk4 - single element runge-kutta done if input matrices are singular
        k1 = constant_val + target_val*target_coeff; % get first k value
        z = target_val + k1*(time_resolution/2); % use this to update z
        
        k2 = constant_val + z*target_coeff; % get 2nd k value
        z = target_val + k2*(time_resolution/2); % update z
        
        k3 = constant_val + z*target_coeff; % get third k value
        z = target_val + k3*time_resolution; % update z
        
        k4 = constant_val + z*target_coeff; % get fourth k value
        
        return_val = (target_val) + (k1 + 2*k2 + 2*k3 + k4)*time_resolution/6; % return the integral approximation for the target value
        
    elseif( isequal(size(target_val),[2,1]) && isequal(size(target_coeff),[2,2]) && isequal(size(constant_val),[2,1]) )
        % rk4 - vector form runge-kutta if input matrices have 2 elements
        % Find k value 
        k1 = [constant_val(1)+target_coeff(1,1)*target_val(1)+target_coeff(1,2)*target_val(2);constant_val(2)+target_coeff(2,1)*target_val(1)+target_coeff(2,2)*target_val(2)];
        z = target_val + k1*(time_resolution/2); % use this to update z
        % Find k value 
        k2 = [constant_val(1)+target_coeff(1,1)*z(1)+target_coeff(1,2)*z(2);constant_val(2)+target_coeff(2,1)*z(1)+target_coeff(2,2)*z(2)];
        z = target_val + k2*(time_resolution/2); % update z
        % Find k value 
        k3 = [constant_val(1)+target_coeff(1,1)*z(1)+target_coeff(1,2)*z(2);constant_val(2)+target_coeff(2,1)*z(1)+target_coeff(2,2)*z(2)];
        z = target_val + k2*time_resolution; % update z
        % Find k value 
        k4 = [constant_val(1)+target_coeff(1,1)*z(1)+target_coeff(1,2)*z(2);constant_val(2)+target_coeff(2,1)*z(1)+target_coeff(2,2)*z(2)];
        
        return_val = (target_val) + (k1 + 2*k2 + 2*k3 + k4)*(time_resolution)/6; %  return the integral approximation for the target value
    else
        error('Incorrect form of input arguments. Must have either 1x1 form for target_val, target_coeff & constant_val OR 2x1 target_val');
    end
end
end
