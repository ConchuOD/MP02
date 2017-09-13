function [V_out,I_out] = Power_Supply_Circuit_Solver(V_in,I_in,Impedances,Sw_1,Sw_2,time_resolution,timesteps,start_time)
% This function solves the circuit problem in M&S 2017 MP02.
% [V_out,I_out] = Power_Supply_Circuit_Solver(V_in,I_in,Impedances,Sw_1,Sw_2,time_resolution,timesteps,time)
%
% V_in takes the form <[V1,V2,V3,input_amplitude]>, where these are the node voltages (1 @ top of the cap,
% 2 @ top of the zener,3 @ the bottom of cap) and input_amplitude is the
% maximum amplitude of the sinusoidal mains voltage.
% I_in takes the form <[I_C,I_Rs,I_L1,R_L2]>, where these are the branchcurrents.
% (I_C through the cap,I_Rs through Rs, I_L1 through resistive load and I_L2 through the inductive load)
% Impedance takes the form <[R_s,C_smoothing]> as these may vary with the tolerance.
% Sw_1 is whether the resistive load is connected. (binary value, 0 open, 1 closed)
% Sw_2 is whether the inductive load is connected. (binary value, 0 open, 1 closed)
% time_resolution is the diffence in time between the voltage/current input
% and the desired output in seconds. Suggested time_resolution is less than 0.0001
% seconds. 
% timesteps is the number of iterations for which the function is to
% determine the next value of the voltages and currents.
% time is the time, in seconds since the mains was turned on.
%
% V_out takes the form <[V1(timesteps);V2(timesteps);V3(timesteps);Input_voltage]>,
% where these are the node voltages (1 @ top of the cap, 2 @ top of the zener,
% 3 @ the bottom of cap, and input_voltages are the historical values of the rectified voltage across the diode bridge.)
% I_out takes the form <[I_C(timesteps);I_Rs(timesteps);I_L1(timesteps);R_L2(timesteps)]>,
% where these are the branchcurrents.
% (I_C through the cap,I_Rs through Rs, I_L1 through resistive load and I_L2 through the inductive load)
%
% Version 2 (07April2017, Andrew Mannion & Conor Dooley)

%% Checking input validity
if nargin ~= 8
    error('Power_Supply_Circuit_Solver() requires 8 input arguments - Please read the help file for more information.')
elseif(~isequal(size(V_in),[1,4]))
    error('V_in must be a row vector of 4 voltages - Please read the help file for more information.')
elseif(~isequal(size(I_in),[1,4]))
    error('I_in must be a row vector of 4 currents - Please read the help file for more information.')
elseif(~(Sw_1 == 1 || Sw_1 == 0) || ~(Sw_2 == 1 ||Sw_2 == 0))
    error('Switches may only be set to either zero or one - Please read the help file for more information.')
elseif(time_resolution > 0.001)
    error('time_resolution is too small to give any accuracy the Runge-Kutta predictions for the next values.')
elseif( timesteps < 1 && rem(timesteps,1) ~= 0)
    error('Please enter a positive integer value for the number of timesteps you wish to calculate for.')
elseif( start_time < 1 && rem(start_time,1) ~= 0)
    error('Please enter a possitive integer value for the number of timesteps that had occured prior to the function call.')
end
%% Determining the state
if(Sw_1 == 0 && Sw_2 == 0)
    config = 1;
elseif(Sw_1 == 1 && Sw_2 == 0)
    config = 2;
elseif(Sw_1 == 0 && Sw_2 == 1)
    config = 3;
elseif(Sw_1 == 1 && Sw_2 == 1)
    config = 4;
else
    error('Could not determine the correct configuration - Please read the help file for more information.');
end
%% declaration of constants
N = 1/20; % Turns ratio of the transformer (assuming ideal transformer)
V_t_Diode = 0.56;  % Forward voltage of the diodes (Volts)
V_t_Zener = 6.5; % Reverse voltage of the Zener (Volts)
R_Z = 5.3; % Reverse resitance of the Zener (ohms)
R_d = 0.2; % Diode resistance in linear model (ohms)
R_L1 = 287; % LHS load resistance (ohms)
R_L2 = 65.7; % RHS load resistance (ohms)
L = 0.1058; % Inductance in RHS load (Henrys)
I_C = zeros(1,timesteps+1); % Array of cap currents
I_Rs = zeros(1,timesteps+1); % Array of smoothing resistor currents
I_L1 = zeros(1,timesteps+1); % Array of LHS load currents
I_L2 = zeros(1,timesteps+1); % Array of RHS load currents
V1 = zeros(1,timesteps+1); % Array of voltages at node 1
V2 = zeros(1,timesteps+1); % Array of voltages at node 2
V3 = zeros(1,timesteps+1); % Array of voltages at node 3
input_voltage = zeros(1,timesteps+1); % Array of input voltages
n = 1; % iterator must start at one
time = start_time; % start at given time
%% Reading in circuit values
% input_voltage(1) = abs(V_in(4))*N;     % Input voltage rectified by bridge and scaled by step down transformer
R_S = Impedances(1); % Value of smoothing resistor (ohms)
C_smoothing = Impedances(2); % Capacitance of electrolytic smoothing capacitor (Farads)
V1(1) = V_in(1);
V2(1) = V_in(2);
V3(1) = V_in(3);
I_C(1)= I_in(1);
I_Rs(1)= I_in(2);
I_L1(1)= I_in(3);
I_L2(1)= I_in(4);
%% If statements determine which state and therefore which differential equation(s) apply.
% target_val contains the variables being refined by the 4th order runge kutta
% target_coeff contains the coefficients relating to these variables
% constant_val are the other terms in each equation
% config 1: Neither load selected
% config 2: Left hand load (resistive) selected
% config 3: Right hand load (inductive) selected
% condig 4: Both loads selected
while(n < timesteps + 1)
    old = n; 
    next = n + 1;
    input_voltage(old) = abs(V_in(4)*sin(100*pi*time))*N; % calculate voltage across diode bridge at this time
    n = n + 1;
    time = time + time_resolution;
    if config == 1 % config 1 is where no load is connected
        if ~(input_voltage(old) > 2*V_t_Diode) && (V1(old) == 0) && (V2(old) == 0) && (V3(old) == 0)
            % Condition 1 - Rectifier diodes not conducting, therefore all internal voltages and currents are at zero (startup only)
            V1(next)=0;
            V2(next)=0;
            V3(next)=0;
            I_Rs(next) = (V1(next)-V2(next))/R_S;
            I_C(next) = 0;
            I_L1(next) = 0;
            I_L2(next) = 0;
        elseif ((input_voltage(old) > 2*V_t_Diode) && ((input_voltage(old)-V1(old)) > V_t_Diode) && ((V2(old)-V3(old)) < 6.5)) % Condition 2
            % Condition 2 - Capacitor is being charged, but Zener is off.
            % rk4_PS is called to solve the diff equation
            target_val = V1(old)-V3(old);
            target_coeff = -1/(2*C_smoothing*R_d);
            constant_val = (input_voltage(old)-(2*V_t_Diode))/(2*C_smoothing*R_d);
            V1_3 = rk4_PowerSupply(target_val, target_coeff, constant_val, time_resolution);
            % Using this to solve for voltages and currents
            V1(next) = input_voltage(old)-2*V_t_Diode-2*(R_d)*((input_voltage(old)-2*V_t_Diode-V1_3)/(2*R_d));
            V2(next) = V1(next);
            V3(next) = 0;
            I_Rs(next) = (V1(next)-V2(next))/R_S;
            I_C(next) = -C_smoothing * ((V1(next)-V3(next))-(V1(old)-V3(old)))*1/time_resolution;
            % Remainder of currents are at zero as load is not connected and Zener is off.
            I_L1(next) = 0;
            I_L2(next) = 0;
        elseif ((input_voltage(old)-V1(old)) > 2*V_t_Diode) && ((V2(old)-V3(old)) > 6.5) % Condition 3
            % Condition 3 - Still charging the capacitor, but the zener diode is now in the reverse conducting mode
            % rk4_PS is called to solve the diff equation
            target_val = V1(old)-V3(old);
            target_coeff = -1/(C_smoothing*(R_S+R_Z)) -1/(C_smoothing*2*R_d);
            constant_val =  V_t_Zener/(C_smoothing*(R_S+R_Z)) + (-2*V_t_Diode+input_voltage(old))/(C_smoothing*2*R_d);
            V1_3 = rk4_PowerSupply(target_val, target_coeff, constant_val, time_resolution);
            % Using this to solve for voltages and currents
            V1(next) = input_voltage(old)-2*V_t_Diode-2*R_d*((-V1_3+input_voltage(old)-2*V_t_Diode)/(2*R_d));
            V2(next) = V1(next)-((V1_3-V_t_Zener)/(R_S+R_Z))*R_S;
            V3(next) = 0;
            I_Rs(next) = (V1(next)-V2(next))/R_S;
            I_C(next) = -C_smoothing * ((V1(next)-V3(next))-(V1(old)-V3(old)))*1/time_resolution;
            % Still no load connected, so load currents are zero
            I_L1(next) = 0;
            I_L2(next) = 0;
        elseif (~((input_voltage(old)-V1(old)) > 2*V_t_Diode) || ~(V3(old) > V_t_Diode)) && ((V2(old)-V3(old)) > V_t_Zener) % Condition 4
            % Condition 4 - The capacitor is being discharged through the zener to ground and the input is no longer connected
            % Both diodes are turned off if one is turned off
            % rk4_PS is called to solve the diff equation
            target_val = V1(old)-V3(old);
            target_coeff = -1/(C_smoothing*(R_S+R_Z));
            constant_val =  V_t_Zener/(C_smoothing*(R_S+R_Z));
            V1_3 = rk4_PowerSupply(target_val, target_coeff, constant_val, time_resolution);
            % Using this to solve for voltages and currents
            V1(next) = V1_3;
            V2(next) = V1_3+R_S*((-V1_3+V_t_Zener)/(R_S+R_Z));
            V3(next) = 0;
            I_Rs(next) = (V1(next)-V2(next))/R_S;
            I_C(next) = -C_smoothing * ((V1(next)-V3(next))-(V1(old)-V3(old)))*1/time_resolution;
            % Still no load connected
            I_L1(next) = 0;
            I_L2(next) = 0;
        else % Condition 5 - Voltage across capacitor no longer enough to exceed Zener voltage
            % voltages no longer change
            V1(next)=V1(old);
            V2(next)=V2(old);
            V3(next)=V3(old);
            % load currents are still at zero as the loads are open circuit
            I_Rs(next) = (V1(next)-V2(next))/R_S;
            I_C(next) = -C_smoothing * ((V1(next)-V3(next))-(V1(old)-V3(old)))*1/time_resolution;
            I_L1(next) = 0;
            I_L2(next) = 0;
        end
    elseif config == 2 % config 2 is where LHS load is connected
        if ~(input_voltage(old) > 2*V_t_Diode) && (V1(old) == 0) && (V2(old) == 0) && (V3(old) == 0) % Condition 1
            % Condition 1 - Rectifier diodes not conducting, therefore all internal voltages and currents are at zero, as this state is only on startup
            V1(next) = 0;
            V2(next) = 0;
            V3(next) = 0;
            I_Rs(next) = (V1(next)-V2(next))/R_S;
            I_C(next) = 0;
            I_L1(next) = 0;
            I_L2(next) = 0;
        elseif ((input_voltage(old) > 2*V_t_Diode) && ((input_voltage(old)-V1(old)) > V_t_Diode) && ((V2(old)-V3(old)) < 6.5)) % Condition 2
            % Condition 2 - Capacitor is being charged, but Zener is off
            % rk4_PS is called to solve the diff equation
            target_val = V1(old)-V3(old);
            target_coeff = -1/(C_smoothing*(R_S+R_L1))-1/(C_smoothing*2*R_d);
            constant_val =  (-2*V_t_Diode+input_voltage(old))/(2*C_smoothing*R_d);
            V1_3 = rk4_PowerSupply(target_val, target_coeff, constant_val, time_resolution);
            % Using this to solve for voltages and currents
            V1(next) = input_voltage(old)-2*V_t_Diode-2*(R_d)*((input_voltage(old)-2*V_t_Diode-V1_3)/(2*R_d));
            V2(next) = V1(next)-R_S*((V1_3)/(R_S+R_L1));
            V3(next) = 0;
            I_Rs(next) = (V1(next)-V2(next))/R_S;
            I_C(next) = -C_smoothing * ((V1(next)-V3(next))-(V1(old)-V3(old)))*1/time_resolution;
            % RH load is not connected but LH is.
            I_L1(next) = I_Rs(next);
            I_L2(next) = 0;
        elseif ((input_voltage(old)-V1(old)) > 2*V_t_Diode) && ((V2(old)-V3(old)) > 6.5) % Condition 3
            % Condition 3 - Still charging the capacitor, but the zener diode is now in the reverse conducting mode
            % rk4_PS is called to solve the diff equation
            target_val = V1(old)-V3(old);
            target_coeff = -1/(2*C_smoothing*R_d) -(R_Z+R_L1)/(C_smoothing*(R_Z*R_L1+R_S*(R_Z+R_L1)));
            constant_val = (input_voltage(old)-2*V_t_Diode)/(2*C_smoothing*R_d)+(V_t_Zener*R_L1)/(C_smoothing*(R_Z*R_L1+R_S*(R_Z+R_L1)));
            V1_3 = rk4_PowerSupply(target_val, target_coeff, constant_val, time_resolution);
            % Using this to solve for voltages and currents
            V1(next) = input_voltage(old)-2*V_t_Diode-2*(R_d)*((input_voltage(old)-2*V_t_Diode-V1_3)/(2*R_d));
            V2(next) = V1(next)-(((V1_3)*(R_Z+R_L1)-V_t_Zener*R_L1)/(R_Z*R_L1 + R_S*(R_Z+R_L1)))*R_S;
            I_Rs(next) = (V1(next)-V2(next))/R_S;
            I_C(next) = -C_smoothing * ((V1(next)-V3(next))-(V1(old)-V3(old)))*1/time_resolution;
            % RH load is not connected but LS is.
            I_L1(next) = (V2(next)-V3(next))/R_L1;
            I_L2(next) = 0;
        elseif (~((input_voltage(old)-V1(old)) > 2*V_t_Diode)) && ((V2(old)-V3(old)) > V_t_Zener) % Condition 4
            % Condition 4 - The capacitor is being discharged through the zener to ground and the input is no longer connected.
            % Both diodes are turned off if one is turned off
            % rk4_PS is called to solve the diff equation
            target_val = V1(old)-V3(old);
            target_coeff = -(R_Z+R_L1)/(C_smoothing*(R_Z*R_L1+R_S*(R_Z+R_L1)));
            constant_val = (V_t_Zener*R_L1)/(C_smoothing*(R_Z*R_L1+R_S*(R_Z+R_L1)));
            V1_3 = rk4_PowerSupply(target_val, target_coeff, constant_val, time_resolution);
            % Using this to solve for voltages and currents
            V1(next) = V1_3;
            V2(next) = V1_3-(((V1_3)*(R_Z+R_L1)-V_t_Zener*R_L1)/(R_Z*R_L1 + R_S*(R_Z+R_L1)))*R_S;
            V3(next) = 0;
            % RH load is not connected but Zener is still on
            I_Rs(next) = (V1(next)-V2(next))/R_S;
            I_C(next) = -C_smoothing * ((V1(next)-V3(next))-(V1(old)-V3(old)))*1/time_resolution;
            I_L1(next) = (V2(next)-V3(next))/R_L1;
            I_L2(next) = 0;
        else % Condition 5 - Voltage across capacitor no longer enough to exceed Zener voltage
            % Voltages no longer change
            V1_3 = rk4_PowerSupply(V1(old)-V3(old),-1/(C_smoothing*(R_S+R_L1)),0,time_resolution);
            V1(next)=V1_3;
            V2(next)=V1_3-(V1_3)/(R_S+R_L1);
            V3(next)=0;
            % Currents through zener are now at zero, path to ground still exists through load
            I_Rs(next) = (V1(next)-V2(next))/R_S;
            I_C(next) = -C_smoothing * ((V1(next)-V3(next))-(V1(old)-V3(old)))*1/time_resolution;
            I_L1(next) = I_Rs(next);
            I_L2(next) = 0;
        end
    elseif config == 3 
        if ~(input_voltage(old) > 2*V_t_Diode) && (V1(old) == 0) && (V2(old) == 0) && (V3(old) == 0) % Condition 1
			% Condition 1 - Rectifier diodes not conducting, therefore all internal voltages and currents are at zero, as this state is only on startup
            V1(next)=0;
            V2(next)=0;
            V3(next)=0;
            I_Rs(next) = (V1(next)-V2(next))/R_S;
            I_C(next) = 0;
            I_L1(next) = 0;
            I_L2(next) = 0;
        elseif (((input_voltage(old)-V1(old)) > 2*V_t_Diode) && ((V2(old)-V3(old)) < 6.5)) % Condition 2
            % Condition 2 - Capacitor is being charged, but Zener is off
            % Set up vectors to pass to rk4_PS
            target_val = [V1(old)-V3(old);I_L2(old)];
            target_coeff = [-1/(2*C_smoothing*R_d),-1/C_smoothing;1/L,(-(R_S+R_L2))/L];
            constant_val = [(input_voltage(old)-2*V_t_Diode)/(2*C_smoothing*R_d);0];
            % rk4_PS is called to solve the diff equation
            Results = rk4_PowerSupply(target_val, target_coeff, constant_val, time_resolution);
            V1(next) = input_voltage(old)-2*V_t_Diode;
            V2(next) = V1(next)-(R_S)*Results(2);
            V3(next) = 0;
            % RHS load is now connected, and its value comes from the evaluation of runge kutta.
            I_Rs(next) = (V1(next)-V2(next))/R_S;
            I_C(next) = (input_voltage(old)-2*V_t_Diode-Results(1))/(2*R_d)-Results(2);
            I_L1(next) = 0;
            I_L2(next) = Results(2);
        elseif ((input_voltage(old)-V1(old)) > 2*V_t_Diode) && ((V2(old)-V3(old)) > V_t_Zener) % Condition 3
            % Still charging the capacitor, but the zener diode is now in the reverse conducting mode
            % Set up vectors to pass to rk4_PS
            target = [V1(old)-V3(old);I_L2(old)];
            target_coeff = [-1/(C_smoothing*(R_Z+R_S))-1/(2*C_smoothing*R_d),R_Z/(C_smoothing*(R_Z+R_S));R_Z/(L*(R_Z+R_S)),(-R_L2*(R_Z+R_S)+R_S*R_Z)/(L*(R_Z+R_S))];
            constant_val = [(V_t_Zener)/(C_smoothing*(R_Z+R_S))+(input_voltage(old)-2*V_t_Diode)/(2*C_smoothing*R_d);(R_S*V_t_Zener)/(L*(R_Z+R_S))];
            % rk4_PS is called to solve the diff equation
            Results = rk4_PowerSupply(target, target_coeff, constant_val, time_resolution);
            V1(next) = input_voltage(old)-2*V_t_Diode;
            V2(next) = V1(next)-R_S*((Results(1)-V_t_Zener-R_Z*Results(2))/(R_Z+R_S));
            V3(next) = 0;
            % I_C calculated from RK4 output
            I_Rs(next) = (V1(next)-V2(next))/R_S;
            I_C(next) = (input_voltage(old)-2*V_t_Diode-Results(1))/(2*R_d)+(-Results(1)+V_t_Zener+R_Z*Results(2))/(R_Z+R_S);
            % RH load connected, LH disconnected still
            I_L1(next) = 0;
            I_L2(next) = Results(2);
        elseif (~((input_voltage(old)-V1(old)) > 2*V_t_Diode)) && ((V2(old)-V3(old)) > 6.5) % Condition 4
            % The capacitor is being discharged through the zener to ground
            % Both diodes are turned off if one is turned off
            % Set up vectors to pass to rk4_PS
            target = [V1(old)-V3(old);I_L2(old)];
            target_coeff = [-1/(C_smoothing*(R_Z+R_S)),R_Z/(C_smoothing*(R_Z+R_S));R_Z/(L*(R_Z+R_S)),(-R_L2*(R_Z+R_S)+R_S*R_Z)/(L*(R_Z+R_S))];
            constant_val = [V_t_Zener/(C_smoothing*(R_Z+R_S));(R_S*V_t_Zener)/(L*(R_Z+R_S))];
            % rk4_PS is called to solve the diff equation
            Results = rk4_PowerSupply(target, target_coeff, constant_val, time_resolution);
            V1(next) = Results(1);
            V2(next) = V1(next)-R_S*((Results(1)-V_t_Zener-R_Z*Results(2))/(R_Z+R_S));
            V3(next) = 0;
            I_Rs(next) = (V1(next)-V2(next))/R_S;
            I_C(next) = -(-Results(1)+V_t_Zener+R_Z*Results(2))/(R_Z+R_S);
            I_L1(next) = 0;
            I_L2(next) = Results(2);
         else % Condition 5
            % Zener not activated, so path to ground exists through the load. Cap discharges
            % Set up vectors to pass to rk4_PS
            target = [V1(old)-V3(old);-I_L2(old)/C_smoothing];
            target_coeff = [0,1;-1/(L*C_smoothing),-(R_S+R_L2)/L];
            constant_val = [0;0];
            % rk4_PS is called to solve the diff equation
            Results = rk4_PowerSupply(target, target_coeff, constant_val, time_resolution);
            % need load current to calculate V2, found in RK4
            I_L2(next) = Results(2)*-C_smoothing;
            V1(next)=Results(1);
            V2(next)=V1(next)-R_S*I_L2(next);
            V3(next)=0;
            I_Rs(next) = (V1(next)-V2(next))/R_S;
            I_C(next) = -C_smoothing*Results(1);
            I_L1(next) = 0;
        end
    elseif config == 4 % config 4 is where LHS and RHS load are connected
        % missing voltage and current equations
        if ~(input_voltage(old) > 2*V_t_Diode) && (V1(old) == 0) && (V2(old) == 0) && (V3(old) == 0 )  % Condition 1
            % Condition 1 - Rectifier diodes not conducting, therefore all internal voltages and currents are at zero, as this state is only on startup
            V1(next)=0;
            V2(next)=0;
            V3(next)=0;
            I_Rs(next) = (V1(next)-V2(next))/R_S;
            I_C(next) = 0;
            I_L1(next) = 0;
            I_L2(next) = 0;
        elseif ((input_voltage(old)-V1(old) > 2*V_t_Diode) && ((V2(old)-V3(old)) < 6.5)) % Condition 2
            % Condition 2 - Capacitor is being charged, but Zener is off
            target_val = [V1(old)-V3(old);I_L2(old)];
            target_coeff = [-1/(2*C_smoothing*R_d)-1/(C_smoothing*(R_L1+R_S)),-R_L1/(C_smoothing*(R_L1+R_S));R_L1/(L*(R_L1*R_S)),-(R_S*R_L1-R_L2*(R_L1+R_S))/(L*(R_L1*R_S))];
            constant_val = [(input_voltage(old)-2*V_t_Diode)/(2*C_smoothing*R_d);0];
            % rk4_PS is called to solve the diff equation
            Results = rk4_PowerSupply(target_val, target_coeff, constant_val, time_resolution);
            % Using this to solve for voltages and currents
            V1(next) = input_voltage(old)-2*V_t_Diode-2*(R_d)*((input_voltage(old)-2*V_t_Diode-Results(1))/(2*R_d));
            V2(next) = V1(next)-R_S*((Results(1)+R_L1*Results(2))/(R_L1+R_S));
            V3(next) = 0;
            % LHS load is an open circuit
            I_Rs(next) = (V1(next)-V2(next))/R_S;
            I_C(next) = 0; % missing
            I_L1(next) = (V2(next)+V3(next))/R_L1;
            % RHS load current comes from the RK4 again
            I_L2(next) = Results(2);
        elseif ((input_voltage(old)-V1(old)) > 2*V_t_Diode) && ((V2(old)-V3(old)) > 6.5 ) % Condition 3
            % Condition 3 - Still charging the capacitor, but the zener diode is now in the reverse conducting mode
            % Setting up the input vectors
            target_val = [V1(old)-V3(old);I_L2(old)];
            target_coeff = [-(R_L1+R_Z)/(C_smoothing*(R_Z*R_L1+R_S*R_L1+R_S*R_Z)),-(R_Z*R_L1)/(C_smoothing*(R_Z*R_L1+R_S*R_L1+R_S*R_Z));1/L-R_S*((R_L1+R_Z)/(L*(R_Z*R_L1+R_S*R_L1+R_S*R_Z))),-R_S*((R_L1*R_Z)/(L*(R_Z*R_L1+R_S*R_L1+R_S*R_Z)))-R_L1/L];
            constant_val = [(V_t_Zener*R_L1)/(C_smoothing*(R_Z*R_L1+R_S*R_L1+R_S*R_Z))+(input_voltage(old)-2*V_t_Diode)/(2*C_smoothing*R_d);-R_S*(-V_t_Zener*R_L1)/(L*(R_Z*R_L1+R_S*R_L1+R_S*R_Z))];
            % rk4_PS is called to solve the diff equation
            Results = rk4_PowerSupply(target_val, target_coeff, constant_val, time_resolution);
            % Using this to solve for voltages and currents
            V1(next) = input_voltage(old)-2*V_t_Diode;%-(R_d)*((input_voltage(old)-2*V_t_Diode-Results(1))/(2*R_d));
            V2(next) = V1(next)-R_S*((Results(1)*(R_L1+R_Z)+Results(2)*(R_Z*R_L1)-V_t_Zener)/(R_Z*R_L1+R_S*R_L1+R_S*R_Z));
            V3(next) = 0;
            % Both loads are now connected. Current through Load 2 comes from the RK4. Load 1 is calculated from voltage drop across resistor
            I_Rs(next) = (V1(next)-V2(next))/R_S;
            I_C(next) = 0; % missing
            I_L1(next) = (V2(next)+V3(next))/R_L1;
            I_L2(next) = Results(2);
        elseif (~((input_voltage(old)-V1(old)) > 2*V_t_Diode)) && ((V2(old)-V3(old)) > 6.5 ) % Condition 4
            % Condition 4 - The capacitor is being discharged through the zener to ground and the input is no longer connected.
            % Setting up the input vectors
            target_val = [V1(old)-V3(old);I_L2(old)];
            target_coeff = [-(R_L1+R_Z)/(C_smoothing*(R_Z*R_L1+R_S*R_L1+R_S*R_Z)),(-(R_Z*R_L1))/(C_smoothing*(R_Z*R_L1+R_S*R_L1+R_S*R_Z));1/L-R_S*((R_L1+R_Z)/(L*(R_Z*R_L1+R_S*R_L1+R_S+R_L1))),-R_S*(R_Z*R_L1)/(L*(R_Z*R_L1+R_S*R_L1+R_S+R_L1))-R_L2];
            constant_val = [(V_t_Zener*R_L1)/(C_smoothing*(R_Z*R_L1+R_S*R_L1+R_S*R_Z));-R_S*(-V_t_Zener*R_L1)/(L*(R_Z*R_L1+R_S*R_L1+R_S+R_L1))];
            % rk4_PS is called to solve the diff equation
            Results = rk4_PowerSupply(target_val, target_coeff, constant_val, time_resolution);
            % Using this to solve for voltages and currents
            V1(next) = Results(1);
            V2(next) = V1(next)-R_S*((Results(1)*(R_L1+R_Z)+Results(2)*(R_Z*R_L1)-V_t_Zener)/(R_Z*R_L1+R_S*R_L1+R_S*R_Z));
            V3(next) = 0;
            I_Rs(next) = (V1(next)-V2(next))/R_S;
            I_C(next) = 0; % missing
            % Both loads are now connected. Current through Load 2 comes from the RK4. Load 1 is calculated from voltage drop across resistor
            I_L1(next) = (V2(next)+V3(next))/R_L1;
            I_L2(next) = Results(2);
        else
            % Condition 5 - Voltage across capacitor no longer enough to exceed Zener voltage
            % Voltages no longer change
            % Setting up the input vectors
            target_val = [V1(old)-V3(old);I_L2(old)];
            target_coeff = [-1/(C_smoothing*(R_L1+R_S)),-1*R_L1/(C_smoothing*(R_L1+R_S));R_L1/(L*(R_L1+R_S)),-((R_L2*(R_L1+R_S)+R_S*R_L1)/(L*(R_L1+R_S)))];
            constant_val = [0;0];
            % rk4_PS is called to solve the diff equation
            Results = rk4_PowerSupply(target_val, target_coeff, constant_val, time_resolution);
            V1(next) = Results(1);
            V2(next) = V1(next)-R_S*((Results(1)+R_L1*Results(2))/(R_L1+R_S));
            V3(next) = 0;
            I_Rs(next) = (V1(next)-V2(next))/R_S;
            I_C(next) = (Results(1)+Results(2)*R_L1)/(R_L1+R_S);
            % Both loads are now connected. Current through Load 2 comes from the RK4. Load 1 is calculated from voltage drop across resistor
            I_L1(next) = (V2(next)+V3(next))/R_L1;
            I_L2(next) = Results(2);
        end
        
    else
        error('This configuration does not exist');
    end
end
V_out = [V1;V2;V3;input_voltage]; % assign array of voltages to return
I_out = [I_C;I_Rs;I_L1;I_L2]; % assign array of currents to return
end
