Mains_RMS = 230; % Root-mean-square value of the mains voltage, the input to the circuit
N = 1/20; % Turns ratio of the transformer (assuming ideal transformer)
Voltage_RHS = (Mains_RMS*(2)^(.5))*N; % Input to circuit to the right of the transformer
C_smoothing = 2200e-6; % Capacitance of electrolytic smoothing capacitor (Farads)
V_t_Diode = 0.56;  % Forward voltage of the diodes (Volts)
V_t_Zener = 6.5; % Reverse voltage of the Zener (Volts)
R_S = 270; % Value of smoothing resistor (ohms)
R_Z = 5.3; % Reverse resitance of the Zener (ohms)
R_d = 0.2; % Diode resistance in linear model (ohms)
R_L1 = 287; % LHS load resistance (ohms)
R_L2 = 65.7; % RHS load resistance (ohms)
L = 0.1058; % Inductance in RHS load (Henrys)
config = 3; % Always start with open circuit load
time = 0; % Time always starts at zero
time_resolution = 0.0001; % Time step
n=1; % Counter of time steps starts at one
Iteration_max = 600000;
I_C = zeros(1,Iteration_max); % Array of cap currents
I_Rs = zeros(1,Iteration_max); % Array of smoothing resistor currents
I_Z = zeros(1,Iteration_max); % Array of Zener currents
I_L1 = zeros(1,Iteration_max); % Array of LHS load currents
I_L2 = zeros(1,Iteration_max); % Array of RHS load currents
I_Total = zeros(1,Iteration_max); % Array of the total current
V1=zeros(1,Iteration_max); % Array of voltages at node 1
V2=zeros(1,Iteration_max); % Array of voltages at node 2
V3=zeros(1,Iteration_max); % Array of voltages at node 3
input_voltage = zeros(1,Iteration_max); % Array of voltages after

%%  If statements determine which state and therefore which differential equation(s) apply.
% target_val contains the variables being refined by the 4th order runge kutta
% target_coeff contains the coefficients relating to these variables
% constant_val are the other terms in each equation
% config 1: Neither load selected
% config 2: Left hand load (resistive) selected
% config 3: Right hand load (inductive) selected
% condig 4: Both loads selected
tic
% for post = 1:3
%     if post == 1
%         R_S = 0.9*R_SMain;
%         C_smoothing = 0.8*C_smoothingMain;
%     elseif post == 3
%        R_S = 1.1*R_SMain;
%        C_smoothing = 1.2*C_smoothingMain;
%     else 
%        R_S = R_SMain;
%        C_smoothing = C_smoothingMain;
%     end
while (n<Iteration_max)
%     if ~(n > 199990)
%         config = 1;
%     elseif (n > 199990) && ~(n > 419990)
%         config = 2;
%     elseif (n > 419990) && ~(n > 546290)
%         config = 4;
%     else
%         config = 3;
%     end
    % iterate voltage, abs represents rectifed value
    input_voltage(n) = abs(Voltage_RHS*sin(100*pi*time));
    if config == 1 % config 1 is where no load is connected
        if ~(input_voltage(n) > 2*V_t_Diode) && (V1(n) == 0) && (V2(n) == 0) && (V3(n) == 0)
            % Condition 1 - Recitifier diodes not conducting, therefore all internal voltages and currents are at zero (startup only)
            V1(n+1)=0;
            V2(n+1)=0;
            V3(n+1)=0;
            I_Rs(n+1) = (V1(n+1)-V2(n+1))/R_S;
            I_C(n+1) = 0;
            I_L1(n+1) = 0;
            I_L2(n+1) = 0;
        elseif ((input_voltage(n) > 2*V_t_Diode) && ((input_voltage(n)-V1(n)) > V_t_Diode) && ((V2(n)-V3(n)) < 6.5)) % Condition 2
            % Condition 2 - Capacitor is being charged, but Zener is off.
            % rk4_PS is called to solve the diff equation
            V1_3 = rk4_PowerSupply(V1(n)-V3(n), -1/(2*C_smoothing*R_d), (input_voltage(n)-(2*V_t_Diode))/(2*C_smoothing*R_d), time_resolution);
            % Using this to solve for voltages and currents
            V1(n+1) = input_voltage(n)-2*V_t_Diode-2*(R_d)*((input_voltage(n)-2*V_t_Diode-V1_3)/(2*R_d));
            V2(n+1) = V1(n+1);
            V3(n+1) = 0;
            I_Rs(n+1) = (V1(n+1)-V2(n+1))/R_S;
            I_C(n+1) = -C_smoothing * ((V1(n+1)-V3(n+1))-(V1(n)-V3(n)))*1/time_resolution;
            % Remainder of currents are at zero as load is not connected and Zener is off.
            I_L1(n+1) = 0;
            I_L2(n+1) = 0;
        elseif ((input_voltage(n)-V1(n)) > 2*V_t_Diode) && ((V2(n)-V3(n)) > 6.5) % Condition 3
            % Condition 3 - Still charging the capacitor, but the zener diode is now in the reverse conducting mode
            % rk4_PS is called to solve the diff equation
            V1_3 = rk4_PowerSupply(V1(n)-V3(n), -1/(C_smoothing*(R_S+R_Z)) -1/(C_smoothing*2*R_d), V_t_Zener/(C_smoothing*(R_S+R_Z)) + (-2*V_t_Diode+input_voltage(n))/(C_smoothing*2*R_d), time_resolution);
            % Using this to solve for voltages and currents
            V1(n+1) = input_voltage(n)-2*V_t_Diode-2*R_d*((-V1_3+input_voltage(n)-2*V_t_Diode)/(2*R_d));
            V2(n+1) = V1(n+1)-((V1_3-V_t_Zener)/(R_S+R_Z))*R_S;
            V3(n+1) = 0;
            I_Rs(n+1) = (V1(n+1)-V2(n+1))/R_S;
            I_C(n+1) = -C_smoothing * ((V1(n+1)-V3(n+1))-(V1(n)-V3(n)))*1/time_resolution;
            % Still no load connected, so load currents are zero
            I_L1(n+1) = 0;
            I_L2(n+1) = 0;
        elseif (~((input_voltage(n)-V1(n)) > 2*V_t_Diode) || ~(V3(n) > V_t_Diode)) && ((V2(n)-V3(n)) > V_t_Zener) % Condition 4
            % Condition 4 - The capacitor is being discharged through the zener to ground and the input is no longer connected
            % Both diodes are turned off if one is turned off
            % rk4_PS is called to solve the diff equation
            V1_3 = rk4_PowerSupply(V1(n), -1/(C_smoothing*(R_S+R_Z)), V_t_Zener/(C_smoothing*(R_S+R_Z)), time_resolution);
            % Using this to solve for voltages and currents
            V1(n+1) = V1_3;
            V2(n+1) = V1_3+R_S*((-V1_3+V_t_Zener)/(R_S+R_Z));
            V3(n+1) = 0;
            I_Rs(n+1) = (V1(n+1)-V2(n+1))/R_S;
            I_C(n+1) = -C_smoothing * ((V1(n+1)-V3(n+1))-(V1(n)-V3(n)))*1/time_resolution;
            % Still no load connected
            I_L1(n+1) = 0;
            I_L2(n+1) = 0;
        else % Condition 5 - Voltage across capacitor no longer enough to exceed Zener voltage
            % voltages no longer change
            V1(n+1)=V1(n);
            V2(n+1)=V2(n);
            V3(n+1)=V3(n);
            % load currents are still at zero as the loads are open circuit
            I_Rs(n+1) = (V1(n+1)-V2(n+1))/R_S;
            I_C(n+1) = -C_smoothing * ((V1(n+1)-V3(n+1))-(V1(n)-V3(n)))*1/time_resolution;
            I_L1(n+1) = 0;
            I_L2(n+1) = 0;
        end
        n = n+1;
        time = time+time_resolution;
    elseif config == 2 % config 2 is where LHS load is connected
        if ~(input_voltage(n) > 2*V_t_Diode) && (V1(n) == 0) && (V2(n) == 0) && (V3(n) == 0) % Condition 0
            % Condition 1 - Recitifier diodes not conducting, therefore all internal voltages and currents are at zero, as this state is only on startup
            V1(n+1) = 0;
            V2(n+1) = 0;
            V3(n+1) = 0;
            I_Rs(n+1) = (V1(n+1)-V2(n+1))/R_S;
            I_C(n+1) = 0;
            I_L1(n+1) = 0;
            I_L2(n+1) = 0;
        elseif ((input_voltage(n) > 2*V_t_Diode) && ((input_voltage(n)-V1(n)) > V_t_Diode) && ((V2(n)-V3(n)) < 6.5)) % Condition 2
            % Condition 2 - Capacitor is being charged, but Zener is off
            % rk4_PS is called to solve the diff equation
            V1_3 = rk4_PowerSupply(V1(n)-V3(n),-1/(C_smoothing*(R_S+R_L1))-1/(C_smoothing*2*R_d), (-2*V_t_Diode+input_voltage(n))/(2*C_smoothing*R_d), time_resolution);
            % Using this to solve for voltages and currents
            V1(n+1) = input_voltage(n)-2*V_t_Diode-2*(R_d)*((input_voltage(n)-2*V_t_Diode-V1_3)/(2*R_d));
            V2(n+1) = V1(n+1)-R_S*((V1_3)/(R_S+R_L1));
            V3(n+1) = 0;
            I_Rs(n+1) = (V1(n+1)-V2(n+1))/R_S;
            I_C(n+1) = -C_smoothing * ((V1(n+1)-V3(n+1))-(V1(n)-V3(n)))*1/time_resolution;
            % RH load is not connected but LH is.
            I_L1(n+1) = I_Rs(n+1);
            I_L2(n+1) = 0;
        elseif ((input_voltage(n)-V1(n)) > 2*V_t_Diode) && ((V2(n)-V3(n)) > 6.5) % Condition 3
            % Condition 3 - Still charging the capacitor, but the zener diode is now in the reverse conducting mode
            % rk4_PS is called to solve the diff equation
            V1_3 = rk4_PowerSupply(V1(n)-V3(n),-1/(2*C_smoothing*R_d) -(R_Z+R_L1)/(C_smoothing*(R_Z*R_L1+R_S*(R_Z+R_L1))),(input_voltage(n)-2*V_t_Diode)/(2*C_smoothing*R_d)+(V_t_Zener*R_L1)/(C_smoothing*(R_Z*R_L1+R_S*(R_Z+R_L1))), time_resolution);
            % Using this to solve for voltages and currents
            V1(n+1) = input_voltage(n)-2*V_t_Diode-2*(R_d)*((input_voltage(n)-2*V_t_Diode-V1_3)/(2*R_d));
            V2(n+1) = V1(n+1)-(((V1_3)*(R_Z+R_L1)-V_t_Zener*R_L1)/(R_Z*R_L1 + R_S*(R_Z+R_L1)))*R_S;
            I_Rs(n+1) = (V1(n+1)-V2(n+1))/R_S;
            I_C(n+1) = -C_smoothing * ((V1(n+1)-V3(n+1))-(V1(n)-V3(n)))*1/time_resolution;
            % RH load is not connected but LS is.
            I_L1(n+1) = (V2(n+1)-V3(n+1))/R_L1;
            I_L2(n+1) = 0;
        elseif (~((input_voltage(n)-V1(n)) > 2*V_t_Diode)) && ((V2(n)-V3(n)) > V_t_Zener) % Condition 4
            % Condition 4 - The capacitor is being discharged through the zener to ground and the input is no longer connected.
            % Both diodes are turned off if one is turned off
            % rk4_PS is called to solve the diff equation
            V1_3 = rk4_PowerSupply(V1(n)-V3(n),-(R_Z+R_L1)/(C_smoothing*(R_Z*R_L1+R_S*(R_Z+R_L1))),(V_t_Zener*R_L1)/(C_smoothing*(R_Z*R_L1+R_S*(R_Z+R_L1))), time_resolution);
            % Using this to solve for voltages and currents
            V1(n+1) = V1_3;
            V2(n+1) = V1_3-(((V1_3)*(R_Z+R_L1)-V_t_Zener*R_L1)/(R_Z*R_L1 + R_S*(R_Z+R_L1)))*R_S;
            V3(n+1) = 0;
            % RH load is not connected but Zener is still on
            I_Rs(n+1) = (V1(n+1)-V2(n+1))/R_S;
            I_C(n+1) = -C_smoothing * ((V1(n+1)-V3(n+1))-(V1(n)-V3(n)))*1/time_resolution;
            I_L1(n+1) = (V2(n+1)-V3(n+1))/R_L1;
            I_L2(n+1) = 0;
        else % Condition 5 - Voltage across capacitor no longer enough to exceed Zener voltage
            % Voltages no longer change
            V1_3 = rk4_PowerSupply(V1(n)-V3(n),-1/(C_smoothing*(R_S+R_L1)),0,time_resolution);
            V1(n+1)=V1_3;
            V2(n+1)=V1_3-(V1_3)/(R_S+R_L1);
            V3(n+1)=0;
            % Currents through zener are now at zero, path to ground still exists through load
            I_Rs(n+1) = (V1(n+1)-V2(n+1))/R_S;
            I_C(n+1) = -C_smoothing * ((V1(n+1)-V3(n+1))-(V1(n)-V3(n)))*1/time_resolution;
            I_L1(n+1) = I_Rs(n+1);
            I_L2(n+1) = 0;
        end
        % iterate timestep
        n = n+1;
        time = time+time_resolution;
        
    elseif config == 3 % missing comments
        if ~(input_voltage(n) > 2*V_t_Diode) && (V1(n) == 0) && (V2(n) == 0) && (V3(n) == 0) % Condition 1
            V1(n+1)=0;
            V2(n+1)=0;
            V3(n+1)=0;
            I_Rs(n+1) = (V1(n+1)-V2(n+1))/R_S;
            I_C(n+1) = 0;
            I_L1(n+1) = 0;
            I_L2(n+1) = 0;
        elseif (((input_voltage(n)-V1(n)) > 2*V_t_Diode) && ((V2(n)-V3(n)) < 6.5)) % Condition 2
            % Condition 2 - Capacitor is being charged, but Zener is off
            % Set up vectors to pass to rk4_PS
            target_val = [V1(n)-V3(n);I_L2(n)];
            target_coeff = [-1/(2*C_smoothing*R_d),-1/C_smoothing;1/L,(-(R_S+R_L2))/L];
            constant_val = [(input_voltage(n)-2*V_t_Diode)/(2*C_smoothing*R_d);0];
            % rk4_PS is called to solve the diff equation
            Results = rk4_PowerSupply(target_val, target_coeff, constant_val, time_resolution);
            V1(n+1) = input_voltage(n)-2*V_t_Diode;
            V2(n+1) = V1(n+1)-(R_S)*Results(2);
            V3(n+1) = 0;
            % RHS load is now connected, and its value comes from the evaluation of runge kutta.
            I_Rs(n+1) = (V1(n+1)-V2(n+1))/R_S;
            I_C(n+1) = (input_voltage(n)-2*V_t_Diode-Results(1))/(2*R_d)-Results(2);
            I_L1(n+1) = 0;
            I_L2(n+1) = Results(2);
        elseif ((input_voltage(n)-V1(n)) > 2*V_t_Diode) && ((V2(n)-V3(n)) > V_t_Zener) % Condition 3
            % Still charging the capacitor, but the zener diode is now in the reverse conducting mode
            % Set up vectors to pass to rk4_PS
            target = [V1(n)-V3(n);I_L2(n)];
            target_coeff = [-1/(C_smoothing*(R_Z+R_S))-1/(2*C_smoothing*R_d),R_Z/(C_smoothing*(R_Z+R_S));R_Z/(L*(R_Z+R_S)),(-R_L2*(R_Z+R_S)+R_S*R_Z)/(L*(R_Z+R_S))];
            constant_val = [(V_t_Zener)/(C_smoothing*(R_Z+R_S))+(input_voltage(n)-2*V_t_Diode)/(2*C_smoothing*R_d);(R_S*V_t_Zener)/(L*(R_Z+R_S))];
            % rk4_PS is called to solve the diff equation
            Results = rk4_PowerSupply(target, target_coeff, constant_val, time_resolution);
            V1(n+1) = input_voltage(n)-2*V_t_Diode;
            V2(n+1) = V1(n+1)-R_S*((Results(1)-V_t_Zener-R_Z*Results(2))/(R_Z+R_S));
            V3(n+1) = 0;
            % I_C calculated from RK4 output
            I_Rs(n+1) = (V1(n+1)-V2(n+1))/R_S;
            I_C(n+1) = (input_voltage(n)-2*V_t_Diode-Results(1))/(2*R_d)+(-Results(1)+V_t_Zener+R_Z*Results(2))/(R_Z+R_S);
            % RH load connected, LH disconnected still
            I_L1(n+1) = 0;
            I_L2(n+1) = Results(2);
        elseif (~((input_voltage(n)-V1(n)) > 2*V_t_Diode)) && ((V2(n)-V3(n)) > 6.5) % Condition 4
            % The capacitor is being discharged through the zener to ground
            % Both diodes are turned off if one is turned off
            % Set up vectors to pass to rk4_PS
            target = [V1(n)-V3(n);I_L2(n)];
            target_coeff = [-1/(C_smoothing*(R_Z+R_S)),R_Z/(C_smoothing*(R_Z+R_S));R_Z/(L*(R_Z+R_S)),(-R_L2*(R_Z+R_S)+R_S*R_Z)/(L*(R_Z+R_S))];
            constant_val = [V_t_Zener/(C_smoothing*(R_Z+R_S));(R_S*V_t_Zener)/(L*(R_Z+R_S))];
            % rk4_PS is called to solve the diff equation
            Results = rk4_PowerSupply(target, target_coeff, constant_val, time_resolution);
            V1(n+1) = Results(1);
            V2(n+1) = V1(n+1)-R_S*((Results(1)-V_t_Zener-R_Z*Results(2))/(R_Z+R_S));
            V3(n+1) = 0;
            I_Rs(n+1) = (V1(n+1)-V2(n+1))/R_S;
            I_C(n+1) = -(-Results(1)+V_t_Zener+R_Z*Results(2))/(R_Z+R_S);
            I_L1(n+1) = 0;
            I_L2(n+1) = Results(2);
        else % Condition 5
            % Zener not activated, so path to ground exists through the load. Cap discharges
            % Set up vectors to pass to rk4_PS
            target = [V1(n)-V3(n);-I_L2(n)/C_smoothing];
            target_coeff = [0,1;-1/(L*C_smoothing),-(R_S+R_L2)/L];
            constant_val = [0;0];
            % rk4_PS is called to solve the diff equation
            Results = rk4_PowerSupply(target, target_coeff, constant_val, time_resolution);
            % need load current to calculate V2, found in RK4
            I_L2(n+1) = Results(2)*-C_smoothing;
            V1(n+1)=Results(1);
            V2(n+1)=V1(n+1)-R_S*I_L2(n+1);
            V3(n+1)=0;
            I_Rs(n+1) = (V1(n+1)-V2(n+1))/R_S;
            I_C(n+1) = -C_smoothing*Results(1);
            I_L1(n+1) = 0;
        end
        % iterate timestep
        n = n+1;
        time = time+time_resolution;
    elseif config == 4 % config 4 is where LHS and RHS load are connected
        % missing voltage and current equations
        if ~(input_voltage(n) > 2*V_t_Diode) && (V1(n) == 0) && (V2(n) == 0) && (V3(n) == 0 )  % Condition 1
            % Condition 1 - Recitifier diodes not conducting, therefore all internal voltages and currents are at zero
            V1(n+1)=0;
            V2(n+1)=0;
            V3(n+1)=0;
            I_Rs(n+1) = (V1(n+1)-V2(n+1))/R_S;
            I_C(n+1) = 0;
            I_L1(n+1) = 0;
            I_L2(n+1) = 0;
        elseif ((input_voltage(n)-V1(n) > 2*V_t_Diode) && ((V2(n)-V3(n)) < 6.5)) % Condition 2
            % Condition 2 - Capacitor is being charged, but Zener is off
            target_val = [V1(n)-V3(n);I_L2(n)];
            target_coeff = [-1/(2*C_smoothing*R_d)-1/(C_smoothing*(R_L1+R_S)),-R_L1/(C_smoothing*(R_L1+R_S));R_L1/(L*(R_L1*R_S)),-(R_S*R_L1-R_L2*(R_L1+R_S))/(L*(R_L1*R_S))];%(-1+R_L1+R_S)/(L*(R_L1+R_S)),(1-R_L2*(R_L1+R_S))/(L*(R_L1+R_S))
            constant_val = [(input_voltage(n)-2*V_t_Diode)/(2*C_smoothing*R_d);0];
            % rk4_PS is called to solve the diff equation
            Results = rk4_PowerSupply(target_val, target_coeff, constant_val, time_resolution);
            % Using this to solve for voltages and currents
            V1(n+1) = input_voltage(n)-2*V_t_Diode-2*(R_d)*((input_voltage(n)-2*V_t_Diode-Results(1))/(2*R_d));
            V2(n+1) = V1(n+1)-R_S*((Results(1)+R_L1*Results(2))/(R_L1+R_S));
            V3(n+1) = 0;
            % LHS load is an open circuit
            I_Rs(n+1) = (V1(n+1)-V2(n+1))/R_S;
            I_C(n+1) = 0; % missing
            I_L1(n+1) = (V2(n+1)+V3(n+1))/R_L1;
            % RHS load current comes from the RK4 again
            I_L2(n+1) = Results(2);
        elseif ((input_voltage(n)-V1(n)) > 2*V_t_Diode) && ((V2(n)-V3(n)) > 6.5 ) % Condition 3
            % Condition 3 - Still charging the capacitor, but the zener diode is now in the reverse conducting mode
            % Setting up the input vectors
            target_val = [V1(n)-V3(n);I_L2(n)];
            target_coeff = [-(R_L1+R_Z)/(C_smoothing*(R_Z*R_L1+R_S*R_L1+R_S*R_Z)),-(R_Z*R_L1)/(C_smoothing*(R_Z*R_L1+R_S*R_L1+R_S*R_Z));1/L-R_S*((R_L1+R_Z)/(L*(R_Z*R_L1+R_S*R_L1+R_S*R_Z))),-R_S*((R_L1*R_Z)/(L*(R_Z*R_L1+R_S*R_L1+R_S*R_Z)))-R_L1/L];
            constant_val = [(V_t_Zener*R_L1)/(C_smoothing*(R_Z*R_L1+R_S*R_L1+R_S*R_Z))+(input_voltage(n)-2*V_t_Diode)/(2*C_smoothing*R_d);-R_S*(-V_t_Zener*R_L1)/(L*(R_Z*R_L1+R_S*R_L1+R_S*R_Z))];
            % rk4_PS is called to solve the diff equation
            Results = rk4_PowerSupply(target_val, target_coeff, constant_val, time_resolution);
            % Using this to solve for voltages and currents
            V1(n+1) = input_voltage(n)-2*V_t_Diode;%-(R_d)*((input_voltage(n)-2*V_t_Diode-Results(1))/(2*R_d));
            V2(n+1) = V1(n+1)-R_S*((Results(1)*(R_L1+R_Z)+Results(2)*(R_Z*R_L1)-V_t_Zener)/(R_Z*R_L1+R_S*R_L1+R_S*R_Z));
            V3(n+1) = 0;
            % Both loads are now connected. Current through Load 2 comes from the RK4. Load 1 is calculated from voltage drop across resistor
            I_Rs(n+1) = (V1(n+1)-V2(n+1))/R_S;
            I_C(n+1) = 0; % missing
            I_L1(n+1) = (V2(n+1)+V3(n+1))/R_L1;
            I_L2(n+1) = Results(2);
        elseif (~((input_voltage(n)-V1(n)) > 2*V_t_Diode)) && ((V2(n)-V3(n)) > 6.5 ) % Condition 4&& isequal(flag,0)
            % Condition 4 - The capacitor is being discharged through the zener to ground and the input is no longer connected.
            % Setting up the input vectors
            target_val = [V1(n)-V3(n);I_L2(n)];
            target_coeff = [-(R_L1+R_Z)/(C_smoothing*(R_Z*R_L1+R_S*R_L1+R_S*R_Z)),(-1*(R_Z*R_L1))/(C_smoothing*(R_Z*R_L1+R_S*R_L1+R_S*R_Z));1/L-R_S*((R_L1+R_Z)/(L*(R_Z*R_L1+R_S*R_L1+R_S+R_L1))),-R_S*(R_Z*R_L1)/(L*(R_Z*R_L1+R_S*R_L1+R_S+R_L1))-R_L2];
            constant_val = [(V_t_Zener*R_L1)/(C_smoothing*(R_Z*R_L1+R_S*R_L1+R_S*R_Z));-R_S*(-V_t_Zener*R_L1)/(L*(R_Z*R_L1+R_S*R_L1+R_S+R_L1))];
            % rk4_PS is called to solve the diff equation
            Results = rk4_PowerSupply(target_val, target_coeff, constant_val, time_resolution);
            % Using this to solve for voltages and currents
            V1(n+1) = Results(1);
            V2(n+1) = V1(n+1)-R_S*((Results(1)*(R_L1+R_Z)+Results(2)*(R_Z*R_L1)-V_t_Zener)/(R_Z*R_L1+R_S*R_L1+R_S*R_Z));
            V3(n+1) = 0;
            I_Rs(n+1) = (V1(n+1)-V2(n+1))/R_S;
            I_C(n+1) = 0; % missing
            % Both loads are now connected. Current through Load 2 comes from the RK4. Load 1 is calculated from voltage drop across resistor
            I_L1(n+1) = (V2(n+1)+V3(n+1))/R_L1;
            I_L2(n+1) = Results(2);
        else
            % Condition 5 - Voltage across capacitor no longer enough to exceed Zener voltage
            % Voltages no longer change
            % Setting up the input vectors
            target_val = [V1(n)-V3(n);I_L2(n)];
            target_coeff = [-1/(C_smoothing*(R_L1+R_S)),-1*R_L1/(C_smoothing*(R_L1+R_S));R_L1/(L*(R_L1+R_S)),-((R_L2*(R_L1+R_S)+R_S*R_L1)/(L*(R_L1+R_S)))];
            constant_val = [0;0];
            % rk4_PS is called to solve the diff equation
            Results = rk4_PowerSupply(target_val, target_coeff, constant_val, time_resolution);
            V1(n+1) = Results(1);
            V2(n+1) = V1(n+1)-R_S*((Results(1)+R_L1*Results(2))/(R_L1+R_S));
            V3(n+1) = 0;
%             if V2(n+1) > 6.499
%                 V2(n+1) = 6.499;
%             end
            I_Rs(n+1) = (V1(n+1)-V2(n+1))/R_S;
            I_C(n+1) = (Results(1)+Results(2)*R_L1)/(R_L1+R_S);
            % Both loads are now connected. Current through Load 2 comes from the RK4. Load 1 is calculated from voltage drop across resistor
            I_L1(n+1) = (V2(n+1)+V3(n+1))/R_L1;
            I_L2(n+1) = Results(2);
        end
        % iterate timestep
        n = n+1;
        time = time+time_resolution;
    else
        error('This config does not exist');
    end

% end
%     figure
%     rs = num2str(R_S);
% plot(xaxis,(I_Rs.^2*R_S))
% title(['Power dissipation in resistor R_S (' rs ' \Omega)'])
% xlabel('Time (seconds)')
% ylabel('Power (watts)')
% 1/(R_S*C_smoothing)
% rs = num2str(R_S);
% cs = num2str(C_smoothing);
% figure
% hold on
% plot(xaxis,V2-V3, 'b')
% plot(xaxis,V1-V3, 'g')
% ylim([0,20])
% title(['Capacitor and Output voltages where R_S = ' rs ' \Omega & C_{smoothing} = ' cs ' Farards'])
% xlabel('Time (seconds)')
% ylabel('Voltage (volts)')
% legend('Voltage across smoothing capacitor.','Voltage across the load.')
% hold off
end
toc
xaxis=[1:Iteration_max]*time_resolution;
figure
hold on
plot(xaxis,input_voltage, 'r')
plot(xaxis,V2-V3, 'b')
plot(xaxis,V1-V3, 'g')
ylim([0,20])
title('Configuration 1')
xlabel('Time (seconds)')
ylabel('Voltage (volts)')
hold off
% figure
% plot(xaxis,(I_Rs.^2*R_S))
% title('Power dissipation in resistor R_S')
% xlabel('Time (seconds)')
% ylabel('Power (watts)')
% % figure
% % plot(xaxis,I_L2)
% max_power = max(I_Rs)^2*R_S
% non_spike_max = 0.04541^2*R_S
