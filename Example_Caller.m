time = 0; % Time always starts at zero
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
Sw_1 = 0;
Sw_2 = 0;
time_resolution = 0.00001;
R_S = 270;
C_smoothing = 2200E-6;
Mains_RMS = 230;
Voltage_RHS = (Mains_RMS*(2)^(.5));
Impedances = [R_S,C_smoothing];
number = 6;
start_time = 0;
while (n<2)
    if ~(n > 199990)
        Sw_1 = 0;
        Sw_2 = 0;
    elseif (n > 199990) && ~(n > 419990)
        Sw_1 = 1;
        Sw_2 = 0;
    elseif (n > 419990) && ~(n > 546290)
        Sw_1 = 1;
        Sw_2 = 1;
    else
        Sw_1 = 0;
        Sw_2 = 1;
    end
    input_voltage(n) = Voltage_RHS*sin(100*pi*time);
    V_in = [V1(n),V2(n),V3(n),230*sqrt(2)];
    I_in = [I_C(n),I_Rs(n),I_L1(n),I_L2(n)];
    
    [V_out,I_out] = Power_Supply_Circuit_Solver(V_in,I_in,Impedances,Sw_1,Sw_2,time_resolution, Iteration_max, start_time) ;
    V1 = V_out(1,:);
    V2 = V_out(2,:);
    V3 = V_out(3,:);
    V_in = V_out(4,:);
    I_C = I_out(1,:);
    I_Rs = I_out(2,:);
    I_L1 = I_out(3,:);
    I_L2 = I_out(4,:);
    
    n = n+1;
    time = time+time_resolution;
end
xaxis=[1:Iteration_max]*time_resolution;
figure
hold on
plot(V_in, 'r')
plot(V2-V3, 'b')
plot(V1-V3, 'g')
ylim([0,20])
hold off
figure
plot(I_L2)
figure
plot(I_Rs)
max_power = max(I_Rs)^2*R_S
non_spike_max = 0.04541^2*R_S
warning('current switch pattern is ***********************NO LONGER the A/A+ question, where the time_resolution is set to 0.0001 and Iteration_max is 60,000');
