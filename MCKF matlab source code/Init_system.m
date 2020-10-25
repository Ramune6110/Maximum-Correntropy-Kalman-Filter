n  = 2;        %number of states
m  = 1;        %number of measurements
N  = 100;     %total number of time steps
F  = [cos(pi/18),-sin(pi/18);sin(pi/18),cos(pi/18)];  %state transition matrix
H  = [1,1];    %observation matrix