clear;
close all;
clc;

%% Variables
d_phi = (pi/180)/10; %rad, differential angle along phi direction in incremenets of 1/100th of a degree
d_theta = d_phi;
l_e_h = [0 1 0]; %unit vector for the electric dipole
l_h_h = [1 0 0]; %unit vector for the magnetic dipole

%% Huygen's Source Analysis

Angle = 1; 

if(Angle == 1)
    theta = 0 * pi/180; %rad
    phi = 0; %rad, looping through phi angle
    F =  zeros(length(0:d_phi:2*pi),1); %initializing the scattered E field as a function of theta
    g = 1; %index counter
    while(phi<=2*pi)
        s_i_h = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)]; %radial unit vector from the center of the dipoles to the field point
        F(g,:) = 0.5*norm(cross(cross(l_e_h, s_i_h), s_i_h)+cross(l_h_h, s_i_h)); %normalized pattern function of the two dipoles
        g = g + 1;
        phi = phi + d_phi;
    end
    

    phi_angle = (0:d_phi:2*pi); %rad, used to set the axis of the plot
    figure;plot(phi_angle, 20*log10((F))); hold all; grid on; ylabel('Magnitude [dB]')
    xlabel('\phi [rad]'); xlim([0 2*pi]); ylim([-120 0]);
    figure;polarplot(F);
    
elseif(Angle == 2)
    phi = 70*pi/180; %rad
    theta = 0;%rad, looping through theta angle
    F =  zeros(length(0:d_theta:2*pi),1); %initializing the scattered E field as a function of theta
    g = 1;%index counter
    while(theta<=2*pi)
        s_i_h = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)]; %radial unit vector from the center of the dipoles to the field point
        F(g,:) = 0.5*norm(cross(cross(l_e_h, s_i_h), s_i_h)+cross(l_h_h, s_i_h)); %normalized pattern function of the two dipoles
        g = g + 1;
        theta = theta + d_theta;
    end
    theta_angle = (0:d_theta:2*pi);
    figure;plot(theta_angle, 20*log10((F))); hold all; grid on; ylabel('Magnitude [dB]')
    xlabel('\theta [rad]'); xlim([0 2*pi]);ylim([-120 0]);
    figure;polarplot(F);
    
else
    disp('Check input')
end


%% Directivity

%Directivity
theta = 0;
range_theta = (0:d_theta:pi);
range_phi = (0:d_phi:2*pi);
arr = zeros(1,length(range_theta).*length(range_phi));
n = 1;
while(theta<=pi)
    phi =0;
    while(phi<=2*pi)
        s_i_h = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)]; %radial unit vector from the center of the dipoles to the field point
        F = 0.5*norm(cross(cross(l_e_h, s_i_h), s_i_h)+cross(l_h_h, s_i_h));%normalized pattern function of the two dipoles
        func = sin(theta)*(abs(F)^2);
        arr(n) = func*d_theta*d_phi;
        phi = phi + d_phi;
        n = n + 1;
    end
    theta = theta + d_theta;
end

solid_angle = sum(arr);
D_dBi = 10*log10(4*pi/solid_angle)
