%% Title Block
%Purpose of code: To plot the co-polarized pattern function of the
%dipole-consisted huygen's source for any constanth theta or phi

clear;
close all;
clc;

%% Variables
d_phi = (pi/180)/10; %rad, differential angle along phi direction in incremenets of 1/100th of a degree
d_theta = d_phi;
l_e_h = [0 1 0]; %unit vector for the electric dipole
l_h_h = [1 0 0]; %unit vector for the magnetic dipole

%% Huygen's Source Analysis

Angle = 1; %1 for theta, 2 for phi

if(Angle == 1)
    theta = 45 * pi/180; %rad, angle specified by user
    phi = 0; %rad, looping through phi angle
    F =  zeros(1,length(0:d_phi:2*pi)); %initializing the scattered E field as a function of theta
    g = 1; %index counter
    while(phi<=2*pi)
        s_i_h = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)]; %radial unit vector from the center of the dipoles to the field point
        E_p = cross(cross(l_e_h, s_i_h), s_i_h);
        H_p = cross(l_h_h, s_i_h);
        F(g) = norm(0.5*H_p+0.5*E_p);        
        g = g + 1; %array counter
        phi = phi + d_phi; %incrementing phi angle by 1/100th of a degree
    end  

    phi_angle = (0:d_phi:2*pi); %rad, used to set the axis of the plot
    %figure;plot(phi_angle, 20*log10(abs(F))); hold all; grid on; ylabel('Magnitude [dB]')
    %xlabel('\phi [rad]'); xlim([0 2*pi]); ylim([-120 0]);
    figure;polarplot(F);title('\theta[deg] = ', theta*180/pi);rlim([0 1])
    
    theta = 135 * pi/180; %rad, angle specified by user
    phi = 0; %rad, looping through phi angle
    F =  zeros(1,length(0:d_phi:2*pi)); %initializing the scattered E field as a function of theta
    g = 1; %index counter
    while(phi<=2*pi)
        s_i_h = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)]; %radial unit vector from the center of the dipoles to the field point
        E_p = cross(cross(l_e_h, s_i_h), s_i_h);
        H_p = cross(l_h_h, s_i_h);
        F(g) = 0.5*norm(E_p+H_p);             
        g = g + 1; %array counter
        phi = phi + d_phi; %incrementing phi angle by 1/100th of a degree
    end  

    %figure;plot(phi_angle, 20*log10(abs(F))); hold all; grid on; ylabel('Magnitude [dB]')
    %xlabel('\phi [rad]'); xlim([0 2*pi]); ylim([-120 0]);
    figure;polarplot(F);title('\theta[deg] = ', theta*180/pi); ;rlim([0 1])
    
elseif(Angle == 2)
    phi = 180*pi/180; %rad, angle specified by user
    theta = 0;%rad, looping through theta angle
    F =  zeros(length(0:d_theta:2*pi),3); %initializing the scattered E field as a function of theta
    g = 1;%index counter
    while(theta<=2*pi)
        s_i_h = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)]; %radial unit vector from the center of the dipoles to the field point
        E_p = cross(cross(l_e_h, s_i_h), s_i_h);
        H_p = cross(l_h_h, s_i_h);
        F(g) = 0.5*norm(E_p+H_p);        
        g = g + 1; %array counter
        theta = theta + d_theta; %incrementing theta angle by 1/100th of a degree
    end
    theta_angle = (0:d_theta:2*pi);
    %figure;plot(theta_angle, 20*log10(abs(F))); hold all; grid on; ylabel('Magnitude [dB]')
    %xlabel('\theta [rad]'); xlim([0 2*pi]);ylim([-120 0]);
    figure;polarplot(F); title('\phi [deg] = ', phi*180/pi);
    
    phi = 270*pi/180; %rad, angle specified by user
    theta = 0;%rad, looping through theta angle
    g = 1;%index counter
    while(theta<=2*pi)
        s_i_h = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)]; %radial unit vector from the center of the dipoles to the field point
        E_p = cross(cross(l_e_h, s_i_h), s_i_h);
        H_p = cross(l_h_h, s_i_h);
        F(g) = 0.5*norm(H_p+E_p);        
        g = g + 1; %array counter
        theta = theta + d_theta; %incrementing theta angle by 1/100th of a degree
    end
    %figure;plot(theta_angle, 20*log10(abs(F))); hold all; grid on; ylabel('Magnitude [dB]')
    %xlabel('\theta [rad]'); xlim([0 2*pi]);ylim([-120 0]);
    figure;polarplot(F); title('\phi [deg] = ', phi*180/pi);
    
else
    disp('Check input')
end


%% Directivity

%Following code is for calculating the directivity of the huygen's source
%via the beam-solid-angle method

theta = 0; %initializing thetas angle
b_S_A = 0; %initializing beam solid angle sum
while(theta<=pi)
    phi =0;
    while(phi<=2*pi)
        s_i_h = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)]; %radial unit vector from the center of the dipoles to the field point
        F = 0.5*norm(cross(cross(l_e_h, s_i_h), s_i_h)+cross(l_h_h, s_i_h));%normalized pattern function of the two dipoles
        b_S_A = b_S_A+(abs(F)^2)*sin(theta)*d_theta*d_phi; %mupltiplying by differential aray, sin(theta)*d_theta*d_phi
        phi = phi + d_phi; %incrementing phi by  1/100th of a degree
    end
    theta = theta + d_theta; %incrementing theta by  1/100th of a degree
end

D_dBi = 10*log10(4*pi/b_S_A)
