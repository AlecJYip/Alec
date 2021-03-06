%% Title Block
%Purpose of code: To determine the E field norm from a differential
%patch on the reflector.
%All equations reference antenna theory and design, stuzmann, third edition

clc;
clear;
close all;

%% Variables
f = 1.5e9; %Hz
c = 3e8; %m/s
lambda = c/f; %m
beta = 2*pi/lambda; %rad/m
D = 18; %m
F = 0.4*D; %m
epsilon = 8.85418782e-12; %F/m
mu = (4*pi)*(10^-7); %H/m
eta = sqrt(mu/epsilon); %ohms
omega = 2*pi*f; %rad/m, radial frequency
field_angle = 2*pi/180; %rad, angle observered from the primary axis of parabolic antenna
d_rho_f = lambda/10; %m, differential distance along x_f axis in incremenets of lambda/10 m
d_phi_f = (pi/180)/1; %rad, differential angle along phi_f direction in incremenets of 1/100th of a degree
y_hat = [0 1 0]; %unit vector in cartisean coordinates
l_e_h = [1 0 0]; %unit vector for the electric dipole
l_h_h = [0 1 0]; %unit vector for the magnetic dipole
d_A = F*d_phi_f*d_rho_f; %m^2, differential patch size for the integration
%0.3*lambda*0.3*lambda
%% Analysis: Calculating scattered E field from parabolic
E_s_mag = zeros(1, length(0:d_phi_f/10:field_angle)); %initializing the scattered E field as a function of theta
g = 1; %counter for array
theta = 0; %rad, angle from reflector axis (referencing global coordinate system)
while(theta <= field_angle)
    
    E_field = [0 0 0];  %V/m, initializing scattered E field for this computiation
    
    phi_f = 0;     %rad, starting phi_f angle
    
    %the following two while loops are the numerical integrator to solve
    %for the scattered E field at a point
    while(phi_f<=2*pi) %rad, iterating from phi_f = 0 to phi = 2*pi
        
        rho = 0; %m, projection of the radial distance onto the aperture plane
        
        while(rho<=D/2) %integrating from the vertix to the rim of the reflector
            
            theta_f = -2*atan(rho/(2*F)); %rad, angular distance in the feed fixed coordinate system. The negative sign accounts for the conversion from feed-fixed to global coordinates
            phi_f = -phi_f; %rad, phi_f is -phi in the global coordinate system
            s_i = F*sec(theta_f/2)^2; %m, distance from feed to reflector, Eqn 9-183 in Stuzmann
            z_f = -s_i*cos(theta_f); %m, z corrdinate in global coordinate system
            s_i_v = [rho*cos(phi_f) rho*sin(phi_f) z_f]; %m, vector pointing from the feed to a point on the reflector in the global coordinate system
            
            %recall that that y_hat = -y_f_hat, z_hat = -z_f_hat, and x_hat
            % = x_f_hat in the feed-fixed coordinate ststem
            
            %Unit Vectors
            s_i_h = s_i_v/norm(s_i_v); %m, unit vector pointing from the origin (focus) to a point on the reflector)
            r_f_h = [sin(theta_f)*cos(phi_f) sin(theta_f)*sin(phi_f) cos(theta_f)]; %r_hat to cartisean
            theta_f_h = [cos(theta_f)*cos(phi_f) cos(theta_f)*sin(phi_f) -sin(theta_f)]; %theta_hat to cartisean
            n_h = -cos(theta_f/2)*r_f_h + sin(theta_f/2)*theta_f_h; %normal vector in the global coordinate system (pointing towards the feed, equation 9-188 in Stuzman
            
            d_a = sqrt(4*F^2+rho^2)/(2*F)*...
                rho*d_phi_f*d_rho_f; %m, differential surface are of the parabolic 16-165 in Stutzman
            
            E_I =  0.5*(cross(cross(l_e_h, s_i_h), s_i_h)+cross(l_h_h, s_i_h))*exp(-j*beta*s_i)/s_i; %normalized pattern function of the two dipoles
                        
            H_I  =  cross((1/eta)*s_i_h, E_I); %A/m, magnetic field normalized without source current Equation 2 in AWPL
            
            J = 2*cross(n_h, H_I); %A/m, equivalent surface current density accross surface of parabolic
            
            phi = -phi_f;%rad, phi_f is -phi in the global coordinate system
            r_hat = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)]; %r_hat to cartisean in the global coordinate system
            phase = beta*dot(r_hat, s_i_v);%m, phase incurred by far-field approximation dot product in Equation 1 in AWPL
            
            integrand = -j*J*omega*mu*...
                exp(j*phase)*...
                (1/(4*pi))*d_a; %V/m, integrand in calculating the scattered e field
            
            E_field = E_field + integrand; %V/m, adding the infantesimal sum of the scattered e field
            rho = rho + d_rho_f; %m, incrementing by the differential
        end
        phi_f = abs(phi_f); %rad, taking the absolute value of the angle for looping
        phi_f = phi_f +  d_phi_f; %incrementating the observation angle (theta) by a 10th of a degree
    end
    E_s_mag(g) = norm(E_field); %V/m, magnitude of the E field
    theta = theta + d_phi_f/10; %rad, incrementating the observation angle (theta) by a 100th of a degree
    g = g + 1;    %incremeneting the array index
end

E_s_squarred = abs(E_s_mag).^2; %magntiude of the scattered e field squarred

%% Analysis: Calculating beam solid angle
% %Used in the gain pattern for recreating figure 3
theta = 0; %initializing thetas angle
b_S_A = 0; %initializing beam solid angle sum
d_phi = d_phi_f;
d_theta = d_phi;
l_e_h = [0 0 1]; %unit vector for the electric dipole
l_h_h = [1 0 0]; %unit vector for the magnetic dipole
while(theta<=pi/2)
    phi =0;
    while(phi<=2*pi)
        s_i_h = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)]; %radial unit vector from the center of the dipoles to the field point
        F = 0.5*norm(cross(cross(l_e_h, s_i_h), s_i_h)+cross(l_h_h, s_i_h));%normalized pattern function of the two dipoles
        b_S_A = b_S_A+(abs(F)^2)*sin(theta)*d_theta*d_phi; %mupltiplying by differential aray, sin(theta)*d_theta*d_phi
        phi = phi + d_phi; %incrementing phi by  1/100th of a degree
    end
    theta = theta + d_theta; %incrementing theta by  1/100th of a degree
end
solid_angle = b_S_A%sr
norm_E = (solid_angle/(4*pi)); %magnitude of E (divided by the source current) from feed squarred multiplied by the beam solid angle and divided by 4*pi

%% Plotting
angle_rim = (0:d_phi_f/10:field_angle); %rad, used in plotting the horizontal axis of the gain plot
D_dBi = 10*log10(E_s_squarred/norm_E); %power pattern of antenna
figure;plot(angle_rim*180/pi,D_dBi);hold all;xlabel('\theta [deg]');xlim([0 2]); %2 degree range on horizontal axis
ylabel('Co-polairzed Pattern [dBi]');title('Figure 3 with Huygen''s Source');grid on;
