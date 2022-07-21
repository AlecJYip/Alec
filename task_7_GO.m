%% Title Block
%Purpose of code: To determine the E field norm from a differential
%patch on the reflector.
%All equations reference antenna theory and design, stuzmann, third edition

clc;
clear;
close all;

%% Variables
f = 1691e6; %Hz
c = 3e8; %m/s
lambda = c/f; %m
beta = 2*pi/lambda; %rad/m
D = 1.8; %m
F = 760/1000; %m
epsilon = 8.85418782e-12; %F/m
mu = (4*pi)*(10^-7); %H/m
eta = sqrt(mu/epsilon); %ohms
omega = 2*pi*f; %rad/m, radial frequency
field_angle = 40*pi/180; %rad, angle observered from the primary axis of parabolic antenna
d_rho_f = lambda/10; %m, differential distance along x_f axis in incremenets of lambda/10 m
d_phi_f = (pi/180)/10; %rad, differential angle along phi_f direction in incremenets of 1/100th of a degree
y_hat = [0 1 0]; %unit vector in cartisean coordinates
y_h = y_hat;%unit vector in cartisean coordinates
x_h = [1 0 0];%unit vector in cartisean coordinates
z_h = [0 0 1];
l_e_h = [0 1 0]; %unit vector for the electric dipole
l_h_h = [-1 0 0]; %unit vector for the magnetic dipole

%% Analysis: Calculating scattered E field from parabolic

E_s_mag = zeros(1, length(0:d_phi_f/10:field_angle)); %initializing the scattered E field as a function of theta
g = 1; %counter for array
theta_c = 0; %rad, angle from reflector axis (referencing global coordinate system)
phi_c = pi/4; %rad, used to set the plane of observation
while(theta_c <= field_angle)
       
    phi_f = 0;     %rad, starting phi_f angle
    
    P = 0; %initializing sum for the integration to obtain the spatial fourier transform of the electrical current
    
    %the following two while loops are the numerical integrator to solve
    %for the scattered E field at a point
    while(phi_f<=2*pi) %rad, iterating from phi_f = 0 to phi = 2*pi
        
        rho = 0; %m, projection of the radial distance onto the aperture plane
        
        while(rho<=D/2) %integrating from the vertix to the rim of the reflector            
            theta_f = -2*atan(rho/(2*F)); %rad, angular distance in the feed fixed coordinate system. The negative sign accounts for the conversion from feed-fixed to global coordinates
            phi_f = -phi_f; %rad, phi_f is -phi in the global coordinate system
            s_i = F*sec(theta_f/2)^2; %m, distance from feed to reflector, Eqn 9-183 in Stuzmann          
            r_f_h = [sin(theta_f)*cos(phi_f) sin(theta_f)*sin(phi_f) cos(theta_f)]; %r_hat to cartisean
            theta_f_h = [cos(theta_f)*cos(phi_f) cos(theta_f)*sin(phi_f) -sin(theta_f)]; %theta_hat to cartisean
            phi_f_h = [-sin(phi_f) cos(phi_f) 0];
            n_h = -cos(theta_f/2)*r_f_h + sin(theta_f/2)*theta_f_h; %normal vector in the global coordinate system (pointing towards the feed, equation 9-188 in Stuzman
            d_a = rho*d_phi_f*d_rho_f; %m^2, differential surface area of circular aperture                 
            z_f = -s_i*cos(theta_f); %m, z corrdinate in global coordinate system
            s_i_v = [rho*cos(phi_f) rho*sin(phi_f) z_f]; %m, vector
            s_i_h = s_i_v/norm(s_i_v); %m, unit vector pointing from the origin (focus) to a point on the reflector)
            
%             s_r_h = z_h;
%             e_perp = cross(n_h,s_r_h)/norm(cross(n_h,s_r_h));
%             if (norm(cross(n_h,s_r_h)) == 0)
%                 e_perp = [0 0 0];
%             end
%             e_p_r = cross(s_r_h,e_perp);
%             e_p_i = cross(s_i_h, e_perp);
%             e_i_h = (cross(cross(l_e_h, s_i_h), s_i_h)+cross(l_h_h, s_i_h))./...
%                 norm(cross(cross(l_e_h, s_i_h), s_i_h)+cross(l_h_h, s_i_h));
%             e_r_h = dot(e_i_h, e_p_i)*e_p_r - dot(e_i_h, e_perp)*e_perp;
%             u_r_h = e_r_h;
            
            E_i = 0.5*(cross(cross(l_e_h, s_i_h), s_i_h)+cross(l_h_h, s_i_h))*exp(-j*beta*s_i)/s_i;
            E_r = 2*dot(n_h,E_i)*n_h-E_i;
            E_a = E_r.*exp(-j*beta*(-z_f))*1;%V/m, aperture E field
            P = P + E_a*d_a*exp(j*beta*rho*sin(theta_c)*cos(phi_c-phi_f)); %unitless, spatial foureier transform of the electric current
            
            rho = rho + d_rho_f; %m, incrementing by the differential
        end
        phi_f = abs(phi_f); %rad, taking the absolute value of the angle for looping
        phi_f = phi_f +  d_phi_f; %incrementating the observation angle (theta) by a 10th of a degree
    end
    Px = dot(x_h,P); %unitless, x comonent of the spatial fourier transform    
    Py = dot(y_h,P); %unitless, y comonent of the spatial fourier transform    
    
    theta_h = [cos(theta_c)*cos(phi_c) cos(theta_c)*sin(phi_c) -sin(theta_c)]; %unitless, theta_hat to cartisean in global coordinate system
    phi_h = [-sin(phi_c) cos(phi_c) 0]; %unitless, phi_hat to cartisean in the global coordinate system    

    E_field = 0.5*beta*(1/(2*pi))*(1+cos(theta_c))*(theta_h*(Px*cos(phi_c)+Py*sin(phi_c))+...
        phi_h*(Py*cos(phi_c)-Px*sin(phi_c))); %V/m, normalized expression for the reflected E field from the parabolic
    
    E_s_mag(g) = norm(E_field); %V/m, magnitude of the E field
    theta_c = theta_c + d_phi_f/10; %rad, incrementating the observation angle (theta) by a 100th of a degree
    g = g + 1;    %incremeneting the array index
end
E_s_squarred = (E_s_mag).^2; %magntiude of the scattered e field squarred

%% Analysis: Calculating beam solid angle
% %Used in the gain pattern for recreating figure 3
theta = 0; %initializing thetas angle
b_S_A = 0; %initializing beam solid angle sum
d_phi = d_phi_f; %rad, differential angle used in integration
d_theta = d_phi; %rad, differential angle used in integration
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
solid_angle = b_S_A; %sr
norm_E = (solid_angle/(4*pi)); %magnitude of E (divided by the source current) from feed squarred multiplied by the beam solid angle and divided by 4*pi

%% Plotting
angle_rim = (0:d_phi_f/10:field_angle); %rad, used in plotting the horizontal axis of the gain plot
D_dBi = 10*log10(E_s_squarred/norm_E); %power pattern of antenna
figure;plot(angle_rim*180/pi,D_dBi);hold all;xlabel('\theta [deg]');xlim([0 field_angle*180/pi]); %2 degree range on horizontal axis
ylabel('Co-polairzed Pattern [dBi]');title('1.8m Results with GOAI');grid on;