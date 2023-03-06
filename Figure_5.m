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
q = 1.14;
D = 18; %m
F = 0.4*D; %m
epsilon = 8.85418782e-12; %F/m
mu = (4*pi)*(10^-7); %H/m
eta = sqrt(mu/epsilon); %ohms
omega = 2*pi*f; %rad/m, radial frequency
field_angle = 2*pi/180; %rad, angle observered from the primary axis of parabolic antenna
d_rho_f = lambda/15; %m, differential distance along x_f axis in incremenets of lambda/10 m
d_phi_f = (pi/180)/5; %rad, differential angle along phi_f direction in incremenets of 1/100th of a degree
y_hat = [0 1 0]; %unit vector in cartisean coordinates
x_h = [1 0 0];
y_h = [0 1 0];
z_h = [0 0 1];
d_rho_f_test = lambda/10;
d_phi_f_test = (pi/180)/5;

d_theta = d_phi_f/4;

d_a_test = sqrt(sqrt(4*F^2+(D/2)^2)/(2*F)*...
    (D/2)*d_phi_f_test*d_rho_f_test)/lambda

d_a = sqrt(sqrt(4*F^2+((D-1)/2)^2)/(2*F)*...
    ((D-1)/2)*d_phi_f*d_rho_f)/lambda
%% Analysis: Determining the Cn values

theta = 1.75*pi/180; %rad, angle from reflector axis (referencing global coordinate system)

phi_f = 0;     %rad, starting phi_f angle

E_field = 0;

%the following two while loops are the numerical integrator to solve
%for the scattered E field at a point
while(phi_f<=2*pi) %rad, iterating from phi_f = 0 to phi = 2*pi

    rho = 0; %m, projection of the radial distance onto the aperture plane

    while(rho<=(D-1)/2) %integrating from the vertix to the rim of the reflector

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

        H_I  =  1i*(1/s_i)*(cross(y_h, s_i_h)/norm(cross(y_h, s_i_h)))...
            *exp(-1i*beta*s_i)*(cos(theta_f)^q); %A/m, magnetic field normalized without source current Equation 2 in AWPL
       
        J = 2*cross(n_h, H_I); %A/m, equivalent surface current density accross surface of parabolic

        phi = -phi_f;%rad, phi_f is -phi in the global coordinate system
        r_hat = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)]; %r_hat to cartisean in the global coordinate system
        phase = beta*dot(r_hat, s_i_v);%m, phase incurred by far-field approximation dot product in Equation 1 in AWPL

        integrand = -1i*J*omega*mu*...
            exp(1i*phase)*...
            (1/(4*pi))*d_a; %V/m, integrand in calculating the scattered e field

        E_field = E_field + integrand; %V/m, adding the infantesimal sum of the scattered e field
        rho = rho + d_rho_f; %m, incrementing by the differential
    end
    phi_f = abs(phi_f); %rad, taking the absolute value of the angle for looping
    phi_f = phi_f +  d_phi_f; %incrementating the observation angle (theta) by a 10th of a degree
end

start_rho = rho;

phi_x = 0*pi/180;

r_f_h = [sin(theta)*cos(phi_x) sin(theta)*sin(phi_x) cos(theta)]; %r_hat to cartisean
e = cross(cross(y_h, r_f_h), r_f_h)...
    /norm(cross(y_h, r_f_h));

C = zeros(1, length(0:d_phi_f_test:2*pi)*length(start_rho:d_rho_f_test:D/2));  %V/m, initializing scattered E field for this computiation
counter = 1;
phi_f = 0;     %rad, starting phi_f angle

%the following two while loops are the numerical integrator to solve
%for the scattered E field at a point
while(phi_f<=2*pi) %rad, iterating from phi_f = 0 to phi = 2*pi

    rho = start_rho; %m, projection of the radial distance onto the aperture plane

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
            rho*d_phi_f_test*d_rho_f_test; %m, differential surface are of the parabolic 16-165 in Stutzman

        H_I  =  1i*(1/s_i)*(cross(y_h, s_i_h)/norm(cross(y_h, s_i_h)))...
            *exp(-1i*beta*s_i)*(cos(theta_f)^q); %A/m, magnetic field normalized without source current Equation 2 in AWPL

        J = 2*cross(n_h, H_I); %A/m, equivalent surface current density accross surface of parabolic

        phi = -phi_f;%rad, phi_f is -phi in the global coordinate system
        r_hat = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)]; %r_hat to cartisean in the global coordinate system
        phase = beta*dot(r_hat, s_i_v);%m, phase incurred by far-field approximation dot product in Equation 1 in AWPL

        E_field_test = -1i*J*omega*mu*...
            exp(1i*phase)*...
            (1/(4*pi))*d_a; %V/m, integrand in calculating the scattered e field

        if(abs(dot(E_field + (-1)*E_field_test,e)) <= abs(dot(E_field+ 1*E_field_test,e)))
            C(counter) = -1;
        else
            C(counter) = 1;
        end
        E_field = E_field+C(counter)*E_field_test;
        counter = counter + 1;
        rho = rho + d_rho_f_test; %m, incrementing by the differential
    end
    phi_f = abs(phi_f); %rad, taking the absolute value of the angle for looping
    phi_f = phi_f +  d_phi_f_test; %incrementating the observation angle (theta) by a 10th of a degree
end


%% Analysis: Calculating scattered E field from parabolic

E_s_mag = zeros(1, length(0*pi/180:d_theta:field_angle)); %initializing the scattered E field as a function of theta
g = 1; %counter for array
theta = 0*pi/180; %rad, angle from reflector axis (referencing global coordinate system)
while(theta <= field_angle)    

    W = [num2str(100*theta/field_angle), '% Complete'];

    disp(W)

    counter = 1;
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
            
            H_I  =  j*(1/s_i)*(cross(y_h, s_i_h)/norm(cross(y_h, s_i_h)))...
                *exp(-j*beta*s_i)*(cos(theta_f)^q); %A/m, magnetic field normalized without source current Equation 2 in AWPL    

            if(rho<((D-1)/2))
                J = 2*cross(n_h, H_I); %A/m, equivalent surface current density accross surface of parabolic
                d_a = sqrt(4*F^2+rho^2)/(2*F)*...
                    rho*d_phi_f*d_rho_f; %m, differential surface are of the parabolic 16-165 in Stutzman
                rho = rho + d_rho_f; %m, incrementing by the differential
            else
                J = C(counter)*2*cross(n_h, H_I); %A/m, equivalent surface current density accross surface of parabolic
                counter = counter + 1;
                d_a = sqrt(4*F^2+rho^2)/(2*F)*...
                    rho*d_phi_f*d_rho_f_test; %m, differential surface are of the parabolic 16-165 in Stutzman
                rho = rho + d_rho_f_test; %m, incrementing by the differential
            end
            phi = -phi_f;%rad, phi_f is -phi in the global coordinate system
            r_hat = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)]; %r_hat to cartisean in the global coordinate system          
            phase = beta*dot(r_hat, s_i_v);%m, phase incurred by far-field approximation dot product in Equation 1 in AWPL   
            
            integrand = -1i*J*omega*mu*...
                exp(1i*phase)*...
                (1/(4*pi))*d_a; %V/m, integrand in calculating the scattered e field       
            
            E_field = E_field + integrand; %V/m, adding the infantesimal sum of the scattered e field
            %rho = rho + d_rho_f; %m, incrementing by the differential            
        end
        phi_f = abs(phi_f); %rad, taking the absolute value of the angle for looping

        if(rho<((D-1)/2))
            phi_f = phi_f +  d_phi_f; %incrementating the observation angle (theta) by a 10th of a degree
        else
            phi_f = phi_f +  d_phi_f_test; %incrementating the observation angle (theta) by a 10th of a degree
        end
    end

    r_f_h = [sin(theta)*cos(phi_x) sin(theta)*sin(phi_x) cos(theta)]; %r_hat to cartisean
    e = cross(cross(y_h, r_f_h), r_f_h)...
        /norm(cross(y_h, r_f_h));

    E_s_mag(g) = norm(dot(E_field, e)); %V/m, magnitude of the E field
    theta = theta + d_theta; %rad, incrementating the observation angle (theta) by a 100th of a degree
    g = g + 1;    %incremeneting the array index
end

E_s_squarred = abs(E_s_mag).^2; %magntiude of the scattered e field squarred



%% Analysis: Calculating beam solid angle  
%Used in the gain pattern for recreating figure 3
theta = 0; %rad, initilazing theta for calculating beam solid angle
integrand = 0; %sr. initilizing the sum (for the integral) for calculating beam solid angle
while(theta <= pi/2)    %calculating the beam solid angle
    F = (cos(theta)^(2*q))*sin(theta)*d_phi_f/10; %normalized pattern function of E field
    integrand = integrand + F;  %adding the infantesimal sum of the scattered e field 
    theta = theta + d_phi_f/10; %incrementating the observation angle (theta) by a 100th of a degree    
end
solid_angle = 2*pi*integrand; %sr, beam solid angle
E_field_max = eta; %V/m, max of the feed E field
norm_E = (solid_angle/(4*pi))*(E_field_max)^2; %magnitude of E (divided by the source current) from feed squarred multiplied by the beam solid angle and divided by 4*pi

%% Plotting
angle_rim = (0*pi/180:d_theta:field_angle); %rad, used in plotting the horizontal axis of the gain plot
D_dBi = 10*log10(E_s_squarred/norm_E); %power pattern of antenna
figure;plot(angle_rim*180/pi,D_dBi, 'LineWidth', 2);hold all;xlabel('\theta [deg]');xlim([0 field_angle*180/pi]); %2 degree range on horizontal axis
ylabel('Co-polairzed Pattern [dBi]');title('Figure 5: Co-Pol Electric Field Magnitude');grid on;set(gca,'FontSize',20);
%ylim([-35 50]);

% x = [0 1 1.75 1.3 0.8 1.1 1.875 0.5 1.5 0.6 1.05 1.7 0.35];
% y = [48.2 -17 -35 23 32 20 2 42 18 40 10 0 45];
% plot(x, y, 'o', 'MarkerFaceColor', 'r')
% 
% legend('Our Code Results', 'AWPL Results')