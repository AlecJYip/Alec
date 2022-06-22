%% Title Block
%Purpose of code: To determine the E field norm from a differential
%patch on the reflector.

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
omega = 2*pi*f; %rad/m
field_angle = 2*pi/180;
d_rho = lambda/10;
d_phi = (pi/180)/10;
y_hat = [0 1 0]; %unit vector in cartisean coordinates

%% Analysis: Numerator
E_field_Magnitude = zeros(1, length(0:d_phi/10:field_angle));
g = 1;
theta = 0;
while(theta <= field_angle)    
    
    E_field = [0 0 0];  %E field array
    
    phi_f = 0;     %starting phi angle
    
    while(phi_f<=2*pi) %iterating from p = 0 to p = D/2                
        
        rho = 0; %m, projection of the radial distance onto the aperture plane       
        
        while(rho<=D/2) 
            
            theta_f = -2*atan(rho/(2*F));
            phi_f = -phi_f;
            s_i = F*sec(theta_f/2)^2; %Eqn 9-183 in Stuzmann
            z_f = -s_i*cos(theta_f); %m
            s_i_vec = [rho*cos(phi_f) rho*sin(phi_f) z_f]; %point on the reflector 
           
            %Unit Vectors
            s_i_hat = s_i_vec/norm(s_i_vec);            
            r_f_hat = [sin(theta_f)*cos(phi_f) sin(theta_f)*sin(phi_f) cos(theta_f)]; %r_hat to cartisean            
            theta_f_hat = [cos(theta_f)*cos(phi_f) cos(theta_f)*sin(phi_f) -sin(theta_f)]; %theta_hat to cartisean              
            n_hat = -cos(theta_f/2)*r_f_hat + sin(theta_f/2)*theta_f_hat; %equation 9-188 in Stuzman
            
            differential_surface_area = sqrt(4*F^2+rho^2)/(2*F)*...
                rho*d_phi*d_rho; %m, 16-165 in Stutzman      
            
            H_over_I  =  j*(1/s_i)*(cross(y_hat, s_i_hat)/norm(cross(y_hat, s_i_hat)))...
                *exp(-j*beta*s_i)*(cos(theta_f)^q); %A/m, Equation 2 in AWPL    
            
            J = 2*cross(n_hat, H_over_I); %A/m, equivalent surface current density  
            
            phi = -phi_f;
            
            r_hat = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)]; %r_hat to cartisean            
            phase = beta*dot(r_hat, s_i_vec);%m, dot product in Equation 1 in AWPL   
            
            integrand = -j*J*omega*mu*...theta_f
                exp(j*phase)*...
                (1/(4*pi))*differential_surface_area;        
            
            E_field = E_field + integrand;
            rho = rho + d_rho;            
        end
        phi_f = phi_f*-1;
        phi_f = phi_f +  d_phi;
    end
    E_field_Magnitude(g) = norm(E_field);
    theta = theta + d_phi/10;    
    g = g + 1;    
end

numerator = abs(E_field_Magnitude).^2;

%% Analysis: Denominator
%Used in the gain pattern for recreating figure 3
theta = 0;
integrand = 0;
while(theta <= pi/2)    
    integrand = integrand + (cos(theta)^(2*q))*sin(theta)*d_phi/10;    
    theta = theta + d_phi/10;    
end
solid_angle = 2*pi*integrand;
E_field_max = eta;
denominator = (solid_angle/(4*pi))*(E_field_max)^2;

%% Plotting
angle_rim = (0:d_phi/10:field_angle);
D_dBi = 10*log10(numerator/denominator); %-17.0809
figure;plot(angle_rim*180/pi,D_dBi);hold all;xlabel('\theta [deg]');xlim([0 2]);
ylabel('Directivity [dBi]');title('Figure 3');grid on;
