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
epsilon = 8.85e-12; %F/m
mu = (4*pi)*(10^-7); %H/m
eta = sqrt(mu/epsilon); %ohms
omega = 2*pi*f; %rad/m
%rim_angle = 2*atan(1/(4*F/D)); %rad
rim_angle = 2*pi/180;
rho_step_size = lambda/10;
angular_step_size = (pi/180)/10;
y_hat = [0 1 0]; %unit vector in cartisean coordinates

%% Analysis: Numerator
theta_f = 0;
E_field_Magnitude = zeros(1, length(0:angular_step_size/10:rim_angle));
k = 1;
while(theta_f <= rim_angle)    
    
    E_field = [0 0 0];  %E field array
    phi_f = 0;     %starting phi angle
    %s_i =  F*sec(theta_f/2)^2; %Eqn 9-183 in Stuzmann    
    
    while(phi_f<=2*pi) %iterating from p = 0 to p = D/2                
        rho = 0; %m, projection of the radial distance onto the aperture plane       
        while(rho<=D/2) 
            z_f = F-(rho^2)/(4*F); %Eqn 9-182 rearranged for z_f
            s_i_vec = [rho*cos(phi_f/2) rho*sin(phi_f/2) z_f]; %point on the reflector 
            s_i = norm(s_i_vec);   
            s_i = F*sec(theta_f/2)^2; %Eqn 9-183 in Stuzmann  
            
            %Unit Vectors
            s_i_hat = s_i_vec/norm(s_i_vec);            
            r_f_hat = [sin(theta_f)*cos(phi_f) sin(theta_f)*sin(phi_f) cos(theta_f)]; %r_hat to cartisean            
            theta_f_hat = [cos(theta_f)*cos(phi_f) cos(theta_f)*sin(phi_f) -sin(theta_f)]; %theta_hat to cartisean              
            n_hat = -cos(theta_f/2)*r_f_hat + sin(theta_f/2)*theta_f_hat; %equation 9-188 in Stuzman
            
            differential_surface_area = rho*sqrt(4*F^2+rho^2)/(2*F)*...
                angular_step_size*rho_step_size; %m, 16-165 in Stutzman      
            
            H_over_I  =  j*(1/s_i)*(cross(y_hat, s_i_hat)/norm(cross(y_hat, s_i_hat)))...
                *exp(-j*beta*s_i)*(cos(theta_f)^q); %A/m, Equation 2 in AWPL    
            
            J = 2*cross(n_hat, H_over_I); %A/m, equivalent surface current density             
            phase = beta*dot(r_f_hat, s_i_vec);%m, dot product in Equation 1 in AWPL   
            
            integrand = -j*J*omega*mu*...
                exp(j*phase)*...
                (1/(4*pi))*differential_surface_area;        
            
            E_field = E_field + integrand;            
            rho = rho + rho_step_size;            
        end        
        phi_f = phi_f + angular_step_size;
    end
    E_field_Magnitude(k) = norm(E_field);
    theta_f = theta_f + angular_step_size/10;    
    k = k + 1;    
end

numerator = abs(E_field_Magnitude).^2;

%% Analysis: Denominator
%Used in the gain pattern for recreating figure 3
theta_f = 0;
integrand = 0;
while(theta_f <= pi/2)    
    integrand = integrand + (cos(theta_f)^(2*q))*sin(theta_f)*angular_step_size/10;    
    theta_f = theta_f + angular_step_size/10;    
end
solid_angle = 2*pi*integrand;
E_field_max = eta/F;
denominator = (solid_angle/(4*pi))*(abs(E_field_max))^2;

%% Plotting
angle_rim = (0:angular_step_size/10:rim_angle);
D_dBi = 10*log10(numerator/denominator);
figure;plot(angle_rim*180/pi,D_dBi);hold all;xlabel('\theta [deg]');xlim([0 2]);
ylabel('Directivity [dBi]');title('Figure 3');grid on;
