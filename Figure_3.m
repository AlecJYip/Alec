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
rim_angle = 2*atan(1/(4*F/D)); %rad
rim_angle = 2*pi/180;
rho_step_size = lambda/10;
angular_step_size = (pi/180)/10;
feed_point = [0 0 F]; %location of the feed in cartisean coordinates
y_hat = [0 1 0]; %unit vector in cartisean coordinates

%% Analysis: Numerator
theta = 0;
E_field_Magnitude = zeros(1, length(0:angular_step_size/10:rim_angle));
k = 1;
while(theta <= rim_angle)    
    
    phi = 0;
    
    E_field = [0 0 0];  %E field array 
    
    rho = 0; %m, projection onto the Z axis   
    
    while(rho<=D/2) %iterating from p = 0 to p = D/2       
        
        phi = 0;       
              
        s_i = F*(sec(theta/2)^2); %m, distance from feed to the reflector        
        
        z = s_i*cos(theta); %m,  %SE simpler/faster to compute z directly from rho. Simpler math, easier to follow.
        
        while(phi<=2*pi)
            
            r_p = [s_i*sin(phi/2) s_i*sin(theta) F-z]; %point on the reflector
            
            %r_p = [ rho*cos(phi) rho*sin(phi) z ];            
            
            s_i_vec = r_p - feed_point; %vector s^i in the paper    
            
            s_i_hat = s_i_vec/norm(s_i_vec);            
            
            differential_surface_area = rho*sqrt(4*F^2+rho^2)/(2*F)*...
                angular_step_size*rho_step_size; %m, 16-165 in Stutzman 
            
            H_over_I  =  1i*(1/s_i)*(cross(y_hat, s_i_hat)/norm(cross(y_hat, s_i_hat)))...
                *exp(-1i*beta*s_i)*(cos(theta)^q); %A/m, Equation 2 in AWPL               
            
            r_f_hat = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)]; %r_hat to cartisean           
            
            theta_hat = [cos(theta)*cos(phi) cos(theta)*sin(phi) -sin(theta)]; %theta_hat to cartisean            
            
            n_hat = -cos(theta/2)*r_f_hat + sin(theta/2)*theta_hat; %equation 9-188 in Stuzman           
            
            J = 2*cross(n_hat, H_over_I); %A/m, equivalent surface current density            
                        
            field_dist = (dot(r_f_hat, s_i_vec));%m, dot product in Equation 1 in AWPL
            %SE This is actually r_hat dot r' where r' is the vector from the origin to the point on the reflector.
            
            integrand = J*omega*mu*(1/(4*pi))*...
                exp(1i*beta*field_dist)*...
                differential_surface_area;
            
            E_field = E_field + integrand;
            
            phi = phi + angular_step_size;
            
        end
        
        rho = rho + rho_step_size;
        
    end
    
    E_field_Magnitude(k) = norm(E_field);
    
    theta = theta + angular_step_size/10;
    
    k = k + 1;
    
end

numerator = E_field_Magnitude.^2;

%% Analysis: Denominator

%Used in the gain pattern for recreating figure 3

theta = 0;

arr2 = 0;

while(theta <= pi/2)
    
    arr2 = arr2 + (cos(theta)^(2*q))*sin(theta)*angular_step_size;
    
    theta = theta + angular_step_size;
    
end

denominator_term_1 = 2*pi*arr2;

denominator_term_2 = (eta^2)/((4*pi)^3);

denominator = denominator_term_1*denominator_term_2;

%% Plotting

angle_rim = (0:angular_step_size/10:rim_angle);

D_dBi = 10*log10(numerator/denominator);

figure;plot(angle_rim*180/pi,D_dBi);hold all;xlabel('\theta [deg]');xlim([0 2]);

ylabel('Directivity [dBi]');title('Figure 3');grid on;


%% Functions
