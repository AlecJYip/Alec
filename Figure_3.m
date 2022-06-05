clc;
clear;
close all;

%% Variables

f = 1.5e9; %Hz

c = 3e8; %m/s

lambda = c/f; %m

beta = 2*pi/lambda; %rad/m

q = 1.14;

Si = 8.5; %m

D = 18; %m

F = 0.4*D;

e = 8.85e-12;

mu = (4*pi)*(10^-7);

eta = sqrt(mu/e);

w = 2*pi*f;

theta_rim_max = 2*pi/180;

rim_step_size = theta_rim_max/200;

power_density = zeros(1,length(0:rim_step_size:theta_rim_max));

theta_rim = 0;

k = 1;

step_size = lambda/200;

range_phi = (0:step_size:2*pi);

s_i_hat = [1 0 0];

%% Analysis: Numerator
while(theta_rim<=theta_rim_max)
    theta = 0;
    n = 1;    
    arr = [0 0 0];
    while(theta<=theta_rim)        
        phi = 0;
        s_i = F*(sec(theta/2)^2);
        p_prime = 2*F*tan(theta/2);
        far_field_dist = sqrt(s_i^2-p_prime^2);
        
        while(phi<=2*pi)
            
            y_hat = [sin(theta)*cos(phi) cos(theta)*cos(phi) -sin(theta)];
                       
            H_over_I  = (1/s_i)*(cross(y_hat, s_i_hat)/magnitude(cross(y_hat, s_i_hat)))...
                *exp(-1i*beta*s_i)*(cos(theta)^q);
            
            n_hat = [-cos(theta/2) sin(theta/2) 0];
            
            J = 2*cross(n_hat, H_over_I);
            
            func = 1i*(1/(4*pi))*w*mu*J...
                *exp(1i*beta*far_field_dist);
            
            arr = arr + (step_size^2)*(func);
                        
            phi = phi + step_size;
                        
            n = n + 1;
        end
        theta = theta + step_size;
    end
    power_density(k)= (1/(2*eta))*(magnitude(arr)^2);
    k = k + 1;
    theta_rim = theta_rim + rim_step_size;
end

%% Analysis: Denominator

theta = 0;

arr2 = 0;

while(theta<=pi)
    phi = 0;
    s_i = F*(sec(theta/2)^2);
    while(phi<=2*pi)
        
        y_hat = [sin(theta)*cos(phi) cos(theta)*cos(phi) -sin(theta)];
                
        H_over_I  = (1/s_i)*(cross(y_hat, s_i_hat)/magnitude(cross(y_hat, s_i_hat)))...
            *exp(-1i*beta*s_i)*(cos(theta)^q);
        
        E = -eta*(H_over_I);
        
        arr2 = arr2 + (step_size^2)*sin(theta)*(1/(2*eta))*magnitude(E)^2;
        
        phi = phi + step_size;
        
        n = n + 1;
    end
    theta = theta + step_size;
end

avg_power_density =  magnitude(arr2);


%% Plotting
angle_rim = (0:rim_step_size:theta_rim_max);
D_dBi = 10*log10(power_density/avg_power_density);
figure;plot(angle_rim*180/pi,D_dBi);hold all;xlabel('\theta [deg]');
ylabel('Directivity [dBi]');title('Figure 3');grid on;

%% Functions
function J = magnitude(Q) 
    n = 1;
    sum_1 = zeros(1, length(Q));
    while(n<=length(Q))        
        sum_1(n) = abs(Q(n))^2;
        n = n + 1;        
    end    
    J = sqrt(sum(sum_1));
end
