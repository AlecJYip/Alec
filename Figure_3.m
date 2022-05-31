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

theta_rim_max = pi/2;

%theta_rim_max = 2*pi/180;

P = zeros(1,length(0:theta_rim_max/200:theta_rim_max));

theta_rim = 0;

k = 1;

step_size = lambda/200;

range_phi = (0:step_size:2*pi);

%% Analysis
while(theta_rim<=theta_rim_max)
    theta = 0;
    n = 1;    
    range_theta = (0:step_size:theta_rim);    
    arr = 0;
    while(theta<=theta_rim)        
        phi =0;
        p = 2*F*tan(theta/2);
        rf = F*(sec(theta/2)^2);
        zf = sqrt(p^2-rf^2);
        while(phi<=2*pi)
            term1 = (cos(theta)^q)/sqrt((cos(theta)^2)*(sin(phi)^2)+(cos(phi)^2));
            term2 = sqrt(abs(...
                (cos(theta)*sin(phi)*cos(theta/2))^2+...
                (cos(phi)*cos(theta/2))^2+...
                (cos(theta)*sin(phi)*sin(theta/2))^2));
            term3 = exp(1i*beta*p*sin(theta)*cos(phi));
            %term3 = exp(-1i*beta*zf);
            func = term1*term2*term3;
            arr = arr + (step_size^2)*func;
            phi = phi + step_size;
            n = n + 1;
        end
        theta = theta + step_size;
    end
    P(k)= arr;
    k = k + 1;
    theta_rim = theta_rim + theta_rim_max/200;
end

%% Plotting
angle_rim = (0:theta_rim_max/200:theta_rim_max);
D_dBi = 10*log10(4*pi./(4*pi*(abs(P).^2)));
figure;plot(angle_rim*180/pi,D_dBi);hold all;xlabel('\theta');ylabel('Directivity [dBi]');title('Figure 3');grid on;
