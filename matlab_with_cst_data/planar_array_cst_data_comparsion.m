
cst_single_dipole_field_data % Data of a signle element

% Single Element E field components
E_theta = re_e_theta_data+ 1j* im_e_theta_data;
E_phi = re_e_phi_data + 1j*im_e_phi_data;

f = 2.4e9;
lambda = 300e6 / f;
k_0 = 2*pi/lambda;


d_theta = 1;
d_phi = 5;

theta = 0:d_theta:180;
phi = 0:d_phi:360;

[THETA ,PHI] = meshgrid(theta,phi);

THETA = deg2rad(THETA);
PHI = deg2rad(PHI); 

k_x = k_0 * sin(THETA) .* cos(PHI);  
k_y = k_0 * sin(THETA) .* sin(PHI);

% E array 

E_array_theta = zeros(length(phi),length(theta));
E_array_phi = zeros(length(phi),length(theta));

% Size of the array
n = 3; % columns
m = 3; % rows

% Space between the elements
d_x = lambda*0.5 ;
d_y = lambda*0.75 ;


% direction
theta_0 = deg2rad(0);  
phi_0 = deg2rad(0);    
 
% steering
beta_x = -k_0 * sin(theta_0) * cos(phi_0);
beta_y = -k_0 * sin(theta_0) * sin(phi_0);



% Theta
for k = 1:n
    for i = 1:m
        E_array_theta = E_array_theta + E_theta.*exp(1j * ((k_x + beta_x)* (i-1) * d_x + (k_y + beta_y) * (k-1) * d_y));
        E_array_phi = E_array_phi + E_phi.*exp(1j * ((k_x + beta_x)* (i-1) * d_x + (k_y + beta_y) * (k-1) * d_y));
    end
end


mag_square_E_array = abs(E_array_theta).^2+abs(E_array_phi).^2;




% E field and intensity of array overall Gain



U = mag_square_E_array;



d_Omega = sin(THETA)*deg2rad(d_theta)*deg2rad(d_phi);
P_rad = sum(U.*d_Omega,"all");

Gain = 4*pi*U/P_rad;


Gain_db = 10*log10(Gain+1e-10);

%Gain_db(Gain_db < -150) = -150;

figure
patternCustom(Gain_db,theta,phi)
title("Gain")

figure
patternCustom(Gain_db, theta, phi, CoordinateSystem="rectangular",Slice="phi", SliceValue=0);
title("Gain")

max_gain = 10*log10(4*pi*(m-1)*(n-1)) 

% Co pol gain calculation

co_pol_theta = E_array_theta.*sin(PHI);
co_pol_phi = E_array_phi.*cos(PHI);

mag_square_E_co = abs(co_pol_theta).^2 + abs(co_pol_phi).^2;

U_co = mag_square_E_co;


Gain_co = 4*pi*U_co/P_rad;


Gain_co_db = 10*log10(Gain_co+1e-10);

figure
patternCustom(Gain_co_db,theta,phi)
title("Gain Co")



% cross pol gain calculation

cross_pol_theta = E_array_theta.*cos(PHI); 
cross_pol_phi = -E_array_phi.*sin(PHI);

mag_square_E_cross = abs(cross_pol_theta).^2 + abs(cross_pol_phi).^2;

U_cross = mag_square_E_cross;

Gain_cross = 4*pi*U_cross/P_rad;


Gain_cross_db = 10*log10(Gain_cross+1e-10);

figure
patternCustom(Gain_cross_db,theta,phi)
title("Gain Cross")




% Cross and Co Gain on the same figure 
figure
patternCustom(Gain_cross_db, theta, phi, CoordinateSystem="rectangular",Slice="phi", SliceValue=0);
hold on
patternCustom(Gain_co_db, theta, phi, CoordinateSystem="rectangular",Slice="phi", SliceValue=0);
hold off
title("Co Gain and Cross Gain")
legend("Cross-Gain","Co-Gain")