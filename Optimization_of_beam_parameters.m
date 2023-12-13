clear all
clc

%%
Plate_thickness = 0.006;
%%
% This section defines the initial values for the optimization loop and the
% step size for the for loop
Start_h_angled = 0.46;
End_h_angled = 0.46;
Stepsize_h_angled = 0.01;
steps_h_angled = round((End_h_angled-Start_h_angled)/Stepsize_h_angled+1);

Start_w_verticle = 0.1;
End_w_verticle = 0.51;
Stepsize_w_verticle = 0.001;
steps_w_verticle = round((End_w_verticle-Start_w_verticle)/Stepsize_w_verticle+1);

Start_h_verticle = 0.46;
End_h_verticle = 0.46;
Stepsize_h_verticle = 0.01;
steps_h_verticle = round((End_h_verticle-Start_h_verticle)/Stepsize_h_verticle+1);

Start_h_notch = 0.005;
End_h_notch = 0.005;
Stepsize_h_notch = 0.005;
steps_h_notch = round((End_h_notch-Start_h_notch)/Stepsize_h_notch+1);

Start_L_notch = 0.06;
End_L_notch = 0.06;
Stepsize_L_notch = 0.01;
steps_L_notch = round((End_L_notch-Start_L_notch)/Stepsize_L_notch+1);

steps = steps_h_angled*steps_w_verticle*steps_h_verticle*steps_h_notch*steps_L_notch;
masses = zeros(1,steps);
ratios = zeros(1,steps);
stiffnesses = zeros(1,steps);
%ratios = zeros(1,steps);
sizes = zeros(1,steps);
w_angled_list = zeros(1,steps);
h_angled_list = zeros(1,steps);
w_verticle_list = zeros(1,steps);
h_verticle_list = zeros(1,steps);
L_notch_list = zeros(1,steps);
h_notch_list = zeros(1,steps);
phi_list = zeros(1,steps);
% 0.495086011065382
% 0.500000000000000
% 0.400000000000000
% 0.500000000000000
% 0.0150000000000000
% 0.00500000000000000
% 0.140314777443616
% w_angled = 0.495086011065382
% h_angled= 0.5;
% w_verticle= Start_w_verticle;
% h_verticle= Start_h_verticle;
% L_notch_input = Start_L_notch;
% h_notch= Start_h_notch;
% phi = 



h_angled= 0;
w_verticle= Start_w_verticle;
h_verticle= Start_h_verticle;
L_notch_input = Start_L_notch;
h_notch= Start_h_notch;

%% here starts the loop which iterates over every dimension possible
for i = 1:steps
h_angled = h_angled+Stepsize_h_angled;
if h_angled>End_h_angled
    h_angled = Start_h_angled;
    w_verticle = w_verticle+Stepsize_w_verticle;
    if w_verticle>End_w_verticle
        w_verticle = Start_w_verticle;
        h_verticle = h_verticle+Stepsize_h_verticle;
        if h_verticle>End_h_verticle
            h_verticle = Start_h_verticle;
            h_notch = h_notch+Stepsize_h_notch;
            if h_notch>End_h_notch
                h_notch = Start_h_notch;
                L_notch_input = L_notch_input+Stepsize_L_notch;
                disp(i/steps)
                disp(max(ratios(1,:)))
                if L_notch_input>End_L_notch
                    L_notch_input = Start_L_notch;
                end
            end
        end
    end
end
D =h_notch/((0.5));
L_notch = L_notch_input + D/2;


%Parameters dependent on above variables
phi = atan((0.5-0.5*h_verticle)/(2-0.5*w_verticle -2*L_notch));
w_angled= h_verticle*cos(phi);



%% C-arm params
g=9.81;
m_carm=500;
Load = m_carm*g;
dist_to_poi = 1;

%% Verticle beam params
% This section defines all kinematics of the beams
t_verticle= Plate_thickness;
E_verticle = 210e9; %Pa
A_verticle= w_verticle*h_verticle-(w_verticle-2*t_verticle)*(h_verticle-2*t_verticle); %m^2
L_verticle= 1.75 + 0.5*h_angled;   %m
Ix_verticle = (1/12)*(w_verticle)*(h_verticle^3) - (1/12)*(w_verticle-2*t_verticle)*(h_verticle-2*t_verticle)^3;
Iy_verticle = (1/12)*(h_verticle)*(w_verticle^3) - (1/12)*(h_verticle-2*t_verticle)*(w_verticle-2*t_verticle)^3;
rho_verticle = 7800;
m_verticle = A_verticle*L_verticle*rho_verticle;

%% Verticle beam stiffness
% Axial stiffness 
Cz_verticle = E_verticle*A_verticle/L_verticle; % N/m
% rotational stiffness
Kx_verticle = E_verticle*Ix_verticle/L_verticle; % Nm/rad
Ky_verticle = E_verticle*Iy_verticle/L_verticle; % Nm/rad
%% Verticle beam forces
% Moment about x axis due to C-arm load
Mx_verticle = Load*dist_to_poi;
%axial load due to poi = Load
My_verticle = Load*0.2;
Fz_verticle = Load;
%% Verticle beam potential energy U

U_verticle = (0.5*My_verticle^2 )/Ky_verticle + (0.5*Mx_verticle^2 )/Kx_verticle +(0.5*Fz_verticle^2)/Cz_verticle;

%% --------------------------------------------------------------------------------
%% Notch flexure params
d1=D;
d2= d1;
t_notch=h_angled; 

%% Ratio of notch holes to width: must be 0.01<beta<0.5

beta = h_notch/D;

%% Material constants for Notch
E_notch= 210e9;
G_notch= 81e9;
rho_notch = 7800;
m_notch = ( w_verticle*h_angled*(L_notch-d1/2) + ...
          w_angled*h_angled*(L_notch-d1/2) + ...
          ((d1+h_notch)*d1-pi*(d1/2)^2)*h_angled ) ...
          *rho_notch;

%% Stiffness calculations for notch (from JPE)

% Point A (axis of rotation) N/m
Cax_notch = 0.48*E_notch*t_notch*sqrt(h_notch/D);
Caz_notch = G_notch*h_notch*t_notch/D;
Cay_notch = 0.56*E_notch*t_notch*sqrt(h_notch/D)*(1/(1.2 +(D/h_notch)));

% Rotational stiffness [Nm/rad]
Kax_notch = (t_notch^2 /12)*Caz_notch;
Kaz_notch = 0.093*E_notch*t_notch*h_notch^2 *sqrt(h_notch/D);
Kay_notch = (t_notch^2 /12)*Cax_notch;

%Total stiffness of the notch (stiffness of end point)
Cx_notch = Cax_notch;
Cz_notch = Caz_notch; 
Cy_notch = Cay_notch; 

%% Forces acting on notch

% lateral load due to verticle beam weight
Fg_verticle = m_verticle*g +Load; 

% torque moment about x axis
Mx_notch = Load*(dist_to_poi+0.5*h_verticle);
%moment due to displaced POI
My_notch = Load*(0.2+0.5*w_verticle) +0.5*w_verticle*m_verticle*g;
%% Potential energy of notch

U_notch = 0.5*(Fg_verticle^2)/(Cz_notch) +0.5*(Mx_notch^2)/(Kax_notch) +0.5*(My_notch^2)/(Kay_notch);
%% ----------------------------------------------------------------------------------------------------
%% Angled beam params
t_angled= Plate_thickness;
E_angled = 210e9; %Pa
A_angled= w_angled*h_angled-(w_angled-2*t_angled)*(h_angled-2*t_angled); %m^2
L_angled= (2-0.5*w_verticle -2*L_notch)/cos(phi);   %m
Iy_angled = (1/12)*w_angled*(h_angled)^3  -(1/12)*(w_angled-t_angled)*(h_angled-t_angled)^3;
Ix_angled = (1/12)*h_angled*(w_angled)^3  -(1/12)*(h_angled-t_angled)*(w_angled-t_angled)^3;
rho_angled = 7800;
G_angled = 81e9;
J_angled =(w_angled*h_angled/12)*(w_angled^2 +h_angled^2) - ((w_angled-2*t_angled)*(h_angled-2*t_angled)/12)*((w_angled-2*t_angled)^2 +(h_angled-2*t_angled)^2);

%% Stiffness of angled beam
%lateral stiffness
Cz_angled = E_angled*A_angled/(L_angled^3); % N/m
% torque stiffness
Kx_angled = G_angled*J_angled/L_angled; % Nm/rad
%bending stiffness
Ky_angled = E_angled*Iy_angled/L_angled ;% Nm/rad

%% Force decomposition

Fz_angled = Load + (m_verticle + m_notch)*g;
Mx_angled = Load*(dist_to_poi+0.5*h_verticle);
My_angled = Load*(0.2+0.5*w_verticle+2*L_notch) +(0.5*w_verticle+2*L_notch)*m_verticle*g + m_notch*L_notch*g;


% New moments acting on
Mx_new = Mx_angled*cos(phi) +My_angled*sin(phi);
My_new = -Mx_angled*sin(phi) +My_angled*cos(phi);


%% Potential energy of Angled beam

U_angled = 0.5*(Fz_angled^2)/(Cz_angled) + 0.5*(My_new^2)/Ky_angled +0.5*(Mx_new^2)/Kx_angled; % add moment and torque


%% ------------------------------------------------------------------------------------------------------------------------------------------------------
%% Sum of Potential energy
U = [U_angled, U_verticle, U_notch];
U1=U;
U=sum(U);
%% Total stiffness

c_total = (0.5*Load^2)/U;

 %% mass calculations
base_plate = 0.46*0.46*0.006*rho_angled;
plates = w_angled*h_angled*Plate_thickness*rho_angled + ...
         w_verticle*h_verticle*Plate_thickness*rho_angled*2 + ...
         base_plate;
m_y_rotation = 0;
m_angled = A_angled*L_angled*rho_angled;
m_verticle = A_verticle*L_verticle*rho_verticle;
m_notch = ( w_verticle*h_angled*(L_notch-d1/2) + ...
          w_angled*h_angled*(L_notch-d1/2) + ...
          ((d1+h_notch)*d1-pi*(d1/2)^2)*h_angled ) ...
          *rho_notch;
Mass = m_notch + m_verticle + m_angled + m_y_rotation + plates;

%% mass over stiffness
% Stores the loop results into arrays
masses(i) = Mass;
Volumes(i) = m_notch/rho_notch + w_verticle*h_verticle*L_verticle +w_angled*h_angled*L_angled;
stiffnesses(i) = c_total;%[c_total, h_angled, w_verticle, h_verticle, L_notch ,h_notch];
ratios(i) = c_total/Mass;%[c_total/Mass; h_angled; w_verticle; h_verticle; L_notch ;h_notch];
sizes(i) = w_verticle;
w_angled_list(i) = w_angled;
h_angled_list(i) = h_angled;
w_verticle_list(i) = w_verticle;
h_verticle_list(i) = h_verticle;
L_notch_list(i) = L_notch;
h_notch_list(i) = h_notch;
phi_list(i) = phi;
L_verticle_list(i)= L_verticle;
L_angled_list(i) = L_angled;
end
%%
X_axis = linspace(1,steps,steps);
%% PLOTS
% close all;
% X_axis2 = X_axis(1,1:100);
% ratios2 = ratios(1,1:100);
% figure()
% hold on
% xlim([0 max(masses)])
% xlabel('Mass (Kg)')
% ylim([0 max(stiffnesses)])
% ylabel('Stiffness (Pa)')
% plot(masses,stiffnesses)
% hold off
% 
% figure()
% hold on
% title("W_{verticle} variation")
% %xlim([0 max(sizes)])
% ylabel('Stiffness/mass')
% %ylim([0 max(ratios)])
% xlabel('W_{verticle} (m)')
% plot(sizes,ratios)
% hold off
% figure()
% hold on
% title('Stiffness/mass variation')
% %xlim([0 max(sizes)])
% ylabel('Stiffness/mass')
% %ylim([0 max(ratios)])
% xlabel('Tests')
% plot(X_axis,ratios)%ratios(1,:))
% hold off
% print("finished")
% 
% figure()
% hold on
% %xlim([0 max(sizes)])
% title('Mass variation')
% ylabel('Mass')
% %ylim([0 max(ratios)])
% xlabel('Tests')
% plot(X_axis,masses)%ratios(1,:))
% hold off
% print("finished")
% 
% figure()
% hold on
% %xlim([0 max(sizes)])
% title('Stiffness variation')
% ylabel('Stiffness')
% %ylim([0 max(ratios)])
% xlabel('Tests')
% plot(X_axis,stiffnesses)%ratios(1,:))
% hold off
% print("finished")
% 
% figure()
% hold on
% %xlim([0 max(sizes)])
% title('Volume variation')
% ylabel('Volume')
% %ylim([0 max(ratios)])
% xlabel('Tests')
% plot(Volumes,ratios)%ratios(1,:))
% hold off
% print("finished")

Dimensions = [w_angled_list;h_angled_list;w_verticle_list;h_verticle_list;L_notch_list;h_notch_list;phi_list;L_verticle_list;L_angled_list];

% best is the index where the stiffness over mass (ratios) is maximum
% This is fed back into the dimension arrays to get out the optimal
% dimensions
best = find(ratios == max(ratios))
best_w_angled = w_angled_list(best)
best_h_angled = h_angled_list(best)
best_w_verticle = w_verticle_list(best)
best_h_verticle = h_verticle_list(best)
best_L_notch = L_notch_list(best)
best_h_notch = h_notch_list(best)
best_phi = phi_list(best)
best_L_verticle= L_verticle_list(best)
best_L_angled = L_angled_list(best)
stiffness_poi = stiffnesses(best)
best_mass_stiffness_ratio = ratios(best)
