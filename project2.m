%%
    % APPM 2360 Project 2 Mariana Trench 
%% 
    % Authors: Noah Abeson, Cameron Hoffman, Rahul Sampangiramiah 
%% 
    % Last Updated: March 23rd, 2020 @ 1:30 AM 
%% Housekeeping    

clear; 
clc; 

%% Importing .csv files 
mariana_depth = importdata('mariana_depth.csv')';
mariana_latitude = importdata('mariana_latitude.csv');
mariana_longitude = importdata('mariana_longitude.csv')'; 

%% Manipulating Data for Incomplete Singular Value Decomposition 
[longitude, latitude] = meshgrid(mariana_longitude, mariana_latitude);
depth = griddata(mariana_longitude, mariana_latitude, mariana_depth,longitude,latitude);

%% Contour Plot of Mariana Trench Depths 
figure(1)
contourf(longitude, latitude, depth) 
title('Mariana Trench Depths at Different Locations')
xlim([138.9 150])
xlabel('Longitude') 
ylabel('Latitude') 
colorbar

%% Function Calls
[Deepest_Point_Sampled, Deepest_Latitude, Deepest_Longitude] = calc_maxdepth(mariana_depth, mariana_latitude, mariana_longitude);
[Mean_Depth_In_Km] = calc_meandepth(mariana_depth) 

%% Eigenvalue Calculations 
unit_vec = randi(10,1440,1);
unit_vec = unit_vec/norm(unit_vec);


AT_A = mariana_depth * mariana_depth'; 


for i = 2:10

unit_vec(:,i) = AT_A * unit_vec(:,i-1) / norm(AT_A * unit_vec(:,i-1));

    if norm(unit_vec(:,i) - unit_vec(:,i-1)) < 0.00001
    
        break; 
    
    end
 
 V1 = unit_vec(:,end); 

end 

N = [1:size(V1)]; 
figure(2)
plot(N,V1)
title('Eigenvector vs. N Location')
xlabel('Location') 
ylabel('Eigenvector Value') 

%% Eigenvector Calculations 
% V50 = zeros(1,50); 
% 
% for i = 2:50 
%     u_1 = randi(10,1440,1);
%     u_1 = u_1 / norm(u_1); 
%     
%     u_1(:,i) = AT_A * u_1(:,i-1); 
%     
%     for j = 1:i-1 
%     u_1(:,i) =  u_1(:,i-1) - symsum(u_1(:,i)' * V50(j)*V50(j),j,1,i-1);
%     end 
%     
%     u_1(:,i) = u_1(:,i)/norm(u_1(:,i))
%     
%     if norm(u_1(:,i) - u_1(:,i-1)) < 10^-3
%         
%         u_1(:,i) = V50(:,i) 
%         break; 
%     end 
%     
% end 

%% Deepest Point Function 
function [Deepest_Point_Sampled, Deepest_Latitude, Deepest_Longitude] = calc_maxdepth(mariana_depth, mariana_latitude, mariana_longitude)
    deepest_point_vec = max(abs(mariana_depth));
    deepest_point_vec = deepest_point_vec';
    [deepest_point,idx_vec] = max(abs(deepest_point_vec));
    
     %% Index for latitude is a little different so this is just quick indexing output 
     for i = 1:1320 
        for j = 1:1440
            if abs(mariana_depth(j,i)) == deepest_point 
                idx_lat = j;
                idx_long = i; 
            end 
        end 
     end 
    
    deepest_lat = mariana_latitude(idx_lat);
    deepest_long = mariana_longitude(idx_long); 
    
    Deepest_Point_Sampled = deepest_point
    Deepest_Latitude = deepest_lat 
    Deepest_Longitude  = deepest_long
end 

%% Mean Depth Function 
function [Mean_Depth_In_Km] = calc_meandepth(mariana_depth) 
r = 1;
    for i = 1:1320 
        for j = 1:1440
            if abs(mariana_depth(j,i)) >= 6000
                nominal_depth_vec(r) = abs(mariana_depth(j,i));
                r = r + 1;
             end 
        end 
    end 
    
     Mean_Depth_In_Km = mean(nominal_depth_vec)/1000; % In km

end 



















