%%Shooting method schrodinger solving
clc; clear all; close all;

h = (6.625 * 10^(-34))/(2*pi); %reduced Planks constant
m0 = 9.109 * 10^-31; %free electron mass in Kg
E_mass = 0.2 *m0; % effective electron mass for GaAs.
q= 1.602e-19; % Conversion factor
delta_E = 1e-5;

well_voltage = 1; %depth of the well in eV
well_width = 10; %well width in nm

%% Defining Mesh
mesh_spacing = 500;
N= mesh_spacing*3;
mesh = ones(1,mesh_spacing*3)*well_voltage*q; %initialising the mesh points voltage to the needed pot
for i = mesh_spacing + 1 : mesh_spacing*2 %setting the well voltage to zero
    mesh(i) = 0;
end

%% wave function definition

Eguess = 1e-5*q; %Initial starting energy guess
kmax = int32(well_voltage/delta_E); %setting the max wave Number value
energies = zeros(1,kmax);
Wave_prefactor = 2 * E_mass * ((well_width*1e-9/mesh_spacing) / h)^2;
matching_point = (N)/3;
error = ones(1,kmax);

for k = 1:kmax
    psi_forward = zeros(1,N);
    psi_backward = zeros(1,N);

    %initialising the first 2 values for wavefunctions
    psi_forward (1) = 0;
    psi_forward (2) = 1e-10;

    psi_backward (N) = 0;
    psi_backward (N-1) = 1e-10;

    % Forward integration
    for i = 2:matching_point
        psi_forward (i+1) = (2 - Wave_prefactor * (Eguess- mesh(i))) * psi_forward(i) - psi_forward(i-1);
    end

    % Reverse integration
    for i = N-1:-1:matching_point+1
        psi_backward (i-1) = (2 - Wave_prefactor * (Eguess- mesh(i))) * psi_backward(i) - psi_backward(i+1);
    end

    %match the forward and backward wave functions by gettinmg the ration
    wave_ratio = psi_forward(matching_point)/psi_backward(matching_point);
    psi_backward_normalised = psi_backward*wave_ratio;

    % after matching we need to find the error in the waves at the matching point+1
    error (k) = abs(psi_backward_normalised (matching_point+1)-psi_forward (matching_point+1))/max(abs(psi_backward_normalised));
    %merging the wave functions
    psi_forward(matching_point+1:N) = psi_backward_normalised(matching_point+1:N);
    %calculating the next energy to simulate the wavevector for
    energies(k) = Eguess/q;
    Eguess = Eguess + delta_E * q;
end

%% Reproducing the wavesvectors for the eigen energies
counter = 1;
for i=2:kmax-2
    if ((error(i) < error(i-1)) && (error(i+1) > error(i)))
        Eigen_energies(counter) = energies(i);
        counter = counter+1;
    end
end

for k = 1:length(Eigen_energies)
    psi_forward = zeros(1,N);
    psi_backward = zeros(1,N);

    ...............
    %initialising the first 2 values for wavefunctions
    psi_forward (1) = 0;
    psi_forward (2) = 1e-10;

    psi_backward (N) = 0;
    psi_backward (N-1) = 1e-10;

    Eguess = Eigen_energies(k)*q;

    % Forward integration
    for i = 2:matching_point
        psi_forward (i+1) = (2 - Wave_prefactor * (Eguess- mesh(i))) * psi_forward(i) - psi_forward(i-1);
    end

    % Reverse integration
    for i = N-1:-1:matching_point+1
        psi_backward (i-1) = (2 - Wave_prefactor * (Eguess- mesh(i))) * psi_backward(i) - psi_backward(i+1);
    end

    %match the forward and backward wave functions by gettinmg the ration
    wave_ratio = psi_forward(matching_point)/psi_backward(matching_point);
    psi_backward_normalised = psi_backward*wave_ratio;
    %merging the wave functions
    psi_forward(matching_point+1:N) = psi_backward_normalised(matching_point+1:N);
    
    Waves(k,:) = psi_forward;

end

%% plotting
figure (1)
hold on
plot(Waves(:,:)'*(q/30))
plot(mesh)
hold off
figure(2)
plot (energies,error)