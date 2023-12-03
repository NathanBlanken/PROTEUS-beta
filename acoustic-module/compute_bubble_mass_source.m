function mass_source = compute_bubble_mass_source(...
    sensed_p,  radii, kgrid, Medium, Microbubble, Transmit)
% =========================================================================
% Compute bubble response and convert to mass source.
% input:    sensed_p: pressure sensed by the bubbles
%           radii:    microbubble radii
%           kgrid:    the k-Wave grid
% output:   mass_source
%           
% MBs vibrate under insonification from transducer, producing pressure
% waves of their own. The responses are estimated with Marmottant model.
% =========================================================================

% Number of microbubbles and signal length:
[N_MB,N] = size(sensed_p);

% Sampling rate for the microbubble module:
fs_MB = Microbubble.SamplingRate;

% Filter settings:
Filter.dt          = 1/fs_MB;
Filter.t_array     = [];
Filter.k_max       = kgrid.k_max;
Filter.TW          = 0.1/(fs_MB*kgrid.dt);
Filter.sound_speed = Medium.SpeedOfSoundMinimum;

% Microbubble driving pulse settings:
pulse.f  = Transmit.CenterFrequency;
pulse.w  = pulse.f * 2 * pi;    
pulse.fs = fs_MB;

% Time vectors for k-Wave signals and microbubble module signals:
t_kwave = (0:(N-1))*kgrid.dt;
M = floor(t_kwave(end)*fs_MB) + 1;
t_MB    = (0:(M-1)) / fs_MB;

% Preallocate the bubble pressure matrix:
mass_source = zeros(N_MB,M,class(sensed_p));

% Resample signal at the sampling rate of the microbubble module:
sensed_p = sinc_interpolation(t_kwave, transpose(sensed_p), t_MB);

% loop over all MBs to compute response 
for i = 1:N_MB
       
    % Microbubble driving pulse:  
    pulse.p = sensed_p(i,:);
    pulse.t = t_MB;
    
    % Get the properties of the liquid, gas, and shell:
    [liquid, gas] = get_microbubble_material_properties(...
        Medium,Microbubble);

    % Shell properties:
    shell = get_microbubble_shell_properties(radii(i),Medium,Microbubble);

    % Bubble properties:
    bubble.R0 = radii(i);          	% Bubble radius (m)
    bubble.dispProgress = false;    % Do not show ODE solver progress  

    % Compute the bubble response and plot the results
    disp('=====================================')
    disp(['Simulating microbubble ' num2str(i) '/' num2str(N_MB) ' ...'])
    disp('=====================================')

    [response, ~] = calcBubbleResponse(liquid, gas, shell, bubble, pulse);
    
    % Compute mass source
    MS = 4*pi*liquid.rho*response.R.^2 .* response.Rdot;
    
    % Filter out unsupported frequencies   
    MS = filterTimeSeries(Filter, Filter, MS,...
        'ZeroPhase',true,'TransitionWidth',Filter.TW,'PPW',2);

    mass_source(i,:) = transpose(MS);
    
end

% Resample signals at the sampling rate of the acoustic module:
mass_source = sinc_interpolation(t_MB, transpose(mass_source), t_kwave);

end


%==========================================================================
% FUNCTIONS
%==========================================================================


function [liquid, gas] = get_microbubble_material_properties(...
    Medium, Microbubble)
% Convert the material properties from the GUI to the format used by the
% microbubble module.

Liquid     = Medium.Vessel;
Gas        = Microbubble.Gas;

% Properties of the liquid:
liquid.k   = Liquid.ThermalConductivity;    % [W/m/K]
liquid.rho = Liquid.Density;                % [kg/m^3]
liquid.cp  = Liquid.SpecificHeat;   	    % [J/kg/K]
liquid.nu  = Liquid.DynamicViscosity;       % [Pa.s]
liquid.c   = Liquid.SpeedOfSound;           % [m/s]

% Environmental conditions:
liquid.T0  = Liquid.Temperature;            % [K]
liquid.P0  = Liquid.Pressure;               % [Pa]

% Thermodynanic model for microbubble oscillations ('Adiabatic', 
% 'Isothermal', or 'Propsperetti'):
liquid.ThermalModel = Microbubble.ThermalModel;

% Properties of the gas:
gas.k      = Gas.ThermalConductivity;     	% [W/m/K]
gas.rho    = Gas.Density;                   % [kg/m^3]
gas.Mg     = Gas.MolarMass;                 % [kg/mol]
gas.gam    = Gas.HeatCapacityRatio;
gas.cp     = Gas.SpecificHeat;              % Constant pressure [J/kg/K]

end


function shell = get_microbubble_shell_properties(R0, Medium, Microbubble)
% Shell properties of a microbubble. According to:
% Marmottant et al., J. Acoust. Soc. Am. 118 6, 2005
% OR
% Segers et al., Soft Matter, 2018, 14, 9550-9561

Shell  = Microbubble.Shell;
Liquid = Medium.Vessel;

shell.sig_0 = Shell.InitialSurfaceTension;

if strcmp(Shell.Model, 'Segers') || strcmp(Shell.Model, 'Custom')
    shell.model = 'SegersTable';
else
    shell.model = Shell.Model;
end


% MODEL CHECK AND MAXIMUM SURFACE TENSION

if R0<0.5e-6 || R0>6e-6
    warning('Microbubble dynamics uncertain for given microbubble radius.')
end

% SURFACE TENSION CURVES

switch shell.model
    case 'Marmottant'
        % MARMOTTANT MODEL    
        
        % Linearised surface tension curve (Marmottant et al., J. Acoust. 
        % Soc. Am. 118 6, December 2005)

        shell.chi   = Shell.Elasticity; % [N/m]

        % Compute the buckling radius (m):
        shell.Rb    = R0/sqrt(1+shell.sig_0/shell.chi);

        % Maximum surface tension [N/m]
        shell.sig_l = Liquid.SurfaceTension;
    
    case {'SegersTable','Custom'}
        % EXPERIMENTAL SURFACE TENSION CURVES (TABLE LOOKUP)

        A_0 = 4*pi*R0^2;         % Initial microbubble area

        % Experimental surface tension curve:
        fit.sig      = Microbubble.Shell.SurfaceTension;            
        fit.A_m_list = Microbubble.Shell.NormalizedArea;

        shell.A_m1 = min(fit.A_m_list);  % Left  domain boundary fit.
        shell.A_m2 = max(fit.A_m_list);  % Right domain boundary fit.

        % Find the normalized area for which the fit equals sig_0:
        A_m0 = interp1(fit.sig, fit.A_m_list, shell.sig_0);

        shell.A_N  = A_0/A_m0; % Reference area surface tension curve.
        
        shell.sig = griddedInterpolant(fit.A_m_list,fit.sig);

        % Surface tension of the surrounding liquid:
        shell.sig_l = max(fit.sig);
    
    otherwise
        error('Unknown shell model')
    
end

% SHELL VISCOSITY
if Microbubble.Advanced
    shell.Ks = Shell.Viscosity; % [N.s/m]    
else
    % Shell viscosity, Segers et al, Soft Matter, 14, 2018
    % Surface dilatational viscosity (N.s/m). Fit to figure 6B:
    c_1=1.5e-9; 
    c_2=8e5; 
    shell.Ks = c_1.*exp(c_2.*R0); 
end


end