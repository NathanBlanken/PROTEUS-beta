function SimulationParameters = reset_simulation_parameters()
% Simulation parameters depend on the transmit frequency f0 (Hz).

CFL = 0.3; % Courant-FriedrichsLewy number
ppwl = 8;  % Points per wavelength

SimulationParameters.IndependentVariable    = 'CFL';
SimulationParameters.CFL                    = CFL;
SimulationParameters.PointsPerWavelength    = ppwl;

SimulationParameters.NumberOfInteractions = 0;
SimulationParameters.HybridSimulation     = true;
SimulationParameters.SensorOnGrid         = false;

SimulationParameters.Solver       = '3DC';
SimulationParameters.DeviceNumber = 0;

end