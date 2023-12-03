function Acquisition = reset_acquisition()

Acquisition.FrameRate = 500;            % (Hz)
Acquisition.NumberOfFrames = 10;    
Acquisition.PulsingScheme = 'Standard';
Acquisition.TimeBetweenPulses = 100e-6; % (s);
Acquisition.NumberOfPulses = 1;

% Use existing ground truth microbubble positions:
Acquisition.Precompute = false;
Acquisition.Folder = '';

end