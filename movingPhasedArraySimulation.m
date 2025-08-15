%% Visualization of the Texas Instruments TIDA-01454 Circular Microphone Board (CMB)

% Define ideal omnidirectional SiSonic SPH1642HT5H-1 mic element
micElement = phased.OmnidirectionalMicrophoneElement('FrequencyRange',[50 10000]); % 50 Hz – 10 kHz per datasheet 
N = 8;   % Number of mic elements
R = 0.024;    % Radius (m) of CMB
c = 343.21; % Speed (m/s) of Sound in Air

% Define UCA with omnidirectional microphones
uca_mics = phased.UCA( ...
    'NumElements', N, ...
    'Radius', R, ...
    'Element', micElement);


viewArray(uca_mics,'ShowNormals',true, ...
    'Title','TIDA-01454 PCM1864-Based Circular Microphone Board')
view(0,90)

%% 2-D Polar Pattern of the Array @ 50 Hz/10 kHz
freq = [50 10000];   % Frequency Response 
pattern(uca_mics,freq,[-180:180],0,'CoordinateSystem','polar', ...
    'Type','powerdb', ...
    'PropagationSpeed',c, ...
    'Normalize',true);

%% 3-D Polar Pattern of the Array @ 50 Hz
fc1 = 50;
pattern(uca_mics,fc1,[-90:0.1:90],[-45:0.1:45],...
    'CoordinateSystem','polar',...
    'Type','efield',...
    'PropagationSpeed',c)
title('3-D Response Pattern of CMB Array @ 50 Hz')

%% 3-D Polar Pattern of the Array @ 10 kHz
fc2 = 10000;
pattern(uca_mics,fc2,[-90:0.1:90],[-45:0.1:45],...
    'CoordinateSystem','polar',...
    'Type','efield', ...
    'PropagationSpeed',c)
title('3-D Response Pattern of CMB Array @ 10 kHz')

%% 8-Channel Microphone Capture with Moving Audio Source
% This section simulates an 8-channel audio capture with a moving speech source.
% The female speech source moves in a circular path around the microphone array,
% simulating a speaker walking around a venue. The movement is implemented by
% dynamically updating the source azimuth angle during simulation.

% Fixed angles for non-moving sources (azimuth; elevation in degrees)
ang_malespeech = [-10; 10]; 
ang_laughter = [20; 0]; 

% Moving source parameters for female speech
initial_angle = -30;   % Starting azimuth in degrees
rotation_speed = 60;   % Degrees per second (full circle in 6 seconds)
radius = 2;            % Distance from array center in meters

% Frequency Resolution (Hz) = fs / frequency_bins  
fs = 48000;
frequency_bins = 6000;

t_duration = 6;  % 6 seconds
t = 0:1/fs:t_duration-1/fs;

% Precompute the moving source angles for each frame
numFrames = floor(t_duration * fs / NSampPerFrame);
t_frame = (0:numFrames-1) * (NSampPerFrame / fs);  % Time at start of each frame
angles_femalespeech = initial_angle + rotation_speed * t_frame;
% Wrap angles to [-180, 180]
angles_femalespeech = mod(angles_femalespeech + 180, 360) - 180;

% Wideband Collector Simulates Channel Phase Shifts/Time Delay at each mic
collector = phased.WidebandCollector('Sensor',uca_mics, ...
    'PropagationSpeed',c, ...
    'SampleRate',fs, ...
    'NumSubbands',frequency_bins, ...
    'ModulatedInput', false);

%% Simulate 8-Channel 24-bit Capture with Moving Source

% Generate White Noise Signal to represent Thermal Noise 
prevS = rng(2008);
noisePwr = 1e-4;

% preallocate
NSampPerFrame = 1000;
NTSample = t_duration*fs;
sigArray = zeros(NTSample,N);
voice_dft = zeros(NTSample,1);
voice_cleanspeech = zeros(NTSample,1);
voice_laugh = zeros(NTSample,1);

% set up audio device writer
player = audioDeviceWriter('SampleRate',fs);

% Female Voice of Interest (moving source)
femaleSpeechFileReader = dsp.AudioFileReader('FemaleSpeech-24bit-48kHz-6secs.wav', ...
    'SamplesPerFrame',NSampPerFrame);
% Male Voice of Interest followed by Music Interference
maleSpeechFileReader = dsp.AudioFileReader('MaleSpeech+Music-24bit-48kHz-20secs.wav', ...
    'SamplesPerFrame',NSampPerFrame);
% Laughter Interference
laughterFileReader = dsp.AudioFileReader('Laughter-24bit-48kHz-7secs.wav', ...
    'SamplesPerFrame',NSampPerFrame);

% Begin Simulation
for m = 1:NSampPerFrame:NTSample
    % Compute current frame index
    frame_index = (m-1)/NSampPerFrame + 1;
    current_ang_femalespeech = [angles_femalespeech(frame_index); 0];
    
    sig_idx = m:m+NSampPerFrame-1;
    x1 = 2*femaleSpeechFileReader();
    x2 = laughterFileReader();
    x3 = maleSpeechFileReader();
    temp = collector([x1 x2 x3], ...
        [current_ang_femalespeech, ang_laughter, ang_malespeech]) + ...
        sqrt(noisePwr)*randn(NSampPerFrame,N);
    player(0.5*temp(:,3));
    sigArray(sig_idx,:) = temp;
    voice_dft(sig_idx) = x1;
    voice_cleanspeech(sig_idx) = x2;
    voice_laugh(sig_idx) = x3;
end

%% Plot Signal
plot(t,sigArray(:,3));
xlabel('Time (sec)'); ylabel ('Amplitude (V)');
title('Signal Received at Channel 3 (Moving Source)'); ylim([-3 3]);

%% Export 8-channel audio simulation to WAV file
% This section extracts the 8 audio channels stored in 'sigArray' variable
% and writes them into a single multi-channel WAV file named 'output_moving_source.wav'.

% Check if sigArray has channels as columns or rows
if size(sigArray, 1) < size(sigArray, 2)
    % If channels are rows, transpose to have channels as columns
    sigArray = sigArray';
end 

% Normalize audio data to [-1, 1] range for WAV compatibility
% This ensures the full dynamic range is used without clipping
maxVal = max(abs(sigArray(:)));
if maxVal > 1
    sigArray = sigArray / maxVal;
end 

% Write the multi-channel audio to a WAV file
filename = 'output_moving_source.wav';
audiowrite(filename, sigArray, fs); 

disp(['8-channel audio with moving source saved to: ' filename]); 