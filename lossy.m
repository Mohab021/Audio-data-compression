clc,clear
close all
%==========================================================================
% Sampling Section - Reading the audio file and sampling
%==========================================================================

% Read the audio file (wav) and returns the sampled data into (sa)
% with a sample rate Fs for that data, Fs = 16KHz.
[sa , Fs] = audioread('Audio_file_test.wav');

% Rotate the array to be suitable for DCT algorithm.
Sa = sa'; %(( Audio signal before compression ))

%==========================================================================
% DCT Section - DCT of the sampled signal
%==========================================================================

% Apply discrete cosine transform to the sampled audio signal.
Sa_dct = dct(Sa);

%==========================================================================
% Storing Section - Storing important coefficients
%==========================================================================

% Sort Sa_dct array in descending order to get indices of sorted values
% from higher to lower.
[sorted_Sa , ind] = sort(abs(Sa_dct) , 'descend');

% Calculate the signal energy.
signal_energy = norm(Sa_dct);

% For loop to check every coefficient of DCT.
imp_cof = 0;
for k = 1 : 1 : length(Sa_dct)
    
    % Add one more coefficient
    % and calculate energy of the range [1 to k].
    cof_energy = norm( Sa_dct( ind(1 : k) ) );
    
    % Check if the selected DCT coefficients
    % still represent 99.5% of the total signal energy.
    if cof_energy < (0.995 * signal_energy)
        
        % Increase the counter by one as
        % the current coefficient is important.
        imp_cof = imp_cof + 1;
    end
end

% No we have imp_cof --> number of important coefficients
% that represent 99.5% of the total signal energy.
% Set the unimportant coefficients to zero (the rest coefficients).
Sa_dct( ind( imp_cof+1 : end ) ) = 0;

%==========================================================================
% IDCT Section - IDCT of the signal
%==========================================================================

% Apply inverse discrete cosine transform to get back
% the sampled audio signal after compression.
Sa_comp = idct(Sa_dct); %(( Audio signal after compression ))

%==========================================================================
% Plotting Section - Original signal, Compressed signal, and the difference
%==========================================================================

% Plotting Audio signal before compression.
subplot(3,1,1);
plot([Sa]' , 'm')
legend('Original Signal')

% Plotting Audio signal after compression.
% Calculate the percentage of important coefficients in the signal.
% Compression percentage.
imp_cof_per_cent = (imp_cof / length(Sa_dct)) * 100;
subplot(3,1,2);
plot([Sa_comp]' , 'b')
legend(['Compressed Signal %' int2str(imp_cof_per_cent)])

% Plotting The difference between the two signals.
subplot(3,1,3);
plot([Sa - Sa_comp]' , 'k')
legend('The difference')

% Plotting all signals together.
figure
plot([Sa]' , 'm'), hold on
plot([Sa_comp]' , 'b'), hold on
plot([Sa - Sa_comp]' , 'k')
legend('Original Signal',...
['Compressed Signal %' int2str(imp_cof_per_cent)],'The difference')

%==========================================================================
% Performance Measures Section - Compression ratio, Retained energy and SNR 
%==========================================================================

% The ratio of the input signal to the output signal.
compression_ratio = imp_cof_per_cent;

% Indicates the amount of energy retained in the compressed signal 
% as a percentage of the energy of the original signal.
retained_signal_energy = norm((Sa_comp) / norm(Sa)) * 100;

% Get the percentage of zero coefficient in the compressed signal.
pr_zero_cof = 100 - compression_ratio;

% Calculate the signal to noise ratio.
SNR = snr(Sa , Sa - Sa_comp);

