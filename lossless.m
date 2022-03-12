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
% Encoding Section - RLE of the sampled signal
%==========================================================================

% Apply Run Length Encoding to the sampled audio signal.

% Set true for values that change
% and set zero for values that repeated.
% This for removing the consecutive repetitions in the audio signal.
A = [true, diff(Sa) ~= 0];

% Store number of consecutive repetitions.
% This for saving the number of consecutive repetitions
% to reconstruct the signal again later.
N = diff(find([A,true]));

% The audio signal after removing the consecutive repetitions.
Sa_comp = Sa(A); %(( Audio signal after compression ))

%==========================================================================
% Decoding Section - IRLE of the sampled signal
%==========================================================================

% Apply Inverse Run Length Encoding to signal.

% The reconstructed signal.
Sa_reconst = [];
% Counter to set the consecutive repetitions in their original index.
c = 0;
% For loop to merge the compressed sampled signal
% with the saved number of consecutive repetitions.
for i = 1 : 1 : length(Sa_comp)
    
    for j = 1 : 1 : N(i) % Loop by the times of the current repetition.
        
        % spread the current element (Sa_comp(i)) by the number of its 
        % saved consecutive repetition to achieve decompression.
        Sa_reconst(j + c) = Sa_comp(i);
    end
    
    % Add to the counter the previous number of repetition to put
    % the next saved values in their right index.
    c = c + N(i);
end

% Compression percentage.
comp_per_cent = length(Sa_comp) / length(Sa) * 100;

%==========================================================================
% Plotting Section - Original signal, Compressed signal, and the difference
%==========================================================================

% Plotting Audio signal before compression.
subplot(3,1,1);
plot([Sa]' , 'm')
legend('Original Signal')

% Plotting Audio signal after compression.
% Use the compressed signal with the consecutive repetitions set as zeros
% to be equally dimensioned with the original signal for plotting.
Sa_comp_wz = ( Sa .* A );
subplot(3,1,2);
plot([Sa_comp_wz]' , 'b')
legend(['Compressed Signal %' int2str(comp_per_cent)])

% Plotting The difference between the two signals.
subplot(3,1,3);
plot([Sa - Sa_comp_wz]' , 'k')
legend('The difference')

% Plotting all signals together.
figure
plot([Sa]' , 'm'), hold on
plot([Sa_comp_wz]' , 'b'), hold on
plot([Sa - Sa_comp_wz]' , 'k')
legend('Original Signal',...
['Compressed Signal %' int2str(comp_per_cent)],'The difference')

%==========================================================================
% Performance Measures Section - Compression ratio, Retained energy and SNR 
%==========================================================================

% the ratio of the input signal to the output signal.
compression_ratio = comp_per_cent;

% indicates the amount of energy retained in the compressed signal 
% as a percentage of the energy of the original signal.
retained_signal_energy = norm((Sa_comp) / norm(Sa)) * 100;

% get the percentage of zero coefficient in the compressed signal.
pr_zero_cof = 0; %The audio signal has been reconstructed without any loss.

% calculate the signal to noise ratio.
SNR = snr(Sa , Sa - (Sa .* A));

