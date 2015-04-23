%% h2 MIMO FRF estimator
% This function estimate a MIMO FRF by using the H2 Algorithm.
%
% input_time_data is an array of matrices; each matrix represents the input
% for a single experimental test, each row is a time-sample, each
% collumn represent a different system input;
%
% output_time_data is an array of matrices; each row is a time-sample, each
% collumn represent an input;
function [freq FRF] = mimo_frf_H2(input_time_data, output_time_data, dt)
    N_measures = length(input_time_data);
    N_inputs = size(input_time_data{1},2);
    N_outputs = size(output_time_data{1},2);
    NFFT = size(input_time_data{1},1);
    for i_inp = 1:N_inputs
        for j_out = 1:N_outputs
            FRF{j_out, i_inp} = H2_FRF(i_inp, j_out, input_time_data, output_time_data);
            FRF{j_out, i_inp} = FRF{j_out, i_inp}(1:NFFT/2+1)
        end        
    end
    Fs = 1/dt;
    freq = Fs/2*linspace(0,1,NFFT/2+1);
end

function FRFji = H2_FRF(i_inp, j_out, input_time_data, output_time_data)
    N_measures = length(input_time_data);
    % cross-power spectra     
     % GUYji = mean (Ui .* conj(Yj))
    % auto-power spectra
     % GYYjj = mean (Yj .* conj(Yj))
     
     GUYji = zeros(size(input_time_data{1},1),1);
     GYYjj = zeros(size(input_time_data{1},1),1);
    for i = 1:N_measures
        GUYji = GUYji + fft(input_time_data{i}(:,i_inp)).*conj(fft(output_time_data{i}(:,j_out)))/N_measures;
        GYYjj = GYYjj + fft(output_time_data{i}(:,j_out)).*conj(fft(output_time_data{i}(:,j_out)))/N_measures;
    end
    FRFji = GYYjj ./ GUYji;
end