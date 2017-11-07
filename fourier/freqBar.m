function [ freq ] = freqBar( abs_FFT, abs_Axis )
    freq = sum(abs_Axis.*abs_FFT.^2)/sum(abs_FFT.^2);
end

%moyenne ponderee : moy_pond_freq = mean(abs_Axis)/sum(abs_FFT.^2);
