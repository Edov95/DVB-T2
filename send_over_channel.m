function [r] = send_over_channel(s, SNR, type, R)
%SEND_OVER_CHANNEL Given the bitstream and the SNR add the noise to the 
% bistream

    SNRpbit = 10.^(SNR/10);   % Eb/No conversion from dB to decimal
    N0_uncoded = 1/SNRpbit;   % since Eb=1
    N0 = N0_uncoded/R;
    sigma = sqrt(N0/2);
    
    if(type = "AWGN")
        r = s + sigma*randn(1,length(s));
    end
    
    if(type == "BSC")
        r = s;
        Pbit = 1 - cdf('Normal',SNRpbit,0,sigma);
        ind = rand(1,length(s));
        r(find(ind > Pbit)) = -s(find(ind > Pbit));
    end
end

