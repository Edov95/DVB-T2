function [ ] = LDPCplot( x, ber_ldpc )


berColor = {'--sb','--+r','--dg','--om'};
leg = {'Simulated Pbit, Rate 1/2','Simulated Pbit, Rate 2/3','Simulated Pbit, Rate 3/4','Simulated Pbit, Rate 5/6'};

f = figure;

    
    EbN0dB = x; 
    % expected BER
    EbN0 = 10.^(EbN0dB/10);
    Q = @(y) 0.5*erfc((y)/sqrt(2));
    Pexp = Q(sqrt(5*EbN0))+4*Q(sqrt(6*EbN0))+12*Q(sqrt(7*EbN0));
        
    semilogy(EbN0dB(1),1,'r');
    hold on;
    for i=1:length(ber_ldpc(1,:))
        h1(i)=plot(EbN0dB,10*log10(ber_ldpc(:,i)),berColor{i}); 
    end
    semilogy(EbN0dB,Pexp);
    hold off; 
    h = [h1];
    h = h(:);    
    leg = leg(1:length(h));
    legend(h,leg);
    xlabel('E_b/N_0 [dB]');
    ylabel('BER');
    grid on;
%     saveas(f,'output/BERvsEbN0');
%     saveas(f,'output/BERvsEbN0','pdf');
end