
close all
% use different signal powers
mc = 1e4;

SNR = 0:0.5:5; %dB
N=100; %length of channel

P1=zeros(1, length(SNR));
P2=zeros(1, length(SNR));
P3=zeros(1, length(SNR));

P31=zeros(1, length(SNR));
P32=zeros(1, length(SNR));
P33=zeros(1, length(SNR));
P34=zeros(1, length(SNR));
P35=zeros(1, length(SNR));

for j=1:mc
    %% Problem 2
    %case 1
    for SNRi = 1:length(SNR);
        h1=(randn(1,1)+randn(1,1)*sqrt(-1))/sqrt(2);
        h2=(randn(1,1)+randn(1,1)*sqrt(-1))/sqrt(2);
        
        h = [h1 h2];

        if (sum(abs(h).^2)*SNR(SNRi)<1)
            P1(SNRi)=P1(SNRi)+1;
        end
    end
    
    %case 2
    for SNRi = 1:length(SNR);
        h1=(randn(1,1)+randn(1,1)*sqrt(-1))/sqrt(2);
        h2=(randn(1,1)+randn(1,1)*sqrt(-1))/sqrt(2)*(randn(1,1)+randn(1,1)*sqrt(-1))/sqrt(2);
        
        h = [h1 h2];

        if (sum(abs(h).^2)*SNR(SNRi)<1)
            P2(SNRi)=P2(SNRi)+1;
        end
    end
    
    %case 3
    for SNRi = 1:length(SNR);
        h1=(randn(1,1)+randn(1,1)*sqrt(-1))/sqrt(2);
        h2=((randn(1,1)+randn(1,1)*sqrt(-1))/sqrt(2)+(randn(1,1)+randn(1,1)*sqrt(-1))/sqrt(2))/2;
        
        h = [h1 h2];

        if (sum(abs(h).^2)*SNR(SNRi)<1)
            P3(SNRi)=P3(SNRi)+1;
        end
    end
    
    %% problem 3
    %case 1
    for SNRi = 1:length(SNR);
        h1=(randn(1,1)+randn(1,1)*sqrt(-1))/sqrt(2);
        h2=(randn(1,1)+randn(1,1)*sqrt(-1))/sqrt(2);
        
        h = [h1 h2].';

        if (sum(abs(h).^2)*SNR(SNRi)<1)
        P31(SNRi)=P31(SNRi)+1;
        end
    end
    
    %case 2 and 3
    for SNRi = 1:length(SNR);
        h1=(randn(1,1)+randn(1,1)*sqrt(-1))/sqrt(2);
        h2=(randn(1,1)+randn(1,1)*sqrt(-1))/sqrt(2);
        h3=(randn(1,1)+randn(1,1)*sqrt(-1))/sqrt(2);
        
        h = [h1 h2 h3].';

        if (sum(abs(h).^2)*SNR(SNRi)<1)
        P32(SNRi)=P32(SNRi)+1;
        end
    end
    
    %case 4
    for SNRi = 1:length(SNR);
        
        h = [];
        for loop=1:10
            h = [h (randn(1,1)+randn(1,1)*sqrt(-1))/sqrt(2)];
        end
        
        if (sum(abs(h).^2)*SNR(SNRi)<1)
        P34(SNRi)=P34(SNRi)+1;
        end
    end
    
    %case 5
    for SNRi = 1:length(SNR);
        
        h1=(randn(1,1)+randn(1,1)*sqrt(-1))/sqrt(2);
        h = [];
        for loop=1:15
            h = [h h1];
        end
        
        if (sum(abs(h).^2)*SNR(SNRi)<1)
        P35(SNRi)=P35(SNRi)+1;
        end
    end
    %case 5
    
end
P1=P1./mc;
P2=P2./mc;
P3=P3./mc;
P31=P31./mc;
P32=P32./mc;
P34=P34./mc;
P35=P35./mc;

figure('NumberTitle','off','Name','Problem 2, case 1,2,3');
semilogy(SNR, P1);
hold on
semilogy(SNR, P2);
hold on
semilogy(SNR, P3);
ylabel('Probability of deep fade')
xlabel('SNR (dB)')
grid on;

figure('NumberTitle','off','Name','Problem 3');
semilogy(SNR, P31);
hold on
semilogy(SNR, P32);
hold on
semilogy(SNR, P33);
hold on
semilogy(SNR, P34);
hold on
semilogy(SNR, P35);
ylabel('Probability of deep fade')
xlabel('SNR (dB)')
grid on;
legend('2x1', '3x1', '1x10', '15x1')

