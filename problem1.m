close all
clear

%% 1.
sigma_w2 = 6;
Mont = 1e4;
SNR_max = 70;
Error_sum = zeros(1,length(0:3:SNR_max));

for SNR = 0:3:SNR_max
    P = 10^(SNR/10)*sigma_w2;
    theta = 1/sqrt((16^2-1)/3)*sqrt(P);
    const = [-15 -13 -11 -9 -7 -5 -3 -1 1 3 5 7 9 11 13 15];
    Err = 0;
    for k = 1:Mont
        H = (randn(1,2)+1i*randn(1,2))/sqrt(2);

        %% 1. 16-PAM
        
        
        s1 = const(randperm(16,1));
        s2 = const(randperm(16,1));

        xtr = theta*[s1,0;0,s2];

        n = ((randn(1,2)+1i*randn(1,2))/sqrt(2))*sqrt(sigma_w2);
        y = H*xtr + n;

        x1_hat = y(1)/H(1,1)/theta;
        x2_hat = y(2)/H(1,2)/theta;

        diff1 = abs(real(x1_hat-const(1)));
        diff2 = abs(real(x2_hat-const(1)));
        idx1 = 1;
        idx2 = 1;
        % We are doing hard decision here. Our assumption is that this is
        % essentially the same as using the ML rule.
        for i = 2:length(const)
            if abs(real(x1_hat-const(i)))<diff1
                diff1 = abs(real(x1_hat-const(i)));
                idx1 = i;
            end
            if abs(real(x2_hat-const(i)))<diff2
                diff2 = abs(real(x2_hat-const(i)));
                idx2 = i;
            end
        end
        x1_hat = const(idx1);
        x2_hat = const(idx2);
        Error_sum(SNR/3+1) = Error_sum(SNR/3+1) + (x1_hat~=s1 || x2_hat~=s2);
    end
end
P_err = Error_sum/Mont;
figure('NumberTitle','off','Name','Problem 1, 1.')
semilogy(0:3:SNR_max,P_err)
%%

%% 2.
clear
sigma_w2 = 6;
Mont = 1e3;
SNR_max = 70;
Error_sum = zeros(1,length(0:3:SNR_max));
R = [1 1; 1 -1];

for SNR = 0:3:SNR_max
    P = 10^(SNR/10)*sigma_w2;
    theta = 1/sqrt(10)*sqrt(P);
    Err = 0;
    for k = 1:Mont
        H = (randn(1,2)+1i*randn(1,2))/sqrt(2);

        %% 1. 16-PAM
        
        %In our opinion, there is no need for rotation in this simulation
        %since there is no repetition in sending rotated codewords and
        %we are only interested in the number of errors.
        s1 = qammod(randperm(16,1)-1,16);
        s2 = qammod(randperm(16,1)-1,16);

        x = [s1 s2]*R;

        n = ((randn(1,2)+1i*randn(1,2))/sqrt(2))*sqrt(sigma_w2);
        y = theta*H*[x(1) 0; 0 x(2)] + n;
        
        x1_hat = y(1)/H(1,1)/theta;
        x2_hat = y(2)/H(1,2)/theta;
        
        s_hat = [x1_hat x2_hat]*pinv(R); 
        
        diff1 = abs(s_hat(1)-qammod(0,16));
        diff2 = abs(s_hat(2)-qammod(0,16));
        idx1 = 0;
        idx2 = 0;
        for i = 1:15
            dist1 = abs(s_hat(1)-qammod(i,16));
            if dist1<diff1
                diff1 = dist1;
                idx1 = i;
            end
            dist2 = abs(s_hat(2)-qammod(i,16));
            if dist2<diff2
                diff2 = dist2;
                idx2 = i;
            end
        end
        s1_hat = qammod(idx1,16);
        s2_hat = qammod(idx2,16);
        Error_sum(SNR/3+1) = Error_sum(SNR/3+1) + (s1_hat~=s1 || s2_hat~=s2);
    end
end
P_err = Error_sum/Mont;
figure('NumberTitle','off','Name','Problem 1, 2.')
semilogy(0:3:SNR_max,P_err)
%%

