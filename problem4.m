clc
close all
clear

%% 4
%1st question
x=0:0.0001:3;
%densGauss = 1/sqrt(2*pi)*exp(-(x).^2);
expApprox=exp(-(x.^2)./2);
figure, plot(x, qfunc(x));
hold on;
plot(x, expApprox);
hold on;
plot(x, abs(expApprox-qfunc(x)), '--');
grid on;
lgd=legend('Q(x)', '$e^{-x^2/2}$', '$|Q(x) - e^{-x^2/2}|$' );
set(lgd,'Interpreter','latex');

%2nd question
smallerThan = 0.01:0.01:2;
%smallerThan = fliplr(smallerThan);
results=zeros(1, length(smallerThan));

rep=5000;
% 
% hArray=[];
% for i=1:rep 
%     hArray = [hArray (randn(1,1)+randn(1,1)*sqrt(-1))];
% end
% 
%     figure, histogram(abs(hArray), 'Normalization','probability');

for i=1:rep
    htest=(randn(1,1)+randn(1,1)*sqrt(-1));
    
    for j=1:length(smallerThan);
        if((abs(htest)^2)<smallerThan(j))
            results(j)=results(j)+1;
           % break;
        end
    end
    
end
        
figure, plot(smallerThan, results./rep);
hold on;
plot(smallerThan, smallerThan);
hold on;
plot(smallerThan, abs(results./rep-smallerThan), '--');
grid on;
lgd=legend('$P(||h||^2<x)$', '$x$', '$|P(||h||^2<x) - x|$' );
set(lgd,'Interpreter','latex');

%3rd question

smallerThan = 0.1:0.1:100;
%smallerThan = fliplr(smallerThan);
results=zeros(1, length(smallerThan));

rep=5000;

for i=1:rep
    htest=(10*randn(1,1)+10*randn(1,1)*sqrt(-1));
    
    for j=1:length(smallerThan);
        if((abs(htest)^2)<smallerThan(j))
            results(j)=results(j)+1;
           % break;
        end
    end
    
end
        
figure, plot(smallerThan, results./rep);
hold on;
plot(smallerThan, smallerThan);
hold on;
plot(smallerThan, abs(results./rep-smallerThan), '--');
grid on;
lgd=legend('$P(||h||^2<x)$', '$x$', '$|P(||h||^2<x) - x|$' );
set(lgd,'Interpreter','latex');

%4th question
k=1:3;
rep=5000;
smallerThan = 0.01:0.01:2;
results=zeros(length(k), length(smallerThan));

figure;
for t=1:length(k)
    
    for i=1:rep
       %htest = chi2rnd(2*k(t)*ones(1,k(t))); %[2k 2k 2k ... as dof] k times
     
        htest = chi2rnd(2*k(t));
        for j=1:length(smallerThan);
            if((sum(htest.^2))<smallerThan(j))
                results(t,j)=results(t,j)+1;
               
            end
        end
    end
    plot(smallerThan, results(t,:)./rep);
    hold on;
    
end
legend('k=1', 'k=2', 'k=3');