close all
clear
%% 4
%1st question
x=0:0.00001:10;
expApprox=exp(-(x.^2)./2);
figure, plot(x, qfunc(x));
hold on;
plot(x, expApprox);
hold on;
plot(x, abs(expApprox-qfunc(x)), 'k--');
grid on;
lgd=legend('$Q(\alpha)$', '$e^{-\alpha^2/2}$', '$|Q(\alpha) - e^{-\alpha^2/2}|$' );
set(lgd,'Interpreter','latex');
xlabel('x')
%2nd question
smallerThan = 0:0.001:2;
results=zeros(1, length(smallerThan));
rep=10000;

for i=1:rep
    htest=(randn(1,1)+randn(1,1)*sqrt(-1));
    for j=1:length(smallerThan)
        if((abs(htest)^2)<smallerThan(j))
            results(j)=results(j)+1;
        end
    end
    
end
        
figure, plot(smallerThan, results./rep);
hold on;
plot(smallerThan, smallerThan);
hold on;
plot(smallerThan, abs(results./rep-smallerThan), 'k--');
grid on;
lgd=legend('$P(||h||^2<\epsilon)$', '$\epsilon$', '$|P(||h||^2<\epsilon) - \epsilon|$' );
xlabel('\epsilon')
ylabel('f(\epsilon)')
set(lgd,'Interpreter','latex');
daspect([1 1 2]);
%3rd question

smallerThan = 0:0.001:100;
results=zeros(1, length(smallerThan));

rep=50000;

for i=1:rep
    htest=(10*randn(1,1)+10*randn(1,1)*sqrt(-1));
    for j=1:length(smallerThan)
        if((abs(htest)^2)<smallerThan(j))
            results(j)=results(j)+1;
        end
    end
    
end
        
figure;
plot(smallerThan, results./rep);
hold on
plot(smallerThan, smallerThan/250)
hold on
plot(smallerThan, abs(results./rep-smallerThan/250), 'k--');
grid on;
lgd=legend('$P(||h||^2<\epsilon)$' ,'$\epsilon/100$','$|P(||h||^2<\epsilon) - \epsilon|$');
xlabel('\epsilon')
ylabel('f(\epsilon)')
set(lgd,'Interpreter','latex');
%daspect([1 1 2]);
%xlim([0 0.5])

%4th question
k=1:3;
rep=50000;
smallerThan = 0:0.01:2;
results=zeros(length(k), length(smallerThan));

figure;
for t=1:length(k)
    for i=1:rep
        hSqtest = chi2rnd(2*k(t));
        for j=1:length(smallerThan)
            if(hSqtest<smallerThan(j))
                results(t,j)=results(t,j)+1;               
            end
        end
    end
    plot(smallerThan, results(t,:)./rep);
    hold on;
end
legend('k=1', 'k=2', 'k=3');

figure, plot(smallerThan, gammainc(smallerThan./2,k(1),'lower'));
hold on;
plot(smallerThan, smallerThan.^k(1))
hold on;
plot(smallerThan, abs(smallerThan.^k(1)-gammainc(smallerThan./2,k(1),'lower')), '--k')
grid on;
lgd=legend('$P(||h||^2<\epsilon)$', '$\epsilon$', '$|P(||h||^2<\epsilon) - \epsilon|$');
xlabel('\epsilon')
ylabel('f(\epsilon)')
set(lgd,'Interpreter','latex');
daspect([1 1 2]);

figure, plot(smallerThan, gammainc(smallerThan./2,k(2),'lower'));
axis([-inf inf 0 0.7]);
hold on;
plot(smallerThan, smallerThan.^k(2))
hold on;
plot(smallerThan, abs(smallerThan.^k(2)-gammainc(smallerThan./2,k(2),'lower')), '--k')
grid on;
lgd=legend('$P(||h||^2<\epsilon)$', '$\epsilon^2$', '$|P(||h||^2<\epsilon) - \epsilon^2|$' );
xlabel('\epsilon')
ylabel('f(\epsilon)')
set(lgd,'Interpreter','latex');

figure, plot(smallerThan, gammainc(smallerThan./2,k(3),'lower'));
axis([-inf inf 0 0.7]);
hold on;
plot(smallerThan, smallerThan.^k(3))
hold on;
plot(smallerThan, abs(smallerThan.^k(3)-gammainc(smallerThan./2,k(3),'lower')), '--k')
grid on;
lgd=legend('$P(||h||^2<\epsilon)$', '$\epsilon^3$', '$|P(||h||^2<\epsilon) - \epsilon^3|$'  );
xlabel('\epsilon')
ylabel('f(\epsilon)')
set(lgd,'Interpreter','latex');