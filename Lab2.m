%% Assignment 1.1a

clear all
clc
close all
dataSet = dlmread('stockdata.tsv');
stockName = {'AstraZenica'; 'Electrolux'; 'Ericsson'; 'Gambio';...
  'Nokia'; 'Swedish Match'; 'Svenska Handelsbanken';};

logReturns = zeros(size(dataSet));
logReturns(:,1) = dataSet(:,1);
logReturns = logReturns(2:end,:);
goodnessOfFitArray = zeros(2,7);
for i=2:size(dataSet,2)
    stockPrice = dataSet(:,i);
    logReturns(:,i) = GetStockReturns(stockPrice);
    
    figure
    subplot(3,1,1)
    plot(logReturns(:,1),logReturns(:,i))
    subplot(3,1,2)
    histogram(logReturns(:,i),15)
    subplot(3,1,3)
    qqplot(logReturns(:,i))
    [hks,pks] = kstest(logReturns(:,i));
    goodnessOfFitArray(:,i-1) = [hks;pks]; 
    
end

%% Assignment 1.1b

numberOfLags = 20;
ACFValues = zeros(numberOfLags,8);
ACFValuesAbs = zeros(numberOfLags,8);
ACFValues(:,1) = 1:numberOfLags;
ACFValuesAbs(:,1) = 1:numberOfLags;
absLogReturns = abs(logReturns);
for i=2:size(dataSet,2)
    for j=1:numberOfLags
        r = ACF(logReturns(:,i),j);
        rlog = ACF(absLogReturns(:,i),j);
        ACFValues(j,i) = r;
        ACFValuesAbs(j,i) = rlog;
    end
figure
hold on
plot(ACFValues(:,1),ACFValues(:,i),'.k','MarkerSize', 20)
title(['Sample ACF for stock: ' num2str(i-1)]);
xlabel('Lag')
ylabel('Sample auto correlation')
axis([0 20 -0.5 1]);
zeroLine = line([0 20],[0 0]);
zeroLine.Color='Black';
% maxLine = line([0 20],[max(ACFValues(:,i)),max(ACFValues(:,i))]);
% maxLine.Color='Black';
% minLine = line([0 20],[min(ACFValues(:,i)),min(ACFValues(:,i))]);
% minLine.Color='Black';
figure
hold on
plot(ACFValues(:,1),ACFValuesAbs(:,i),'.k','MarkerSize', 20)
title(['Sample ACF for absolute value of log-returns for stock: ' num2str(i-1)])
xlabel('Lag')
ylabel('Sample auto correlation')
axis([0 20 -0.5 1]);
zeroLine = line([0 20],[0 0]);
zeroLine.Color='Black';
% maxLine = line([0 20],[max(ACFValuesAbs(:,i)),max(ACFValuesAbs(:,i))]);
% maxLine.Color='Black';
% minLine = line([0 20],[min(ACFValuesAbs(:,i)),min(ACFValuesAbs(:,i))]);
% minLine.Color='Black';
end

%% Assignment 1.1c
meanLogReturns = zeros(2,7);
varLogReturns = zeros(2,7);
meanLogReturns = mean(logReturns(:,2:8));
varLogReturns = var(logReturns(:,2:8));
for i=2:size(dataSet,2)
    cov(logReturns(:,i), logReturns(:,2));
    R = corrcoef(logReturns(:,i),logReturns(:,2))
end

%% Assignment 1.2

delta = linspace(-2,2,10);
utility = zeros(1,10);
n=0;
for i=-3:1:3
    n=n+1;
    for j=1:10
        if(i>0)
            utility(n,j) = 1-exp(-i*delta(1,j));
        elseif (i~=0)
            utility(n,j) = exp(-i*delta(1,j))-1;
        else
            utility(n,j) = delta(1,j);
        end
    end
    hold on
    %plot(delta(1,:),utility(n,:))
    %axis([-2 2 -8 8])
end
k=0.5:0.1:15;
expectedUtil = zeros(1,7);
for i=1:length(k)
    for stock=1:7
        expectedUtil(i,stock)=1-exp(-k(i)*(meanLogReturns(stock)'-(k(i)*varLogReturns(stock))/2));
    end
end

for stock=1:7
    hold on
    plot(k',expectedUtil(:,stock))
end
%% Assignment 3

ericssonStock = logReturns(:,4);
gambioStock = logReturns(:,5);
meanEricssonLogReturns = mean(ericssonStock);
meanGambioLogReturns = mean(gambioStock);
meanVector = [meanEricssonLogReturns, meanGambioLogReturns]';
covMatrix = cov(ericssonStock, gambioStock);
w1 = 0:0.01:1;
w2 = 1-w1;
k=4;
for i=1:length(w1)
    w = [w1(i) w2(i)]';
    expectedUtil(i)=1-exp(-k*(meanVector'*w-k/2*w'*covMatrix*w));
end

plot(w1,expectedUtil);

fun = @(w)-(meanVector'*w-k/2*w'*covMatrix*w);
x0=[0.4 0.6]';
A=[];
b=[];
Aeq=[1 1];
beq=1;
lb = [0 0];
ub = [1 1];
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub);

%% Assignment 4

meanVector = mean(logReturns(:,2:8))';
covMatrix = cov(logReturns(:,2:8));
w1 = 0:0.01:1;
w2 = 1-w1;
K=0:0.1:10;

for i=1:length(K)
    k=K(i);
    fun = @(w)-(meanVector'*w-k/2*w'*covMatrix*w);
    expUtil = @(w)1-exp(-k*(meanVector'*w-k/2*w'*covMatrix*w));
    x0=[0 0 0 0 0 0 0]';
    A=[];
    b=[];
    Aeq=[1 1 1 1 1 1 1];
    beq=1;
    lb = [0 0 0 0 0 0 0];
    ub = [1 1 1 1 1 1 1];
    x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub);
    expectedUtil(i) = expUtil(x);
    expectedUtilNaive(i) = expUtil(ones(7,1).*1/7);
end
hold on
plot(K, expectedUtil)
plot(K,expectedUtilNaive)





