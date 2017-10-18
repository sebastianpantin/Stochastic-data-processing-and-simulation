%% Lab 1, assignment 1.1a
clc;
clf;
sigma = 1;
mu = 0;
normalDist = normrnd(mu,sigma,[100,1]);
figure
histogram(normalDist,15);
title('Histogram')
figure
probplot('normal', normalDist);
title('pp-plot')
figure
qqplot(normalDist);
title('qq-plot');
[hchi,pchi] = chi2gof(normalDist);
[hks,pks] = kstest(normalDist);

%% Assignment 1.1b
clc;
clf;
a = 2;
b = 2;
gammaDist = gamrnd(a,b,[100,1]);
figure
histogram(gammaDist,15);
title('Histogram')
figure
probplot('normal', gammaDist);
title('pp-plot')
figure
qqplot(gammaDist);
title('qq-plot');
[hchi,pchi] = chi2gof(gammaDist);
[hks,pks] = kstest(gammaDist);

%% Assignment 1.2a
clear all;
sigma = 1;
mu = 0;
alpha = 0.05;
N = 1000;
isInRange = zeros(N,1);
for i = 1:N
    normalDist = normrnd(mu,sigma,[100,1]);
    [CI, intervalW] = sigmaCI(normalDist,alpha);
    stdVector(i) = 1^2;
    intervalWVector(i) = intervalW;
    if stdVector(i)> CI(1) && stdVector(i)<CI(2)
        isInRange(i) = 1;
    end
end
sum(isInRange)

%% Assignment 1.2b
clear all
clc;
clf;
a = 2;
b = 2;
N = 1000;
alpha = 0.05;
isInRange = zeros(N,1);
stdVector = a*b^2;
for i = 1:N
    gammaDist = gamrnd(a,b,[100,1]);
    [CI, intervalW] = sigmaCI(gammaDist,alpha);
    intervalWVector(i) = intervalW;
    if stdVector> CI(1) && stdVector<CI(2)
        isInRange(i) = 1;
    end
end
sum(isInRange)
%% Assignment 2.1
clear all
clc
clf
n = 100;
r = rand(n,1);
Y = normrnd(0,1, [n,1]);
Z = trnd(1,n,1);
eps = 0.05;

X=Y;
X(r < eps) = Z(r < eps);
figure
histogram(X,15)
title('Histogram')
figure
probplot('normal', X)
title('pp-plot')
figure
qqplot(X)
title('qq-plot')

[hchi,pchi] = chi2gof(X);

%% Assignment 2.2
clear all
clc
clf
n = 100;
N = 1000;
meanArray = zeros(N,1);
for i=1:N
    r = rand(n,1);
    Y = normrnd(0,1, [n,1]);
    Z = trnd(1,n,1);
    eps = 0.05;
    
    X=Y;
    X(r < eps) = Z(r < eps);
    meanArray(i) = mean(X);
end
sortedMeanArray = sort(meanArray);
mean25 = sortedMeanArray(25);
mean975 = sortedMeanArray(975);
figure
histogram(sortedMeanArray)
title('sorted means')
% se uppgift 1.1
%% Assignment 2.3
clear all
clc
clf
n = 100;
N = 1000;
meanArray = zeros(N,1);
medianArray = zeros(N,1);
k=0.1*n;
for i=1:N
    r = rand(n,1);
    Y = normrnd(0,1, [n,1]);
    Z = trnd(1,n,1);
    eps = 0.05;
    
    X=Y;
    X(r < eps) = Z(r < eps);
    X = sort(X);
    
    medianArray(i) = median(X);
    meanArray(i) = sum(X([k+1 n-k],1)/(n-2*k));
end
sortedMedianArray = sort(medianArray);
sortedMeanArray = sort(meanArray);
mean25 = sortedMeanArray(25);
mean975 = sortedMeanArray(975);
figure
histogram(sortedMeanArray)
title('sorted means')
figure
histogram(sortedMedianArray)
title('sorted median')




