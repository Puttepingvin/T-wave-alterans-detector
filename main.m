%Implements the alterans detector described in the attatched article
load('AlternansData.mat')
data = ecg_h(1,:);
time = t_h;

facit_h = [128:137 370:698];
facit_p1 = 110:332;
facit_p2 = [16:30 45:60 80:90 110:120 130:135 145:155 160:170	185:200	215:230	245:255	270:280	295:310];
facit = facit_p1;
close all
resall = {};
sensitivity_all = zeros(12,1);
specificity_all = zeros(12,1);
PPV_all = zeros(12,1);
NPV_all = zeros(12,1);
accuracy_all = zeros(12,1);
for m = 3:3
data = ecg_p1(m,:);
time = t_p1;
data_pre = data/mean(abs(data));

numpsegs = 7; %m in literature
Fs = 1000;
d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',46,'HalfPowerFrequency2',50, ...
               'DesignMethod','butter','SampleRate',Fs);

data = data/mean(abs(data));
data = filtfilt(d,data);
data = lowpass(data,100,Fs);

X = {};
lastR = 1;
len = find(time > length(data), 1)-2;
if isempty(len)
    len = length(time);
end
Q = zeros(len-1,1);
R = zeros(len-1,1);
knots = zeros(len-1,1);




for i = 1:len-1
    R(i) = time(i)+100 + find(data(time(i)+100:time(i+1)) == max(data(time(i)+100:time(i+1))));
    X{i} = data(lastR:R(i));
    lastR = R(i);
end
for i = 1:len-1
    Q(i) = length(X{i}) - 101 + find(X{i}(end-100:end) == min(X{i}(end-100:end)),1);
    ti = (Q(i)-50):Q(i)-1;
    slopes =X{i}(ti - 3) + X{i}(ti - 1) - X{i}(ti + 1) - X{i}(ti + 3);
    minslope = 0;
    for j = 1:length(slopes)
        if slopes(j) < 0 || slopes(j) > minslope
            minslope = slopes(j);
        end
        
    end
    slopes =X{i}(ti - 3) + X{i}(ti - 1) - X{i}(ti + 1) - X{i}(ti + 3);
    knots(i) = time(i) + Q(i) - length(ti) - 1 + find(slopes == minslope,1) - 66;
end


 spln = spline(knots, data(knots), 1:length(data));
  data = data-spln;


for i = 1:len-2
    X{i} = data(R(i):R(i+1));
end
dtau = zeros(len,1);
for i = 2:len-1
    dtau(i) = round((12*sqrt(length(X{i})) - (Q(i-1) - length(X{i-1})))/numpsegs);
    
end
windowsize = 10;
P = zeros(windowsize,7);
stepsize = 1;
alldists = [];
thresh = 0.7;
skips = 20;
result = R*0;
pred = result;
for i = (skips+windowsize/2):stepsize:(len-windowsize/2-1-skips)
    Xall = [];
    for j = -(windowsize/2):(windowsize/2)
        if length(X{i+j}) > 7*dtau(i+j) && length(X{i+j-1}) > 7*dtau(i+j-1) && dtau(i+j) > 0  && dtau(i+j-1) > 0
            P(j+1+windowsize/2,:) = [X{i+j}((1:7)*dtau(i+j)) - X{i+j-1}((1:7)*dtau(i+j-1))];
        end
    end
    odd = P(1:2:windowsize,:);
    even = P(2:2:windowsize,:);
	minst = min([size(odd,1) size(even,1)] );
    c1 = [mean(odd(1:minst-1,:),2), mean(even(1:minst-1,:),2)];
    c2 = [mean(even(1:minst-1,:),2),mean(odd(2:minst,:),2)];
    dist = norm(mean(c1) - mean(c2));
    if 1
        figure
        hold on
        plot(mean(c1(:,1)), mean(c1(:,2)), '*r', 'markersize', 20)
        plot(mean(c2(:,1)), mean(c2(:,2)), '*b', 'markersize', 20)
        for j = 1:2:windowsize
            plot(P(j),P(j+1), '.r');
            plot(P(j+1),P(j+2), '.b');
        end
        legend('Ceter of mass for odd beats','Center of mass for even beats', 'Odd beats', 'Even beats')
        xlabel('x_{i+1} - x_{i}')
        ylabel('x_{i+2} - x_{i+1}')
        grid
        axis equal
        hold off
        pause
    end
    if dist > thresh
        pred(i) = 1;
        if isempty(find(facit == i, 1))
            result(i) = -1;
        else
            result(i) = 1;
        end
    else
        pred(i) = -1;
        if isempty(find(facit == i, 1))
            result(i) = 1;
        else
            result(i) = -1;
        end
    end
    alldists = [alldists dist];
end
t = (1:length(alldists))*stepsize + (skips+windowsize/2) - (stepsize/2);
figure
plot(t,alldists)
hold on
plot(t,t./t*thresh)
hold off
xlabel('beat number')
ylabel('index')
title('TWA index over time for the human')
accuracy_all(m) = sum(result==1)/sum(result~=0);
sensitivity_all(m) = sum(result==1 & pred==1)/(sum(result==1 & pred==1) + sum(result==-1 & pred==-1));
specificity_all(m) = sum(result==1 & pred==-1)/(sum(result==-1 & pred==1) + sum(result==1 & pred==-1));
PPV_all(m) = sum(result==1 & pred==1)/(sum(result==1 & pred==1) + sum(result==-1 & pred==1));
NPV_all(m) = sum(result==1 & pred==-1)/(sum(result==-1 & pred==-1) + sum(result==1 & pred==-1));
resall{m} = result;
end
mean(accuracy_all)
mean(sensitivity_all)
mean(rmmissing(specificity_all))
mean(rmmissing(PPV_all))
mean(rmmissing(NPV_all))