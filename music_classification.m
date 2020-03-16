clear all
close all
clc

files = dir(fullfile('Adele','*.wav'));
cd Adele
for i = 1:length(files)
    [y,Fs] = audioread(files(i).name);
    y = mean(y,2);
    y = y(1:10:end);
    [s,w,t] = spectrogram(y,500,250);
    Adele(:,i) = s(:);
end
cd ../

files = dir(fullfile('Linkinp','*.wav'));
cd Linkinp
for i = 1:length(files)
    [y,Fs] = audioread(files(i).name);
    y = mean(y,2);
    y = y(1:10:end);
    [s,w,t] = spectrogram(y,500,250);
    Linkp(:,i) = s(:);
end
cd ../

files = dir(fullfile('Bach','*.wav'));
cd Bach
for i = 1:length(files)
    [y,Fs] = audioread(files(i).name);
    y = mean(y,2);
    y = y(1:10:end);
    [s,w,t] = spectrogram(y,500,250);
    Bach(:,i) = s(:);
end
cd ../

musclass = [Adele Linkp Bach];

%%
clear all 
close all
clc

files = dir(fullfile('test1','*.wav'));
cd test1
for i = 1:length(files)
    [y,Fs] = audioread(files(i).name);
    y = mean(y,2);
    y = y(1:10:end);
    [s,w,t] = spectrogram(y,500,250);
    test1(:,i) = s(:);
end
cd ../
save('test1.mat','test1')


%%
clear all
close all
clc

load('Adele.mat')
load('Linkp.mat')
load('Bach.mat')

feature = 10;
[U,S,V,w,sort1,sort2,sort3] = music_trainer(abs(Adele),abs(Linkp),abs(Bach),feature);

plot(sort1,zeros(1,length(sort1)),'bo','LineWidth',2)
hold on
% histogram(sort2,20)
plot(sort2,zeros(1,length(sort1)),'ro','LineWidth',2)
% histogram(sort3,20)
plot(sort3,zeros(1,length(sort1)),'ko','LineWidth',2)
legend('Adele','Linkin Park','Bach')

title('Projection of data')
threshold1 = 45;
threshold2 = 100 ;
y = -5:0.01:5;
plot(threshold1*ones(1,length(y)),y,'r-')
plot(threshold2*ones(1,length(y)),y,'m-')

%%
load('test1.mat')

testmat = U'*abs(test1);
pval = w'*testmat;

answer = [1 1 1 3 3 3 2 2 2];

classify = zeros(1,length(pval));

classify(pval<=threshold1) = 3;

classify((threshold1<=pval)&(pval<=threshold2)) = 2;
classify(pval>=threshold2) = 1;


%%
function [U,S,V,w,sort1,sort2,sort3] = music_trainer(mus1,mus2,mus3,feature)
    nd = size(mus1,2);
    
    [U,S,V] = svd([mus1 mus2 mus3],'econ');
    
    musics = S*V'; % projection onto principal components
    U = U(:,1:feature);
    class1 = musics(1:feature,1:nd);
    class2 = musics(1:feature,nd+1:nd+nd);
    class3 = musics(1:feature,2*nd+1:end);
    
    m1 = mean(class1,2);
    m2 = mean(class2,2);
    m3 = mean(class3,2);
    
    Sw = 0; % within class variances
  
    for k=1:nd
        Sw = Sw + (class1(:,k)-m1)*(class1(:,k)-m1)';
    end
 
    for k=1:nd
        Sw = Sw + (class2(:,k)-m2)*(class2(:,k)-m2)';
    end
    
    for k=1:nd
        Sw = Sw + (class3(:,k)-m3)*(class3(:,k)-m3)';
    end
    x_bar = (nd*m1+nd*m2+nd*m3)/(3*nd);
    Sb = nd*(m1-x_bar)*(m1-x_bar)'; % between class 
    Sb = Sb + nd*(m2-x_bar)*(m2-x_bar)';
    Sb = Sb + nd*(m3-x_bar)*(m3-x_bar)';
    
    [V2,D] = eig(Sb,Sw); % linear discriminant analysis
    [~,ind] = max(abs(diag(D)));
    w = V2(:,ind); w = w/norm(w,2);
    
    v1 = w'*class1; 
    v2 = w'*class2;
    v3 = w'*class3;
    
%     if mean(vdog)>mean(vcat)
%         w = -w;
%         vdog = -vdog;
%         vcat = -vcat;
%     end
    % dog < threshold < cat
    
    sort1 = sort(v1);
    sort2 = sort(v2);
    sort3 = sort(v3);
    
%     t1 = length(sortdog);
%     t2 = 1;
%     while sortdog(t1)>sortcat(t2)
%         t1 = t1-1;
%         t2 = t2+1;
%     end
%     threshold = (sortdog(t1)+sortcat(t2))/2;
end