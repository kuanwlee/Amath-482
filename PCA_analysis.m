clear all
close all
clc

load('cam1_1.mat')
load('cam2_1.mat')
load('cam3_1.mat')

% implay(vidFrames1_1)

I = rgb2gray(vidFrames1_1(:,:,:,1));

objectRegion = [315,215,20,20];
objectRegion2 = [265,285,20,20];
objectRegion3 = [330,280,20,20];
 
position1 = trackobj(objectRegion,vidFrames1_1);
position2 = trackobj(objectRegion2,vidFrames2_1);
position3 = trackobj(objectRegion3,vidFrames3_1);
%%
l = length(position1(1,:));
start2 = 20;
start3 = 6;

X1 = [position1(1,:);position1(2,:);...
    position2(1,start2:start2+l-1);position2(2,start2:start2+l-1);...
    position3(1,start3:start3+l-1);position3(2,start3:start3+l-1)];

[U1,S1,V1] = svd(X1,'econ');

plot(log(diag(S1)))
title('SIngular values')
ylabel('log(S)')
xlabel('i')

figure(2)
hold on
x = 1:1:length(V1(:,1));
plot(x,V1(:,1),x,V1(:,2),x,V1(:,3),x,V1(:,4),'LineWidth',2)
legend('First mode','Second mode','Third mode','Fourth mode')

title('temporal mode')


%% Shaking camera
clear all
close all
clc

load('cam1_2.mat')
load('cam2_2.mat')
load('cam3_2.mat')

% implay(vidFrames1_2)

start1 = 15;
start2 = 1;
start3 = 19;
objectRegion1 = [325,315,30,30];
objectRegion2 = [280,325,50,50];
objectRegion3 = [370,250,30,30];
 
position1 = trackobj(objectRegion1,vidFrames1_2);
position2 = trackobj(objectRegion2,vidFrames2_2(:,:,:,5:end));
position3 = trackobj(objectRegion3,vidFrames3_2);

%%
l = 150;
X2 = [position1(1,start1:start1+l-1);position1(2,start1:start1+l-1);...
    position2(1,start2:start2+l-1);position2(2,start2:start2+l-1);...
    position3(1,start3:start3+l-1);position3(2,start3:start3+l-1)];

[U2,S2,V2] = svd(X2,'econ');

plot(log(diag(S2)),'LineWidth',2)
title('Singular values')
ylabel('log(S)')
xlabel('i')

figure(2)
hold on
x = 1:1:length(V2(:,1));
plot(x,V2(:,1),x,V2(:,2),x,V2(:,3),x,V2(:,4),'LineWidth',2)
legend('First mode','Second mode','Third mode','Fourth mode')

title('temporal mode')
%%
clear all
close all
clc

load('cam1_3.mat')
load('cam2_3.mat')
load('cam3_3.mat')

start1 = 19;
start2 = 46;
start3 = 10;

objectRegion1 = [325,300,30,30];
objectRegion2 = [260,320,20,20];
objectRegion3 = [360,220,20,20];

position1 = trackobj(objectRegion1,vidFrames1_3);
position2 = trackobj(objectRegion2,vidFrames2_3(:,:,:,1:end));
position3 = trackobj(objectRegion3,vidFrames3_3);
 
% implay(vidFrames1_3)
%%
start1 = 19;
start2 = 46;
start3 = 10;
l = 150;

X3 = [position1(1,start1:start1+l-1);position1(2,start1:start1+l-1);...
    position2(1,start2:start2+l-1);position2(2,start2:start2+l-1);...
    position3(1,start3:start3+l-1);position3(2,start3:start3+l-1)];

[U3,S3,V3] = svd(X3,'econ');

plot(log(diag(S3)),'LineWidth',2)
title('Singular values')
ylabel('log(S)')
xlabel('i')

figure(2)
hold on
x = 1:1:length(V3(:,1));
plot(x,V3(:,1),x,V3(:,2),x,V3(:,3),x,V3(:,4),'LineWidth',2)
legend('First mode','Second mode','Third mode','Fourth mode')

title('temporal mode')

%%
clear all
close all
clc

load('cam1_4.mat')
load('cam2_4.mat')
load('cam3_4.mat')

start1 = 19;
start2 = 46;
start3 = 10;

objectRegion1 = [380,290,30,30];
objectRegion2 = [260,320,20,20];
objectRegion3 = [360,220,20,20];

position1 = trackobj(objectRegion1,vidFrames1_3);
position2 = trackobj(objectRegion2,vidFrames2_3(:,:,:,1:end));
position3 = trackobj(objectRegion3,vidFrames3_3);
%%
objectRegion2 = [260,270,30,30];

trackobj(objectRegion2,vidFrames2_4(:,:,:,1:end));


 %%
function [position] = trackobj(objectRegion,frames)

v = VideoWriter('motion.avi');
open(v)
writeVideo(v,frames)
close(v)
videoFileReader = vision.VideoFileReader('motion.avi');
videoPlayer = vision.VideoPlayer('Position',[100,100,680,520]);
objectFrame = videoFileReader();

objectImage = insertShape(objectFrame,'Rectangle',objectRegion,'Color','red');
figure;
imshow(objectImage);
title('Red box shows object region');

% points = detectMinEigenFeatures(rgb2gray(objectFrame),'ROI',objectRegion);
x = objectRegion(1):1:objectRegion(1)+objectRegion(3);
y = objectRegion(2):1:objectRegion(2)+objectRegion(4);
[X,Y] = meshgrid(x,y);
points(1,:) = X(:);
points(2,:) = Y(:);
pointImage = insertMarker(objectFrame,points','+','Color','blue');
figure;
imshow(pointImage);
title('Detected interest points');
tracker = vision.PointTracker('MaxBidirectionalError',6);
blkmatcher = vision.BlockMatcher();
initialize(tracker,points',objectFrame);
k = 1;
while ~isDone(videoFileReader)
      frame = videoFileReader();
      [points,validity] = tracker(frame);
      posi = points;
      position(1,k) = mean(posi(:,1));
      position(2,k) = mean(posi(:,2));
      out = insertMarker(frame,points(validity, :),'+');
      videoPlayer(out);
      pause(0.1)
      k = k+1;
end

end
