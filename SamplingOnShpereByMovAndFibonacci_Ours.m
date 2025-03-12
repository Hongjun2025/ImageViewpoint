% https://extremelearning.com.au/how-to-evenly-distribute-points-on-a-sphere-more-effectively-than-the-canonical-fibonacci-lattice/
% 另一种比较“美观”的做法是构造斐波那契网格，
%本质是上区域上的斐波那契网格经过保面积变换投影到球面上的样子
% SamplingOnShpereByMovAndFibonacci_Ours.m
%补充根据距离移动

%parameters
mvStep = 0.0017;
itNum = 3773; %the number of itrating
n = 900; %
%---------
tic

goldenRatio = (1 + 5^0.5)/2;
Li = linspace(0, n-1,n);
Li = Li';
theta = 2*pi * Li / goldenRatio;
phi = acos(1 - 2*(Li+0.5)/n ); %acos(1 - 2*(Li)/n); %

x = cos(theta) .* sin(phi) ; 
y = sin(theta) .* sin(phi) ; 
z = cos(phi) ;
% x = cos(theta) .* cos(phi) ; 
% y = sin(theta) .* cos(phi) ; 
% z = sin(phi) ;
Pts=[x(:),y(:),z(:)];   
PtsInit = Pts;   
ptNum = size(Pts(:,1),1);


AdjN = 13;%the number of neighbours
% %移动重回的点
% ChongHeN = 1;
% while ChongHeN > 0
%     [NeiID, dists] = knnsearch(Pts,Pts,'K',2 );%cosine
%     idx0 = find(dists(:,2)<0.000001);
%     ChongHeN = length(idx0)
%     if(ChongHeN<1)
%         break;
%     end
%     Pts(idx0,:) = Pts(idx0,:) + (rand(ChongHeN,3)-0.5)*mvStep;
% end
%
flag = 1;
while itNum > 0  %moving
    itNum = itNum-1;
    [NeiID, dists] = knnsearch(Pts,Pts,'K',2 );%cosine
    mvStep = 0.3*0.5*(max(dists(:,2))-min(dists(:,2)));
    if flag == 1
        pq = (Pts - Pts(NeiID(:,2),:));%最近邻反方向 
        pq0 = pq./( repmat(dists(:,2),1,3) );%每行单位向量 
        Pts = Pts + mvStep*pq0;%Moving
    else
        pq = rand(ptNum,3)-0.5;%random directions
        pnLen = sqrt(sum(pq.^2,2));
        pq0 = pq./( repmat(pnLen,1,3) );%每行单位向量 
        Pts = Pts + mvStep*pq0;%Moving
    end
    flag = 1-flag;
    NewL = sqrt(sum(Pts.^2,2));
    Pts = Pts./( repmat(NewL,1,3) );%New position
%     Dens= EvaluationByNeighbourDens(Pts, AdjN);
     %h = 
 
end
TmOurs = toc

%vis-------------------
fig=figure('Color','w');
[X,Y,Z] = sphere;
% S = surf(X,Y,Z,'EdgeColor', 'none','FaceColor',[0.0,0.96,0.85]);
S = surf(X,Y,Z,'EdgeColor', [0.,0.0,0.9],'FaceColor',[0.,0.9,0.0]);
hold on
axis equal
lighting('gouraud') % 使用Phong光照模型（高亮反射）
lightangle(90,30);
view(-45,30);
set(fig, 'Alpha', 0.5);
plot3(PtsInit(:,1), PtsInit(:,2),PtsInit(:,3),'y.');
plot3(Pts(:,1), Pts(:,2),Pts(:,3),'r.');



%Evaluation, by adjacent points density analysis-----
% AdjN = 15;%the number of neighbours
alpha = 0.05;
AveNeibDistOld = EvaluationByNeighbourDens(PtsInit, AdjN);
[UniHksOld,UniPOld] = MyKStestforUniformDistribution(AveNeibDistOld,alpha);
[NrmHksOld,NrmPOld] = MyKStestforNormalDistribution(AveNeibDistOld,alpha);

AveNeibDist= EvaluationByNeighbourDens(Pts, AdjN);
%test: K-S
[UniHks,UniP] = MyKStestforUniformDistribution(AveNeibDist,alpha);
[NrmHks,NrmP] = MyKStestforNormalDistribution(AveNeibDist,alpha);
%%%%-------------


figHist=figure('Color','w');
subplot(1,2,1); hist(AveNeibDistOld); title('DensOld');
subplot(1,2,2); hist(AveNeibDist); title('DensNew');

%save to file
SVTF = 1;
if SVTF>0
    svFlder=['Results',datestr(now,'yyyymmdd')];
    mkdir(svFlder);
    svFlder=[svFlder,'\'];
    svFile0=['SamplingOnSphere',num2str(ptNum),'OursTm',num2str(TmOurs,6)];
    svFile1=[svFlder,svFile0,'.xls']
    xlswrite(svFile1,{'Fx','Fy','Fz','OurX','OurY','OurZ'},1,'A1:F1' );
    xlswrite(svFile1,[PtsInit,Pts],1,['A2:F',num2str(ptNum+1)] );
    IdxInfoStr ={'minAveD','maxAveD','stdAveD','UniHks',...
        'UniP','NrmHks','NrmP','alpha'};
    IdxInfoValOld =[min(AveNeibDistOld), max(AveNeibDistOld),std(AveNeibDistOld),UniHksOld,...
        UniPOld, NrmHksOld,NrmPOld,alpha];
    IdxInfoVal =[min(AveNeibDist), max(AveNeibDist),std(AveNeibDist),UniHks,...
        UniP, NrmHks,NrmP,alpha];
    xlswrite(svFile1,IdxInfoStr,2,'A1:G1' );
    xlswrite(svFile1,IdxInfoValOld,2,'A2:G2' );
    xlswrite(svFile1,IdxInfoVal,2,'A3:G3' );
    xlswrite(svFile1,{'Fibo';'Ours'},2,'H2:H3' );
%svFile2=[svFile0,'.png'], saveas(gcf,svFile2 );
%saveas(fig1,svFile2 );
end