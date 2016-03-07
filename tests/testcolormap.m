close all; 

nvals = 8; 
N = 100; 
x = linspace(0,1,N); 
cmp1 = parula(10); cmp2 = hot(10); cmp3 = copper(10); cmp4 = colorcube(30); cmp5 = bone(10); 

mycolormap = [cmp1(1,:); cmp2(3,:) ; cmp1(6,:); cmp3(5,:); cmp1(3,:);cmp4(21,:);cmp1(9,:);cmp5(2,:)];

%6 <-> 3%
mycolormap = [cmp1(1,:); cmp2(3,:) ; cmp4(21,:); cmp3(5,:); cmp1(3,:);cmp1(6,:);cmp1(9,:);cmp5(2,:)];
%7 <-> 5%
mycolormap = [cmp1(1,:); cmp2(3,:) ; cmp4(21,:); cmp3(5,:); cmp1(9,:);cmp1(6,:);cmp1(3,:);cmp5(2,:)];
%6 <-> 8%
mycolormap = [cmp1(1,:); cmp2(3,:) ; cmp4(21,:); cmp3(5,:); cmp1(9,:);cmp5(2,:);cmp1(3,:);cmp1(6,:)];
%1 <-> 8%
mycolormap = [cmp1(6,:); cmp2(3,:) ; cmp4(21,:); cmp3(5,:); cmp1(9,:);cmp5(2,:);cmp1(3,:);cmp1(1,:)];
%3 <-> 8%
mycolormap = [cmp1(6,:); cmp2(3,:) ; cmp1(1,:); cmp3(5,:); cmp1(9,:);cmp5(2,:);cmp1(3,:);cmp4(21,:)];


figure; 
cmp = mycolormap;
set(gcf,'Name',['mycolormap ' num2str(nvals)]);
set(gcf,'DefaultAxesColorOrder',cmp);


set(gcf,'Color','w');
set(gcf,'DefaultLineLineWidth',1.2);
set(gcf,'DefaultAxesLineWidth',1.0);
set(gcf,'DefaultAxesFontName','Arial');
set(gcf,'DefaultTextFontName','Arial');

hold on;


for i = 1:nvals 
    
    plot(x,i*ones(N,1));
   
end
ylim([0,nvals+1]);
hold off; 


