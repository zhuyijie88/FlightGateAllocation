filename_xlsx='InputData.xlsx';
filename_txt='InputData.txt';

load ../data_backup_set9_200_add_generate_random/each_20_max_plane_choose_by_random/best242/data250.mat;

[~, ~]=read_data_p1(filename_txt, filename_xlsx);
load plane_struct.mat;
load gate_struct.mat;
load passen_struct.mat;

x1d=x(1,:); % 对x1d进行分析

%% a. 给出成功分配到登机口的航班数量和比例，按宽、窄体机分别画线状图。

num_of_success_plane=numel(find(x1d<num_of_gate+1))

% input x1d
figure(1)
% wide plane
for i=1:1:size(x1d,2)
    if(x1d(i)<70 & plane_struct.u(i)==1)
    xstart=plane_struct.A(i);
    xend=plane_struct.D(i);
    
    xx=linspace(xstart,xend,floor(abs(xstart-xend)));
    yy=x1d(i)*ones(size(xx));
    mycolor=[rand rand rand];
    while(mean(mycolor)<0.5)
        mycolor=[rand rand rand];
    end
    
    plot(xx,yy,'.','MarkerSize', 25, 'Color', mycolor);
    drawnow; hold on;
    text((xstart+xend)/2,x1d(i),plane_struct.puck_id(i),'fontsize',10);
    end
end


set(gca,'linewidth',2,'FontSize',20,'FontName','Times New Roman');
set(gcf,'position',[60 0 1800 900]);
xlabel("Time/hour");
ylabel("Gate");
xlim([-45, 24*60+45]);
box on;

xtic=[0:60*2:60*24];
xlab=string(xtic/60);
xlab=strcat(xlab, ":00");
xlab(1)="20/Jan/18";
xlab(end)="21/Jan/18";

set(gca,'xtick',xtic);
set(gca,'xticklabel',xlab);   %style 3

ytic=[1:6:69];
ylab=gate_struct.gate_id(1:6:end)';
set(gca,'ytick',ytic);
set(gca,'yticklabel',ylab);   %style 3

figure(2);
% narrow plane
% input x1d
% wide plane
for i=1:1:size(x1d,2)
    if(x1d(i)<70 & plane_struct.u(i)==0)
    xstart=plane_struct.A(i);
    xend=plane_struct.D(i);
    
    xx=linspace(xstart,xend,floor(abs(xstart-xend)));
    yy=x1d(i)*ones(size(xx));
    mycolor=[rand rand rand];
    while(mean(mycolor)<0.5)
        mycolor=[rand rand rand];
    end
    
    plot(xx,yy,'.','MarkerSize', 25, 'Color', mycolor);
    drawnow; hold on;
    text((xstart+xend)/2,x1d(i),plane_struct.puck_id(i),'fontsize',10);
    end
end


set(gca,'linewidth',2,'FontSize',20,'FontName','Times New Roman');
set(gcf,'position',[60 0 1800 900]);
xlabel("Time/hour");
ylabel("Gate");
xlim([-45, 24*60+45]);
box on;

xtic=[0:60*2:60*24];
xlab=string(xtic/60);
xlab=strcat(xlab, ":00");
xlab(1)="20/Jan/18";
xlab(end)="21/Jan/18";

set(gca,'xtick',xtic);
set(gca,'xticklabel',xlab);   %style 3

ytic=[1:6:69];
ylab=gate_struct.gate_id(1:6:end)';
set(gca,'ytick',ytic);
set(gca,'yticklabel',ylab);   %style 3



%% b. 给出T和S登机口的使用数目和被使用登机口的平均使用率（登机口占用时间比率），要求画线状图。
% calculate the T S total number;
figure(2)
num_of_T=numel(unique(x1d(find(x1d<=28))));
num_of_S=numel(unique(x1d(find(x1d>=29&x1d<=69))));

time_used_of_gate=zeros(69,1);

for i=1:1:size(x1d,2)
    if(x1d(i)<70)
    xstart=plane_struct.A(i);
    xend=plane_struct.D(i);
    
    time_used_of_gate(x1d(i))=time_used_of_gate(x1d(i))+min(xend,24*60)-max(0,xstart);
    
%     xx=linspace(xstart,xend,floor(abs(xstart-xend)));
%     yy=x1d(i)*ones(size(xx));
%     hold on;
%     mycolor=[rand rand rand];
%     while(mean(mycolor)<0.5)
%         mycolor=[rand rand rand];
%     end
%     
%     plot(xx,yy,'.','MarkerSize', 25, 'Color', mycolor);
%     drawnow;    
%     text((xstart+xend)/2,x1d(i),plane_struct.puck_id(i));
    end
end
rate_used_of_gate=time_used_of_gate/(24*60);
for i=1:69
    mycolor=[rand rand rand];
    while(mean(mycolor)<0.5)
        mycolor=[rand rand rand];
    end
    
    x_rate=linspace(0,rate_used_of_gate(i),floor(1000*rate_used_of_gate(i)));
    y_rate=i*ones(size(x_rate));
    plot(x_rate,y_rate,'.','MarkerSize', 25, 'Color', mycolor);
    drawnow; hold on;
    text(rate_used_of_gate(i)+0.02,i,num2str(rate_used_of_gate(i),3),'fontsize',10);
end


set(gca,'linewidth',2,'FontSize',20,'FontName','Times New Roman');
set(gcf,'position',[60 0 1800 900]);
xlabel("Utilization rate");
ylabel("Gate");
xlim([0 1]);
box on;

ytic=[1:6:69];
ylab=gate_struct.gate_id(1:6:end)';
set(gca,'ytick',ytic);
set(gca,'yticklabel',ylab);   %sty
% 
% text(1,69,[num2str(num_of_T), "gate T"]) 
% text(1,60,[num2str(num_of_S), "gate S"]) 

%% c. 仅限问题二、三：给出换乘失败旅客数量和比率。
[time tensity time_of_passen each_of_passen time_gap_plane]=calculate_passenger_time(x1d,plane_struct,passen_struct,gate_struct,num_of_gate,num_of_plane, 3);
% 换乘失败的乘客数量和比率
fail_num=sum(each_of_passen(find((time_gap_plane-time_of_passen)<0)));
fail_rate=fail_num/sum(each_of_passen);
%% d. 仅限问题二、三：总体旅客换乘时间分布图：换乘时间在5分钟以内的中转旅客比率，10分钟以内的中转旅客比率，15分钟以内的中转旅客比率，……
figure(3)
time_transf=[5:5:100];
rate_time_transf=[];
for i_time=time_transf
    rate_time_transf=[rate_time_transf sum(each_of_passen(find(time_of_passen<=i_time & time_of_passen>i_time-5)))];
end
rate_time_transf=rate_time_transf/sum(rate_time_transf);    
list=find(rate_time_transf~=0);
time_transf=time_transf(list);
rate_time_transf=rate_time_transf(list);



for i=1:size(time_transf,2)
    mycolor=[rand rand rand];
    while(mean(mycolor)<0.5)
        mycolor=[rand rand rand];
    end
    
    x_rate=linspace(0.005,rate_time_transf(i),floor(1000*rate_time_transf(i)));
    y_rate=time_transf(i)*ones(size(x_rate));
    plot(x_rate,y_rate,'.','MarkerSize', 55, 'Color', mycolor);
    drawnow; hold on;
    text(rate_time_transf(i)+0.02,time_transf(i),num2str(rate_time_transf(i),3),'fontsize',20);
end




set(gca,'linewidth',2,'FontSize',20,'FontName','Times New Roman');
set(gcf,'position',[60 0 1800 900]);
xlabel("Rate");
ylabel("Transfer time");
xlim([0 1]);
box on;

ytic=time_transf;
ylab=strcat(string(ytic-5), "~", string(ytic));
set(gca,'ytick',ytic);
set(gca,'yticklabel',ylab);   %sty

ylim([min(time_transf)-5, max(time_transf)+5]);



%% e. 仅限问题二、三：总体旅客换乘紧张度分布图：紧张度在0.1以内的中转旅客比率，0.2以内的中转旅客比率，0.3以内的中转旅客比率，……
figure(4)
rate_of_time=time_of_passen./time_gap_plane;
time_transf=[0.1:0.1:1];
rate_time_transf=[];
for i_time=time_transf
    rate_time_transf=[rate_time_transf sum(each_of_passen(find(rate_of_time<=i_time & rate_of_time>i_time-0.1)))];
end
rate_time_transf=rate_time_transf/sum(rate_time_transf);    
list=find(rate_time_transf~=0);
time_transf=time_transf(list);
rate_time_transf=rate_time_transf(list);



for i=1:size(time_transf,2)
    mycolor=[rand rand rand];
    while(mean(mycolor)<0.5)
        mycolor=[rand rand rand];
    end
    
    x_rate=linspace(0.005,rate_time_transf(i),floor(1000*rate_time_transf(i)));
    y_rate=time_transf(i)*ones(size(x_rate));
    plot(x_rate,y_rate,'.','MarkerSize', 55, 'Color', mycolor);
    drawnow; hold on;
    text(rate_time_transf(i)+0.02,time_transf(i),num2str(rate_time_transf(i),3),'fontsize',20);
end




set(gca,'linewidth',2,'FontSize',20,'FontName','Times New Roman');
set(gcf,'position',[60 0 1800 900]);
xlabel("Rate");
ylabel("Time tension");
xlim([0 1]);
box on;

ytic=time_transf;
ylab=strcat(string(ytic-0.1), "~", string(ytic));
set(gca,'ytick',ytic);
set(gca,'yticklabel',ylab);   %sty

ylim([min(time_transf)-0.1, max(time_transf)+0.1]);

%% 5 generate a pair for Excel CSV
[~,Pucks]=xlsread(filename_xlsx,'Pucks');
New_Pucks=cell(size(Pucks,1),size(Pucks,2)+1);
New_Pucks(:,1:end-1)=Pucks;
New_Pucks(1,end)={'问题二登机口'};
puck_id=string(Pucks(2:end,1));
for i=1:size(x1d,2)
    if(x1d(i)<num_of_gate+1) % this is feasible
        gate_name=gate_struct.gate_id(x1d(i));
        plane_pk=plane_struct.puck_id(i);
        id_change=find(puck_id==plane_pk);
        New_Pucks(id_change+1,end)={gate_name};
    end
end




