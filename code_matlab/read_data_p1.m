function [plane_struct gate_struct]=read_data_p1(filename_txt, filename_xlsx)
% INPUT: filename_txt is the name of txt file
% INPUT: filename_xlsx is the name of xlsx file
% OUTPUT: plane_struct stores all informations of planes
% OUTPUT: gate_struct stores all informations of gates
date1=datetime('18/1/19');
date2=datetime('18/1/20');
date3=datetime('18/1/21');
%% read the Pcuks information
% [pucks{1:12}]=textread('InputData.txt','%q %q %q %q %q %q %q %q %q %q %q %q');
[pucks{1:12}]=textread(filename_txt,'%q %q %q %q %q %q %q %q %q %q %q %q');
% load pucks.mat;
% [~,Pucks]=xlsread('InputData.xlsx','Pucks');
[~,Pucks]=xlsread(filename_xlsx,'Pucks');
head=Pucks(1,:);
arrive_date=datetime(Pucks(2:end,2));
depart_date=datetime(Pucks(2:end,7));
arrive_time=string(pucks{1,3});
depart_time=string(pucks{1,8});
num_of_planes=size(arrive_time,1);

puck_id=string(Pucks(2:end,1));
flight_id_arrive=string(Pucks(2:end,4));
flight_id_depart=string(Pucks(2:end,9));

for i=1:num_of_planes
    A=strsplit(arrive_time(i,1),':');
    B=strsplit(depart_time(i,1),':');
    hour1(i)=double(A(1,1));
    mint1(i)=double(A(1,2));
    hour2(i)=double(B(1,1));
    mint2(i)=double(B(1,2));
end
arrive_time=(hour1*60+mint1)';
depart_time=(hour2*60+mint2)';

[list1,~]=find((arrive_date==date1 & depart_date==date2));
[list2,~]=find((arrive_date==date2 & depart_date==date2));
[list3,~]=find((arrive_date==date2 & depart_date==date3));
% for list1 2 & 3
Atime1=ones(size(list1))*-45;
Dtime1=depart_time(list1);
Atime2=arrive_time(list2);
Dtime2=depart_time(list2);
Atime3=arrive_time(list3);
Dtime3=ones(size(list1))*(24*60+45);

Atime=[Atime1; Atime2; Atime3];
Dtime=[Dtime1; Dtime2; Dtime3];

% size of the plane narrow or wide, wide==1
size_of_plane=string(pucks{1,6});
wide_ref=string({"332" "333" "33E" "33H" "33L" "773"})';
narrow_ref=string({"319" "320" "321" "323" "325" "738" "73A" "73E" "73H" "73L"})';
size_of_plane=ismember(size_of_plane,wide_ref);
size_of_plane=[size_of_plane(list1); size_of_plane(list2); size_of_plane(list3)];

% arrive and depart type of the plane
arrive_type=string(pucks{1,5});
depart_type=string(pucks{1,10});
ad_type=strcat(arrive_type,depart_type); % II ID DI DD
AD_type=zeros(num_of_planes,4);
AD_type(:,1)=ad_type(:,1)=="II";
AD_type(:,2)=ad_type(:,1)=="ID";
AD_type(:,3)=ad_type(:,1)=="DI";
AD_type(:,4)=ad_type(:,1)=="DD";
AD_type=[AD_type(list1,:); AD_type(list2,:); AD_type(list3,:)];

% the id of the flight
puck_id=[puck_id(list1); puck_id(list2); puck_id(list3)];
flight_id_arrive=[flight_id_arrive(list1); flight_id_arrive(list2); flight_id_arrive(list3)];
flight_id_depart=[flight_id_depart(list1); flight_id_depart(list2); flight_id_depart(list3)];
% sort the data by arrive time
% [~,list]=sort(Atime);
% AD_type(:,:)=AD_type(list,:);
% size_of_plane=size_of_plane(list);
% Atime=Atime(list);
% Dtime=Dtime(list);

% create struct for plane
plane_struct=struct(...
    'A_type',double(AD_type), ...
    'u',double(size_of_plane), ...
    'A',double(Atime), ...
    'D',double(Dtime),...
    'flight_id_arrive',flight_id_arrive, ...
    'flight_id_depart',flight_id_depart, ...
    'puck_id', puck_id);

%% read the Passenger information
[~,Tickets]=xlsread(filename_xlsx,'Tickets');
head=Tickets(1,:);
num_of_tick=size(Tickets,1)-1;
num_of_pasen=load('num_of_pasen.txt');
pasen_id=string(Tickets(2:end,1));
f_id_arrive=string(Tickets(2:end,3));
f_id_depart=string(Tickets(2:end,5));

arrive_date=datetime(Tickets(2:end,4));
depart_date=datetime(Tickets(2:end,6));
% only consider the passengers arrive at date1
date1=datetime('18/1/20');
[list,~]=find((arrive_date==date1 & depart_date==date1));



passen_struct=struct(...
    'pasen_id', pasen_id(list), ...
    'pasen_num', num_of_pasen(list), ...
    'f_id_arrive', f_id_arrive(list), ...
    'f_id_depart', f_id_depart(list));

%% read the Gates information
% [~,Gates]=xlsread('InputData.xlsx','Gates');
[~,Gates]=xlsread(filename_xlsx,'Gates');
head=Gates(1,:);
num_of_gates=size(Gates,1)-1;
% size of the gate narrow or wide, wide==1
size_of_gate=string(Gates(2:end,6));
size_of_gate=ismember(size_of_gate,"W");

% arrive and depart type of the plane
arrive_type_gate=string(Gates(2:end,4));
depart_type_gate=string(Gates(2:end,5));
ad_type_gate=strcat(arrive_type_gate,depart_type_gate);
AD_type_gate=zeros(num_of_gates,4);
AD_type_gate(:,1)=ad_type_gate(:,1)=="II";
AD_type_gate(:,2)=ad_type_gate(:,1)=="ID";
AD_type_gate(:,3)=ad_type_gate(:,1)=="DI";
AD_type_gate(:,4)=ad_type_gate(:,1)=="DD";
for i=1:num_of_gates
    if(ad_type_gate(i)=="D, ID")
        AD_type_gate(i,:)=[0 1 0 1];
    else if(ad_type_gate(i)=="D, ID, I")
            AD_type_gate(i,:)=[1 1 1 1];
        else if(ad_type_gate(i)=="D, II")
                AD_type_gate(i,:)=[1 0 1 0];
            else if(ad_type_gate(i)=="ID, I")
                    AD_type_gate(i,:)=[1 1 0 0];
                else if(ad_type_gate(i)=="DD, I")
                        AD_type_gate(i,:)=[0 0 1 1];
                    end
                end
            end
        end
    end
end

% add gate type for TN TC TS SD SC SS SE
TS=string(Gates(2:end,2));
NS=string(Gates(2:end,3));
NS(NS=="North")="N";
NS(NS=="Center")="C";
NS(NS=="South")="S";
NS(NS=="East")="E";

TSNS=strcat(TS,NS);

% add the id of gate
gate_id=string(Gates(2:end,1));

% create struct for gates
gate_struct=struct(...
    'T_type',double(AD_type_gate), ...
    'v',double(size_of_gate), ...
    'TS', TS, ...
    'NS', NS, ...
    'TSNS', TSNS, ...
    'gate_id', gate_id);

%% save data
save('plane_struct.mat', 'plane_struct');
save('passen_struct.mat', 'passen_struct');
save('gate_struct.mat', 'gate_struct');
end


