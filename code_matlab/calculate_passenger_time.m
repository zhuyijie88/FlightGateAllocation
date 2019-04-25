function [time tensity time_of_passen each_of_passen time_gap_plane]=calculate_passenger_time(x1d,plane_struct,passen_struct,gate_struct,num_of_gate,num_of_plane, num_of_p)
% INPUT: x_old is the solution
% INPUT: num_of_mut is the max number of mutation times
% INPUT: plane_struct stores the all information of planes
% INPUT: passen_struct stores the all information of passengers
% INPUT: gate_struct stores the all information of gates
% INPUT: num_of_gate is number of gates
% INPUT: num_of_plane is number of planes
% INPUT: num_of_p is the number of problem, == 2 or 3
% OUTPUT: time is the mean time for passenger transfering
% OUTPUT: tensity is the time tensity for passenger transfering
% OUTPUT: time_of_passen is the time of each passenger
% OUTPUT: each_of_passen is the number of passengers in each passenger rank 
% OUTPUT: time_gap_plane is the gap time between two flights for each passenger rank

num_of_passen=size(passen_struct.pasen_id,1);
type=[...
    "DTDT", "DTDS", "DTIT", "DTIS", ...
    "DSDT", "DSDS", "DSIT", "DSIS", ...
    "ITDT", "ITDS", "ITIT", "ITIS", ...
    "ISDT", "ISDS", "ISIT", "ISIS"];
time1=[ ...
    15, 20, 35, 40, ...
    20, 15, 40, 35, ...
    35, 40, 20, 30, ...
    40, 45, 30, 20];
num=[ ...
    0, 1, 0, 1, ...
    1, 0, 1, 0, ...
    0, 1, 0, 1, ...
    1, 2, 1, 0];

% store the data in table 1 of word
p2_dict=struct(...
    'type', type, ...
    'time', time1, ...
    'num', num ...
    );

type=[...
    "TNTN", "TNTC", "TNTS", "TNSN", "TNSC", "TNSS", "TNSE", ...
    "TCTN", "TCTC", "TCTS", "TCSN", "TCSC", "TCSS", "TCSE", ...
    "TSTN", "TSTC", "TSTS", "TSSN", "TSSC", "TSSS", "TSSE", ...
    "SNTN", "SNTC", "SNTS", "SNSN", "SNSC", "SNSS", "SNSE", ...
    "SCTN", "SCTC", "SCTS", "SCSN", "SCSC", "SCSS", "SCSE", ...
    "SSTN", "SSTC", "SSTS", "SSSN", "SSSC", "SSSS", "SSSE", ...
    "SETN", "SETC", "SETS", "SESN", "SESC", "SESS", "SESE", ...
    ];
time1=[...
    10, 15, 20, 25, 20, 25, 25, ...
    15, 10, 15, 20, 15, 20, 20, ...
    20, 15, 10, 25, 20, 25, 25, ...
    25, 20, 25, 10, 15, 20, 20, ...
    20, 15, 20, 15, 10, 15, 15, ...
    25, 20, 25, 20, 15, 10, 20, ...
    25, 20, 25, 20, 15, 20, 10
    ];
% store the data in table 2 of word
p3_dict=struct(...
    'type', type, ...
    'time', time1 ...
    );

%% calculate the second problem
% 1st, rank_of_passen ==>> flight_id_arrive
% 2nd, rank_of_passen ==>> flight_id_depart
% 3rd, flight_id_arrive ==>> gate_type_arrive (gate_struct.TS)
% 4th, flight_id_depart ==>> gate_type_depart (gate_struct.TS)
% 5th, flight_id ==> flight_property_I_or_D (International or Domestic)
% 5th, gate_type_arrive (gate_struct.TS)  }
%                                         } ==> the Flow Time (流程时间)
%      gate_type_depart (gate_struct.TS)  }
id1_date2=55;
id2_date2=249;

if(num_of_p==2)
    time_of_passen=[];
    each_of_passen=[];
    time_gap_plane=[];
    for i_pas=1:num_of_passen
        flight_id_arrive=find(passen_struct.f_id_arrive(i_pas)==plane_struct.flight_id_arrive(id1_date2:id2_date2));
        flight_id_depart=find(passen_struct.f_id_arrive(i_pas)==plane_struct.flight_id_arrive(id1_date2:id2_date2));
        if(x1d(flight_id_arrive)>0 & x1d(flight_id_arrive)<num_of_gate+1) & ...
                (x1d(flight_id_depart)>0 & x1d(flight_id_depart)<num_of_gate+1)
            % if true, the passenger will arrive the airport and leave the
            % airport successfully, calculate the time.
            type_DI=find(plane_struct.A_type(flight_id_arrive,:)==1);
            switch type_DI
                case 1
                    arrive_type_DI="I";
                    depart_type_DI="I";
                case 2
                    arrive_type_DI="I";
                    depart_type_DI="D";
                case 3
                    arrive_type_DI="D";
                    depart_type_DI="I";
                case 4
                    arrive_type_DI="D";
                    depart_type_DI="D";
            end
            
            arrive_gate_type_TS=gate_struct.TS(x1d(flight_id_arrive));
            depart_gate_type_TS=gate_struct.TS(x1d(flight_id_depart));
            
            tot_type=strcat(arrive_type_DI,arrive_gate_type_TS,depart_type_DI,depart_gate_type_TS);
            
            time_of_passen=[time_of_passen p2_dict.time(find(p2_dict.type==tot_type))];
            each_of_passen=[each_of_passen passen_struct.pasen_num(i_pas)];        % store each number of passengers for each guy in sheet2
            
            time_gap_plane=[time_gap_plane plane_struct.D(flight_id_depart)-plane_struct.A(flight_id_arrive)];
            
        end
        
    end
    time=sum(time_of_passen.*each_of_passen)/sum(each_of_passen);
    tensity=mean(time_of_passen./time_gap_plane);
end



%% calculate the third problem
% 1st, rank_of_passen ==>> flight_id_arrive
% 2nd, rank_of_passen ==>> flight_id_depart
% 3rd, flight_id_arrive ==>> gate_location_arrive (gate_struct.TSNS)
% 4th, flight_id_depart ==>> gate_location_depart (gate_struct.TSNS)
% 5th, gate_location ==> gate_property_I_or_D (International or Domestic)
% 5th, gate_location_arrive (gate_struct.TSNS)  }
%                                               } ==> the Flow Time (流程时间)
%      gate_location_depart (gate_struct.TSNS)  }
if(num_of_p==3) % only need to add the walk time
    time_of_passen=[];
    each_of_passen=[];
    time_gap_plane=[];
    for i_pas=1:num_of_passen
        flight_id_arrive=find(passen_struct.f_id_arrive(i_pas)==plane_struct.flight_id_arrive(id1_date2:id2_date2));
        flight_id_depart=find(passen_struct.f_id_arrive(i_pas)==plane_struct.flight_id_arrive(id1_date2:id2_date2));
        if(x1d(flight_id_arrive)>0 & x1d(flight_id_arrive)<num_of_gate+1) & ...
                (x1d(flight_id_depart)>0 & x1d(flight_id_depart)<num_of_gate+1)
            % if true, the passenger will arrive the airport and leave the
            % airport successfully, calculate the time.
            type_DI=find(plane_struct.A_type(flight_id_arrive,:)==1);
            switch type_DI
                case 1
                    arrive_type_DI="I";
                    depart_type_DI="I";
                case 2
                    arrive_type_DI="I";
                    depart_type_DI="D";
                case 3
                    arrive_type_DI="D";
                    depart_type_DI="I";
                case 4
                    arrive_type_DI="D";
                    depart_type_DI="D";
            end
            
            arrive_gate_type_TS=gate_struct.TS(x1d(flight_id_arrive));
            depart_gate_type_TS=gate_struct.TS(x1d(flight_id_depart));         
            tot_type=strcat(arrive_type_DI,arrive_gate_type_TS,depart_type_DI,depart_gate_type_TS);
            
            % calculate the gate type of TN TS TC SN SC SS SE
            arrive_gate_type_NS=gate_struct.NS(x1d(flight_id_arrive));
            depart_gate_type_NS=gate_struct.NS(x1d(flight_id_depart));
            tot_type_TNSC=strcat(arrive_gate_type_TS,arrive_gate_type_NS,depart_gate_type_TS,depart_gate_type_NS);
            
            if p2_dict.num(find(p2_dict.type==tot_type))==2 % walk time need to multiple 2
                walk_time=2*20;
                if tot_type_TNSC=="TCSC"
                    walk_time=2*15;
                end
            else                       
                walk_time=p3_dict.time(find(p3_dict.type==tot_type_TNSC));
            end      
            
            time_of_passen=[time_of_passen p2_dict.time(find(p2_dict.type==tot_type))+ ...
                p2_dict.num(find(p2_dict.type==tot_type))*8+ ...  % added the shuttle time, 8 min for each time
                walk_time]; % added the walk time
            each_of_passen=[each_of_passen passen_struct.pasen_num(i_pas)];        % store each number of passengers for each guy in sheet2
            
            time_gap_plane=[time_gap_plane plane_struct.D(flight_id_depart)-plane_struct.A(flight_id_arrive)];
            
        end
        
    end
    time=sum(time_of_passen.*each_of_passen)/sum(each_of_passen);
    tensity=mean(time_of_passen./time_gap_plane);    
    
end


end