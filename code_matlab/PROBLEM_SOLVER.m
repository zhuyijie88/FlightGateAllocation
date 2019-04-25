%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Annotation: the code generate a set of feasible solutions stored in
% a matrix x(num_of_samples,num_of_planes);
% num_of_samples is the number of solutions;
% num_of_planes is the number of planes;
% each element in x is the rank of the gate from 1~70,
% where 1~28 is the T gate, 29~69 is the S gate, 
% and the 70th ghost gate is used to place the unscheduled planes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
clc;
%% choose the problem you want to solve
% problem_choice=1, solve the problem 1;
% problem_choice=2, solve the problem 2;
% problem_choice=3, solve the problem 3;
problem_choice=1;
%% read data
% load data
[plane_struct gate_struct]=read_data_p1('InputData.txt', 'InputData.xlsx');
load gate_struct.mat;
load passen_struct.mat;
load plane_struct.mat;
% make the time in date=19 or 21 different with each other
[listdate19,~]=find(plane_struct.A<0);
plane_struct.A(listdate19)=0-listdate19*50;
[listdate21,~]=find(plane_struct.D>1440);
plane_struct.D(listdate21)=1440+listdate21*50;

% shift the time to be postive
shift_time=100-min(plane_struct.A);
plane_struct.A=plane_struct.A+shift_time;
plane_struct.D=plane_struct.D+shift_time;

%% input parameter
num_of_plane=size(plane_struct.A,1); % number of plane
num_of_gate=size(gate_struct.v,1); % number of gate
num_of_samples=20; % number of samples will be generated in each generation
lamda=0.8; % the weight parameter used in fitness(target) function
MAXGEN=2500; % the maximum generation that will be generated
NUM_PARENT=4; % the maximum number of parents will be chosed in each crossover
NUM_MUT=50; % the maximum attempt times in each mutation
MAX_ITER=floor(num_of_samples*10); % the maximum iteration times in each generation
myseed=sum(clock); % the seed for usage in random functions

%% program main
iter_times=0;
x=generate2(plane_struct,gate_struct,num_of_gate,num_of_plane,myseed,num_of_samples); % generate the first x samples
for i_gen=1:MAXGEN
    save(['data', num2str(i_gen)], '*'); % save the data in each generation
    clear x_sons;
    x_sons=[]; % to store the new solutions
    iter_times=0;
    [x_sorted, grade,~,~]=sort_by_target(x,plane_struct,passen_struct,gate_struct,num_of_gate,num_of_plane,lamda,problem_choice); % sort the samples from good to bad
    x=unique(x,'rows'); % delete the same rows in x
    
    while(size(x_sons,1)<num_of_samples & iter_times<MAX_ITER) % iteration until num_of_samples new solution are generated or the itertation times > MAX_ITER
        if(mod(iter_times,1)==0)
            disp([size(x_sons,1), iter_times, MAX_ITER, i_gen])
        end
        iter_times=iter_times+1;
        
        for i=1:NUM_PARENT
            % choose parent genes and generate new genes
            [id_male, id_female]=choose_parents(grade);
            
            % i==even, use method corssover to generate new solution
            % i==odd, use method crossover_replace to generate new solution
            if(mod(i,2)==0)
                [x1_son, x2_son]=crossover(x(id_male,:), x(id_female,:),plane_struct,gate_struct,num_of_gate,num_of_plane);
            else
                [x1_son, x2_son]=crossover_replace(x(id_male,:), x(id_female,:),plane_struct,gate_struct,num_of_gate,num_of_plane);
            end
            
            [~, okay]=check_1d(x1_son,plane_struct,gate_struct,num_of_gate,num_of_plane); % check the x1_son is feasible or not
            if  (okay==1)  % remove duplicate solution and keep the new one
                x_tot=[x; x1_son];
                dimen1=size(x_tot,1);
                dimen2=size(unique(x_tot,'rows'),1);
                if dimen1==dimen2
                    x_sons=[x_sons; x1_son];
                end

            else
                [x_new, success]=mutation(x1_son,NUM_MUT,plane_struct,gate_struct,num_of_gate,num_of_plane); % x1_son is infeasible, x_new is feasible if success==1 by mutation
                if(success==1) % remove duplicate solution and keep the new one
                    x_tot=[x; x_new];
                    dimen1=size(x_tot,1);
                    dimen2=size(unique(x_tot,'rows'),1);
                    if dimen1==dimen2
                        x_sons=[x_sons; x_new];
                    end
                    
                end
            end
            
            [~, okay]=check_1d(x2_son,plane_struct,gate_struct,num_of_gate,num_of_plane);
            if  (okay==1) % remove duplicate solution and keep the new one
                x_tot=[x; x2_son];
                dimen1=size(x_tot,1);
                dimen2=size(unique(x_tot,'rows'),1);
                if dimen1==dimen2
                    x_sons=[x_sons; x2_son];
                end
                
            else
                [x_new, success]=mutation(x2_son,NUM_MUT,plane_struct,gate_struct,num_of_gate,num_of_plane); % x2_son is infeasible, x_new is feasible if success==1 by mutation
                if(success==1) % remove duplicate solution and keep the new one
                    x_tot=[x; x_new];
                    dimen1=size(x_tot,1);
                    dimen2=size(unique(x_tot,'rows'),1);
                    if dimen1==dimen2
                        x_sons=[x_sons; x_new];
                    end
                    
                end
            end
        end
        
        
        % mutation the each gene ranked by iter_times in x by reducing the number of ghost gate directly
        if(iter_times<size(x,1))
            % mutation the iter_times ranked gene in x
            x_old=x(iter_times,:);
            if(randi([1,2])==1)
                [x_new, success]=mutation_ghost(x_old,NUM_MUT,plane_struct,gate_struct,num_of_gate,num_of_plane);
            else
                [x_new, success]=mutation_ghost_random(x_old,NUM_MUT,plane_struct,gate_struct,num_of_gate,num_of_plane);
            end
            
            if(success==1) % remove duplicate solution and keep the new one
                x_tot=[x; x_new];
                dimen1=size(x_tot,1);
                dimen2=size(unique(x_tot,'rows'),1);
                if dimen1==dimen2
                    x_sons=[x_sons; x_new];
                end
            end
        end
        
        

        % generate random x by order of random list of plane
        [x_new]=generate2_by_random(plane_struct,gate_struct,num_of_gate,num_of_plane,myseed,1);
        if(mod(iter_times,2)==0)
            % mutation the feasible sample.
            [x_new, success]=mutation_good(x_new,NUM_MUT,plane_struct,gate_struct,num_of_gate,num_of_plane);
        end
        x_tot=[x; x_new]; % remove duplicate solution and keep the new one
        dimen1=size(x_tot,1);
        dimen2=size(unique(x_tot,'rows'),1);
        if dimen1==dimen2
            x_sons=[x_sons; x_new];
        end
        x_sons=unique(x_sons,'row');
    end
    
    x_tot=[x; x_sons];
    [x_sorted, grade,~,~]=sort_by_target(x_tot,plane_struct,passen_struct,gate_struct,num_of_gate,num_of_plane,lamda,problem_choice); % sort the x from good to bad
    grade=(grade-min(grade))/(max(grade)-min(grade));
    
    % choose the new result by random delete the worse x
    while(size(x_sorted,1)>num_of_samples)
        grade_rand=rand();
        i_rand=randi([1,size(x_sorted,1)]);
        if(grade_rand>grade(i_rand))
            x_sorted(i_rand,:)=[];
        end
    end
    % update the x
    x=x_sorted(1:num_of_samples,:);
end



%% choose parents samples
% INPUT: grade is the fitness of the solution set x
% OUTPUT: id_male and id_female are the ranks of parents
function [id_male, id_female]=choose_parents(grade)
n=length(grade);
prob=grade/sum(grade);
weight(1)=prob(1);
num_id=0;
id_male=0;
id_female=0;
for i=2:n
    weight(i)=weight(i-1)+prob(i);
end
value=rand();
for i=1:n
    value=value-weight(i);
    if(value<0)
        id_male=i;
        id_female=id_male;
        break;
    end
    
end
while(id_female==id_male)
    value=rand();
    for i=1:n
        value=value-weight(i);
        if(value<0)
            id_female=i;
            break;
        end
    end
    
end
end

%% generate random solution set x (may infeasible)
% INPUT: plane_struct stores the all information of planes
% INPUT: gate_struct stores the all information of gates
% INPUT: num_of_plane is number of planes
% INPUT: gate_struct stores the all information of gates
% INPUT: myseed the seed for rand function
% INPUT: num_of_samples is the number of generated solutions
% OUTPUT: x is the generated solutions, x maybe not feasible
function [x]=generate(plane_struct,gate_struct,num_of_gate,num_of_plane,myseed,num_of_samples)
x_temp=zeros(num_of_samples,1);
x=zeros(num_of_samples,num_of_plane);  %generate n*m samples
for i=1:num_of_plane
    A=find(gate_struct.v==plane_struct.u(i));  % find the list id of gate-plane W/N type are the same
    B=find(gate_struct.T_type*plane_struct.A_type(i,:)'==1); % find the list id of gate-plane I/D type are the same
    
    C=intersect(A,B); % choose the shared elements in A & B
    C=[C; num_of_gate+1];
    
    rand('seed',sum(100*clock)+myseed*i);
    x_temp=C(randi([1,length(C)],num_of_samples,1));
    
    x(:,i)=x_temp;
end

end

%% generate random solution set x (feasible)
% INPUT: plane_struct stores the all information of planes
% INPUT: gate_struct stores the all information of gates
% INPUT: num_of_plane is number of planes
% INPUT: gate_struct stores the all information of gates
% INPUT: myseed the seed for rand function
% INPUT: num_of_samples is the number of generated solutions
% OUTPUT: x is the generated solutions, x are all feasible
function [x]=generate2(plane_struct,gate_struct,num_of_gate,num_of_plane,myseed,num_of_samples)
x_temp=zeros(num_of_samples,1);
x=zeros(num_of_samples,num_of_plane);  %generate n*m samples
for i=1:num_of_plane
    i_in_generate2=i
    A=find(gate_struct.v==plane_struct.u(i));  % find the list id of gate-plane W/N type are the same
    B=find(gate_struct.T_type*plane_struct.A_type(i,:)'==1); % find the list id of gate-plane I/D type are the same
    
    C=intersect(A,B); % choose the shared elements in A & B
    C=[C; num_of_gate+1];
    C_2d=repmat(C,1,num_of_samples);
    for j=1:i-1
        is_j_in_c=ismember(x(:,j),C);                       % air[j].gate.id is a member of air[i].gate.id==C
        for k=1:length(is_j_in_c)
            
            if(is_j_in_c(k)==1 && x(k,j)~=num_of_gate+1)  % if for kth choice for j has the same gate in C, delete this gate in C for this k column
                is_conflict=min([...
                    abs(plane_struct.D(i)-plane_struct.D(j)), ...
                    abs(plane_struct.D(i)-plane_struct.A(j)), ...
                    abs(plane_struct.A(i)-plane_struct.D(j)), ...
                    abs(plane_struct.A(i)-plane_struct.A(j)) ...
                    ])<45;  % less than 45 means conflict
                is_conflict=max(is_conflict, ...
                    (plane_struct.D(j)-plane_struct.A(i))* ...
                    (plane_struct.D(i)-plane_struct.A(j))>0);
                
                if(is_conflict==1)
                    C_2d(find(C_2d(:,k)==x(k,j)),k)=-1; % make the same gate as the ghost gate id
                end
                
            end
        end
        
    end
    rand('seed',sum(100*clock)+myseed*i);
    for k=1:num_of_samples
        C=C_2d(:,k);
        C(find(C==-1))=[];
        x_temp(k)=C(randi([1,length(C)],1));
    end
    
    x(:,i)=x_temp;
end

end


%% generate random x by random list of plane
% INPUT: plane_struct stores the all information of planes
% INPUT: gate_struct stores the all information of gates
% INPUT: num_of_plane is number of planes
% INPUT: gate_struct stores the all information of gates
% INPUT: myseed the seed for rand function
% INPUT: num_of_samples is the number of generated solutions
% OUTPUT: x is the generated solutions, x are all feasible
function [x]=generate2_by_random(plane_struct,gate_struct,num_of_gate,num_of_plane,myseed,num_of_samples)
x_temp=zeros(num_of_samples,1);
x=zeros(num_of_samples,num_of_plane);  %generate n*m samples

list_plane=1:num_of_plane;
list_plane=list_plane(randperm(numel(list_plane)));

for id1=1:num_of_plane
    i=list_plane(id1);  % choose the i by random list
    i_in_generate2=i;
    A=find(gate_struct.v==plane_struct.u(i));  % find the list id of gate-plane W/N type are the same
    B=find(gate_struct.T_type*plane_struct.A_type(i,:)'==1); % find the list id of gate-plane I/D type are the same
    
    C=intersect(A,B); % choose the shared elements in A & B
    C=[C; num_of_gate+1];
    C_2d=repmat(C,1,num_of_samples);
    for jd1=1:id1-1
        j=list_plane(jd1);  % choose the j by random list
        is_j_in_c=ismember(x(:,j),C);                       % air[j].gate.id is a member of air[i].gate.id==C
        for k=1:length(is_j_in_c)
            
            if(is_j_in_c(k)==1 && x(k,j)~=num_of_gate+1)  % if for kth choice for j has the same gate in C, delete this gate in C for this k column
                is_conflict=min([...
                    abs(plane_struct.D(i)-plane_struct.D(j)), ...
                    abs(plane_struct.D(i)-plane_struct.A(j)), ...
                    abs(plane_struct.A(i)-plane_struct.D(j)), ...
                    abs(plane_struct.A(i)-plane_struct.A(j)) ...
                    ])<45;  % less than 45 means conflict
                is_conflict=max(is_conflict, ...
                    (plane_struct.D(j)-plane_struct.A(i))* ...
                    (plane_struct.D(i)-plane_struct.A(j))>0);
                
                if(is_conflict==1)
                    C_2d(find(C_2d(:,k)==x(k,j)),k)=-1; % make the same gate as the ghost gate id
                end
                
            end
        end
        
    end
    rand('seed',sum(100*clock)+myseed*i);
    for k=1:num_of_samples
        C=C_2d(:,k);
        C(find(C==-1))=[];
        x_temp(k)=C(randi([1,length(C)],1));
    end
    
    x(:,i)=x_temp;
end

end

%% generate feasible gates set for index'th element :compare 1:index-1 and index+1:end sections in sub_x
% INPUT: sub_x is the sub solution that need to be improved
% INPUT: index is the id of bug gene in sub_x
% INPUT: plane_struct stores the all information of planes
% INPUT: gate_struct stores the all information of gates
% INPUT: gate_struct stores the all information of gates
% INPUT: num_of_plane is number of planes
% OUTPUT: C is the set of feasible gates
function [C]=generate_feasible_gates(sub_x,index,plane_struct,gate_struct,num_of_gate,num_of_plane)
sub_x_length=length(sub_x);
A=find(gate_struct.v==plane_struct.u(index));  % find the list id of gate-plane W/N type are the same
B=find(gate_struct.T_type*plane_struct.A_type(index,:)'==1); % find the list id of gate-plane I/D type are the same

C=intersect(A,B); % choose the shared elements in A & B
for i=1:sub_x_length
    if(i==index)
        continue;
    end
    is_i_in_c=ismember(sub_x(i),C);
    if(is_i_in_c==1)
        
        is_conflict=min([...
            abs(plane_struct.D(i)-plane_struct.D(index)), ...
            abs(plane_struct.D(i)-plane_struct.A(index)), ...
            abs(plane_struct.A(i)-plane_struct.D(index)), ...
            abs(plane_struct.A(i)-plane_struct.A(index)) ...
            ])<45;  % less than 45 means conflict
        is_conflict=max(is_conflict, ...
            (plane_struct.D(index)-plane_struct.A(i))* ...
            (plane_struct.D(i)-plane_struct.A(index))>0);
        
        if(is_conflict==1)
            C(find(C==sub_x(i)))=[];
        end
        
    end
end
C=[C; num_of_gate+1];
end

%% check time information for the x solution, return the feasible set of solutions of x
% INPUT: x is the solutions set that need to be checked
% INPUT: plane_struct stores the all information of planes
% INPUT: gate_struct stores the all information of gates
% INPUT: num_of_gate is number of gates
% INPUT: num_of_plane is number of planes
% OUTPUT: x_okay is the feasible solutions set
function [x_okay]=check(x,plane_struct,gate_struct,num_of_gate,num_of_plane)
[m,n]=size(x);
x_okay=[];

if isempty(find(x<1|x>num_of_gate+1))~=1
    pause;
end

for ii=1:m  % each sample
    x_temp=x(ii,:);
    judge=1;
    for i=1:n
        for j=i+1:n
            if(x_temp(i)==x_temp(j) && x_temp(i)<num_of_gate+1)
                judge=min([...
                    abs(plane_struct.D(i)-plane_struct.D(j)), ...
                    abs(plane_struct.D(i)-plane_struct.A(j)), ...
                    abs(plane_struct.A(i)-plane_struct.D(j)), ...
                    abs(plane_struct.A(i)-plane_struct.A(j)) ...
                    ])>=45;
                judge=min(judge, ...
                    (plane_struct.D(j)-plane_struct.A(i))* ...
                    (plane_struct.D(i)-plane_struct.A(j))<=0);
            end
            if(judge==0)
                break;
            end
        end
        if(judge==0)
            break;
            %             disp("x warning, find x time constraint not satisfied !!!!!!");
        end
    end
    if(judge==1)
        x_okay=[x_okay; x_temp];
    end
end
end

%% check time information for the x_1d solution
% INPUT: x_1d is the solution that need to be checked
% INPUT: plane_struct stores the all information of planes
% INPUT: gate_struct stores the all information of gates
% INPUT: num_of_gate is number of gates
% INPUT: num_of_plane is number of planes
% OUTPUT: x1d_okay is the feasible solutions
% OUTPUT: okay is the check result
function [x1d_okay, okay]=check_1d(x_1d,plane_struct,gate_struct,num_of_gate,num_of_plane)
n=length(x_1d);
x1d_okay=[];
okay=0;

if isempty(find(x_1d<1|x_1d>num_of_gate+1))~=1
    pause;
end
x_temp=x_1d;
judge=1;
for i=1:n
    for j=i+1:n
        if(x_temp(i)==x_temp(j) && x_temp(i)<num_of_gate+1)
            judge=min([...
                abs(plane_struct.D(i)-plane_struct.D(j)), ...
                abs(plane_struct.D(i)-plane_struct.A(j)), ...
                abs(plane_struct.A(i)-plane_struct.D(j)), ...
                abs(plane_struct.A(i)-plane_struct.A(j)) ...
                ])>=45;
            judge=min(judge, ...
                (plane_struct.D(j)-plane_struct.A(i))* ...
                (plane_struct.D(i)-plane_struct.A(j))<=0);
        end
        if(judge==0)
            break;
        end
    end
    if(judge==0)
        break;
    end
end
if(judge==1)
    x1d_okay=x_temp;
    okay=1;
end
end

%% check time information for the x_1d solution
% INPUT: x_1d is the solution that need to be checked
% INPUT: plane_struct stores the all information of planes
% INPUT: gate_struct stores the all information of gates
% INPUT: num_of_gate is number of gates
% INPUT: num_of_plane is number of planes
% OUTPUT: x1d_okay is the feasible solutions
% OUTPUT: okay is the check result
function [x1d_okay, okay]=check_1d_last(x_1d,plane_struct,gate_struct,num_of_gate,num_of_plane)
n=length(x_1d);
x1d_okay=[];
okay=0;

if isempty(find(x_1d<1|x_1d>num_of_gate+1))~=1
    pause;
end
x_temp=x_1d;
judge=1;
j=n;
for i=1:n-1
    if(x_temp(i)==x_temp(j) && x_temp(i)<num_of_gate+1)
        judge=min([...
            abs(plane_struct.D(i)-plane_struct.D(j)), ...
            abs(plane_struct.D(i)-plane_struct.A(j)), ...
            abs(plane_struct.A(i)-plane_struct.D(j)), ...
            abs(plane_struct.A(i)-plane_struct.A(j)) ...
            ])>=45;
        judge=min(judge, ...
            (plane_struct.D(j)-plane_struct.A(i))* ...
            (plane_struct.D(i)-plane_struct.A(j))<=0);
    end
    if(judge==0)
        return;
    end

end
if(judge==1)
    x1d_okay=x_temp;
    okay=1;
end
end

%% check time information for the x_1d solution
% INPUT: x_1d is the solution that need to be checked
% INPUT: plane_struct stores the all information of planes
% INPUT: gate_struct stores the all information of gates
% INPUT: num_of_gate is number of gates
% INPUT: num_of_plane is number of planes
% OUTPUT: x1d_okay is the feasible solutions
% OUTPUT: okay is the check result
function [x1d_okay, okay]=check_1d_first(x_1d,plane_struct,gate_struct,num_of_gate,num_of_plane)
n=length(x_1d);
x1d_okay=[];
okay=0;

if isempty(find(x_1d<1|x_1d>num_of_gate+1))~=1
    disp("x_1d error, find x_1d<0 or x_1d>num_of_gate+1 !!!!!!");
    pause;
end
x_temp=x_1d;
judge=1;
j=1;
for i=2:n
    %    for j=i+1:n
    if(x_temp(i)==x_temp(j) && x_temp(i)<num_of_gate+1)
        judge=min([...
            abs(plane_struct.D(i)-plane_struct.D(j)), ...
            abs(plane_struct.D(i)-plane_struct.A(j)), ...
            abs(plane_struct.A(i)-plane_struct.D(j)), ...
            abs(plane_struct.A(i)-plane_struct.A(j)) ...
            ])>=45;
        judge=min(judge, ...
            (plane_struct.D(j)-plane_struct.A(i))* ...
            (plane_struct.D(i)-plane_struct.A(j))<=0);
    end
    if(judge==0)
        return;
    end
end
if(judge==1)
    x1d_okay=x_temp;
    okay=1;
end
end

%% sort solution x by target
% INPUT: x_1d is the solution that need to be checked
% INPUT: plane_struct stores the all information of planes
% INPUT: passen_struct stores the all information of passengers
% INPUT: gate_struct stores the all information of gates
% INPUT: num_of_gate is number of gates
% INPUT: num_of_plane is number of planes
% INPUT: lamda is the weight parameter in fitness fuunctions
% INPUT: problem_choice is the choice of problem need to be solved
% OUTPUT: x_sorted the sorted solutions set
% OUTPUT: grade is the value of fitness, larger is better
% OUTPUT: num_used_planes is the number of used planes
function [x_sorted, grade,num_open_gates,num_used_planes]=sort_by_target(x,plane_struct,passen_struct,gate_struct,num_of_gate,num_of_plane,lamda,problem_choice)
[m,n]=size(x);  % m is the num of samples, n is the num of planes
grade=zeros(m,1); % to store the values of target functions

% calculate the grade for each sample
for i=1:m
    num_open_gates=numel(unique(x(i,:)));  % calculate the num of open gates
    if(find(x(i,:)==num_of_gate+1)) % it has the ghost gate used
        num_open_gates=num_open_gates-1;
    end
    
    % default solve the problem 1
    num_used_planes=num_of_plane-numel(find(x(i,:)==num_of_gate+1));  % calculate the num of used planes
    grade(i,1)=lamda*(1-num_open_gates/num_of_gate)+(1-lamda)*num_used_planes/num_of_plane; % the grade less is better
    
    %solve the problem 2
    if(problem_choice==2)
    [time tensity time_of_passen each_of_passen time_gap_plane]=calculate_passenger_time(x(i,:),plane_struct,passen_struct,gate_struct,num_of_gate,num_of_plane, problem_choice);
    grade(i,1)=grade(i,1)*0.2+0.8*(1-time/45);
    end
    
    %solve the problem 3
    if(problem_choice==3)
    [time tensity time_of_passen each_of_passen time_gap_plane]=calculate_passenger_time(x(i,:),plane_struct,passen_struct,gate_struct,num_of_gate,num_of_plane, problem_choice);
    grade(i,1)=grade(i,1)*0.2+0.8*(1-tensity);
    end
    
end


% sort the x from good to bad by grade, return the sorted x_sorted
[grade index]=sort(grade,'descend');
x_sorted=x(index,:);

% print the best result
num_open_gates=numel(unique(x_sorted(1,:)));  % calculate the num of open gates
if(find(x_sorted(1,:)==num_of_gate+1)) % it has the ghost gate used
    num_open_gates=num_open_gates-1;
end
num_used_planes=num_of_plane-numel(find(x_sorted(1,:)==num_of_gate+1));  % calculate the num of used planes
disp("num_open_gates num_used_planes=");
disp([num_open_gates, num_used_planes]);
end

%% intersect two family samples x1_old+x2_old=x1_new+x2_new
% INPUT: x1_old & x2_old are parents
% INPUT: plane_struct stores the all information of planes
% INPUT: gate_struct stores the all information of gates
% INPUT: num_of_gate is number of gates
% INPUT: num_of_plane is number of planes
% OUTpUT: x1_new & x2_new are the sons

% old 1: 11111111111111111111111111111111111111111111111111
% old 2: 22222222222222222222222222222222222222222222222222
% cut_point:                               ^
% new 1: 11111111111111111111111111111111112222222222222222
% new 2: 22222222222222222222222222222222221111111111111111
function [x1_new, x2_new]=crossover(x1_old, x2_old,plane_struct,gate_struct,num_of_gate,num_of_plane)
m=length(x1_old);
cut_point=randi([1,m-1]);
x1_new=[x1_old(1:cut_point) x2_old(cut_point+1:end)];
x2_new=[x2_old(1:cut_point) x1_old(cut_point+1:end)];
end

%% replace one section between two family samples x1_old+x2_old=x1_new+x2_new
% INPUT: x1_old & x2_old are parents
% INPUT: plane_struct stores the all information of planes
% INPUT: gate_struct stores the all information of gates
% INPUT: num_of_gate is number of gates
% INPUT: num_of_plane is number of planes
% OUTpUT: x1_new & x2_new are the sons

% old 1: 11111111111111111111111111111111111111111111111111
% old 2: 22222222222222222222222222222222222222222222222222
% cut_point:               ^                ^
% new 1: 11111111111111111112222222222222222211111111111111
% new 2: 22222222222222222221111111111111111122222222222222
function [x1_new, x2_new]=crossover_replace(x1_old, x2_old,plane_struct,gate_struct,num_of_gate,num_of_plane)
m=length(x1_old);
cut_point1=randi([2,m-2]);
cut_point2=randi([cut_point1+1, m-1]);
x1_new=[x1_old(1:cut_point1) x2_old(cut_point1+1:cut_point2) x1_old(cut_point2+1:end)];
x2_new=[x2_old(1:cut_point1) x1_old(cut_point1+1:cut_point2) x2_old(cut_point2+1:end)];
end

%% genovariation a infeasible gene into a feasible gene, try most k times
% INPUT: x_old is parent
% INPUT: num_of_mut is the max number of mutation times
% INPUT: plane_struct stores the all information of planes
% INPUT: gate_struct stores the all information of gates
% INPUT: num_of_gate is number of gates
% INPUT: num_of_plane is number of planes
% OUTPUT: x_new is the son
% OUTPUT: success means the mutation success or not

% old: 1111111111111111111111E11111111111111111E11111111111
% bug:                       ^                 ^
% new: 1111111111111111111111111111111111111111111111111111
% input x_old is 1D, num_of_var is the maximum times of the genovariation
function [x_new, success]=mutation(x_old,num_of_mut,plane_struct,gate_struct,num_of_gate,num_of_plane) % x_old is infeasible, x_new is feasible if success==1
irank=0;
success=0;
% [~, okay]=check_1d(x_old,plane_struct,gate_struct,num_of_gate,num_of_plane);
n=length(x_old);
for i=1:n
    [~, okay_last]=check_1d_last(x_old(1:i),plane_struct,gate_struct,num_of_gate,num_of_plane);
    [~, okay_first]=check_1d_first(x_old(i:end),plane_struct,gate_struct,num_of_gate,num_of_plane);
    
    if((okay_last+okay_first)~=2 & irank<num_of_mut) % mean the ith plane is conflict with the previous section

        C=generate_feasible_gates(x_old,i,plane_struct,gate_struct,num_of_gate,num_of_plane);
        C=[C; num_of_gate+1];
        

        rand('seed',sum(100*clock)*i);
        x_old(i)=C(randi([1,length(C)],1));
        irank=irank+1;
    end
end

if(randi([1,2])==1)
    [x_old, ~]=mutation_ghost(x_old,num_of_mut,plane_struct,gate_struct,num_of_gate,num_of_plane);
else
    [x_old, ~]=mutation_ghost_random(x_old,num_of_mut,plane_struct,gate_struct,num_of_gate,num_of_plane);
end

[x_new, success]=check_1d(x_old,plane_struct,gate_struct,num_of_gate,num_of_plane);

num_70=numel(find(x_old==num_of_gate+1));
if(num_70>80)
    success=0;
end
end

%% genovariation a infeasible gene into a feasible gene, try most k times
% INPUT: x_old is parent
% INPUT: num_of_mut is the max number of mutation times
% INPUT: plane_struct stores the all information of planes
% INPUT: gate_struct stores the all information of gates
% INPUT: num_of_gate is number of gates
% INPUT: num_of_plane is number of planes
% OUTPUT: x_new is the son
% OUTPUT: success means the mutation success or not

% old: 1111111111111111111111E11111111111111111E11111111111
% bug:                       ^                 ^
% new: 1111111111111111111111111111111111111111111111111111
% input x_old is 1D, num_of_var is the maximum times of the genovariation
function [x_new, success]=mutation_good(x_old,num_of_mut,plane_struct,gate_struct,num_of_gate,num_of_plane) % x_old is feasible, x_new is feasible if success==1
irank=0;
success=0;
% [~, okay]=check_1d(x_old,plane_struct,gate_struct,num_of_gate,num_of_plane);
n=length(x_old);
rand_num=randi([1,n],1);
rand_list=randi([1,n],1,rand_num);
rand_list=unique(rand_list);
for i=rand_list
    [~, okay_last]=check_1d_last(x_old(1:i),plane_struct,gate_struct,num_of_gate,num_of_plane);
    [~, okay_first]=check_1d_first(x_old(i:end),plane_struct,gate_struct,num_of_gate,num_of_plane);
    
    if((okay_last+okay_first)==2 & irank<num_of_mut) % mean the ith plane is conflict with the previous section

        C=generate_feasible_gates(x_old,i,plane_struct,gate_struct,num_of_gate,num_of_plane);
        C=[C; num_of_gate+1];
        
        rand('seed',sum(100*clock)*i);
        x_old(i)=C(randi([1,length(C)],1));
        irank=irank+1;
    end
end

if(randi([1,2])==1)
    [x_old, ~]=mutation_ghost(x_old,num_of_mut,plane_struct,gate_struct,num_of_gate,num_of_plane);
else
    [x_old, ~]=mutation_ghost_random(x_old,num_of_mut,plane_struct,gate_struct,num_of_gate,num_of_plane);
end

[x_new, success]=check_1d(x_old,plane_struct,gate_struct,num_of_gate,num_of_plane);

num_70=numel(find(x_old==num_of_gate+1));
if(num_70>80)
    success=0;
end
end



%% genovariation a infeasible gene into a feasible gene, try most k times
% INPUT: x_old is parent
% INPUT: num_of_mut is the max number of mutation times
% INPUT: plane_struct stores the all information of planes
% INPUT: gate_struct stores the all information of gates
% INPUT: num_of_gate is number of gates
% INPUT: num_of_plane is number of planes
% OUTPUT: x_new is the son
% OUTPUT: success means the mutation success or not

% old: 111111111111111111111170111111111111111170E11111111111
% bug:                       ^                 ^
% new: 111111111111111111111111111111111111111111111111111111
% input x_old is 1D, num_of_var is the maximum times of the genovariation
function [x_new, success]=mutation_ghost(x_old,num_of_mut,plane_struct,gate_struct,num_of_gate,num_of_plane) % mutation the plane for that in ghost gates
irank=0;
success=0;
% [~, okay]=check_1d(x_old,plane_struct,gate_struct,num_of_gate,num_of_plane);
n=length(x_old);

for i=1:n

    if(x_old(i)==(num_of_gate+1) & irank<num_of_mut) % mean the ith plane is conflict with the previous section
        C=generate_feasible_gates(x_old,i,plane_struct,gate_struct,num_of_gate,num_of_plane);
        
        if(isempty(C)==1)
            irank=num_of_mut; % if C=[]; break out
            return;
        else
            rand('seed',sum(100*clock)*i);
            x_old(i)=C(randi([1,length(C)],1));
            irank=irank+1;
        end
    end
end

[x_new, success]=check_1d(x_old,plane_struct,gate_struct,num_of_gate,num_of_plane);
end



%% genovariation a infeasible gene into a feasible gene, try most k times
% INPUT: x_old is parent
% INPUT: num_of_mut is the max number of mutation times
% INPUT: plane_struct stores the all information of planes
% INPUT: gate_struct stores the all information of gates
% INPUT: num_of_gate is number of gates
% INPUT: num_of_plane is number of planes
% OUTPUT: x_new is the son
% OUTPUT: success means the mutation success or not

% old: 111111111111111111111170111111111111111170E11111111111
% bug:                       ^                 ^
% new: 111111111111111111111111111111111111111111111111111111
% input x_old is 1D, num_of_var is the maximum times of the genovariation
function [x_new, success]=mutation_ghost_random(x_old,num_of_mut,plane_struct,gate_struct,num_of_gate,num_of_plane) % mutation the plane for that in ghost gates by random
irank=0;
success=0;
n=length(x_old);

list=find(x_old==num_of_gate+1); % find the list of x(i)=70
% for i=1:n
for i=list
    
    if(x_old(i)==(num_of_gate+1) & irank<num_of_mut) % mean the ith plane is conflict with the previous section
        C=generate_feasible_gates(x_old,i,plane_struct,gate_struct,num_of_gate,num_of_plane);
        
        if(isempty(C)==1)
            irank=num_of_mut; % if C=[]; break out
            return;
        else
            rand('seed',sum(100*clock)*i);
            x_old(i)=C(randi([1,length(C)],1));
            %             [~, okay]=check_1d(x_old(1:i),plane_struct,gate_struct,num_of_gate,num_of_plane);
            irank=irank+1;
        end
    end
end

[x_new, success]=check_1d(x_old,plane_struct,gate_struct,num_of_gate,num_of_plane);
end
