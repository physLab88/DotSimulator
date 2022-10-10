%% Run Python SET Simulation and plot results with selected SET data
%Michael McConnell
%9/29/2014

%% Instructions
%loadSETData.m MUST be run first for this script to work properly.

%Python with SciPy must also be installed. I would recommend the Anaconda
%distribution found here: http://continuum.io/downloads

%This script asks for the SET parameters for the simulation, executes the
%simulation, and then plots the simulation with a selected column from
%Gds_data, which is created by loadSETData.m

%Please note, that the first time you run this script you must specify all
%SET parameters. However, on subsequent runs, since you are mostly likely
%going to be chaning only one parameter at a time, leaving any field blank
%will use the value from the previous simulation. This can save you a lot
%of time. For example, if you run your first simulation at 3K and discover
%you need to increase your temperature but want to keep all other
%parameters the same, change temperature to whatever value you want and
%just hit enter for all the other parameters.

%Maximum number of island electrons: this is the maximum number of
%electrons the simulation will consider. If the simulation results do not
%make sense, try increasing this number. Including more electrons in the
%simulation will make it run slower but is necessary at high Vds values and
%higher temperatures as more electrons have enough energy to tunnel onto
%the island under these conditions.

%The simulation is currently set up to only sweep Vds and leave Vg at 0mV,
%but it can sweep both. If you need that data, let me know.



G0=7.7480917346E-5; %Conductance quantum in Siemens;

temp=input('Temperature (K): ','s');
if ~isempty(temp)
    T=temp;
end
temp=input('Minimum Vds (mV): ','s');
if ~isempty(temp)
    vds_start=temp;
end
temp=input('Maximum Vds (mV): ','s');
if ~isempty(temp)
    vds_end=temp;
end
temp=input('Number of points in Vds array: ','s');
if ~isempty(temp)
    numVdspoints=temp;
end
temp=input('Minimum Vg (mV): ','s');
if ~isempty(temp)
    vg_start=temp;
end
temp=input('Maximum Vg (mV): ','s');
if ~isempty(temp)
    vg_end=temp;
end
temp=input('Number of points in Vg array: ','s');
if ~isempty(temp)
    numVgpoints=temp;
end
temp=input('Source-island capacitance (F): ','s');
if ~isempty(temp)
    Cs=temp;
end
temp=input('Drain-island capacitance (F): ','s');
if ~isempty(temp)
    Cd=temp;
end
temp=input('Gate-island capacitance (F): ','s');
if ~isempty(temp)
    Cg=temp;
end
temp=input('Source-island conductance (S): ');
if ~isempty(temp)
    Gs=temp;
    Gs=Gs/G0;
    Gs_string=num2str(Gs);
end
temp=input('Drain-island conductance (S): ');
if ~isempty(temp)
    Gd=temp;
    Gd=Gd/G0;
    Gd_string=num2str(Gd);
end
temp=input('Maximum number of island electrons: ','s');
if ~isempty(temp)
    num_e=temp;
end
% temp=input('Compare with which data set (index)? ');
% if ~isempty(temp)
%     index=temp;
% end

Vds_sim=linspace(str2num(vds_start),str2num(vds_end),str2num(numVdspoints));
%command=['python AllInOne_mcconnell.py ' T ' ' vds_start ' ' vds_end ' '...
%    numVdspoints ' ' Cs ' ' Cd ' ' Gs_string ' ' Gd_string ' ' num_e ' '...
%    vg_start ' ' vg_end ' ' numVgpoints ' ' Cg];
command=['C:\Python27_32\python.exe guidiamonds.py ' T ' ' vds_start ' ' vds_end ' '...
    numVdspoints ' ' Cs ' ' Cd ' ' Gs_string ' ' Gd_string ' ' num_e ' '...
    vg_start ' ' vg_end ' ' numVgpoints ' ' Cg];
[status,result]=system(command);
disp(['Simulation Output: ' result])
Gds_sim=load('simData.dat');

% figure(1)
% plot(Vds_sim(1:length(Vds_sim)-1), Gds_sim, Vds_data*1e3, Gds_data(:,index))
% legend('Simulation','Data','Location','SouthEast');
% xlabel('V_d_s (V)')
% ylabel('Differential Transconductance (S)')
% lims=ylim;
% ylim([0 lims(2)*1.1])



