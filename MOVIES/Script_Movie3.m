%Script reproducing Movie3
[times,particles]=Data_collect("Dataset1",1);

for i=1:1:176
name="p"+num2str(i);
particles.(name)(3,:)=-particles.(name)(3,:); %Fix the angle convention
end

t=simulate_tracer(times,particles,[0.75;50;35],30,0.25,0.215);
tracer.t1=t; 
clear t;

v=movie_traj(particles,times,tracer,0.5,1,30);
