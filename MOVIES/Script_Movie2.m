%Script reproducing Movie2
[times,particles]=Data_collect("Dataset1",0.3);

for i=1:1:176
name="p"+num2str(i);
particles.(name)(3,:)=-particles.(name)(3,:); %Fix the angle convention
end

t=simulate_tracer(times,particles,[0.75;50;30],30,0.25,0.215);
tracer.t1=t; 
clear t;

v=movie_traj(particles,times,tracer,0.5,1,30);