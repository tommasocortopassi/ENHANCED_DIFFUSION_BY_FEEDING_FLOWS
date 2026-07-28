%Script reproducing Movie1
[times,particles]=Data_collect("Dataset1",0.3);

for i=1:1:176
name="p"+num2str(i);
particles.(name)(3,:)=-particles.(name)(3,:); %Fix the angle convention
end
 
p.p1=particles.p105;
%Treat the Stylonychia as a tracer 
% ([z;x;y] convention and rescale units of measurement
% to show the trail in the movie
p.p1=[ones(1,4889);p.p1(1,:)/200;p.p1(2,:)/200];

%Fit linearly missing NaN values
p=fit_trajectories(p,0);

%Treat the Stylonychia as a tracer to show the trail in the movie
v=movie_traj([],times,p,0.1,1,30);