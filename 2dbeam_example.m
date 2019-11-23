nel = 64;
tel = 3; % 1:linear, 2:lagrangian 2nd order, 3:serendipity 2nd order
E = 2.18820e9;
poisson = 0.2;
P = 1000;
t = 0.3;
ro = 2400;
g = 9.81;
l = 6;
h = 0.5;
state = 1; %1: plane stress, 2: plane strain
dofs = 2;

beam = beam_sqr(E,t,poisson,ro,state,l,h);
square = beam.conect_sqr(nel,tel);
Nodes = square.Nodes;
Elem = square.Elem;

%cond de cotorno
r = find(Nodes(:,1)==0 & (Nodes(:,2)==0.25));
desp = [r ones(size(r)) ones(size(r)) zeros(size(r)) zeros(size(r))];
r = find(Nodes(:,1)==0 & Nodes(:,2)~=0.25);
desp = [desp; r ones(size(r)) zeros(size(r)) zeros(size(r)) zeros(size(r))];

dof_list = beam.dof_list(dofs,desp,nel,tel).dof_list;
total_dof = beam.dof_list(dofs,desp,nel,tel).total_dof;

%forces
fp.node = find(Nodes(:,1)==l & (Nodes(:,2)==0.25));
fp.fx = 0;
fp.fy = -P;

K = beam.stiffness(nel,tel,desp,dofs);
F = beam.forces(nel,tel,desp,dofs,fp);
U = K\F;

U(end,1)
newNodes = Nodes;

 for i = 1:size(Nodes,1)
     if dof_list(i,2) < total_dof
         newNodes(i,1) = U(dof_list(i,2),1)*10+Nodes(i,1);
     end
     if dof_list(i,3) < total_dof
         newNodes(i,2) = U(dof_list(i,3),1)*10+Nodes(i,2);
     end
 end

dt = delaunayTriangulation(Nodes(:,1),Nodes(:,2));
triplot(dt)
hold on
triplot(dt.ConnectivityList,newNodes(:,1), newNodes(:,2),'r')
