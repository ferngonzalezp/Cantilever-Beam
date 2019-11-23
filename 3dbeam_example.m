nel = 1728;
E = 2.18820e9;
poisson = 0.2;
P = 1000;
t = 0.3;
ro = 2400;
g = 9.81;
l = 6;
h = 0.5;
dofs = 3;
P = 1000;

b1 = beam3d(E,t,poisson,ro,l,h);
C = b1.conect_hex(nel);
Nodes = C.Nodes;
Elem = C.Elem;

%cond de contorno
r = find(Nodes(:,1)==0 & (Nodes(:,2)==0.25));
desp = [r ones(size(r)) ones(size(r)) ones(size(r)) zeros(size(r)) zeros(size(r)) zeros(size(r))];
r = find(Nodes(:,1)==0 & Nodes(:,2)~=0.25);
desp = [desp; r ones(size(r)) zeros(size(r)) ones(size(r)) zeros(size(r)) zeros(size(r)) zeros(size(r))];

dof_list = b1.dof_list(dofs,desp,nel).dof_list;
total_dof = b1.dof_list(dofs,desp,nel).total_dof;
dof_free = b1.dof_list(dofs,desp,nel).dof_free;

%forces
fp.nodes = find(Nodes(:,1)==l & (Nodes(:,3)==0.5));
fp.fx = 0;
fp.fy = 0;
fp.fz = -P/length(fp.nodes);

K = b1.stiffness(nel,desp,dofs);
F = b1.forces(nel,desp,dofs,fp);

U = K\F;
U(end,1)

f = [1 2 4 3 1 5 7 3 4 8 7 5 6 8 4 2 6 5 1 2];
nodesz_free = find(dof_list(:,4)<=dof_free);
Nodes(nodesz_free,3) = Nodes(nodesz_free,3) + 10.*U(dof_list(nodesz_free,4));

for i = 1:nel
    v = [Nodes(Elem(i,:),1), Nodes(Elem(i,:),2), Nodes(Elem(i,:),3)];
    patch('Faces',f,'Vertices',v,'EdgeColor','r','FaceColor','w');
    hold on
end
view(40,10)

