classdef beam3d
    properties
        t
        E
        poisson
        ro
        g
        l
        h
    end
    methods
        function V = conect_hex(obj,nel)
            if mod(nel,nel^(1/3))==0
                nx = int64(nel^(1/3));
                ny = nx;
                nz = nx;
            elseif mod(nel,2)==0
                div = divisors(nel);
                d1 = div(floor(end/2));
                div = divisors(nel/d1);
                d2 = div(floor(end/2));
                d3 = nel/(d1*d2);
                nx = max([d1 d2 d3]);
                ny = min([d1 d2 d3]);
                nz = nel/(nx*ny);
            end
            x = linspace(0,obj.l,nx+1);
            y = linspace(0,obj.t,ny+1);
            z = linspace(0,obj.h,nz+1);
            [X,Y,Z]=meshgrid(x,y,z);
            X = X(:);
            Y = Y(:);
            Z = Z(:);
            Nodes = [X Y Z];
            ni = 1;
            n1 = ni;
            k = 1;
            step = int64((min([nx ny nz])));
            for i = 1:nel
                n2 = n1+1;
                n3 = n1 + step + 1;
                n4 = n3 +1;
                n5 = (n1 + (nx+1)*(ny+1));
                n6 = n5 + 1;
                n7 = n5 + step + 1;
                n8 = n7 + 1;
                V.Elem(i,:) = [n1 n2 n3 n4 n5 n6 n7 n8];
                if Nodes((n3),1) < obj.l
                    n1 = n3;
                else
                    if Nodes((n5),3) < obj.h
                        n1 = k + (nx+1)*(ny+1);
                        k = n1;
                    else
                        n1 = ni+1;
                        ni = n1;
                        k = n1;
                    end
                end
            end
            V.Nodes = Nodes;
        end
        function V = dof_list(obj,dofs,desp,nel)
            Nodes = obj.conect_hex(nel).Nodes;
            dof_fixed = size(Nodes,1)*dofs*10;
            dof_free = 1;
            V.total_dof = size(Nodes,1)*dofs*10;
            V.dof_list = zeros(size(Nodes,1),dofs+1);
            V.dof_list(:,1)=(1:1:size(Nodes,1));
            for i=1:size(desp,1)
            if desp(i,2) == 1
                a = dof_fixed;
                dof_fixed = dof_fixed + 1;
            else
                a = 0;
            end
            if desp(i,3) == 1
                b = dof_fixed;
                dof_fixed = dof_fixed + 1;
            else
                b = 0;
            end
            if desp(i,4) == 1
                c = dof_fixed;
                dof_fixed = dof_fixed + 1;
            else
                c = 0;
            end
            V.dof_list(desp(i,1),:) = [desp(i,1), a, b, c];
            end
        for i = 1:size(V.dof_list,1)

                if V.dof_list(i,2) == 0
                    a = dof_free;
                    dof_free = dof_free+1;
                    V.dof_list(i,2) = a;
                end
                if V.dof_list(i,3) == 0
                    b = dof_free;
                    dof_free = dof_free+1;
                    V.dof_list(i,3) =b;
                end
                if V.dof_list(i,4) == 0
                    c = dof_free;
                    dof_free = dof_free+1;
                    V.dof_list(i,4) =c;
                end
        end
        V.dof_free = dof_free;
        end
        function D = material(obj)
            D = eye(6);
            D(4,4) = (1-2*obj.poisson)/(2*(1-obj.poisson));
            D(5,5) = D(4,4);
            D(6,6) = D(4,4);
            D(1,2) = obj.poisson/(1-obj.poisson);
            D(1,3) = D(1,2);
            D(2,1) = D(1,2);
            D(2,3) = D(1,2);
            D(3,1) = D(1,2);
            D(3,2) = D(1,2);
            
            D = (obj.E*(1-obj.poisson))/((1+obj.poisson)*(1-2*obj.poisson)).*D;
        end
        function v = shapefun(obj,x,y,z)
            v.f = zeros(1,8);
            v.df = zeros(3,8);
            xi = -1; yi = -1; zi = -1;
            v.f(1) = 1/8*(1+xi*x)*(1+yi*y)*(1+zi*z);
            v.df(1,1) = 1/8*(xi)*(1+yi*y)*(1+zi*z);
            v.df(2,1) = 1/8*(1+xi*x)*(yi)*(1+zi*z);
            v.df(3,1) = 1/8*(1+xi*x)*(1+yi*y)*(zi);
            xi = -1; yi = 1; zi = -1;
            v.f(2) = 1/8*(1+xi*x)*(1+yi*y)*(1+zi*z);
            v.df(1,2) = 1/8*(xi)*(1+yi*y)*(1+zi*z);
            v.df(2,2) = 1/8*(1+xi*x)*(yi)*(1+zi*z);
            v.df(3,2) = 1/8*(1+xi*x)*(1+yi*y)*(zi);
            xi = 1; yi = -1; zi = -1;
            v.f(3) = 1/8*(1+xi*x)*(1+yi*y)*(1+zi*z);
            v.df(1,3) = 1/8*(xi)*(1+yi*y)*(1+zi*z);
            v.df(2,3) = 1/8*(1+xi*x)*(yi)*(1+zi*z);
            v.df(3,3) = 1/8*(1+xi*x)*(1+yi*y)*(zi);
            xi = 1; yi = 1; zi = -1;
            v.f(4) = 1/8*(1+xi*x)*(1+yi*y)*(1+zi*z);
            v.df(1,4) = 1/8*(xi)*(1+yi*y)*(1+zi*z);
            v.df(2,4) = 1/8*(1+xi*x)*(yi)*(1+zi*z);
            v.df(3,4) = 1/8*(1+xi*x)*(1+yi*y)*(zi);
            xi = -1; yi = -1; zi = 1;
            v.f(5) = 1/8*(1+xi*x)*(1+yi*y)*(1+zi*z);
            v.df(1,5) = 1/8*(xi)*(1+yi*y)*(1+zi*z);
            v.df(2,5) = 1/8*(1+xi*x)*(yi)*(1+zi*z);
            v.df(3,5) = 1/8*(1+xi*x)*(1+yi*y)*(zi);
            xi = -1; yi = 1; zi = 1;
            v.f(6) = 1/8*(1+xi*x)*(1+yi*y)*(1+zi*z);
            v.df(1,6) = 1/8*(xi)*(1+yi*y)*(1+zi*z);
            v.df(2,6) = 1/8*(1+xi*x)*(yi)*(1+zi*z);
            v.df(3,6) = 1/8*(1+xi*x)*(1+yi*y)*(zi);
            xi = 1; yi = -1; zi = 1;
            v.f(7) = 1/8*(1+xi*x)*(1+yi*y)*(1+zi*z);
            v.df(1,7) = 1/8*(xi)*(1+yi*y)*(1+zi*z);
            v.df(2,7) = 1/8*(1+xi*x)*(yi)*(1+zi*z);
            v.df(3,7) = 1/8*(1+xi*x)*(1+yi*y)*(zi);
            xi = 1; yi = 1; zi = 1;
            v.f(8) = 1/8*(1+xi*x)*(1+yi*y)*(1+zi*z);
            v.df(1,8) = 1/8*(xi)*(1+yi*y)*(1+zi*z);
            v.df(2,8) = 1/8*(1+xi*x)*(yi)*(1+zi*z);
            v.df(3,8) = 1/8*(1+xi*x)*(1+yi*y)*(zi);
        end
        function v = strain(obj,x,y,z,i,nel)
            C = obj.conect_hex(nel);
            Nodes = C.Nodes;
            Elem = C.Elem;
            N = obj.shapefun(x,y,z);
            Je = zeros(3);
            for j = 1:8
                Je = Je + [N.df(1,j)*Nodes(Elem(i,j),1) N.df(1,j)*Nodes(Elem(i,j),2) N.df(1,j)*Nodes(Elem(i,j),3); N.df(2,j)*Nodes(Elem(i,j),1) N.df(2,j)*Nodes(Elem(i,j),2) N.df(2,j)*Nodes(Elem(i,j),3); N.df(3,j)*Nodes(Elem(i,j),1) N.df(3,j)*Nodes(Elem(i,j),2) N.df(3,j)*Nodes(Elem(i,j),3)];
            end
            Be = [];
            for j = 1:8
                Jei = inv(Je);
                bcd = N.df(1,j).*Jei(:,1);
                bcd = bcd + N.df(2,j).*Jei(:,2);
                bcd = bcd + N.df(3,j).*Jei(:,3);
                b = bcd(1); c = bcd(2); d = bcd(3);
                Be = [ Be [b 0 0; 0 c 0; 0 0 d; c b 0; d 0 b; 0 d c]];
            end
            v.Be = Be;
            v.Je = Je;
        end
        function v = integrate(obj,f)
            v = f(-sqrt(1/3),-sqrt(1/3),-sqrt(1/3))+f(-sqrt(1/3),sqrt(1/3),-sqrt(1/3))+f(sqrt(1/3),-sqrt(1/3),-sqrt(1/3))+f(sqrt(1/3),sqrt(1/3),-sqrt(1/3));
            v = v + f(-sqrt(1/3),-sqrt(1/3),sqrt(1/3))+f(-sqrt(1/3),sqrt(1/3),sqrt(1/3))+f(sqrt(1/3),-sqrt(1/3),sqrt(1/3))+f(sqrt(1/3),sqrt(1/3),sqrt(1/3));
        end
        function v = stiffness(obj,nel,desp,dofs)
            dof = obj.dof_list(dofs,desp,nel);
            dof_free = dof.dof_free;
            D = obj.material();
            K = sparse(dof_free-1, dof_free-1);
            Elem = obj.conect_hex(nel).Elem;
            for i = 1:size(Elem,1)
                conect = index(Elem(i,:),3,dof.dof_list);
                f = @(x,y,z) obj.strain(x,y,z,i,nel).Be'*D*obj.strain(x,y,z,i,nel).Be*obj.t*det(obj.strain(x,y,z,i,nel).Je);
                Ke = obj.integrate(f);
                
                for k = 1: size(Ke,1)
                                if conect(k) < dof.total_dof
                                    for j = 1: size(Ke,2)
                                        if conect(j) < dof.total_dof
                                        K(conect(k),conect(j)) = Ke(k,j)+ K(conect(k),conect(j));
                                        end
                                    end
                                end
                end
            end
            v = K;
        end
        function v = forces(obj,nel,desp,dofs,fp)
            dof = obj.dof_list(dofs,desp,nel);
            dof_free = dof.dof_free;
            F = sparse(dof_free-1,1);
            Elem = obj.conect_hex(nel).Elem;
            F(dof.dof_list(fp.nodes,2),1) = F(dof.dof_list(fp.nodes,2),1) + fp.fx;
            F(dof.dof_list(fp.nodes,3),1) = F(dof.dof_list(fp.nodes,3),1) + fp.fy;
            F(dof.dof_list(fp.nodes,4),1) = F(dof.dof_list(fp.nodes,4),1) + fp.fz;
            for i = 1:size(Elem,1)
                conect = index(Elem(i,:),3,dof.dof_list);
                Fe = [];
            for j = 1:length(Elem(i,:))
                        fbe = @(x,y,z) obj.shapefun(x,y,z).f(j).*[0;0; -obj.ro*obj.g].*obj.t.*det(obj.strain(x,y,z,i,nel).Je);
                        fbe = obj.integrate(fbe);
                        Fe = [Fe;fbe];
            end
            for k = 1: size(Fe,1)
                     if conect(k) < dof.total_dof
                       F(conect(k),1) = Fe(k,1)+ F(conect(k),1);
                     end
            end
            end
            v = F;    
        end
        function obj = beam3d(E,t,poisson,ro,l,h)
            obj.t=t;
            obj.E=E;
            obj.poisson=poisson;
            obj.ro=ro;
            obj.g=9.81;
            obj.l=l;
            obj.h=h;
        end
         
     end
end