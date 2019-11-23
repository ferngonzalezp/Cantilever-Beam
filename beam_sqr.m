classdef beam_sqr
    properties
        t
        E
        poisson
        ro
        g
        state
        l
        h
    end
    methods
        function V = conect_sqr(obj,nel,tel)
            if tel~=1
                nel = nel*4;
            end
            if mod(nel,sqrt(nel))==0
                nx = sqrt(nel);
                ny = nx;
                [x,y] = meshgrid(linspace(0,obj.l,sqrt(nel)+1),linspace(0,obj.h,sqrt(nel)+1));
                x = x(:);
                y = y(:);
                Nodes = [x y];
            elseif mod(nel,2)==0
                d = divisors(nel);
                nx = d(length(d)/2+1);
                ny = d(length(d)/2);
                [x,y] = meshgrid(linspace(0,obj.l,nx+1),linspace(0,obj.h,ny+1));
                x = x(:);
                y = y(:);
                Nodes = [x y];
            end
            if tel==1
                Elem = zeros(nel,4);
                n = size(Nodes,1);
                j = 1;
                k = 1;
                for i = 1:size(Elem,1)
                    Elem(i,1) = j;
                    Elem(i,2) = j + (ny+1);
                    Elem(i,3) = Elem(i,2) + 1;
                    Elem(i,4) = Elem(i,3) - (ny+1);
                    j = j + ny + 1;
                    if j >= n - min([nx,ny]) - 1
                        j = k+1;
                        k = j;
                    end
                end
            elseif tel == 2
                Elem = zeros(nel/4,9);
                n = size(Nodes,1);
                j = 1;
                k = 1;
                for i = 1:size(Elem,1)
                    Elem(i,1) = j;
                    Elem(i,2) = j + (ny+1);
                    Elem(i,3) = Elem(i,2) + (ny+1);
                    Elem(i,4) = Elem(i,3) +1;
                    Elem(i,5) = Elem(i,4) + 1;
                    Elem(i,6) = Elem(i,5) - (ny+1);
                    Elem(i,7) = Elem(i,6) - (ny+1);
                    Elem(i,8) = Elem(i,7) - 1;
                    Elem(i,9) = Elem(i,8) + (ny+1);
                    j = j + 2*ny + 2;
                    if j >= n - min([nx,ny]) - 1
                     j = k+2;
                     k = j;
                    end
                end
            elseif tel == 3
                freenodes = zeros(1,nel/4);
                j = 1;
                k=1;
                n = size(Nodes,1);
                for i = 1:nel/4
                    freenodes(i) = j + ny + 2;
                    j = j + 2*ny+2;
                    if j >= n-ny-2
                        j = k+2;
                        k = j;
                    end
                end
                Nodes(freenodes',:)=[];
                Elem = zeros(nel/4,8);
                n = size(Nodes,1);
                j = 1;
                k = 1;
                for i = 1:size(Elem,1)
                    Elem(i,1) = j;
                    Elem(i,2) = find(Nodes(:,2) == Nodes(Elem(i,1),2) & Nodes(:,1) > Nodes(Elem(i,1),1) ,1);
                    Elem(i,3) = find(Nodes(:,2) == Nodes(Elem(i,2),2) & Nodes(:,1) > Nodes(Elem(i,2),1) ,1);
                    Elem(i,4) = Elem(i,3) + 1;
                    Elem(i,5) = Elem(i,4) + 1;
                    Elem(i,6) = Elem(i,2) + 1;
                    Elem(i,7) = Elem(i,1) + 2;
                    Elem(i,8) = Elem(i,1) + 1;
                    j = Elem(i,3);
                    if j >= n - min([nx,ny]) - 1
                     j = k+2;
                     k = j;
                    end
                end
            end
            V.Nodes = Nodes;
            V.Elem = Elem;
        end
        function V = dof_list(obj,dofs,desp,nel,tel)
            Nodes = obj.conect_sqr(nel,tel).Nodes;
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
            V.dof_list(desp(i,1),:) = [desp(i,1), a, b];
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
        end
        V.dof_free = dof_free;
        end
        function D = material(obj)
                if obj.state == 1
                    D = 1/(1-obj.poisson^2).*[obj.E obj.E*obj.poisson 0; obj.E*obj.poisson obj.E 0; 0 0 obj.E*(1-obj.poisson^2)/(2*(1+obj.poisson))];
                else
                    a = 1-obj.poisson^2; b=obj.poisson+obj.poisson^2; b=c; d=a;
                    D = 1/(a*d-b*c).*[a*obj.E b*obj.E 0; obj.E*c obj.E*d 0; 0 0 (a*d-b*c)*obj.E/(2*(1+obj.poisson))];
                end
        end
        function v = shapefun1(obj,chi,ep)
            v.f = [1/4*(1-chi)*(1-ep) 1/4*(1+chi)*(1-ep) 1/4*(1+chi)*(1+ep) 1/4*(1-chi)*(1+ep)];
            v.df = [  ep/4 - 1/4,    1/4 - ep/4,  ep/4 + 1/4, - ep/4 - 1/4; chi/4 - 1/4, - chi/4 - 1/4, chi/4 + 1/4,  1/4 - chi/4];
        end
        function v = shapefun2(obj,chi,ep)
            f = zeros(1,9);
            df = zeros(2,9);
            chii = [-1 0 1 1 1 0 -1 -1 0];
            epi = [-1 -1 -1 0 1 1 1 0 0];
            for i=1:9
                if i==1 || i==3 || i==5 || i==7
                    f(i) = 1/4*(chi^2+chi*chii(i))*(ep^2+ep*epi(i));
                    df(1,i) = 1/4*(2*chi+chii(i))*(ep^2+ep*epi(i));
                    df(2,i) = 1/4*(chi^2+chi*chii(i))*(2*ep+epi(i));
                elseif i==9
                    f(i) = (1-chi^2)*(1-ep^2);
                    df(1,i) = (-2*chi)*(1-ep^2);
                    df(2,i) = (1-chi^2)*(-2*ep);
                else
                    f(i) = 1/2*epi(i)^2*(ep^2+ep*epi(i))*(1-chi^2) + 1/2*chii(i)^2*(chi^2+chi*chii(i))*(1-ep^2);
                    df(1,i) = -chi*epi(i)^2*(ep^2+ep*epi(i)) + 1/2*chii(i)^2*(2*chi+chii(i))*(1-ep^2);
                    df(2,i) = 1/2*epi(i)^2*(2*ep+epi(i))*(1-chi^2) - ep*chii(i)^2*(chi^2+chi*chii(i));
                end
                v.f = f;
                v.df = df;
            end
        end
        function v = shapefun3(obj,chi,ep)
            f = zeros(1,8);
            df = zeros(2,8);
            chii = [-1 0 1 1 1 0 -1 -1];
            epi = [-1 -1 -1 0 1 1 1 0];
            for i=1:8
                if mod(i,2) == 0
                    if i==4 || i==8
                        f(i) = 1/2*(1+chi*chii(i))*(1-ep^2);
                        df(1,i) = 1/2*(chii(i))*(1-ep^2);
                        df(2,i) = 1/2*(1+chi*chii(i))*(-2*ep);
                    else
                        f(i) = 1/2*(1+ep*epi(i))*(1-chi^2);
                        df(1,i) =  1/2*(1+ep*epi(i))*(-2*chi);
                        df(2,i) = 1/2*(epi(i))*(1-chi^2);
                    end
                else
                    f(i) = 1/4*(1+chi*chii(i))*(1+ep*epi(i))*(chi*chii(i)+ep*epi(i)-1);  
                    df(1,i) = 1/4*(chii(i))*(1+ep*epi(i))*(chi*chii(i)+ep*epi(i)-1) + 1/4*(1+chi*chii(i))*(1+ep*epi(i))*(chii(i));
                    df(2,i) = 1/4*(1+chi*chii(i))*(epi(i))*(chi*chii(i)+ep*epi(i)-1) + 1/4*(1+chi*chii(i))*(1+ep*epi(i))*(epi(i));
                end
            end
            v.f = f;
            v.df = df;
        end
        function v = strain(obj,chi,ep,i,nel,tel)
            Nodes = obj.conect_sqr(nel,tel).Nodes;
            Elem = obj.conect_sqr(nel,tel).Elem;
            if tel == 1
                N = obj.shapefun1(chi,ep);
            elseif tel==2
                N = obj.shapefun2(chi,ep);
            elseif tel==3
                N = obj.shapefun3(chi,ep);
            end
            Je = [0 0; 0 0];
            for j = 1:length(N.f)
                Je = Je + [N.df(1,j)*Nodes(Elem(i,j),1) N.df(1,j)*Nodes(Elem(i,j),2); N.df(2,j)*Nodes(Elem(i,j),1) N.df(2,j)*Nodes(Elem(i,j),2)];
            end
            Be = [];
            for j = 1:length(N.f)
                b = Je(2,2)*N.df(1,j)-Je(1,2)*N.df(2,j);
                c = Je(1,1)*N.df(2,j)-Je(2,1)*N.df(1,j);
                Be = [Be  1/det(Je).*[b 0;0 c; c b]];
            end
            v.Be = Be;
            v.Je = Je;
        end
        function v = integrate2(obj,f)
            v = f(-sqrt(1/3),-sqrt(1/3))+f(-sqrt(1/3),sqrt(1/3))+f(sqrt(1/3),-sqrt(1/3))+f(sqrt(1/3),sqrt(1/3));
        end
        function v = integrate3(obj,f)
            w1 = 5/9; w2 = 8/9;
            q1 = -sqrt(3/5); q2 = 0; q3 = -q1;
            v = w1^2.*f(q1,q1)+w2.*w1.*f(q2,q1)+w1^2.*f(q3,q1)+w1.*w2.*f(q1,q2)+q2^2.*f(q2,q2)+w1.*w2.*f(q3,q2)+w1^2.*f(q1,q3)+w1.*w2.*f(q2,q3)+w1^2.*f(q3,q3);
        end
        function v = stiffness(obj,nel,tel,desp,dofs)
            dof = obj.dof_list(dofs,desp,nel,tel);
            dof_free = dof.dof_free;
            D = obj.material();
            K = sparse(dof_free-1, dof_free-1);
            Elem = obj.conect_sqr(nel,tel).Elem;
            for i = 1:size(Elem,1)
                conect = index(Elem(i,:),2,dof.dof_list);
                f = @(chi,ep) obj.strain(chi,ep,i,nel,tel).Be'*D*obj.strain(chi,ep,i,nel,tel).Be*obj.t*det(obj.strain(chi,ep,i,nel,tel).Je);
                if tel==1 || tel==3
                    Ke = obj.integrate2(f);
                else
                    Ke = obj.integrate3(f);
                end
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
        function v = forces(obj,nel,tel,desp,dofs,fp)
            dof = obj.dof_list(dofs,desp,nel,tel);
            dof_free = dof.dof_free;
            F = sparse(dof_free-1,1);
            Elem = obj.conect_sqr(nel,tel).Elem;
            for i = 1:size(Elem,1)
                conect = index(Elem(i,:),2,dof.dof_list);
                Fe = [];
                for j = 1:length(Elem(i,:))
                    if tel==1
                        fbe = @(chi,ep) obj.shapefun1(chi,ep).f(j).*[0; -obj.ro*obj.g].*obj.t.*det(obj.strain(chi,ep,i,nel,tel).Je);
                        fbe = obj.integrate2(fbe);
                        Fe = [Fe;fbe];
                    elseif tel==2
                        fbe = @(chi,ep) obj.shapefun2(chi,ep).f(j).*[0; -obj.ro*obj.g].*obj.t.*det(obj.strain(chi,ep,i,nel,tel).Je);
                        fbe = obj.integrate3(fbe);
                        Fe = [Fe;fbe];
                    elseif tel==3
                        fbe = @(chi,ep) obj.shapefun3(chi,ep).f(j).*[0; -obj.ro*obj.g].*obj.t.*det(obj.strain(chi,ep,i,nel,tel).Je);
                        fbe = obj.integrate2(fbe);
                        Fe = [Fe;fbe];
                    end
                end
                for k = 1: size(Fe,1)
                     if conect(k) < dof.total_dof
                       F(conect(k),1) = Fe(k,1)+ F(conect(k),1);
                     end
                end
            end
            F(dof.dof_list(fp.node,2),1) = F(dof.dof_list(fp.node,2),1) + fp.fx;
            F(dof.dof_list(fp.node,3),1) = F(dof.dof_list(fp.node,3),1) + fp.fy;
            v = F;
        end
        function obj = beam_sqr(E,t,poisson,ro,state,l,h)
            obj.t=t;
            obj.E=E;
            obj.poisson=poisson;
            obj.ro=ro;
            obj.g=9.81;
            obj.state=state;
            obj.l=l;
            obj.h=h;
        end
    end
end