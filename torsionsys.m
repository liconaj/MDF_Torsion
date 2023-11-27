classdef torsionsys
    %TORSIONSYS Summary of this class goes here
    %   Detailed explanation goes here
    properties (GetAccess = public)
        mesh
        meshindex
        sketch
        matA
        vecF
        Phi
        imgPhi
        theta
        tauxz
        tauyz
        tau
        J, W, L, G, h, T
    end

    properties (GetAccess = private)
        targetT
        profile_data
        n % numero de nodos
        m % tamaño matriz
    end

    methods
        function obj = torsionsys(profile_data, params)
            %TORSIONSYS Construct an instance of this class
            %   Detailed explanation goes here
            obj.profile_data = profile_data;
            obj.W = params.W;
            obj.L = params.L;
            obj.G = params.G;
            obj.h = params.h;
            obj.targetT = params.T;
            obj = obj.setupemesh();
            obj.J = obj.calcJ();
        end
        function obj = setupsystem(obj)
            obj.matA = -4*speye(obj.n);
            for y = 1:obj.m
                for x = 1:obj.m
                    ii = obj.meshindex(y, x);
                    if isnan(ii) || ii == 0
                        continue
                    end
                    ncase = obj.mesh(y-1:y+1, x-1:x+1);
                    adjcoefs = obj.getadjcoefs(ncase);
                    jj(1) = obj.meshindex(y, x-1);
                    jj(2) = obj.meshindex(y-1, x);
                    jj(3) = obj.meshindex(y, x+1);
                    jj(4) = obj.meshindex(y+1, x);
                    for idx = 1:length(adjcoefs)
                        c = adjcoefs(idx);
                        if c > 0
                            obj.matA(ii, jj(idx)) = c;
                        end
                    end
                end
            end
        end
        function obj = solvesystem(obj)
            eq = @(T) obj.calcT(T) - obj.targetT;
            Ttheta = fsolve(eq, obj.targetT, optimoptions(@fsolve, "Display", "off"));
            obj.theta = Ttheta / (obj.J * obj.G);
            obj.Phi = obj.calcPhi(obj.theta);
            obj.T = 2*sum(obj.Phi) * obj.h^2;
            obj.imgPhi = obj.mesh;
            for y = 1:obj.m
                for x = 1:obj.m
                    idx = obj.meshindex(y, x);
                    if isnan(idx) || idx == 0
                        continue
                    end
                    obj.imgPhi(y, x) = obj.Phi(idx);
                end
            end

            obj.tauxz = nan(obj.m-2);
            obj.tauyz = nan(obj.m-2);
            for y = 2:obj.m-1
                for x = 2:obj.m-1
                    if isnan(obj.imgPhi(y,x))
                        continue
                    end
                    dPhix = obj.imgPhi(y,x+1)-obj.imgPhi(y,x-1);
                    obj.tauxz(y-1,x-1) = dPhix/(2*obj.h);
                    dPhiy = obj.imgPhi(y+1,x)-obj.imgPhi(y-1,x);
                    obj.tauyz(y-1,x-1) = -dPhiy/(2*obj.h);
                end
            end
            obj.tau = (obj.tauxz.^2+obj.tauyz.^2).^0.5;
        end
        function Phi = calcPhi(obj, theta)
            F = -2*obj.G*theta;
            obj.vecF = F * ones(obj.n, 1);
            Phi = obj.matA\obj.vecF;
        end
        function T = calcT(obj, Ttheta)
            theta_guess = Ttheta / (obj.J * obj.G);
            T = 2*sum(obj.calcPhi(theta_guess)) * obj.h^2;
        end
        function showPhi(obj)
            figure
            xt = (1:obj.m)*obj.h;
            yt = (1:obj.m)*obj.h;
            contourf(xt,yt,obj.imgPhi(end:-1:1,:), 20, 'LineColor', 'flat')
            title("Distribución \Phi")
            colorbar
            daspect([1 1 1])
            % set(gca,'Color','k')
        end
        function showtau(obj)
            figure
            xt = (2:obj.m-1)*obj.h;
            yt = (2:obj.m-1)*obj.h;
            contourf(xt,yt,obj.tau(end:-1:1,:), 20, 'LineColor', 'flat')
            title("Distribución magnitud \tau")
            hc = colorbar;
            title(hc, "GPa")
            colormap turbo
            daspect([1 1 1])
        end
        function showtaucomps(obj)
            figure
            xt = (2:obj.m-1)*obj.h;
            yt = (2:obj.m-1)*obj.h;
            
            subplot(1,2,1)
            contourf(xt,yt,obj.tauxz(end:-1:1,:), 20, 'LineColor', 'flat')
            title("\tau_xz")
            hc = colorbar;
            title(hc, "GPa")
            colormap turbo
            daspect([1 1 1])
            
            subplot(1,2,2)
            contourf(xt,yt,obj.tauyz(end:-1:1,:), 20, 'LineColor', 'flat')
            title("\tau_yz")
            hc = colorbar;
            title(hc, "GPa")
            colormap turbo
            daspect([1 1 1])
        end
    end
    methods (Access = private)
        function obj = setupemesh(obj)
            size = length(obj.profile_data);
            cellsize = obj.W/size;
            scale = cellsize/obj.h;
            obj.sketch = imresize(obj.profile_data, scale, "nearest");
            obj.m = length(obj.sketch) + 1;
            rawmesh = nan(obj.m);
            for y = 1:obj.m
                for x = 1:obj.m
                    if x == 1
                        xs = x;
                    elseif x == obj.m
                        xs = x - 1;
                    else
                        xs = [x-1, x];
                    end
                    if y == 1
                        ys = y;
                    elseif y == obj.m
                        ys = y - 1;
                    else
                        ys = [y-1, y];
                    end
                    neighbors = obj.sketch(ys, xs);
                    val = sum(neighbors, "all");
                    if val > 2 % && ~any(isnan(neighbors(:)))
                        rawmesh(y,x) = 1;
                    elseif val > 0
                        rawmesh(y,x) = 0;
                    end
                end
            end
            obj.mesh = rawmesh;
            obj.meshindex = zeros(obj.m);
            obj.n = 0;
            idx = 0;
            for y = 1:obj.m
                for x = 1:obj.m
                    if ~isequal(rawmesh(y,x), 1)
                        continue
                    end
                    idx = idx + 1;
                    obj.n = obj.n + 1;
                    sumy = sum(rawmesh(y,[x-1, x+1]));
                    sumx = sum(rawmesh([y-1, y+1], x));
                    obj.mesh(y,x) = sumy + sumx;
                    obj.meshindex(y,x) = idx;
                end
            end
        end

        function J = calcJ(obj)
            J = 0;
            dA = obj.h^2;
            cx = obj.W/2;
            cy = obj.W/2;
            [sm, sn] = size(obj.sketch);
            for jj = 1:sm
                for ii = 1:sn
                    if ~obj.sketch(jj,ii)
                        continue
                    end
                    x = ii*obj.h;
                    y = jj*obj.h;
                    r2 = (x-cx)^2+(y-cy)^2;
                    J = J + r2 * dA;
                end
            end
        end
    end
    methods (Static)
        function coefs = getadjcoefs(ncase)
            ntype = ncase(2,2);
            switch ntype
                case 4
                    coefs = [1 1 1 1];
                case 3
                    coefs = torsionsys.shiftcoefs(ncase, [0 1 1 1], 0, 0);
                case 2
                    coefs = torsionsys.shiftcoefs(ncase, [1 1 0 0], 4, 1);
            end
        end
        function coefs = shiftcoefs(nodescase, defcoefs, nodeseek, seekcorner)
            % adjnodes => [1 2 3 4]
            % "side"  . 2 .   "corner"   1 . 2
            %         1 . 3              . . .
            %         . 4 .              4 . 3
            adjnodes = zeros(1, 4);
            y = 2; x = 2;
            if seekcorner
                adjnodes(1) = nodescase(y-1, x-1);
                adjnodes(2) = nodescase(y-1, x+1);
                adjnodes(3) = nodescase(y+1, x+1);
                adjnodes(4) = nodescase(y+1, x-1);
            else
                adjnodes(1) = nodescase(y, x-1);
                adjnodes(2) = nodescase(y-1, x);
                adjnodes(3) = nodescase(y, x+1);
                adjnodes(4) = nodescase(y+1, x);
            end
            shiftsteps = find(adjnodes == nodeseek) - 1;
            coefs = circshift(defcoefs, shiftsteps);
        end
    end
end