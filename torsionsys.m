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
    end

    properties (GetAccess = private)
        W, L, G, h, T
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
            obj.T = params.T;
            obj = obj.setupemesh();
        end
        function obj = setupsystem(obj)
            F = -2*obj.G*pi/4;
            obj.vecF = F * ones(obj.n, 1);
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
            obj.matA = obj.matA / obj.h^2;
        end
        function obj = solvesystem(obj)
            obj.Phi = obj.matA\obj.vecF;
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
        end
        function showPhi(obj)
            figure
            xt = (1:obj.m)*obj.h;
            yt = (1:obj.m)*obj.h;
            contourf(xt,yt,obj.imgPhi, 20, 'LineColor', 'flat')
            title("Distribución \Phi")
            colorbar
            daspect([1 1 1])
            % set(gca,'Color','k')
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
                    val = sum(obj.sketch(ys, xs), "all");
                    if val > 2
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